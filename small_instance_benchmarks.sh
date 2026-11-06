#!/usr/bin/env bash
# Pairs from data_benchmarks/small_instances that can fit in A0.
#
# "Small" in that archive means small for mpfree/AIDA. Matching distance
# builds one dual line per distinct bigrade. cube2_density_05000 has
# ~14000 grades per file and exhausts RAM. A file is feasible when it has
# at most MAX_GRADES distinct bigrades (default 200), so a pair has at
# most ~330 A0 lines (~5e4 vertices). Pairs stay inside one suffix:
# .scc_k with .scc_k, .minpres_k with .minpres_k, plain .scc with .scc.
# The trailing index is reflected homology (.scc_2 is H_0), not --dim.
#
#   ./small_instance_benchmarks.sh
#   ./small_instance_benchmarks.sh --no-aida
#   ./small_instance_benchmarks.sh --max-grades 400
#   ./small_instance_benchmarks.sh --fresh
#
# play currently builds A0 and nothing more (no matching distance). Each
# row stores the A0 line count, the vertex/edge/face counts of the
# arrangement, and the wall time (AIDA preprocessing included unless
# --no-aida). Results go to tests/small_instance_results_a0.txt. The old
# distance logs (small_instance_results.txt, *_refine.txt) are not touched.
# Ctrl+N skips the current pair and records ??? for its measurement
# columns. Ctrl+C keeps completed rows.

set -u
set -o pipefail
# Each solver is its own process group, so Ctrl+N can stop that pair
# without stopping this script.
set -m

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

DATA_DIR="$ROOT/data_benchmarks/small_instances"
RESULTS=""
RESULTS_WAS_SET=0
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"

MAX_GRADES=200
fresh=0
USE_AIDA=1

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Unordered pairs of feasible files in data_benchmarks/small_instances.
A file is feasible when it has at most --max-grades distinct bigrades
(default 200). Pairs are only within the same suffix (.scc_k,
.minpres_k, or plain .scc). Ordered by the sum of those grade counts,
smallest first.

Options:
  --max-grades N  Per-file unique-grade cap (default: 200)
  --no-aida       Skip AIDA (presentations are already minimal)
  --aida          Run AIDA (default)
  --fresh         Overwrite the results file instead of resuming
  --results PATH  Results TSV (default: tests/small_instance_results_a0.txt)
  --play PATH     play binary (default: build/play)
  -h, --help      Show this help

Ctrl+N skips the pair that is running and writes ??? in its measurement
columns. That row is kept, so a later run does not retry it. Ctrl+C stops
the suite and keeps completed rows. Re-run without --fresh to skip pairs
already recorded under the same max-grades and flags.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --max-grades)
            [[ $# -ge 2 ]] || { echo "error: --max-grades needs an integer" >&2; exit 1; }
            MAX_GRADES="$2"
            shift 2
            ;;
        --no-aida) USE_AIDA=0; shift ;;
        --aida) USE_AIDA=1; shift ;;
        --fresh) fresh=1; shift ;;
        --results)
            [[ $# -ge 2 ]] || { echo "error: --results needs a path" >&2; exit 1; }
            RESULTS="$2"
            RESULTS_WAS_SET=1
            shift 2
            ;;
        --play)
            [[ $# -ge 2 ]] || { echo "error: --play needs a path" >&2; exit 1; }
            PLAY_BIN="$2"
            shift 2
            ;;
        -h|--help) usage; exit 0 ;;
        *)
            echo "error: unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

case "$MAX_GRADES" in
    ''|*[!0-9]*) echo "error: --max-grades must be a non-negative integer" >&2; exit 1 ;;
esac

if [[ "$USE_AIDA" -eq 1 ]]; then
    play_flags="A0 only, AIDA"
else
    play_flags="A0 only, --no-aida"
fi

if [[ "$RESULTS_WAS_SET" -eq 0 ]]; then
    RESULTS="$ROOT/tests/small_instance_results_a0.txt"
fi

if [[ ! -x "$PLAY_BIN" ]]; then
    echo "error: play not found or not executable: $PLAY_BIN" >&2
    echo "Build it with: cmake -S . -B build && cmake --build build" >&2
    exit 1
fi
if [[ ! -d "$DATA_DIR" ]]; then
    echo "error: missing directory $DATA_DIR" >&2
    exit 1
fi

# name, group, unique grades. Group is the suffix class so H_0 is not
# paired with H_1. Sorted by grade count.
mapfile -t files < <(python3 - "$DATA_DIR" "$MAX_GRADES" << 'PY'
import sys
from pathlib import Path

data = Path(sys.argv[1])
cap = int(sys.argv[2])

def group_of(name):
    if ".minpres_" in name:
        return "minpres_" + name.rsplit("_", 1)[-1]
    if ".scc_" in name:
        return "scc_" + name.rsplit("_", 1)[-1]
    if name.endswith(".scc"):
        return "scc"
    return "other"

def unique_grades(path):
    lines = path.read_text(errors="replace").splitlines()
    try:
        i = next(k for k, line in enumerate(lines) if line.strip() == "scc2020")
    except StopIteration:
        return None
    grades = set()
    for line in lines[i + 3:]:
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if line == "scc2020":
            break
        parts = line.split(";", 1)[0].split()
        if len(parts) >= 2:
            grades.add((parts[0], parts[1]))
    return len(grades)

rows = []
for path in sorted(data.iterdir()):
    if not path.is_file():
        continue
    n = unique_grades(path)
    if n is None or n > cap:
        continue
    rows.append((n, group_of(path.name), path.name))
rows.sort()
for n, group, name in rows:
    print(f"{name}\t{group}\t{n}")
PY
)

if [[ ${#files[@]} -eq 0 ]]; then
    echo "error: no feasible files under $DATA_DIR with <= $MAX_GRADES grades" >&2
    exit 1
fi

# Unordered pairs within a group, smallest grade-sum first.
mapfile -t pairs < <(printf '%s\n' "${files[@]}" | python3 -c '
import sys
from collections import defaultdict
files = []
for line in sys.stdin:
    line = line.rstrip("\n")
    if not line:
        continue
    name, group, n = line.split("\t")
    files.append((group, name, int(n)))
by = defaultdict(list)
for group, name, n in files:
    by[group].append((name, n))
pairs = []
for group, items in by.items():
    items.sort()
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, na = items[i]
            b, nb = items[j]
            pairs.append((na + nb, a, b, na, nb, group))
pairs.sort()
for total, a, b, na, nb, group in pairs:
    print(f"{a}\t{b}\t{na}\t{nb}\t{group}")
')

n_pairs=${#pairs[@]}
if [[ "$n_pairs" -eq 0 ]]; then
    echo "error: feasible files do not form a same-suffix pair" >&2
    exit 1
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/md-small-inst.XXXXXX")"
child_pid=""
play_sec=""
a0_lines=""
a0_vertices=""
a0_edges=""
a0_faces=""
key_reader_pid=""
old_stty=""

now() {
    if [[ -n "${EPOCHREALTIME:-}" ]]; then
        printf '%s' "$EPOCHREALTIME"
    else
        date +%s.%N
    fi
}

elapsed() {
    awk -v s="$1" -v e="$2" 'BEGIN { printf "%.3f", e - s }'
}

# A0 line count and arrangement size from play's stdout. Empty if absent.
parse_play_stdout() {
    a0_lines=""
    a0_vertices=""
    a0_edges=""
    a0_faces=""
    if [[ ! -f "$WORKDIR/stdout" ]]; then
        return
    fi
    read -r a0_lines a0_vertices a0_edges a0_faces < <(awk '
        /^A0 lines:/ { l = $3 }
        /^Arrangement has/ { v = $3; e = $5; f = $8 }
        END { print l, v, e, f }' "$WORKDIR/stdout")
}

write_row() {
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$1" "$2" "$3" "$4" "$5" "$6" "$7" >> "$RESULTS"
}

# Stop play and anything it spawned (AIDA, mpfree).
kill_child() {
    local pid="$1"
    [[ -n "$pid" ]] || return 0
    kill -TERM -- -"$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null || true
}

restore_tty() {
    if [[ -n "$old_stty" ]]; then
        stty "$old_stty" < /dev/tty 2>/dev/null || true
        old_stty=""
    fi
}

stop_key_reader() {
    if [[ -n "$key_reader_pid" ]]; then
        kill -TERM "$key_reader_pid" 2>/dev/null || true
        wait "$key_reader_pid" 2>/dev/null || true
        key_reader_pid=""
    fi
    restore_tty
}

# Ctrl+N is delivered as a key, not a signal. Read it from the terminal
# while play runs. A background read would otherwise get SIGTTIN.
start_key_reader() {
    [[ -r /dev/tty && -w /dev/tty ]] || return 0
    old_stty="$(stty -g < /dev/tty 2>/dev/null || true)"
    [[ -n "$old_stty" ]] || return 0
    stty -icanon -echo min 1 time 0 < /dev/tty || { old_stty=""; return 0; }
    # Same process group as this script. A background job of its own would
    # be stopped with SIGTTIN as soon as it read the terminal.
    set +m
    (
        trap '' INT TTIN
        while IFS= read -r -n 1 key < /dev/tty; do
            [[ "$key" == $'\x0e' ]] || continue
            : > "$WORKDIR/skip"
            if [[ -f "$WORKDIR/child_pid" ]]; then
                pid="$(cat "$WORKDIR/child_pid" 2>/dev/null || true)"
                [[ -n "$pid" ]] || continue
                kill -TERM -- -"$pid" 2>/dev/null || kill -TERM "$pid" 2>/dev/null || true
            fi
        done
    ) &
    key_reader_pid=$!
    set -m
}

on_exit() {
    local status=$?
    if [[ -n "$child_pid" ]]; then
        kill_child "$child_pid"
        wait "$child_pid" 2>/dev/null || true
        child_pid=""
    fi
    stop_key_reader
    rm -rf "$WORKDIR"
    return "$status"
}

on_signal() {
    echo >&2
    echo "Interrupted. Results so far are in $RESULTS" >&2
    exit 130
}

trap on_exit EXIT
trap on_signal INT TERM

write_header() {
    local started
    started="$(date -Iseconds 2>/dev/null || date)"
    cat > "$RESULTS" <<EOF
# A0 small_instances benchmarks
# source: data_benchmarks/small_instances
# rule: same suffix, each file has <= max_grades distinct bigrades
# started: ${started}
# max_grades: ${MAX_GRADES}
# play_flags: ${play_flags}
# play: ${PLAY_BIN}
# files: ${#files[@]}  pairs: ${n_pairs}
# Ctrl+N records ??? for the current pair and continues
# Ctrl+C keeps every completed row
file_a	file_b	a0_lines	vertices	edges	faces	seconds
EOF
}

header_matches() {
    local got_max got_flags
    got_max="$(awk -F': ' '/^# max_grades:/{print $2; exit}' "$RESULTS")"
    got_flags="$(awk -F': ' '/^# play_flags:/{print $2; exit}' "$RESULTS")"
    [[ "$got_max" == "$MAX_GRADES" && "$got_flags" == "$play_flags" ]]
}

declare -A done_pairs=()
if [[ "$fresh" -eq 1 || ! -f "$RESULTS" ]]; then
    write_header
else
    if ! header_matches; then
        echo "error: $RESULTS was produced with different max-grades or flags." >&2
        echo "Pass --fresh to overwrite, or reuse the same options." >&2
        exit 1
    fi
    while IFS=$'\t' read -r a b a0 _rest; do
        [[ -z "${a:-}" || "$a" == \#* || "$a" == file_a ]] && continue
        if [[ -n "${a0:-}" ]]; then
            done_pairs["$a"$'\t'"$b"]=1
        fi
    done < "$RESULTS"
fi

run_timed() {
    local start end
    : > "$WORKDIR/stdout"
    : > "$WORKDIR/stderr"
    rm -f "$WORKDIR/skip" "$WORKDIR/child_pid"
    start="$(now)"
    "$@" > "$WORKDIR/stdout" 2> "$WORKDIR/stderr" &
    child_pid=$!
    printf '%s\n' "$child_pid" > "$WORKDIR/child_pid"
    wait "$child_pid"
    run_status=$?
    if [[ -f "$WORKDIR/skip" ]]; then
        kill_child "$child_pid"
        wait "$child_pid" 2>/dev/null || true
        skip_pair=1
    fi
    rm -f "$WORKDIR/child_pid"
    child_pid=""
    end="$(now)"
    run_sec="$(elapsed "$start" "$end")"
}

echo "Feasible files (<= $MAX_GRADES distinct bigrades): ${#files[@]}" >&2
for entry in "${files[@]}"; do
    IFS=$'\t' read -r name group n <<<"$entry"
    echo "  $n grades  $group  $name" >&2
done
echo "Pairs: $n_pairs   flags: $play_flags" >&2
echo "Results: $RESULTS" >&2
echo "Ctrl+N skips the current pair (recorded as ???). Ctrl+C stops the suite." >&2
start_key_reader
if [[ ${#done_pairs[@]} -gt 0 ]]; then
    echo "Resuming: ${#done_pairs[@]} pair(s) already complete" >&2
fi

pair_index=0
skipped=0
user_skipped=0
skip_pair=0
for entry in "${pairs[@]}"; do
    pair_index=$((pair_index + 1))
    IFS=$'\t' read -r a_base b_base ga gb group <<<"$entry"
    key="$a_base"$'\t'"$b_base"
    if [[ -n "${done_pairs[$key]+x}" ]]; then
        skipped=$((skipped + 1))
        continue
    fi

    echo "[$pair_index/$n_pairs] $group  $a_base  vs  $b_base  (grades $ga + $gb)" >&2

    play_args=( "$PLAY_BIN" "$DATA_DIR/$a_base" "$DATA_DIR/$b_base" )
    if [[ "$USE_AIDA" -eq 0 ]]; then
        play_args+=( --no-aida )
    fi
    skip_pair=0
    run_timed "${play_args[@]}"
    play_sec="$run_sec"
    if [[ "$skip_pair" -eq 1 ]]; then
        echo "    skipped  ???" >&2
        write_row "$a_base" "$b_base" '???' '???' '???' '???' '???'
        skip_pair=0
        user_skipped=$((user_skipped + 1))
        continue
    fi
    parse_play_stdout
    if [[ "$run_status" -ne 0 || -z "$a0_lines" || -z "$a0_vertices" ]]; then
        echo "    error  (${play_sec}s, exit $run_status)" >&2
        if [[ -s "$WORKDIR/stderr" ]]; then
            tail -n 8 "$WORKDIR/stderr" >&2 || true
        fi
        write_row "$a_base" "$b_base" error "" "" "" "$play_sec"
    else
        echo "    $a0_lines lines, $a0_vertices vertices, $a0_edges edges, $a0_faces faces  (${play_sec}s)" >&2
        write_row "$a_base" "$b_base" \
            "$a0_lines" "$a0_vertices" "$a0_edges" "$a0_faces" "$play_sec"
    fi
done

echo "Done. $pair_index pairs considered, $skipped skipped (already complete), $user_skipped skipped with Ctrl+N." >&2
echo "Results: $RESULTS" >&2
