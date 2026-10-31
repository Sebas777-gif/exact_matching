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
#   ./small_instance_benchmarks.sh --no-refine
#   ./small_instance_benchmarks.sh --max-grades 400
#   ./small_instance_benchmarks.sh --fresh
#
# Default is AIDA + incremental refine (play's default). Rows append to
# tests/small_instance_results.txt. Ctrl+C keeps completed rows and, if
# A0 was already printed, a partial row for the pair in progress.

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

DATA_DIR="$ROOT/data_benchmarks/small_instances"
RESULTS="$ROOT/tests/small_instance_results.txt"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"

MAX_GRADES=200
fresh=0
USE_AIDA=1
NO_REFINE=0

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
  --no-refine     A0 sample only
  --refine        Incremental comparison-line refine (default)
  --fresh         Overwrite the results file instead of resuming
  --results PATH  Results TSV (default: tests/small_instance_results.txt)
  --play PATH     Exact solver binary (default: build/play)
  -h, --help      Show this help

Ctrl+C stops the current solver and keeps completed rows. Re-run without
--fresh to skip pairs whose distance was already recorded under the same
max-grades and flags.
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
        --no-refine) NO_REFINE=1; shift ;;
        --refine) NO_REFINE=0; shift ;;
        --fresh) fresh=1; shift ;;
        --results)
            [[ $# -ge 2 ]] || { echo "error: --results needs a path" >&2; exit 1; }
            RESULTS="$2"
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

if [[ "$USE_AIDA" -eq 0 && "$NO_REFINE" -eq 1 ]]; then
    exact_flags="--no-aida --no-refine"
elif [[ "$USE_AIDA" -eq 1 && "$NO_REFINE" -eq 0 ]]; then
    exact_flags="AIDA (refine)"
elif [[ "$USE_AIDA" -eq 1 && "$NO_REFINE" -eq 1 ]]; then
    exact_flags="AIDA --no-refine"
else
    exact_flags="--no-aida (refine)"
fi

if [[ ! -x "$PLAY_BIN" ]]; then
    echo "error: exact solver not found or not executable: $PLAY_BIN" >&2
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
current_a=""
current_b=""
exact_dist=""
exact_sec=""
a0_lines=""

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

normalize_dist() {
    local raw
    raw="$(printf '%s' "$1" | tr -d '[:space:]')"
    case "$raw" in
        inf|Inf|INF|infinity|Infinity|+inf|+Inf|+INF)
            printf 'infinity'
            ;;
        *)
            printf '%s' "$raw"
            ;;
    esac
}

parse_play_stdout() {
    local line
    exact_dist=""
    a0_lines=""
    if [[ ! -f "$WORKDIR/stdout" ]]; then
        return
    fi
    line="$(grep -E '^Matching distance:' "$WORKDIR/stdout" | tail -n 1 || true)"
    line="${line#Matching distance: }"
    line="${line#"${line%%[![:space:]]*}"}"
    if [[ -n "$line" ]]; then
        exact_dist="$(normalize_dist "$line")"
    fi
    line="$(grep -E '^A0 lines:' "$WORKDIR/stdout" | tail -n 1 || true)"
    line="${line#A0 lines: }"
    line="${line#"${line%%[![:space:]]*}"}"
    a0_lines="$line"
}

write_row() {
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$1" "$2" "$3" "$4" "$5" >> "$RESULTS"
}

flush_partial() {
    if [[ -z "$current_a" ]]; then
        return
    fi
    parse_play_stdout
    if [[ -z "$exact_dist" && -z "$a0_lines" ]]; then
        return
    fi
    write_row "$current_a" "$current_b" \
        "${exact_dist:-}" "${exact_sec:-}" "${a0_lines:-}"
    current_a=""
}

on_exit() {
    local status=$?
    if [[ -n "$child_pid" ]]; then
        kill -TERM "$child_pid" 2>/dev/null || true
        wait "$child_pid" 2>/dev/null || true
        child_pid=""
    fi
    flush_partial
    rm -rf "$WORKDIR"
    return "$status"
}

on_signal() {
    interrupted=1
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
# matching-distance small_instances benchmarks
# source: data_benchmarks/small_instances
# rule: same suffix, each file has <= max_grades distinct bigrades
# started: ${started}
# max_grades: ${MAX_GRADES}
# exact: ${exact_flags}
# play: ${PLAY_BIN}
# files: ${#files[@]}  pairs: ${n_pairs}
# Ctrl+C keeps every completed row (and a partial row if A0 lines were printed)
file_a	file_b	exact_distance	exact_seconds	a0_lines
EOF
}

header_matches() {
    local got_max got_exact
    got_max="$(awk -F': ' '/^# max_grades:/{print $2; exit}' "$RESULTS")"
    got_exact="$(awk -F': ' '/^# exact:/{print $2; exit}' "$RESULTS")"
    [[ "$got_max" == "$MAX_GRADES" && "$got_exact" == "$exact_flags" ]]
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
    while IFS=$'\t' read -r a b e_dist e_sec a0; do
        [[ -z "${a:-}" || "$a" == \#* || "$a" == file_a ]] && continue
        if [[ -n "${e_dist:-}" ]]; then
            done_pairs["$a"$'\t'"$b"]=1
        fi
    done < "$RESULTS"
fi

run_timed() {
    local start end
    : > "$WORKDIR/stdout"
    : > "$WORKDIR/stderr"
    start="$(now)"
    "$@" > "$WORKDIR/stdout" 2> "$WORKDIR/stderr" &
    child_pid=$!
    wait "$child_pid"
    run_status=$?
    child_pid=""
    end="$(now)"
    run_sec="$(elapsed "$start" "$end")"
}

echo "Feasible files (<= $MAX_GRADES distinct bigrades): ${#files[@]}" >&2
for entry in "${files[@]}"; do
    IFS=$'\t' read -r name group n <<<"$entry"
    echo "  $n grades  $group  $name" >&2
done
echo "Pairs: $n_pairs   flags: $exact_flags" >&2
echo "Results: $RESULTS" >&2
if [[ ${#done_pairs[@]} -gt 0 ]]; then
    echo "Resuming: ${#done_pairs[@]} pair(s) already complete" >&2
fi

pair_index=0
skipped=0
for entry in "${pairs[@]}"; do
    pair_index=$((pair_index + 1))
    IFS=$'\t' read -r a_base b_base ga gb group <<<"$entry"
    key="$a_base"$'\t'"$b_base"
    if [[ -n "${done_pairs[$key]+x}" ]]; then
        skipped=$((skipped + 1))
        continue
    fi

    echo "[$pair_index/$n_pairs] $group  $a_base  vs  $b_base  (grades $ga + $gb)" >&2

    current_a="$a_base"
    current_b="$b_base"
    exact_dist=""
    exact_sec=""
    a0_lines=""

    play_args=( "$PLAY_BIN" "$DATA_DIR/$a_base" "$DATA_DIR/$b_base" )
    if [[ "$USE_AIDA" -eq 0 ]]; then
        play_args+=( --no-aida )
    fi
    if [[ "$NO_REFINE" -eq 1 ]]; then
        play_args+=( --no-refine )
    fi
    run_timed "${play_args[@]}"
    exact_sec="$run_sec"
    parse_play_stdout
    if [[ "$run_status" -ne 0 || -z "$exact_dist" ]]; then
        if [[ -z "$exact_dist" ]]; then
            exact_dist="error"
        fi
        echo "    $exact_dist  a0=${a0_lines:-?}  (${exact_sec}s, exit $run_status)" >&2
        if [[ -s "$WORKDIR/stderr" ]]; then
            tail -n 8 "$WORKDIR/stderr" >&2 || true
        fi
    else
        echo "    $exact_dist  a0=$a0_lines  (${exact_sec}s)" >&2
    fi

    write_row "$current_a" "$current_b" \
        "$exact_dist" "$exact_sec" "${a0_lines:-}"
    current_a=""
    current_b=""
    exact_dist=""
    exact_sec=""
    a0_lines=""
done

echo "Done. $pair_index pairs considered, $skipped skipped (already complete)." >&2
echo "Results: $RESULTS" >&2
