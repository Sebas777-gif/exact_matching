#!/usr/bin/env bash
# Benchmark exact matching distance (./build/play) against Hera (matching_dist)
# on every unordered pair of tests/random_bifiltrations/*.bif.
#
# Each completed pair is appended to tests/results.txt immediately, so Ctrl+C
# keeps every row that finished (and any partial row for the pair in progress).
#
# Usage:
#   ./benchmarks.sh [--no-refine|--refine] [--dim N] [--fresh]
#                   [--play PATH] [--hera PATH]
#
# If neither --no-refine nor --refine is given and stdin is a terminal, the
# script asks. Otherwise it defaults to refinement (the same default as play).

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

BIF_DIR="$ROOT/tests/random_bifiltrations"
RESULTS="$ROOT/tests/results.txt"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"
HERA_BIN="${HERA_BIN:-}"

dim=0
fresh=0
no_refine=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Compute matching distance for every unordered pair of .bif files in
tests/random_bifiltrations, once with the exact solver (play) and once
with Hera. Rows are flushed to tests/results.txt after each pair.

Options:
  --no-refine     Pass --no-refine to play (faster; may miss the true max)
  --refine        Full per-face refinement (play default)
  --dim N         Homology dimension for both solvers (default: 0)
  --fresh         Overwrite tests/results.txt instead of resuming
  --play PATH     Exact solver binary (default: build/play)
  --hera PATH     Hera matching_dist binary
  -h, --help      Show this help

If neither --no-refine nor --refine is given and the script is running
interactively, you will be asked which play mode to use.

Ctrl+C stops the current solver, writes any finished measurements for
the pair in progress, and exits. Re-run without --fresh to skip pairs
already in the results file (same dim and refine setting).
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-refine) no_refine=1; shift ;;
        --refine)    no_refine=0; shift ;;
        --dim)
            [[ $# -ge 2 ]] || { echo "error: --dim needs an integer" >&2; exit 1; }
            dim="$2"
            shift 2
            ;;
        --fresh) fresh=1; shift ;;
        --play)
            [[ $# -ge 2 ]] || { echo "error: --play needs a path" >&2; exit 1; }
            PLAY_BIN="$2"
            shift 2
            ;;
        --hera)
            [[ $# -ge 2 ]] || { echo "error: --hera needs a path" >&2; exit 1; }
            HERA_BIN="$2"
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

if [[ -z "$no_refine" ]]; then
    if [[ -t 0 ]]; then
        echo "Exact solver refinement:"
        echo "  n / --no-refine : maximize F only on A0 vertices (faster)"
        echo "  y / --refine    : per-face comparison lines (play default)"
        read -r -p "Use --no-refine? [y/N] " ans
        case "$ans" in
            y|Y|yes|YES) no_refine=1 ;;
            *)           no_refine=0 ;;
        esac
    else
        no_refine=0
    fi
fi

if [[ -z "$HERA_BIN" ]]; then
    for candidate in \
        "$HOME/hera_aida_matching/extern/hera/build/matching/matching_dist" \
        "$ROOT/../hera_aida_matching/extern/hera/build/matching/matching_dist" \
        "$ROOT/../hera/build/matching/matching_dist"
    do
        if [[ -x "$candidate" ]]; then
            HERA_BIN="$candidate"
            break
        fi
    done
    if [[ -z "$HERA_BIN" ]] && command -v matching_dist >/dev/null 2>&1; then
        HERA_BIN="$(command -v matching_dist)"
    fi
fi

if [[ ! -x "$PLAY_BIN" ]]; then
    echo "error: exact solver not found or not executable: $PLAY_BIN" >&2
    echo "Build it with: cmake -S . -B build && cmake --build build" >&2
    exit 1
fi
if [[ -z "$HERA_BIN" || ! -x "$HERA_BIN" ]]; then
    echo "error: Hera matching_dist not found. Pass --hera PATH." >&2
    exit 1
fi
if [[ ! -d "$BIF_DIR" ]]; then
    echo "error: missing directory $BIF_DIR" >&2
    exit 1
fi

shopt -s nullglob
files=( "$BIF_DIR"/*.bif )
shopt -u nullglob
n_files=${#files[@]}
if [[ "$n_files" -lt 2 ]]; then
    echo "error: need at least two .bif files in $BIF_DIR" >&2
    exit 1
fi
n_pairs=$(( n_files * (n_files - 1) / 2 ))

if [[ "$no_refine" -eq 1 ]]; then
    exact_flags="--no-refine"
else
    exact_flags="(refine)"
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/md-bench.XXXXXX")"
child_pid=""
current_a=""
current_b=""
exact_dist=""
exact_sec=""
hera_dist=""
hera_sec=""
interrupted=0

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

write_row() {
    # Always append; the file is closed after the write, so a crash or
    # Ctrl+C cannot lose earlier rows.
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$1" "$2" "$3" "$4" "$5" "$6" >> "$RESULTS"
}

flush_partial() {
    if [[ -z "$current_a" ]]; then
        return
    fi
    if [[ -z "$exact_dist" && -z "$hera_dist" ]]; then
        return
    fi
    write_row "$current_a" "$current_b" \
        "${exact_dist:-}" "${exact_sec:-}" \
        "${hera_dist:-}" "${hera_sec:-}"
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
# matching-distance benchmarks
# started: ${started}
# dim: ${dim}
# exact: ${exact_flags}
# play: ${PLAY_BIN}
# hera: ${HERA_BIN}
# files: ${n_files}  pairs: ${n_pairs}
# Ctrl+C keeps every completed row (and a partial row for the pair in progress)
file_a	file_b	exact_distance	exact_seconds	hera_distance	hera_seconds
EOF
}

header_matches() {
    local got_dim got_exact
    got_dim="$(awk -F': ' '/^# dim:/{print $2; exit}' "$RESULTS")"
    got_exact="$(awk -F': ' '/^# exact:/{print $2; exit}' "$RESULTS")"
    [[ "$got_dim" == "$dim" && "$got_exact" == "$exact_flags" ]]
}

declare -A done_pairs=()
if [[ "$fresh" -eq 1 || ! -f "$RESULTS" ]]; then
    write_header
else
    if ! header_matches; then
        echo "error: $RESULTS was produced with different dim/refine settings." >&2
        echo "Pass --fresh to overwrite, or reuse the same --dim / --no-refine." >&2
        exit 1
    fi
    while IFS=$'\t' read -r a b e_dist e_sec h_dist h_sec; do
        [[ -z "${a:-}" || "$a" == \#* || "$a" == file_a ]] && continue
        if [[ -n "${e_dist:-}" && -n "${h_dist:-}" ]]; then
            done_pairs["$a"$'\t'"$b"]=1
        fi
    done < "$RESULTS"
fi

run_timed() {
    # Run "$@" in the background so Ctrl+C is delivered to the solver.
    # Sets run_status, run_sec; stdout/stderr in $WORKDIR/{stdout,stderr}.
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

echo "Benchmarking $n_files files ($n_pairs unordered pairs)" >&2
echo "  play: $PLAY_BIN  flags: --dim $dim ${exact_flags}" >&2
echo "  hera: $HERA_BIN  flags: -d $dim" >&2
echo "  results: $RESULTS" >&2
if [[ ${#done_pairs[@]} -gt 0 ]]; then
    echo "  resuming: ${#done_pairs[@]} pair(s) already complete" >&2
fi

pair_index=0
skipped=0
for ((i = 0; i < n_files; i++)); do
    for ((j = i + 1; j < n_files; j++)); do
        pair_index=$((pair_index + 1))
        a_path="${files[i]}"
        b_path="${files[j]}"
        a_base="$(basename "$a_path")"
        b_base="$(basename "$b_path")"
        key="$a_base"$'\t'"$b_base"
        if [[ -n "${done_pairs[$key]+x}" ]]; then
            skipped=$((skipped + 1))
            continue
        fi

        echo "[$pair_index/$n_pairs] $a_base  vs  $b_base" >&2

        current_a="$a_base"
        current_b="$b_base"
        exact_dist=""
        exact_sec=""
        hera_dist=""
        hera_sec=""

        echo "  hera..." >&2
        run_timed "$HERA_BIN" -d "$dim" "$a_path" "$b_path"
        hera_sec="$run_sec"
        if [[ "$run_status" -eq 0 ]]; then
            hera_line="$(tail -n 1 "$WORKDIR/stdout" 2>/dev/null || true)"
            hera_dist="$(normalize_dist "$hera_line")"
            if [[ -z "$hera_dist" ]]; then
                hera_dist="error"
            fi
        else
            hera_dist="error"
        fi
        echo "    $hera_dist  (${hera_sec}s)" >&2

        echo "  exact..." >&2
        play_args=( "$PLAY_BIN" "$a_path" "$b_path" --dim "$dim" )
        if [[ "$no_refine" -eq 1 ]]; then
            play_args+=( --no-refine )
        fi
        run_timed "${play_args[@]}"
        exact_sec="$run_sec"
        if [[ "$run_status" -eq 0 ]]; then
            exact_line="$(grep -E '^Matching distance:' "$WORKDIR/stdout" | tail -n 1 || true)"
            exact_line="${exact_line#Matching distance: }"
            exact_dist="$(normalize_dist "$exact_line")"
            if [[ -z "$exact_dist" ]]; then
                exact_dist="error"
            fi
        else
            exact_dist="error"
        fi
        echo "    $exact_dist  (${exact_sec}s)" >&2

        write_row "$current_a" "$current_b" \
            "$exact_dist" "$exact_sec" "$hera_dist" "$hera_sec"
        current_a=""
        current_b=""
        exact_dist=""
        exact_sec=""
        hera_dist=""
        hera_sec=""
    done
done

echo "Done. $pair_index pairs considered, $skipped skipped (already complete)." >&2
echo "Results: $RESULTS" >&2
