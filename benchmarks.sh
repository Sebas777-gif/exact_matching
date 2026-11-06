#!/usr/bin/env bash
# Build A0 (./build/play) for every unordered pair of
# tests/random_bifiltrations/*.bif, with Hera's matching_dist as a
# reference distance. play currently prints no matching distance.
#
# Each completed pair is appended to tests/results_a0.txt immediately, so
# Ctrl+C keeps every row that finished (and any partial row for the pair in
# progress). The old distance log tests/results.txt is not touched.
#
# Usage:
#   ./benchmarks.sh [--dim N] [--no-aida] [--fresh]
#                   [--play PATH] [--hera PATH] [--no-hera]

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

BIF_DIR="$ROOT/tests/random_bifiltrations"
RESULTS="$ROOT/tests/results_a0.txt"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"
HERA_BIN="${HERA_BIN:-}"

dim=0
fresh=0
USE_AIDA=1
USE_HERA=1

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Build A0 for every unordered pair of .bif files in
tests/random_bifiltrations (play) and compute Hera's matching distance
as a reference. Rows are flushed to tests/results_a0.txt after each pair.

Options:
  --dim N         Homology dimension for both tools (default: 0)
  --no-aida       Pass --no-aida to play (skip mpfree/AIDA)
  --fresh         Overwrite tests/results_a0.txt instead of resuming
  --play PATH     play binary (default: build/play)
  --hera PATH     Hera matching_dist binary
  --no-hera       Do not run Hera
  -h, --help      Show this help

Ctrl+C stops the current solver, writes any finished measurements for
the pair in progress, and exits. Re-run without --fresh to skip pairs
already in the results file (same dim and play flags).
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-aida) USE_AIDA=0; shift ;;
        --no-hera) USE_HERA=0; shift ;;
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

if [[ "$USE_HERA" -eq 1 && -z "$HERA_BIN" ]]; then
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
    echo "error: play not found or not executable: $PLAY_BIN" >&2
    echo "Build it with: cmake -S . -B build && cmake --build build" >&2
    exit 1
fi
if [[ "$USE_HERA" -eq 1 && ( -z "$HERA_BIN" || ! -x "$HERA_BIN" ) ]]; then
    echo "warning: Hera matching_dist not found; reference column skipped. Pass --hera PATH." >&2
    USE_HERA=0
fi
if [[ "$USE_HERA" -eq 0 ]]; then
    HERA_BIN="(skipped)"
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

if [[ "$USE_AIDA" -eq 1 ]]; then
    play_flags="A0 only, AIDA"
else
    play_flags="A0 only, --no-aida"
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/md-bench.XXXXXX")"
child_pid=""
current_a=""
current_b=""
a0_lines=""
a0_vertices=""
a0_edges=""
a0_faces=""
play_sec=""
hera_dist=""
hera_sec=""

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

# A0 line count and arrangement size from play's stdout.
parse_play_stdout() {
    read -r a0_lines a0_vertices a0_edges a0_faces < <(awk '
        /^A0 lines:/ { l = $3 }
        /^Arrangement has/ { v = $3; e = $5; f = $8 }
        END { print l, v, e, f }' "$WORKDIR/stdout")
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
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" >> "$RESULTS"
}

flush_partial() {
    if [[ -z "$current_a" ]]; then
        return
    fi
    if [[ -z "$a0_lines" && -z "$hera_dist" ]]; then
        return
    fi
    write_row "$current_a" "$current_b" \
        "${a0_lines:-}" "${a0_vertices:-}" "${a0_edges:-}" "${a0_faces:-}" \
        "${play_sec:-}" "${hera_dist:-}" "${hera_sec:-}"
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
# A0 benchmarks
# started: ${started}
# dim: ${dim}
# play_flags: ${play_flags}
# play: ${PLAY_BIN}
# hera: ${HERA_BIN}
# files: ${n_files}  pairs: ${n_pairs}
# Ctrl+C keeps every completed row (and a partial row for the pair in progress)
file_a	file_b	a0_lines	vertices	edges	faces	seconds	hera_distance	hera_seconds
EOF
}

header_matches() {
    local got_dim got_flags
    got_dim="$(awk -F': ' '/^# dim:/{print $2; exit}' "$RESULTS")"
    got_flags="$(awk -F': ' '/^# play_flags:/{print $2; exit}' "$RESULTS")"
    [[ "$got_dim" == "$dim" && "$got_flags" == "$play_flags" ]]
}

declare -A done_pairs=()
if [[ "$fresh" -eq 1 || ! -f "$RESULTS" ]]; then
    write_header
else
    if ! header_matches; then
        echo "error: $RESULTS was produced with different dim/play settings." >&2
        echo "Pass --fresh to overwrite, or reuse the same --dim / --no-aida." >&2
        exit 1
    fi
    while IFS=$'\t' read -r a b a0 _v _e _f _s h_dist _hs; do
        [[ -z "${a:-}" || "$a" == \#* || "$a" == file_a ]] && continue
        if [[ -n "${a0:-}" && -n "${h_dist:-}" ]]; then
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
echo "  play: $PLAY_BIN  flags: --dim $dim  ${play_flags}" >&2
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
        a0_lines=""
        a0_vertices=""
        a0_edges=""
        a0_faces=""
        play_sec=""
        hera_dist=""
        hera_sec=""

        if [[ "$USE_HERA" -eq 1 ]]; then
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
        else
            hera_dist="skip"
        fi

        echo "  A0..." >&2
        play_args=( "$PLAY_BIN" "$a_path" "$b_path" --dim "$dim" )
        if [[ "$USE_AIDA" -eq 0 ]]; then
            play_args+=( --no-aida )
        fi
        run_timed "${play_args[@]}"
        play_sec="$run_sec"
        parse_play_stdout
        if [[ "$run_status" -ne 0 || -z "$a0_lines" || -z "$a0_vertices" ]]; then
            a0_lines="error"
            echo "    error  (${play_sec}s, exit $run_status)" >&2
            tail -n 8 "$WORKDIR/stderr" >&2 || true
        else
            echo "    $a0_lines lines, $a0_vertices vertices, $a0_edges edges, $a0_faces faces  (${play_sec}s)" >&2
        fi

        write_row "$current_a" "$current_b" \
            "$a0_lines" "$a0_vertices" "$a0_edges" "$a0_faces" \
            "$play_sec" "$hera_dist" "$hera_sec"
        current_a=""
        current_b=""
        a0_lines=""
        hera_dist=""
    done
done

echo "Done. $pair_index pairs considered, $skipped skipped (already complete)." >&2
echo "Results: $RESULTS" >&2
