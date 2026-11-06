#!/usr/bin/env bash
# The 18 pairs from tests/results.txt whose old (--no-refine) distance was
# finite (the rest of random_bifiltrations are infinity in ~1-3s). Old
# times were 3-104 minutes.
#
# play currently builds A0 and nothing more (no matching distance). Each
# row stores the A0 line count, the vertex/edge/face counts, and the wall
# time, plus Hera's distance as a reference if Hera is found.
#
# Usage:
#   ./finite_benchmarks.sh                # AIDA (play default)
#   ./finite_benchmarks.sh --no-aida
#   ./finite_benchmarks.sh --fresh
#
# Results go to tests/finite_results_a0_<mode>.txt unless --results PATH
# is given. Ctrl+C keeps completed rows and any partial measurements for
# the pair in progress.

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

BIF_DIR="$ROOT/tests/random_bifiltrations"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"
HERA_BIN="${HERA_BIN:-}"

dim=0
fresh=0
USE_AIDA=1
USE_HERA=1
RESULTS=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Build A0 for the 18 finite-distance pairs from tests/results.txt (dim 0,
old --no-refine). Those were the only pairs in the first 1577 random
bifiltration comparisons that did not return infinity. play prints no
distance yet; Hera's distance is recorded as a reference.

Options:
  --no-aida       skip mpfree/AIDA
  --aida          use play's default preprocess (mpfree then AIDA, default)
  --dim N         Homology dimension (default: 0, matching results.txt)
  --fresh         Overwrite the results file instead of resuming
  --results PATH  Results TSV (default: tests/finite_results_a0_<mode>.txt)
  --play PATH     play binary (default: build/play)
  --hera PATH     Hera matching_dist binary
  --no-hera       Do not run Hera
  -h, --help      Show this help

Ctrl+C stops the current solver, writes any finished measurements for
the pair in progress, and exits. Re-run without --fresh to skip pairs
already complete in the results file (same dim and flags).
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-aida) USE_AIDA=0; shift ;;
        --aida) USE_AIDA=1; shift ;;
        --no-hera) USE_HERA=0; shift ;;
        --dim)
            [[ $# -ge 2 ]] || { echo "error: --dim needs an integer" >&2; exit 1; }
            dim="$2"
            shift 2
            ;;
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

if [[ "$USE_AIDA" -eq 1 ]]; then
    mode_tag="aida"
    play_flags="A0 only, AIDA"
else
    mode_tag="noaida"
    play_flags="A0 only, --no-aida"
fi

if [[ -z "$RESULTS" ]]; then
    RESULTS="$ROOT/tests/finite_results_a0_${mode_tag}.txt"
fi

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

# Finite dim-0 pairs from tests/results.txt (old --no-refine, started
# 2026-10-21). Ordered by that run's exact_seconds, shortest first.
# Columns: file_a, file_b, old seconds, old --no-refine distance, hera.
# The old distances come from an early version of F on A0 (before the
# unbounded-face fixes); several are far above Hera. Trust Hera's column.
pairs=(
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_100_dim_5_n_top_simplices_100_sample_9.bif	185.498	85	34.4853"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_100_dim_5_n_top_simplices_100_sample_5.bif	223.486	39.5	39.2855"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_10.bif	r_bif_v_100_dim_5_n_top_simplices_50_sample_9.bif	378.730	44.5	44.4961"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_50_dim_3_n_top_simplices_50_sample_3.bif	651.718	34.5	34.4945"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_50_dim_2_n_top_simplices_50_sample_10.bif	1578.122	35	20.9998"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_100_dim_3_n_top_simplices_100_sample_2.bif	1602.201	44.5	44.4919"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_1.bif	r_bif_v_50_dim_2_n_top_simplices_50_sample_5.bif	1636.060	35	34.987"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_50_dim_2_n_top_simplices_50_sample_4.bif	1721.101	29	28.9896"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_2.bif	r_bif_v_50_dim_2_n_top_simplices_50_sample_5.bif	1779.539	30	25.4931"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_50_dim_2_n_top_simplices_50_sample_6.bif	1867.189	34	25.9974"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_100_dim_3_n_top_simplices_100_sample_1.bif	2213.119	41	34.4977"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_100_dim_3_n_top_simplices_100_sample_8.bif	2237.961	41.5	30.4884"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_1.bif	r_bif_v_100_dim_3_n_top_simplices_100_sample_3.bif	2354.656	41.5	41.4983"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_2.bif	r_bif_v_100_dim_3_n_top_simplices_100_sample_3.bif	2582.402	38	37.9846"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_100_dim_3_n_top_simplices_100_sample_9.bif	2584.192	23.5	23.4991"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_100_dim_3_n_top_simplices_100_sample_4.bif	2959.544	31	30.9924"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_3.bif	r_bif_v_100_dim_2_n_top_simplices_150_sample_10.bif	4549.712	34.5	24.9894"
    "r_bif_v_100_dim_2_n_top_simplices_100_sample_1.bif	r_bif_v_100_dim_2_n_top_simplices_100_sample_2.bif	6243.116	46	45.9965"
)
n_pairs=${#pairs[@]}

missing=0
for entry in "${pairs[@]}"; do
    IFS=$'\t' read -r a_base b_base _rest <<<"$entry"
    for f in "$a_base" "$b_base"; do
        if [[ ! -f "$BIF_DIR/$f" ]]; then
            echo "error: missing $BIF_DIR/$f" >&2
            missing=1
        fi
    done
done
if [[ "$missing" -ne 0 ]]; then
    exit 1
fi

WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/md-finite-bench.XXXXXX")"
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
# A0 finite-pair benchmarks
# source: 18 finite rows of tests/results.txt (dim 0, old --no-refine)
# started: ${started}
# dim: ${dim}
# play_flags: ${play_flags}
# play: ${PLAY_BIN}
# hera: ${HERA_BIN}
# pairs: ${n_pairs}
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
        echo "Pass --fresh to overwrite, or reuse the same flags / --results." >&2
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

echo "Finite-pair benchmark: $n_pairs pairs from tests/results.txt" >&2
echo "  play: $PLAY_BIN  flags: --dim $dim  ${play_flags}" >&2
echo "  hera: $HERA_BIN  flags: -d $dim" >&2
echo "  results: $RESULTS" >&2
if [[ ${#done_pairs[@]} -gt 0 ]]; then
    echo "  resuming: ${#done_pairs[@]} pair(s) already complete" >&2
fi

pair_index=0
skipped=0
for entry in "${pairs[@]}"; do
    pair_index=$((pair_index + 1))
    IFS=$'\t' read -r a_base b_base old_sec old_exact old_hera <<<"$entry"
    a_path="$BIF_DIR/$a_base"
    b_path="$BIF_DIR/$b_base"
    key="$a_base"$'\t'"$b_base"
    if [[ -n "${done_pairs[$key]+x}" ]]; then
        skipped=$((skipped + 1))
        continue
    fi

    echo "[$pair_index/$n_pairs] $a_base  vs  $b_base" >&2
    echo "  old --no-refine: $old_exact in ${old_sec}s  (hera $old_hera)" >&2

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

echo "Done. $pair_index pairs considered, $skipped skipped (already complete)." >&2
echo "Results: $RESULTS" >&2
