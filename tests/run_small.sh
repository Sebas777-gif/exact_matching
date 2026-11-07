#!/usr/bin/env bash
# Regression tests for the A0 arrangement.
#
# play currently builds A0 and nothing more: no evaluation of F and no
# matching distance. Each case checks the number of A0 lines and the
# vertex/edge/face counts of the CGAL arrangement. The matching distance
# each pair was built for is kept in the comments, for when F comes back.
#
# Usage:
#   ./tests/run_small.sh                 # all cases
#   ./tests/run_small.sh --fast          # tiny cases only
#   PLAY_BIN=./build/play ./tests/run_small.sh
#
# Always --no-aida so AIDA/mpfree are not in the loop.

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"
TIMEOUT_SEC="${TIMEOUT_SEC:-30}"
FAST=0
failed=0
passed=0

for arg in "$@"; do
    case "$arg" in
        --fast) FAST=1 ;;
        -h|--help)
            echo "Usage: $0 [--fast]"
            echo "  --fast   skip the medium cases (~140 A0 lines)"
            exit 0
            ;;
        *)
            echo "error: unknown argument: $arg" >&2
            exit 1
            ;;
    esac
done

if [[ ! -x "$PLAY_BIN" ]]; then
    echo "error: play not executable: $PLAY_BIN" >&2
    echo "Build it with: cmake -S . -B build && cmake --build build" >&2
    exit 1
fi

PLAY_OUT=""
cleanup() {
    if [[ -n "${PLAY_OUT}" && -f "${PLAY_OUT}" ]]; then
        rm -f "${PLAY_OUT}"
    fi
}
interrupted() {
    cleanup
    echo
    echo "interrupted" >&2
    exit 130
}
trap cleanup EXIT
trap interrupted INT TERM

# Foreground, not $(timeout ...): bash defers INT inside command
# substitution, so Ctrl+C would otherwise finish the case and continue.
run_play() {
    PLAY_OUT="$(mktemp)"
    if command -v timeout >/dev/null 2>&1; then
        timeout --foreground --signal=KILL "${TIMEOUT_SEC}s" \
            "$PLAY_BIN" "$@" >"$PLAY_OUT" 2>&1
    else
        "$PLAY_BIN" "$@" >"$PLAY_OUT" 2>&1
    fi
}

# "lines vertices edges faces" from play's output.
a0_stats() {
    awk '/^A0 lines:/ { l = $3 }
         /^Arrangement has/ { v = $3; e = $5; f = $8 }
         END { if (l != "" && v != "") print l, v, e, f }'
}

# run_case NAME "LINES VERTICES EDGES FACES" play-args...
run_case() {
    local name="$1" expected="$2"
    shift 2
    local start end elapsed status out got
    printf '%-42s' "$name"
    start="$(date +%s.%N 2>/dev/null || date +%s)"
    run_play "$@"
    status=$?
    end="$(date +%s.%N 2>/dev/null || date +%s)"
    elapsed="$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.3f", e - s }')"
    out="$(cat "$PLAY_OUT")"
    cleanup
    PLAY_OUT=""
    case "$status" in
        130|143)
            echo "INTERRUPTED  (${elapsed}s)"
            exit 130
            ;;
        124|137)
            echo "TIMEOUT  (${elapsed}s)"
            failed=$((failed + 1))
            return
            ;;
        0) ;;
        *)
            echo "FAIL exit $status  (${elapsed}s)"
            printf '%s\n' "$out" | tail -n 20
            failed=$((failed + 1))
            return
            ;;
    esac
    got="$(printf '%s\n' "$out" | a0_stats)"
    if [[ "$got" == "$expected" ]]; then
        echo "ok  $got  (${elapsed}s)"
        passed=$((passed + 1))
    else
        echo "FAIL got '$got' expected '$expected'  (${elapsed}s)"
        failed=$((failed + 1))
    fi
}

cd "$ROOT"

echo "play: $PLAY_BIN" >&2
echo "flags: --no-aida" >&2
echo "timeout: ${TIMEOUT_SEC}s" >&2
echo "columns: A0 lines, vertices, edges, faces" >&2
echo

# Rank-1 intervals (d_match = 1). Grades (0,0), (1,1), (2,2): all three
# duals and s = 1 meet at (s, t) = (1, 0).
run_case "interval scc" "5 4 15 12" \
    tests/interval_long.scc tests/interval_short.scc --no-aida
run_case "interval bif" "5 4 15 12" \
    tests/interval_long.bif tests/interval_short.bif --no-aida

# Two origin vertices merging at (51,0) vs a single infinite bar
# (d_match = 25.5). The .scc version adds unused grades at x, y in
# {10, 20, 30} that only cut A0; the .bif version does not.
run_case "far_merge vs A0 anchors" "10 24 62 39" \
    tests/far_merge.scc tests/a0_anchors.scc --no-aida
run_case "far_merge bif" "5 4 15 12" \
    tests/far_merge.bif tests/a0_anchors.bif --no-aida

# Horizontal merge at (50,0) vs vertical merge at (0,50) (d_match = 25).
run_case "horizontal vs vertical merge" "5 6 18 13" \
    tests/merge_horizontal.scc tests/merge_vertical.scc --no-aida

# Different numbers of infinite bars (d_match = infinity).
run_case "unmatched infinite bars" "4 3 11 9" \
    tests/free_one.scc tests/free_two.scc --no-aida

# Raw H_2 presentations of two prism 2-skeleta (d_match = infinity).
run_case "prism H_2" "7 12 34 23" \
    tests/prism_1.bif tests/prism_2.bif --dim 2 --no-aida

# L of two intervals vs one spanning interval (d_match = 1.5). The
# maximum of F sits at (s, t) = (1, 1/2), on an A0 edge, not an A0 vertex.
run_case "ell vs span" "7 9 29 21" \
    tests/ell.scc tests/span.scc --no-aida

if [[ "$FAST" -eq 0 ]]; then
    echo
    echo "Medium cases:"

    # far_merge with 140 unused relation grades (d_match = 25.5).
    run_case "far_merge vs 140 A0 anchors" "144 5183 10581 5399" \
        tests/far_merge.scc tests/a0_anchors_medium.scc --no-aida

    # Merge pair with ~140 distinct grades (d_match = 25).
    run_case "horizontal vs vertical merge medium" "145 5493 11222 5730" \
        tests/merge_horizontal_medium.scc tests/merge_vertical_medium.scc \
        --no-aida

    # Real 10-vertex 1-skeleta, H_0 (d_match = 3 and 2.5).
    BIF=tests/random_bifiltrations
    run_case "v10 bif sample 1 vs 2" "74 1966 4221 2256" \
        "$BIF/r_bif_v_10_dim_2_n_top_simplices_50_sample_1.bif" \
        "$BIF/r_bif_v_10_dim_2_n_top_simplices_50_sample_2.bif" \
        --dim 0 --no-aida
    run_case "v10 bif sample 1 vs 5" "71 1806 3883 2078" \
        "$BIF/r_bif_v_10_dim_2_n_top_simplices_50_sample_1.bif" \
        "$BIF/r_bif_v_10_dim_2_n_top_simplices_50_sample_5.bif" \
        --dim 0 --no-aida
fi

echo
if [[ "$failed" -eq 0 ]]; then
    echo "All $passed tests passed."
    exit 0
fi
echo "$passed passed, $failed failed."
exit 1
