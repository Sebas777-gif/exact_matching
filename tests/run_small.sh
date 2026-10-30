#!/usr/bin/env bash
# Regression tests for maximize_F pitfalls that showed up on the large
# random .bif pairs (type-III / s→+∞ overestimates, unmatched infinite
# bars).
#
# Usage:
#   ./tests/run_small.sh                 # incremental refine (default)
#   ./tests/run_small.sh --no-refine     # A0 sample only
#   ./tests/run_small.sh --fast          # tiny cases only (~0.2s each)
#   ./tests/run_small.sh --fast --no-refine
#   ./tests/run_medium.sh                # larger timing set
#   PLAY_BIN=./build/play ./tests/run_small.sh
#
# Always --no-aida so AIDA/mpfree are not in the loop.

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"
FAST=0
NO_REFINE=0
failed=0
passed=0

TIMEOUT_WAS_SET=0
if [[ -n "${TIMEOUT_SEC+x}" ]]; then
    TIMEOUT_WAS_SET=1
fi

for arg in "$@"; do
    case "$arg" in
        --fast) FAST=1 ;;
        --no-refine) NO_REFINE=1 ;;
        -h|--help)
            echo "Usage: $0 [--fast] [--no-refine]"
            echo "  --fast        skip the 10–30s medium cases"
            echo "  --no-refine   A0 sample only (default is incremental refine)"
            exit 0
            ;;
        *)
            echo "error: unknown argument: $arg" >&2
            exit 1
            ;;
    esac
done

if [[ "$TIMEOUT_WAS_SET" -eq 0 ]]; then
    if [[ "$NO_REFINE" -eq 1 ]]; then
        TIMEOUT_SEC=30
    else
        TIMEOUT_SEC=300
    fi
fi
TIMEOUT_SEC="${TIMEOUT_SEC:-300}"

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

abort_if_interrupted() {
    local status="$1" elapsed="$2"
    case "$status" in
        130|143)
            echo "INTERRUPTED  (${elapsed}s)"
            exit 130
            ;;
    esac
}

# Compare a play distance string to an expected value. Integers and
# halves (25, 25.5) must match exactly; infinity is a literal.
dist_ok() {
    local got="$1" expected="$2"
    if [[ "$expected" == "infinity" ]]; then
        [[ "$got" == "infinity" ]]
        return
    fi
    awk -v g="$got" -v e="$expected" 'BEGIN {
        if (g == "" || g + 0 != g) exit 1
        d = g - e
        if (d < 0) d = -d
        exit (d > 1e-6)
    }'
}

run_case() {
    local name="$1" expected="$2"
    shift 2
    local start end elapsed status out dist
    printf '%-42s' "$name"
    start="$(date +%s.%N 2>/dev/null || date +%s)"
    set +e
    run_play "$@"
    status=$?
    set -e
    end="$(date +%s.%N 2>/dev/null || date +%s)"
    elapsed="$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.3f", e - s }')"
    abort_if_interrupted "$status" "$elapsed"
    out="$(cat "$PLAY_OUT")"
    cleanup
    PLAY_OUT=""
    if [[ "$status" -eq 124 || "$status" -eq 137 ]]; then
        echo "TIMEOUT  (${elapsed}s)"
        failed=$((failed + 1))
        return
    fi
    if [[ "$status" -ne 0 ]]; then
        echo "FAIL exit $status  (${elapsed}s)"
        printf '%s\n' "$out" | tail -n 20
        failed=$((failed + 1))
        return
    fi
    dist="$(printf '%s\n' "$out" | awk '/^Matching distance:/{d=$3} END{print d}')"
    if dist_ok "$dist" "$expected"; then
        echo "ok  $dist  (${elapsed}s)"
        passed=$((passed + 1))
    else
        echo "FAIL got '$dist' expected $expected  (${elapsed}s)"
        failed=$((failed + 1))
    fi
}

PLAY_FLAGS=(--no-aida)
if [[ "$NO_REFINE" -eq 1 ]]; then
    PLAY_FLAGS+=(--no-refine)
fi

cd "$ROOT"

echo "play: $PLAY_BIN" >&2
if [[ "$NO_REFINE" -eq 1 ]]; then
    echo "flags: --no-aida --no-refine" >&2
else
    echo "flags: --no-aida (incremental refine)" >&2
fi
echo "timeout: ${TIMEOUT_SEC}s" >&2
echo

# Rank-1 intervals, matching distance = L_inf of the death grades = 1.
run_case "interval scc" 1 \
    tests/interval_long.scc tests/interval_short.scc "${PLAY_FLAGS[@]}"
run_case "interval bif" 1 \
    tests/interval_long.bif tests/interval_short.bif "${PLAY_FLAGS[@]}"

# Two origin vertices merging at (51,0) vs a single infinite bar, with
# extra A0 dual lines at x=10,20,30 and y=10,20,30. Those lines make
# unbounded s→+∞ wedges; a sample that leaves the face mixes those
# partitions with r=-30 and used to report 30 (or 51) instead of 25.5.
run_case "far_merge vs A0 anchors" 25.5 \
    tests/far_merge.scc tests/a0_anchors.scc "${PLAY_FLAGS[@]}"
run_case "far_merge bif" 25.5 \
    tests/far_merge.bif tests/a0_anchors.bif "${PLAY_FLAGS[@]}"

# Horizontal merge at (50,0) vs vertical merge at (0,50). Type-III
# slices; the closed-form c|r+q_x| on an A0 outer ray is 50, the true
# limit (and d_match) is 25.
run_case "horizontal vs vertical merge" 25 \
    tests/merge_horizontal.scc tests/merge_vertical.scc "${PLAY_FLAGS[@]}"

# Different numbers of infinite bars.
run_case "unmatched infinite bars" infinity \
    tests/free_one.scc tests/free_two.scc "${PLAY_FLAGS[@]}"

# Prism surfaces: H_2 of the 2-skeleton is rank 1 vs 2 in the raw
# d-simplex presentation (no 3-simplices). Same infinity as Hera --dim 2.
run_case "prism H_2" infinity \
    tests/prism_1.bif tests/prism_2.bif --dim 2 "${PLAY_FLAGS[@]}"

# Two intervals in an L vs one spanning interval. True max 1.5 sits at
# (s, t) = (1, 1/2) on an A0 edge, not at an A0 vertex. Incremental refine
# finds it; --no-refine samples vertices (F ≤ 1) and one open-face point
# (~1.057).
run_ell_no_refine_miss() {
    printf '%-42s' "ell vs span --no-refine < 1.5"
    start="$(date +%s.%N 2>/dev/null || date +%s)"
    set +e
    run_play tests/ell.scc tests/span.scc --no-aida --no-refine
    status=$?
    set -e
    end="$(date +%s.%N 2>/dev/null || date +%s)"
    elapsed="$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.3f", e - s }')"
    abort_if_interrupted "$status" "$elapsed"
    out="$(cat "$PLAY_OUT")"
    cleanup
    PLAY_OUT=""
    dist="$(printf '%s\n' "$out" | awk '/^Matching distance:/{d=$3} END{print d}')"
    if [[ "$status" -eq 0 ]] && awk -v g="$dist" 'BEGIN { exit (g+0==g && g < 1.5 - 1e-6 && g >= 0) ? 0 : 1 }'; then
        echo "ok  $dist  (${elapsed}s)"
        passed=$((passed + 1))
    else
        echo "FAIL got '$dist' expected < 1.5  (${elapsed}s)"
        failed=$((failed + 1))
    fi
}

if [[ "$NO_REFINE" -eq 1 ]]; then
    run_ell_no_refine_miss
else
    run_case "ell vs span" 1.5 \
        tests/ell.scc tests/span.scc "${PLAY_FLAGS[@]}"
    run_ell_no_refine_miss
fi

if [[ "$FAST" -eq 0 ]]; then
    echo
    echo "Medium cases:"

    # Same 25.5 as far_merge vs a0_anchors, with 140 unused relation
    # grades so A0 has enough dual lines for the s→+∞ interior-point
    # bug to have room to fire. ~19s with --no-refine.
    run_case "far_merge vs 140 A0 anchors" 25.5 \
        tests/far_merge.scc tests/a0_anchors_medium.scc "${PLAY_FLAGS[@]}"

    # Type-III pair with ~140 distinct grades (70 dummy duals on each
    # side). Closed form c|r+q_x| is 50; true distance 25. ~19s with
    # --no-refine.
    run_case "horizontal vs vertical merge medium" 25 \
        tests/merge_horizontal_medium.scc tests/merge_vertical_medium.scc \
        "${PLAY_FLAGS[@]}"

    # Real 10-vertex 1-skeleta, same H_0 pipeline as the slow random
    # benchmark pairs. ~13–17s with --no-refine.
    run_case "v10 bif sample 1 vs 2" 3 \
        tests/random_bifiltrations/r_bif_v_10_dim_2_n_top_simplices_50_sample_1.bif \
        tests/random_bifiltrations/r_bif_v_10_dim_2_n_top_simplices_50_sample_2.bif \
        --dim 0 "${PLAY_FLAGS[@]}"
    run_case "v10 bif sample 1 vs 5" 2.5 \
        tests/random_bifiltrations/r_bif_v_10_dim_2_n_top_simplices_50_sample_1.bif \
        tests/random_bifiltrations/r_bif_v_10_dim_2_n_top_simplices_50_sample_5.bif \
        --dim 0 "${PLAY_FLAGS[@]}"
fi

echo
if [[ "$failed" -eq 0 ]]; then
    echo "All $passed tests passed."
    exit 0
fi
echo "$passed passed, $failed failed."
exit 1
