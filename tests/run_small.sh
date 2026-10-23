#!/usr/bin/env bash
# Small regression tests for maximize_F pitfalls that showed up on the
# large random .bif pairs (type-III / s→+∞ overestimates, unmatched
# infinite bars). Each case is a handful of grades and should finish in
# well under a minute.
#
# Usage:
#   ./tests/run_small.sh
#   PLAY_BIN=./build/play ./tests/run_small.sh
#
# Requires --no-aida so AIDA/mpfree are not in the loop.

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"
TIMEOUT_SEC="${TIMEOUT_SEC:-60}"
failed=0
passed=0

if [[ ! -x "$PLAY_BIN" ]]; then
    echo "error: play not executable: $PLAY_BIN" >&2
    echo "Build it with: cmake -S . -B build && cmake --build build" >&2
    exit 1
fi

run_timeout() {
    if command -v timeout >/dev/null 2>&1; then
        timeout --signal=KILL "${TIMEOUT_SEC}s" "$@"
    else
        "$@"
    fi
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
    out="$(run_timeout "$PLAY_BIN" "$@" 2>&1)"
    status=$?
    set -e
    end="$(date +%s.%N 2>/dev/null || date +%s)"
    elapsed="$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.3f", e - s }')"
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

cd "$ROOT"

# Rank-1 intervals, matching distance = L_inf of the death grades = 1.
run_case "interval scc --no-refine" 1 \
    tests/interval_long.scc tests/interval_short.scc --no-aida --no-refine
run_case "interval scc refine" 1 \
    tests/interval_long.scc tests/interval_short.scc --no-aida
run_case "interval bif --no-refine" 1 \
    tests/interval_long.bif tests/interval_short.bif --no-aida --no-refine

# Two origin vertices merging at (51,0) vs a single infinite bar, with
# extra A0 dual lines at x=10,20,30 and y=10,20,30. Those lines make
# unbounded s→+∞ wedges; a sample that leaves the face mixes those
# partitions with r=-30 and used to report 30 (or 51) instead of 25.5.
run_case "far_merge vs A0 anchors --no-refine" 25.5 \
    tests/far_merge.scc tests/a0_anchors.scc --no-aida --no-refine
run_case "far_merge vs A0 anchors refine" 25.5 \
    tests/far_merge.scc tests/a0_anchors.scc --no-aida
run_case "far_merge bif --no-refine" 25.5 \
    tests/far_merge.bif tests/a0_anchors.bif --no-aida --no-refine

# Horizontal merge at (50,0) vs vertical merge at (0,50). Type-III
# slices; the closed-form c|r+q_x| on an A0 outer ray is 50, the true
# limit (and d_match) is 25.
run_case "horizontal vs vertical merge --no-refine" 25 \
    tests/merge_horizontal.scc tests/merge_vertical.scc --no-aida --no-refine
run_case "horizontal vs vertical merge refine" 25 \
    tests/merge_horizontal.scc tests/merge_vertical.scc --no-aida

# Different numbers of infinite bars.
run_case "unmatched infinite bars" infinity \
    tests/free_one.scc tests/free_two.scc --no-aida --no-refine

# Prism surfaces: H_2 of the 2-skeleton is rank 1 vs 2 in the raw
# d-simplex presentation (no 3-simplices). Same infinity as Hera --dim 2.
run_case "prism H_2 --no-aida" infinity \
    tests/prism_1.bif tests/prism_2.bif --dim 2 --no-aida --no-refine

echo
if [[ "$failed" -eq 0 ]]; then
    echo "All $passed tests passed."
    exit 0
fi
echo "$passed passed, $failed failed."
exit 1
