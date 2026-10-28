#!/usr/bin/env bash
# Medium maximizer cases for manual code tweaks. Same pitfalls as the
# large random .bif pairs (unbounded A0, type-III limits, real H_0
# 1-skeleta). Default is incremental refine. Sized so --no-refine
# finishes in tens of seconds. --aida without --no-refine switches the
# synthetics to a 25-pair band calibrated for AIDA+refine at ~30-60s
# each (the 140-anchor / 8x8 / merge-70 .bif files take minutes once
# comparison lines are inserted).
#
# Usage:
#   ./tests/run_medium.sh                 # --no-aida, incremental refine
#   ./tests/run_medium.sh --no-refine     # --no-aida, A0 sample only
#   ./tests/run_medium.sh --aida          # mpfree+AIDA, incremental refine (~30-60s)
#   ./tests/run_medium.sh --aida --no-refine
#   ./tests/run_medium.sh --no-aida       # explicit default preprocess skip
#   PLAY_BIN=./build/play ./tests/run_medium.sh --aida
#   TIMEOUT_SEC=90 ./tests/run_medium.sh --aida
#
# Prints one line per case and a TSV summary.
#
# Default --no-aida uses dummy-grade .scc synthetics (empty relation
# columns expand A0 only). --aida cannot keep those columns, so the
# synthetics switch to .bif 1-skeleta with short extra H_0 bars at the
# same grades; mpfree+AIDA keep the dual lines and the constructed
# distances. Extra bars are a star into the origin: isolated pairs
# add extra infinite bars and send d_match to infinity.

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"
USE_AIDA=0
NO_REFINE=0
failed=0
passed=0
rows=()

TIMEOUT_WAS_SET=0
if [[ -n "${TIMEOUT_SEC+x}" ]]; then
    TIMEOUT_WAS_SET=1
fi

usage() {
    cat <<EOF
Usage: $0 [--aida|--no-aida] [--no-refine]

  --no-aida     skip mpfree/AIDA (default). Synthetics use dummy-grade .scc
  --aida        run play's default preprocess (mpfree then AIDA). Synthetics
                use .bif 1-skeleta so extra A0 grades survive minimization.
                With refine (default), 25 pairs sized for ~30-60s each replace
                the 140-anchor / merge-70 files, which take minutes to refine.
  --no-refine   A0 sample only. Default is incremental comparison-line refine

Environment: PLAY_BIN (default $ROOT/build/play), TIMEOUT_SEC
(default 60, or 30 with --no-refine, or 90 with --aida refine).
EOF
}

for arg in "$@"; do
    case "$arg" in
        --aida) USE_AIDA=1 ;;
        --no-aida) USE_AIDA=0 ;;
        --no-refine) NO_REFINE=1 ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "error: unknown argument: $arg" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ "$TIMEOUT_WAS_SET" -eq 0 ]]; then
    if [[ "$NO_REFINE" -eq 1 ]]; then
        TIMEOUT_SEC=30
    elif [[ "$USE_AIDA" -eq 1 ]]; then
        # 30-60s cases (merge K=8 ~54s) need slack under load.
        TIMEOUT_SEC=90
    else
        TIMEOUT_SEC=60
    fi
fi
TIMEOUT_SEC="${TIMEOUT_SEC:-60}"

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
    if [[ ${#rows[@]} -gt 0 ]]; then
        echo
        printf '%s\n' $'name\texpected\tgot\tseconds\tstatus'
        printf '%s\n' "${rows[@]}"
    fi
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
            interrupted
            ;;
    esac
}

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
    local start end elapsed status out dist mark
    printf '%-44s' "$name"
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
        rows+=("$name	$expected	TIMEOUT	$elapsed	timeout")
        failed=$((failed + 1))
        return
    fi
    if [[ "$status" -ne 0 ]]; then
        echo "FAIL exit $status  (${elapsed}s)"
        printf '%s\n' "$out" | tail -n 8
        rows+=("$name	$expected	exit$status	$elapsed	fail")
        failed=$((failed + 1))
        return
    fi
    dist="$(printf '%s\n' "$out" | awk '/^Matching distance:/{d=$3} END{print d}')"
    if dist_ok "$dist" "$expected"; then
        echo "ok  $dist  (${elapsed}s)"
        mark="ok"
        passed=$((passed + 1))
    else
        echo "FAIL got '$dist' expected $expected  (${elapsed}s)"
        mark="fail"
        failed=$((failed + 1))
    fi
    rows+=("$name	$expected	$dist	$elapsed	$mark")
}

cd "$ROOT"
BIF=tests/random_bifiltrations

if [[ "$USE_AIDA" -eq 1 ]]; then
    FM=tests/far_merge.bif
    FM41=tests/far_merge_41.bif
    A0M=tests/a0_anchors_medium.bif
    A0G=tests/a0_anchors_grid.bif
    MH=tests/merge_horizontal_medium.bif
    MV=tests/merge_vertical_medium.bif
    MH40=tests/merge_h40_medium.bif
    MV40=tests/merge_v40_medium.bif
    MH60=tests/merge_h60_medium.bif
    MV60=tests/merge_v60_medium.bif
    FLAGS=(--dim 0)
    FLAG_DESC="AIDA (mpfree+AIDA)"
else
    FM=tests/far_merge.scc
    FM41=tests/far_merge_41.scc
    A0M=tests/a0_anchors_medium.scc
    A0G=tests/a0_anchors_grid.scc
    MH=tests/merge_horizontal_medium.scc
    MV=tests/merge_vertical_medium.scc
    MH40=tests/merge_h40_medium.scc
    MV40=tests/merge_v40_medium.scc
    MH60=tests/merge_h60_medium.scc
    MV60=tests/merge_v60_medium.scc
    FLAGS=(--no-aida)
    FLAG_DESC="--no-aida"
fi
if [[ "$NO_REFINE" -eq 1 ]]; then
    FLAGS+=(--no-refine)
    FLAG_DESC="$FLAG_DESC --no-refine"
else
    FLAG_DESC="$FLAG_DESC (incremental refine)"
fi

echo "play: $PLAY_BIN" >&2
echo "flags: $FLAG_DESC" >&2
echo "timeout: ${TIMEOUT_SEC}s" >&2
echo

# --- synthetics, known distances ---------------------------------------
if [[ "$USE_AIDA" -eq 1 && "$NO_REFINE" -eq 0 ]]; then
    # AIDA + incremental refine, ~30-60s each. Extra bars are a star
    # into the origin (offset 200, step 3). Isolated pairs would add
    # extra infinite bars and send d_match to infinity. Changing L
    # keeps runtime and changes the constructed distance (L/2).
    far_for() {
        case "$1" in
            51) echo tests/far_merge.bif ;;
            *) echo "tests/far_merge_$1.bif" ;;
        esac
    }
    for L in 41 51 61 71 81; do
        expected="$(awk -v L="$L" 'BEGIN { print L / 2 }')"
        far="$(far_for "$L")"
        for K in 33 34 35; do
            run_case "far_merge_$L vs xanc_$K" "$expected" \
                "$far" "tests/aida_xanc_${K}.bif" "${FLAGS[@]}"
        done
    done
    for L in 41 51 61; do
        expected="$(awk -v L="$L" 'BEGIN { print L / 2 }')"
        far="$(far_for "$L")"
        for K in 21 22; do
            run_case "far_merge_$L vs yanc_$K" "$expected" \
                "$far" "tests/aida_yanc_${K}.bif" "${FLAGS[@]}"
        done
    done
    for L in 40 50 60 70; do
        expected="$(awk -v L="$L" 'BEGIN { print L / 2 }')"
        run_case "merge K=8 L=$L" "$expected" \
            "tests/aida_mh8_${L}.bif" "tests/aida_mv8_${L}.bif" "${FLAGS[@]}"
    done
else
    run_case "far_merge vs 140 A0 anchors" 25.5 \
        "$FM" "$A0M" "${FLAGS[@]}"
    run_case "far_merge_41 vs 140 A0 anchors" 20.5 \
        "$FM41" "$A0M" "${FLAGS[@]}"
    run_case "far_merge vs 10x10 A0 grid" 25.5 \
        "$FM" "$A0G" "${FLAGS[@]}"

    run_case "horizontal vs vertical merge L=50" 25 \
        "$MH" "$MV" "${FLAGS[@]}"
    run_case "horizontal vs vertical merge L=40" 20 \
        "$MH40" "$MV40" "${FLAGS[@]}"
    run_case "horizontal vs vertical merge L=60" 30 \
        "$MH60" "$MV60" "${FLAGS[@]}"
fi

# --- real 10-vertex 1-skeleta (H_0) ------------------------------------
v10() {
    local i="$1" j="$2" expected="$3"
    run_case "v10 n50 sample $i vs $j" "$expected" \
        "$BIF/r_bif_v_10_dim_2_n_top_simplices_50_sample_${i}.bif" \
        "$BIF/r_bif_v_10_dim_2_n_top_simplices_50_sample_${j}.bif" \
        --dim 0 "${FLAGS[@]}"
}

v10 1 2 3
v10 1 5 2.5
v10 2 8 3
v10 3 7 4.5
v10 4 9 6
v10 7 10 3
v10 8 10 3

echo
printf '%s\n' $'name\texpected\tgot\tseconds\tstatus'
printf '%s\n' "${rows[@]}"
echo
if [[ "$failed" -eq 0 ]]; then
    echo "All $passed medium cases passed."
    exit 0
fi
echo "$passed passed, $failed failed."
exit 1
