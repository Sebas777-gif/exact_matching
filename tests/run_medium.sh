#!/usr/bin/env bash
# Medium A0 cases for timing code tweaks.
#
# play currently builds A0 and nothing more (no F, no matching distance).
# Each case reports the A0 line count, the vertex/edge/face counts, and
# the wall time. The known matching distance of each pair and, if Hera's
# matching_dist is found, Hera's value are printed next to it for
# reference only; they are not checked.
#
# Usage:
#   ./tests/run_medium.sh                 # --no-aida (default)
#   ./tests/run_medium.sh --aida          # mpfree+AIDA, .bif synthetics
#   PLAY_BIN=./build/play ./tests/run_medium.sh --aida
#   TIMEOUT_SEC=300 ./tests/run_medium.sh --aida
#   HERA_BIN=/path/to/matching_dist ./tests/run_medium.sh
#
# Default --no-aida uses dummy-grade .scc synthetics (empty relation
# columns expand A0 only). --aida cannot keep those columns, so the
# synthetics switch to .bif 1-skeleta with short extra H_0 bars at the
# same grades; mpfree+AIDA keep the dual lines and the constructed
# distances.

set -u
set -o pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PLAY_BIN="${PLAY_BIN:-$ROOT/build/play}"
HERA_BIN="${HERA_BIN:-}"
TIMEOUT_SEC="${TIMEOUT_SEC:-60}"
HERA_TIMEOUT_SEC="${HERA_TIMEOUT_SEC:-30}"
USE_AIDA=0
USE_HERA=1
failed=0
passed=0
rows=()
HEADER=$'name\tknown_dist\ta0_lines\tvertices\tedges\tfaces\tseconds\tstatus\thera\thera_seconds'

usage() {
    cat <<EOF
Usage: $0 [--aida|--no-aida] [--no-hera]

  --no-aida     skip mpfree/AIDA (default). Synthetics use dummy-grade .scc
  --aida        run play's default preprocess (mpfree then AIDA). Synthetics
                use .bif 1-skeleta so extra A0 grades survive minimization.
  --no-hera     do not run Hera for the reference distance

Environment: PLAY_BIN (default $ROOT/build/play), HERA_BIN
(matching_dist; auto-detected if unset, skipped if not found),
TIMEOUT_SEC (default 60), HERA_TIMEOUT_SEC (default 30).
EOF
}

for arg in "$@"; do
    case "$arg" in
        --aida) USE_AIDA=1 ;;
        --no-aida) USE_AIDA=0 ;;
        --no-hera) USE_HERA=0 ;;
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

if [[ ! -x "$PLAY_BIN" ]]; then
    echo "error: play not executable: $PLAY_BIN" >&2
    echo "Build it with: cmake -S . -B build && cmake --build build" >&2
    exit 1
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
if [[ "$USE_HERA" -eq 1 && ( -z "$HERA_BIN" || ! -x "$HERA_BIN" ) ]]; then
    echo "warning: Hera matching_dist not found; reference column skipped. Set HERA_BIN." >&2
    USE_HERA=0
fi

PLAY_OUT=""
HERA_OUT=""
cleanup() {
    if [[ -n "${PLAY_OUT}" && -f "${PLAY_OUT}" ]]; then
        rm -f "${PLAY_OUT}"
    fi
    if [[ -n "${HERA_OUT}" && -f "${HERA_OUT}" ]]; then
        rm -f "${HERA_OUT}"
    fi
}
print_rows() {
    echo
    printf '%s\n' "$HEADER"
    printf '%s\n' "${rows[@]}"
}
interrupted() {
    cleanup
    echo
    echo "interrupted" >&2
    if [[ ${#rows[@]} -gt 0 ]]; then
        print_rows
    fi
    exit 130
}
trap cleanup EXIT
trap interrupted INT TERM

# Foreground, not $(timeout ...): bash defers INT inside command
# substitution, so Ctrl+C would otherwise finish the case and continue.
run_cmd() {
    local outvar="$1" tsec="$2"
    shift 2
    local tmp
    tmp="$(mktemp)"
    printf -v "$outvar" '%s' "$tmp"
    if command -v timeout >/dev/null 2>&1; then
        timeout --foreground --signal=KILL "${tsec}s" \
            "$@" >"$tmp" 2>&1
    else
        "$@" >"$tmp" 2>&1
    fi
}

run_play() {
    run_cmd PLAY_OUT "$TIMEOUT_SEC" "$PLAY_BIN" "$@"
}

run_hera() {
    run_cmd HERA_OUT "$HERA_TIMEOUT_SEC" "$HERA_BIN" "$@"
}

# Hera reads 1-critical .bif only. Dummy-grade .scc synthetics have a
# matching .bif 1-skeleton; use that.
bif_for_hera() {
    local f="$1"
    local b
    case "$f" in
        *.scc)
            b="${f%.scc}.bif"
            ;;
        *.bif)
            b="$f"
            ;;
        *)
            return 1
            ;;
    esac
    [[ -f "$b" ]] || return 1
    printf '%s' "$b"
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

abort_if_interrupted() {
    local status="$1" elapsed="$2"
    case "$status" in
        130|143)
            echo "INTERRUPTED  (${elapsed}s)"
            interrupted
            ;;
    esac
}

# "lines vertices edges faces" from play's output.
a0_stats() {
    awk '/^A0 lines:/ { l = $3 }
         /^Arrangement has/ { v = $3; e = $5; f = $8 }
         END { if (l != "" && v != "") print l, v, e, f }'
}

# run_case NAME KNOWN_DIST play-args...
run_case() {
    local name="$1" known="$2"
    shift 2
    local play_args=("$@")
    local start end elapsed status stats
    local hera_dim=0 hera_a="" hera_b="" hera_dist="skip" hera_sec=""
    local hera_status=0 hera_raw="" pos=() arg i
    local lines="" verts="" edges="" faces=""

    i=0
    while [[ $i -lt ${#play_args[@]} ]]; do
        arg="${play_args[$i]}"
        case "$arg" in
            --dim)
                i=$((i + 1))
                hera_dim="${play_args[$i]}"
                ;;
            --*)
                ;;
            *)
                pos+=("$arg")
                ;;
        esac
        i=$((i + 1))
    done
    if [[ "$USE_HERA" -eq 1 && ${#pos[@]} -ge 2 ]]; then
        hera_a="$(bif_for_hera "${pos[0]}" || true)"
        hera_b="$(bif_for_hera "${pos[1]}" || true)"
    fi

    printf '%-40s' "$name"

    if [[ -n "$hera_a" && -n "$hera_b" ]]; then
        start="$(date +%s.%N 2>/dev/null || date +%s)"
        run_hera -d "$hera_dim" "$hera_a" "$hera_b"
        hera_status=$?
        end="$(date +%s.%N 2>/dev/null || date +%s)"
        hera_sec="$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.3f", e - s }')"
        abort_if_interrupted "$hera_status" "$hera_sec"
        if [[ "$hera_status" -eq 124 || "$hera_status" -eq 137 ]]; then
            hera_dist="TIMEOUT"
        elif [[ "$hera_status" -ne 0 ]]; then
            hera_dist="error"
        else
            hera_raw="$(tail -n 1 "$HERA_OUT" 2>/dev/null || true)"
            hera_dist="$(normalize_dist "$hera_raw")"
            [[ -n "$hera_dist" ]] || hera_dist="error"
        fi
        rm -f "$HERA_OUT"
        HERA_OUT=""
    fi

    start="$(date +%s.%N 2>/dev/null || date +%s)"
    run_play "${play_args[@]}"
    status=$?
    end="$(date +%s.%N 2>/dev/null || date +%s)"
    elapsed="$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.3f", e - s }')"
    abort_if_interrupted "$status" "$elapsed"
    stats="$(a0_stats < "$PLAY_OUT")"
    read -r lines verts edges faces <<<"$stats"

    local ref="known $known"
    if [[ "$hera_dist" != "skip" ]]; then
        ref="$ref, hera $hera_dist (${hera_sec}s)"
    fi

    local mark
    if [[ "$status" -eq 124 || "$status" -eq 137 ]]; then
        echo "TIMEOUT  (${elapsed}s)  $ref"
        mark="timeout"
        failed=$((failed + 1))
    elif [[ "$status" -ne 0 || -z "$stats" ]]; then
        echo "FAIL exit $status  (${elapsed}s)  $ref"
        tail -n 8 "$PLAY_OUT"
        mark="fail"
        failed=$((failed + 1))
    else
        printf 'ok  %s lines, %s vertices  (%ss)  %s\n' "$lines" "$verts" "$elapsed" "$ref"
        mark="ok"
        passed=$((passed + 1))
    fi
    rm -f "$PLAY_OUT"
    PLAY_OUT=""
    rows+=("$name	$known	$lines	$verts	$edges	$faces	$elapsed	$mark	$hera_dist	$hera_sec")
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

echo "play: $PLAY_BIN" >&2
if [[ "$USE_HERA" -eq 1 ]]; then
    echo "hera: $HERA_BIN (reference only)" >&2
else
    echo "hera: skipped" >&2
fi
echo "flags: $FLAG_DESC" >&2
echo "timeout: ${TIMEOUT_SEC}s (hera ${HERA_TIMEOUT_SEC}s)" >&2
echo

# --- synthetics, known distances ---------------------------------------
run_case "far_merge vs 140 A0 anchors" 25.5 "$FM" "$A0M" "${FLAGS[@]}"
run_case "far_merge_41 vs 140 A0 anchors" 20.5 "$FM41" "$A0M" "${FLAGS[@]}"
run_case "far_merge vs 10x10 A0 grid" 25.5 "$FM" "$A0G" "${FLAGS[@]}"

run_case "horizontal vs vertical merge L=50" 25 "$MH" "$MV" "${FLAGS[@]}"
run_case "horizontal vs vertical merge L=40" 20 "$MH40" "$MV40" "${FLAGS[@]}"
run_case "horizontal vs vertical merge L=60" 30 "$MH60" "$MV60" "${FLAGS[@]}"

# --- real 10-vertex 1-skeleta (H_0) ------------------------------------
v10() {
    local i="$1" j="$2" known="$3"
    run_case "v10 n50 sample $i vs $j" "$known" \
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

print_rows
echo
if [[ "$failed" -eq 0 ]]; then
    echo "All $passed medium cases built A0."
    exit 0
fi
echo "$passed passed, $failed failed."
exit 1
