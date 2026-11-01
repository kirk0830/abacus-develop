#!/usr/bin/env bash
set -uo pipefail

R=/root/deepmodeling/abacus-develop
TEMPLATE=/root/abacus-refactor/temp/repro_case
RUN_CASES=/root/d5_verify/run_cases.sh
COMPARE=/root/d2verify2/compare5.py
COMPARE_NEW=/root/d2verify2/cases
MPI_BUILD="$R/build-stage3-mpi"
SERIAL_BUILD="$R/build-stage3-serial"
MPI_BIN="$MPI_BUILD/abacus_basic_para"
CTEST_REGEX='MODULE_BASE_(matrix2d|ndarray|block_cyclic|container|blas_connector)'

export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMP_NUM_THREADS=1 HWLOC_COMPONENTS=-gl
WITH_D0=0
LOG_DIR=/root/abacus-refactor/logs
CASES_DIR=$COMPARE_NEW
RUN_VARIANT=

usage() {
    cat <<'EOF'
Usage: scripts/stage3_double_build_gate.sh [--with-d0]
       [--log-dir <dir>] [--cases-dir <dir>]

Without --with-d0, run the MPI/serial builds, related ctest subsets, and
the D0 prerequisite checks. With --with-d0, also run all four D0 cases and
require SUMMARY ok=36 diff=0 missing=0.
EOF
}

while (( $# > 0 )); do
    case "$1" in
        --with-d0)
            WITH_D0=1
            shift
            ;;
        --log-dir)
            if (( $# < 2 )) || [[ -z "$2" ]]; then
                echo "ERROR: --log-dir requires a non-empty directory" >&2
                exit 2
            fi
            LOG_DIR=$2
            shift 2
            ;;
        --log-dir=*)
            LOG_DIR=${1#*=}
            if [[ -z "$LOG_DIR" ]]; then
                echo "ERROR: --log-dir requires a non-empty directory" >&2
                exit 2
            fi
            shift
            ;;
        --cases-dir)
            if (( $# < 2 )) || [[ -z "$2" ]]; then
                echo "ERROR: --cases-dir requires a non-empty directory" >&2
                exit 2
            fi
            CASES_DIR=$2
            shift 2
            ;;
        --cases-dir=*)
            CASES_DIR=${1#*=}
            if [[ -z "$CASES_DIR" ]]; then
                echo "ERROR: --cases-dir requires a non-empty directory" >&2
                exit 2
            fi
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ "$CASES_DIR" != "$COMPARE_NEW" ]]; then
    echo "ERROR: --cases-dir must match compare5.py NEW=$COMPARE_NEW" >&2
    exit 2
fi

if ! mkdir -p -- "$LOG_DIR"; then
    echo "ERROR: cannot create log directory: $LOG_DIR" >&2
    exit 2
fi

JOBS="$(nproc 2>/dev/null)"
if [[ ! "$JOBS" =~ ^[1-9][0-9]*$ ]]; then
    JOBS=1
fi

STAGE_IDS=()
STAGE_STATES=()
STAGE_NOTES=()

record_stage() {
    local stage_id=$1
    local state=$2
    local note=$3
    STAGE_IDS+=("$stage_id")
    STAGE_STATES+=("$state")
    STAGE_NOTES+=("$note")
    printf '[%s] %s: %s\n' "$stage_id" "$state" "$note"
}

run_logged() {
    local log_file=$1
    shift
    "$@" 2>&1 | tee "$log_file"
    return "$?"
}

run_logged_append() {
    local log_file=$1
    shift
    "$@" 2>&1 | tee -a "$log_file"
    return "$?"
}

log_has_error() {
    local log_file=$1
    grep -Eiq '(^|[^[:alnum:]_])(fatal[[:space:]]+)?error([^[:alnum:]_]|$)' "$log_file"
}

cleanup_run_variant() {
    if [[ -n "$RUN_VARIANT" && -f "$RUN_VARIANT" ]]; then
        rm -f -- "$RUN_VARIANT"
    fi
}
trap cleanup_run_variant EXIT

echo "[S1] START: configure and build MPI"
S1_LOG="$LOG_DIR/stage3-S1-mpi-configure-build.log"
run_logged "$S1_LOG" \
    cmake -S "$R" -B "$MPI_BUILD" \
    -DENABLE_MPI=ON \
    -DENABLE_LIBXC=ON \
    -DBUILD_TESTING=ON \
    -DCMAKE_EXE_LINKER_FLAGS="-Wl,--no-undefined"
S1_CONFIGURE_RC=$?
S1_BUILD_RC=1
if (( S1_CONFIGURE_RC == 0 )); then
    run_logged_append "$S1_LOG" cmake --build "$MPI_BUILD" -j"$JOBS"
    S1_BUILD_RC=$?
else
    printf 'Build skipped because configure failed with exit code %d.\n' \
        "$S1_CONFIGURE_RC" | tee -a "$S1_LOG"
fi
S1_LOG_ERROR=0
if log_has_error "$S1_LOG"; then
    S1_LOG_ERROR=1
fi
if (( S1_CONFIGURE_RC == 0 && S1_BUILD_RC == 0 && S1_LOG_ERROR == 0 )); then
    S1_OK=1
    record_stage "S1" "PASS" \
        "MPI configure/build succeeded (configure=0, build=0)"
else
    S1_OK=0
    record_stage "S1" "FAIL" \
        "MPI configure/build failed (configure=$S1_CONFIGURE_RC, build=$S1_BUILD_RC, error_text=$S1_LOG_ERROR)"
fi

echo
echo "[S2] START: configure and build serial"
S2_LOG="$LOG_DIR/stage3-S2-serial-configure-build.log"
run_logged "$S2_LOG" \
    cmake -S "$R" -B "$SERIAL_BUILD" \
    -DENABLE_MPI=OFF \
    -DENABLE_LIBXC=OFF \
    -DBUILD_TESTING=ON \
    -DCMAKE_EXE_LINKER_FLAGS="-Wl,--no-undefined"
S2_CONFIGURE_RC=$?
S2_BUILD_RC=1
if (( S2_CONFIGURE_RC == 0 )); then
    run_logged_append "$S2_LOG" cmake --build "$SERIAL_BUILD" -j"$JOBS"
    S2_BUILD_RC=$?
else
    printf 'Build skipped because configure failed with exit code %d.\n' \
        "$S2_CONFIGURE_RC" | tee -a "$S2_LOG"
fi
S2_LOG_ERROR=0
if log_has_error "$S2_LOG"; then
    S2_LOG_ERROR=1
fi
if (( S2_CONFIGURE_RC == 0 && S2_BUILD_RC == 0 && S2_LOG_ERROR == 0 )); then
    S2_OK=1
    record_stage "S2" "PASS" \
        "serial configure/build succeeded (configure=0, build=0)"
else
    S2_OK=0
    record_stage "S2" "FAIL" \
        "serial configure/build failed (configure=$S2_CONFIGURE_RC, build=$S2_BUILD_RC, error_text=$S2_LOG_ERROR)"
fi

echo
echo "[S3] START: related MPI ctest subset"
if (( S1_OK == 0 )); then
    record_stage "S3" "FAIL" "not run because the MPI build failed"
else
    S3_LOG="$LOG_DIR/stage3-S3-mpi-ctest.log"
    run_logged "$S3_LOG" \
        env OMP_NUM_THREADS=1 ctest --test-dir "$MPI_BUILD" \
        -N -R "$CTEST_REGEX"
    S3_DISCOVERY_RC=$?
    S3_TEST_COUNT="$(
        awk '/^Total Tests:/ { print $3; exit }' "$S3_LOG"
    )"
    if [[ ! "$S3_TEST_COUNT" =~ ^[1-9][0-9]*$ ]]; then
        S3_TEST_COUNT=0
    fi
    S3_CTEST_RC=1
    if (( S3_DISCOVERY_RC == 0 && S3_TEST_COUNT > 0 )); then
        run_logged_append "$S3_LOG" \
            env OMP_NUM_THREADS=1 ctest --test-dir "$MPI_BUILD" \
            -R "$CTEST_REGEX" --output-on-failure
        S3_CTEST_RC=$?
    fi
    if (( S3_DISCOVERY_RC == 0 && S3_CTEST_RC == 0 && S3_TEST_COUNT > 0 )); then
        record_stage "S3" "PASS" \
            "MPI related ctest passed ($S3_TEST_COUNT tests, OMP_NUM_THREADS=1)"
    else
        record_stage "S3" "FAIL" \
            "MPI related ctest failed (discovered=$S3_TEST_COUNT, discovery_rc=$S3_DISCOVERY_RC, ctest_rc=$S3_CTEST_RC)"
    fi
fi

echo
echo "[S4] START: related serial ctest subset"
if (( S2_OK == 0 )); then
    record_stage "S4" "FAIL" "not run because the serial build failed"
else
    S4_LOG="$LOG_DIR/stage3-S4-serial-ctest.log"
    run_logged "$S4_LOG" \
        env OMP_NUM_THREADS=1 ctest --test-dir "$SERIAL_BUILD" \
        -N -R "$CTEST_REGEX"
    S4_DISCOVERY_RC=$?
    S4_TEST_COUNT="$(
        awk '/^Total Tests:/ { print $3; exit }' "$S4_LOG"
    )"
    if [[ ! "$S4_TEST_COUNT" =~ ^[1-9][0-9]*$ ]]; then
        S4_TEST_COUNT=0
    fi
    S4_CTEST_RC=1
    if (( S4_DISCOVERY_RC == 0 && S4_TEST_COUNT > 0 )); then
        run_logged_append "$S4_LOG" \
            env OMP_NUM_THREADS=1 ctest --test-dir "$SERIAL_BUILD" \
            -R "$CTEST_REGEX" --output-on-failure
        S4_CTEST_RC=$?
    fi
    if (( S4_DISCOVERY_RC == 0 && S4_CTEST_RC == 0 && S4_TEST_COUNT > 0 )); then
        record_stage "S4" "PASS" \
            "serial related ctest passed ($S4_TEST_COUNT tests)"
    else
        record_stage "S4" "FAIL" \
            "serial related ctest failed (discovered=$S4_TEST_COUNT, discovery_rc=$S4_DISCOVERY_RC, ctest_rc=$S4_CTEST_RC)"
    fi
fi

echo
echo "[S5] START: D0 prerequisite checks"
S5_LOG="$LOG_DIR/stage3-S5-d0-precheck.log"
: | tee "$S5_LOG" >/dev/null
S5_OK=1

record_s5_check() {
    local state=$1
    local description=$2
    printf '  %s %s\n' "$state" "$description" | tee -a "$S5_LOG"
}

if [[ -d "$TEMPLATE" ]]; then
    record_s5_check "PASS" "template directory exists: $TEMPLATE"
else
    record_s5_check "FAIL" "template directory missing: $TEMPLATE"
    S5_OK=0
fi

if [[ -f "$RUN_CASES" && -x "$RUN_CASES" ]]; then
    record_s5_check "PASS" "run_cases.sh exists and is executable: $RUN_CASES"
else
    record_s5_check "FAIL" "run_cases.sh missing or not executable: $RUN_CASES"
    S5_OK=0
fi

if [[ -f "$COMPARE" ]]; then
    record_s5_check "PASS" "compare5.py exists: $COMPARE"
else
    record_s5_check "FAIL" "compare5.py missing: $COMPARE"
    S5_OK=0
fi

if [[ -f "$MPI_BIN" ]]; then
    record_s5_check "PASS" "MPI test binary exists: $MPI_BIN"
else
    record_s5_check "FAIL" "MPI test binary missing: $MPI_BIN"
    S5_OK=0
fi

if [[ -f "$COMPARE" ]] &&
    grep -Fq "NEW = \"$CASES_DIR\"" "$COMPARE"; then
    record_s5_check "PASS" "cases directory matches compare5.py NEW: $CASES_DIR"
else
    record_s5_check "FAIL" "cases directory does not match compare5.py NEW: $CASES_DIR"
    S5_OK=0
fi

if (( S5_OK == 1 )); then
    record_stage "S5" "PASS" "all D0 prerequisites are present"
else
    record_stage "S5" "FAIL" "one or more D0 prerequisites are missing or inconsistent"
fi

if (( WITH_D0 == 0 )); then
    record_stage "S6" "SKIP" "full D0 run not requested; use --with-d0"
elif (( S1_OK == 0 || S5_OK == 0 )); then
    record_stage "S6" "FAIL" "not run because S1 or S5 failed"
else
    echo
    echo "[S6] START: full D0 four-case comparison"
    S6_LOG="$LOG_DIR/stage3-S6-d0-run.log"
    S6_COMPARE_LOG="$LOG_DIR/stage3-S6-d0-compare.log"
    RUN_VARIANT="$(mktemp "${TMPDIR:-/tmp}/stage3_run_cases.XXXXXX")"
    S6_TEMP_RC=$?
    S6_SED_RC=1
    S6_RUN_RC=1
    S6_COMPARE_RC=1
    S6_OK=0
    S6_DIFF=
    S6_MISSING=

    if (( S6_TEMP_RC == 0 )); then
        {
            printf 'READ %s\n' "$RUN_CASES"
            printf 'WRITE %s\n' "$RUN_VARIANT"
            printf 'CASE OUTPUT DIR %s\n' "$CASES_DIR"
        } | tee "$S6_LOG"
        sed "s#^BIN=.*#BIN=$MPI_BIN#" \
            "$RUN_CASES" >"$RUN_VARIANT" 2>>"$S6_LOG"
        S6_SED_RC=$?
        if (( S6_SED_RC == 0 )) &&
            ! grep -Fq "BIN=$MPI_BIN" "$RUN_VARIANT"; then
            S6_SED_RC=1
            printf 'FAIL: temporary run_cases variant has no BIN=%s\n' \
                "$MPI_BIN" | tee -a "$S6_LOG"
        fi
    fi

    if (( S6_TEMP_RC == 0 && S6_SED_RC == 0 )); then
        run_logged_append "$S6_LOG" \
            env OMP_NUM_THREADS=1 bash "$RUN_VARIANT"
        S6_RUN_RC=$?
        run_logged "$S6_COMPARE_LOG" env OMP_NUM_THREADS=1 python3 "$COMPARE"
        S6_COMPARE_RC=$?
    fi

    S6_SUMMARY="$(
        grep -E '^SUMMARY ok=[0-9]+ diff=[0-9]+ missing=[0-9]+' \
            "$S6_COMPARE_LOG" 2>/dev/null | tail -n 1
    )"
    if [[ "$S6_SUMMARY" =~ ^SUMMARY\ ok=([0-9]+)\ diff=([0-9]+)\ missing=([0-9]+) ]]; then
        S6_OK_COUNT=${BASH_REMATCH[1]}
        S6_DIFF=${BASH_REMATCH[2]}
        S6_MISSING=${BASH_REMATCH[3]}
        if (( S6_RUN_RC == 0 && S6_COMPARE_RC == 0 &&
            S6_OK_COUNT == 36 && S6_DIFF == 0 && S6_MISSING == 0 )); then
            record_stage "S6" "PASS" \
                "D0 full comparison matched (ok=36 diff=0 missing=0)"
        else
            record_stage "S6" "FAIL" \
                "D0 comparison failed (ok=$S6_OK_COUNT diff=$S6_DIFF missing=$S6_MISSING, run_rc=$S6_RUN_RC, compare_rc=$S6_COMPARE_RC)"
        fi
    else
        record_stage "S6" "FAIL" \
            "D0 comparison produced no complete SUMMARY (run_rc=$S6_RUN_RC, compare_rc=$S6_COMPARE_RC)"
    fi
fi

echo
echo "Stage summary"
printf '%-6s %-6s %s\n' "Stage" "State" "Note"
OVERALL_RC=0
for i in "${!STAGE_IDS[@]}"; do
    printf '%-6s %-6s %s\n' \
        "${STAGE_IDS[$i]}" "${STAGE_STATES[$i]}" "${STAGE_NOTES[$i]}"
    if [[ "${STAGE_STATES[$i]}" == "FAIL" ]]; then
        OVERALL_RC=1
    fi
done

if (( OVERALL_RC == 0 )); then
    echo "Overall: PASS"
else
    echo "Overall: FAIL"
fi
exit "$OVERALL_RC"
