#!/bin/bash

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 MATRIX_NAME [float|double] [hot|cold] [1th|2th|8th]"
    exit 1
fi

# get UZP_PROJECT_ROOT from environment variable if set
UZP_PROJECT_ROOT=${UZP_PROJECT_ROOT:-"/uzp-artifact"}

MATRIX_NAME=$1
FLOAT_DATA_TYPE=$2
CACHE_MODE=$3
EXECUTION_THREADS=$4

# Colors
BYellow='\033[1;33m'      # Yellow
Color_Off='\033[0m'       # Text Reset
BackRed='\033[0;41m'      # Back Red
SCRIPT_TAG="${BYellow}[SCRIPT]${Color_Off}"
ERROR_TAG="${BackRed}[ERROR]${Color_Off}"

# Validate FLOAT_DATA_TYPE
if [[ "$FLOAT_DATA_TYPE" != "float" && "$FLOAT_DATA_TYPE" != "double" ]]; then
    echo -e "${ERROR_TAG} Invalid FLOAT_DATA_TYPE: $FLOAT_DATA_TYPE. Must be 'float' or 'double'."
    exit 1
else
    echo -e "${SCRIPT_TAG} FLOAT_DATA_TYPE: $FLOAT_DATA_TYPE"
fi

# Validate CACHE_MODE
if [[ "$CACHE_MODE" != "hot" && "$CACHE_MODE" != "cold" ]]; then
    echo -e "${ERROR_TAG} Invalid CACHE_MODE: $CACHE_MODE. Must be 'hot' or 'cold'."
    exit 1
else
    echo -e "${SCRIPT_TAG} CACHE_MODE: $CACHE_MODE"
fi

# Validate EXECUTION_THREADS
if [[ "${EXECUTION_THREADS}" != "1th" && "${EXECUTION_THREADS}" != "2th" && "${EXECUTION_THREADS}" != "8th" ]];
then
    echo -e "${ERROR_TAG} Invalid EXECUTION_THREADS: ${EXECUTION_THREADS}. Must be '1th', '2th', or '8th'."
    exit 1
else
    echo -e "${SCRIPT_TAG} EXECUTION_THREADS: ${EXECUTION_THREADS}"
fi

TEMP_DIRECTORY="/tmp"
MATRIX_FILE="${TEMP_DIRECTORY}/${MATRIX_NAME}.tar.gz"

# Run the pldi25-singularity/scripts/z_polyhedrator_manager.py
if [ ! -f "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.1d.uzp" ]; then
    python3 ../lib/scripts/z_polyhedrator_manager.py "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.mtx" "$TEMP_DIRECTORY/${MATRIX_NAME}";
fi

# Tune the UZP file
if [ -f "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.1d.uzp" ]; then
    if [ ! -f "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.tuned.uzp" ]; then
            # Build if necessary
            echo -e "${SCRIPT_TAG} Building spf_aggregator..."
	    make -C ../uzp-tuners spf_aggregator

            echo -e "${SCRIPT_TAG} Running spf_aggregator..."
            $PROJECT_ROOT/uzp-tuners/spf_aggregator "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.1d.uzp" "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.tuned.uzp"
    fi
else
    echo -e "${ERROR_TAG} File \"$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.1d.uzp\" not found."
    exit 1
fi

# Generate headers for SPFv3
python3 uzp-genex/automatic_spf_loop_generator.py "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.tuned.uzp" uzp-genex/generated_shapes.h

echo -e "${SCRIPT_TAG} Building executable..."

if [ "$FLOAT_DATA_TYPE" == "float" ]; then
    CMD_DT="-float"
else
    CMD_DT=""
fi

if [ "${CACHE_MODE}" == "hot" ]; then
    KERNEL_REPS=100;
else
    KERNEL_REPS=1;
fi

if [ "${EXECUTION_THREADS}" == "1th" ];
then
    CMD="spf_matvec"
    OMP_NUM_THREADS="1"
    PIN_CORES="10"
elif [ "${EXECUTION_THREADS}" == "2th" ];
then
    CMD="spf_matvec_omp"
    OMP_NUM_THREADS="2"
    PIN_CORES="10,12"
else
    CMD="spf_matvec_omp"
    OMP_NUM_THREADS="8"
    PIN_CORES="10,12,14,1,2,4,6,8"
fi

# Recompile and launch
make -C uzp-genex -B ${CMD} OUTPUT_TYPE=GFLOPS KERNEL_REPS=${KERNEL_REPS}

# RUNCOUNT=10
# i=0; while [[ $i -lt RUNCOUNT ]]; do
#     OMP_NUM_THREADS=${OMP_NUM_THREADS} numactl -C${PIN_CORES} hugectl --heap ./uzp-genex/${CMD} $CMD_DT "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.tuned.uzp"
#     i=$((i+1));
# done
