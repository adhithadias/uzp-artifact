#!/bin/bash

# ----------------------------------------------------------------------------------------------------
# Where is Intel OneAPI?
export ONEAPIROOT=/opt/intel/oneapi/
export MKLROOT=${ONEAPIROOT}/mkl/latest/
# ----------------------------------------------------------------------------------------------------

if [ "$#" -lt 5 ]; then
    echo "Usage: $0 GROUP MATRIX_NAME [float|double] [hot|cold] [1th|2th|8th]"
    exit 1
fi

GROUP=$1
MATRIX_NAME=$2
FLOAT_DATA_TYPE=$3
CACHE_MODE=$4
EXECUTION_THREADS=$5

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

# Check if the matrix file already exists in the temporary directory
if [ ! -f "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.rb" ]; then
    echo -e "${SCRIPT_TAG} Matrix file \"$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.rb\" not found. Downloading..."
    # Download the matrix file using wget
    wget -O "$MATRIX_FILE" "https://suitesparse-collection-website.herokuapp.com/RB/${GROUP}/${MATRIX_NAME}.tar.gz"

    echo -e "${SCRIPT_TAG} Extracting file \"$MATRIX_FILE\". Contents:"
    # Extract the matrix file
    tar xzfv "$MATRIX_FILE" -C "$TEMP_DIRECTORY"
else
    echo -e "${SCRIPT_TAG} Matrix file already downloaded in $TEMP_DIRECTORY"
fi

export LD_LIBRARY_PATH=${MKLROOT}/lib

# Run the spmv
echo -e "${SCRIPT_TAG} Building executable..."

if [ "$FLOAT_DATA_TYPE" == "float" ]; then
    CMD="spmv_mkl"
else
    CMD="spmvd_mkl"
fi

if [ "${CACHE_MODE}" == "hot" ]; then
    KERNEL_REPS=100;
else
    KERNEL_REPS=1;
fi

if [ "${EXECUTION_THREADS}" == "1th" ];
then
    OMP_NUM_THREADS="1"
    PIN_CORES="10"
elif [ "${EXECUTION_THREADS}" == "2th" ];
then
    CMD="${CMD}.omp"
    OMP_NUM_THREADS="2"
    PIN_CORES="10,12"
else
    CMD="${CMD}.omp"
    OMP_NUM_THREADS="8"
    PIN_CORES="10,12,14,1,2,4,6,8"
fi

# Recompile and launch
make -C spmv-mkl -B ${CMD} KERNEL_REPS=${KERNEL_REPS} MKLROOT=${MKLROOT} ONEAPIROOT=${ONEAPIROOT}

RUNCOUNT=10
i=0; while [[ $i -lt RUNCOUNT ]]; do
    OMP_NUM_THREADS=${OMP_NUM_THREADS} numactl -C${PIN_CORES} hugectl --heap ./spmv-mkl/${CMD} "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}.rb"
    i=$((i+1));
done

cd - > /dev/null
