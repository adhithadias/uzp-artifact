#!/bin/bash

if [ "$#" -lt 5 ]; then
    echo "Usage: $0 GROUP MATRIX_NAME [float] [hot|cold] [1th]"
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
if [[ "$FLOAT_DATA_TYPE" != "float" ]]; then
    echo -e "${ERROR_TAG} Invalid FLOAT_DATA_TYPE: $FLOAT_DATA_TYPE. Must be 'float' ('double' not supported by MACVETH artifact)."
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
if [[ "${EXECUTION_THREADS}" != "1th" ]];
then
    echo -e "${ERROR_TAG} Invalid EXECUTION_THREADS: ${EXECUTION_THREADS}. Must be '1th' (multithreaded execution not supported by MACVETH artifact)."
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

# Get the specific data-specific codes from the DS-SpMV repository
MACVETH_REPO_URL="https://github.com/gabriel-rodriguez/DS-SpMV/raw/refs/heads/main/pact22/"
if [ ! -f "$TEMP_DIRECTORY/${MATRIX_NAME}/${MATRIX_NAME}-d1g4.ast.c" ]; then
	# Download the skeleton file
	wget -P ${TEMP_DIRECTORY}/${MATRIX_NAME}/ ${MACVETH_REPO_URL}/${MATRIX_NAME}-d1g4.ast.c.bz2

	# Download the fragment files
	fragment_count=0
	fragment_count_str="00"
	while wget -P ${TEMP_DIRECTORY}/${MATRIX_NAME}/ ${MACVETH_REPO_URL}/${MATRIX_NAME}-d1g4.fragments-${fragment_count_str}.mv.c.bz2;
	do
		fragment_count=$((fragment_count+1))
		printf -v fragment_count_str "%02d" ${fragment_count}
	done

	# Bunzip everything in parallel
	bunzip2 ${TEMP_DIRECTORY}/${MATRIX_NAME}/*.c.bz2;
fi

echo -e "${SCRIPT_TAG} Building Macveth code (this may take a long time, consider modifying this macveth_spmv.sh script to launch in parallel by adding the appropriate -j knob)..."

if [ "${CACHE_MODE}" == "hot" ]; then
    KERNEL_REPS=100;
else
    KERNEL_REPS=1;
fi

# Copy build files to temp directory and rebuild as needed
cp spmv-macveth/* ${TEMP_DIRECTORY}/${MATRIX_NAME}/
make -C ${TEMP_DIRECTORY}/${MATRIX_NAME} -j1 ${MATRIX_NAME}-d1g4.mv KERNEL_REPS=${KERNEL_REPS} # Change -j1 as needed to parallelize build

RUNCOUNT=10
i=0; while [[ $i -lt RUNCOUNT ]]; do
    numactl -C10 hugectl --heap ${TEMP_DIRECTORY}/${MATRIX_NAME}/${MATRIX_NAME}-d1g4.mv ${TEMP_DIRECTORY}/${MATRIX_NAME}/${MATRIX_NAME}.rb
    i=$((i+1));
done
