
"""uzp.py

Minimal entry-point guard added so this module can be executed as a script
or imported without running top-level code.
"""

import sys
import os

UZP_SPMV_EXECUTORS_DIR = "spmv-executors"
MATRIX_NAME = "sparse_matrix"

# get PROJECT_ROOT from environment variable if set, else error message
UZP_PROJECT_ROOT = os.getenv("UZP_PROJECT_ROOT", None)
if UZP_PROJECT_ROOT is None:
    print("ERROR: UZP_PROJECT_ROOT environment variable is not set.")
    sys.exit(1)

def main(argv=None):
	# create 3x4 sparse matrix
    coordinates = [(0, 0), (1, 2), (2, 3)]
    values = [1.0, 2.0, 3.0]
    shape = (3, 4)
    from scipy.sparse import coo_matrix
    sparse_matrix = coo_matrix((values, zip(*coordinates)), shape=shape)
    
    # save sparse matrix to .mtx file in tmp directory
    # directory here is /tmp/sparse_matrix/sparse_matrix.mtx
    import os
    import tempfile
    tmp_dir = tempfile.gettempdir()
    tmp_dir = os.path.join(tmp_dir, MATRIX_NAME)
    os.makedirs(tmp_dir, exist_ok=True)
    mtx_file_path = os.path.join(tmp_dir, f"{MATRIX_NAME}.mtx")
    from scipy.io import mmwrite
    mmwrite(mtx_file_path, sparse_matrix)
    print(f"Sparse matrix saved to {mtx_file_path}")
    
    # Run the uzp_spmv.sh script to generate tuned .uzp file
    # before running it cd UZP_PROJECT_ROOT/spmv-executors
    # this script uses z_polyhedrator to generate the uzp file
    uzp_spmv_script = os.path.join(UZP_PROJECT_ROOT, UZP_SPMV_EXECUTORS_DIR, "uzp_spmv_modified.sh")
    import subprocess
    try:
        subprocess.run(
            [uzp_spmv_script, MATRIX_NAME, "float", "hot", "1th"],
            check=True,
            cwd=os.path.join(UZP_PROJECT_ROOT, UZP_SPMV_EXECUTORS_DIR)
        )
        print(f"UZP generation completed successfully for {MATRIX_NAME}")
    except subprocess.CalledProcessError as e:
        print(f"UZP generation failed: {e}")
    
    # generate a file similar to spf_matvect.c here that reads the generated uzp file and performs SpMV
    ######## ALL YOUR COMPILER LOGIC TO CREATE spf_matvect.c GOES HERE ########

    # compile and run the spf_matvect executor on the generated uzp file
    spf_matvect_c = os.path.join(UZP_PROJECT_ROOT, "spf_matvect.c")
    spf_matvect_o = os.path.join(tempfile.gettempdir(), "spf_matvect.o")
    spf_matvec_executable = os.path.join(tempfile.gettempdir(), "spf_matvec")
    uzp_file_path = os.path.join(tmp_dir, f"{MATRIX_NAME}.tuned.uzp")
    
    try:
        subprocess.run(
            [
                "gcc-11", "-O2", "-march=native", "-ffast-math", "-ftree-vectorize",
                "-floop-unroll-and-jam", "-DREPS=1", "-I", os.path.join(UZP_PROJECT_ROOT, UZP_SPMV_EXECUTORS_DIR, "uzp-genex"),
                "-I", ".", "-DGEN_EXECUTOR_SPMV_V2", "-DPOLYBENCH_GFLOPS",
                "-c", spf_matvect_c, "-o", spf_matvect_o
            ],
            check=True
        )
        subprocess.run(
            [
                "gcc-11", "-O2", "-march=native", "-ffast-math", "-ftree-vectorize",
                "-floop-unroll-and-jam", "-DREPS=1", "-I", os.path.join(UZP_PROJECT_ROOT, UZP_SPMV_EXECUTORS_DIR, "uzp-genex"),
                "-I", ".", "-DGEN_EXECUTOR_SPMV_V2", "-DPOLYBENCH_GFLOPS",
                os.path.join(UZP_PROJECT_ROOT, UZP_SPMV_EXECUTORS_DIR, "uzp-genex", "polybench.o"),
                os.path.join(UZP_PROJECT_ROOT, UZP_SPMV_EXECUTORS_DIR, "uzp-genex", "spf_structure.o"),
                os.path.join(UZP_PROJECT_ROOT, UZP_SPMV_EXECUTORS_DIR, "uzp-genex", "spf_executors.o"),
                os.path.join(UZP_PROJECT_ROOT, UZP_SPMV_EXECUTORS_DIR, "uzp-genex", "spf_executors_uninc.o"),
                spf_matvect_o,
                "-o", spf_matvec_executable,
                "-lm"
            ],
            check=True
        )
        subprocess.run(
            [
                "numactl", "-C10", "hugectl", "--heap", spf_matvec_executable,
                "-float", uzp_file_path
            ],
            check=True,
            env={**os.environ, "OMP_NUM_THREADS": "1"}
        )
        print(f"SPF SpMV execution completed successfully for {MATRIX_NAME}")
    except subprocess.CalledProcessError as e:
        print(f"SPF SpMV execution failed: {e}")    

if __name__ == "__main__":
	main()
    
    