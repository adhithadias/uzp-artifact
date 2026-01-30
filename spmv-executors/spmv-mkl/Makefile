ifndef POLYTEST
POLYTEST=-DPOLYBENCH_GFLOPS
endif
POLYFLAGS=${POLYTEST} -DPOLYBENCH_CACHE_SIZE="(65*1024)"

CC=gcc-11
CFLAGS=-Ofast -march=native -mtune=native -ftree-vectorize -ftree-loop-vectorize -fvect-cost-model=unlimited -fsimd-cost-model=unlimited -fprefetch-loop-arrays -D_POSIX_C_SOURCE=200112L -Wno-incompatible-pointer-types -I $(MKLROOT)/include -DREPS=$(KERNEL_REPS)

BINARIES=spmv_mkl spmvd_mkl spmv_mkl.omp spmvd_mkl.omp

all: ${BINARIES}

%: %.c polybench.c
	${CC} ${CFLAGS} ${POLYFLAGS} -DMKL_ILP64 -I . -o $@ $^ -lm -lrbio -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread # Single threaded

%.omp: %.c polybench.c
#	${CC} ${CFLAGS} ${POLYFLAGS} -DMKL_ILP64 -I . -o $@ $^ -lrbio -lm -L${MKLROOT}/lib/intel64 -L/home/gabriel/workspace/intel/oneapi/compiler/2022.0.2/linux/compiler/lib/intel64_lin/ -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl # Multithreaded MKL
	${CC} ${CFLAGS} ${POLYFLAGS} -DMKL_ILP64 -I . -o $@ $^ -lrbio -lm -L$(MKLROOT)/lib/intel64 -L$(ONEAPIROOT)/compiler/2024.1/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread # Multithreaded MKL

clean:
	rm -f ${BINARIES}
