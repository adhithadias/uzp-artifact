ifndef POLYTEST
POLYTEST=-DPOLYBENCH_GFLOPS
endif
POLYFLAGS=${POLYTEST} -DPOLYBENCH_CACHE_SIZE="(65*1024)"

CC=gcc
CFLAGS=-Ofast -march=native -mtune=native -ftree-vectorize -ftree-loop-vectorize -fvect-cost-model=unlimited -fsimd-cost-model=unlimited -fprefetch-loop-arrays -D_POSIX_C_SOURCE=200112L -Wno-incompatible-pointer-types -DREPS=$(KERNEL_REPS)
CFLAGS_OMP=-fopenmp ${CFLAGS}

BINARIES=spmv_csr spmvd_csr spmv_csr.omp spmvd_csr.omp

all: ${BINARIES}

%: %.c polybench.c
	${CC} ${CFLAGS} ${POLYFLAGS} -I . -o $@ $^ -lm -lrbio

%.omp: %.c polybench.c
	${CC} ${CFLAGS_OMP} ${POLYFLAGS} -I . -o $@ $^ -lrbio -lm

clean:
	rm -f ${BINARIES}
