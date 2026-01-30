ifndef POLYTEST
POLYTEST=-DPOLYBENCH_GFLOPS
endif
POLYFLAGS=${POLYTEST} -DPOLYBENCH_CACHE_SIZE="(65*1024)"

CC=gcc
CFLAGS=-Ofast -march=native -mtune=native -ftree-vectorize -ftree-loop-vectorize -fvect-cost-model=unlimited -fsimd-cost-model=unlimited -fprefetch-loop-arrays -D_POSIX_C_SOURCE=200112L -Wno-incompatible-pointer-types -DREPS=$(KERNEL_REPS)

MACVETH=/home/gabriel/workspace/MACVETH-test4/build/macveth
MACVETH_OPTS=--mvmarch=alderlake --misa=avx2 --simd-cost-model=unlimited --min-redux-size=4 --redux-win-size=-1

.PRECIOUS: %.S %.o %.mv.c
.SUFFIXES:

%.S: %.c
	${CC} ${CFLAGS} ${POLYFLAGS} -I . -S $< -o $@

%.o: %.S
	${CC} ${CFLAGS} ${POLYFLAGS} -c -o $@ $<
	objdump -d $@ > $@.dump
	PYTHONOPTIMIZE=2 python smart_fix_asm.py $< $@.dump $<.tmpfix.S
	${CC} ${CFLAGS} ${POLYFLAGS} -c -o $@ $<.tmpfix.S
	rm $<.tmpfix.S $@.dump

.SECONDEXPANSION:
%.mv: %.ast.c spmv_poly.c polybench.c $$(subst .mv.c,.mv.o,$$(wildcard %.fragments-*.mv.c))
	${CC} ${CFLAGS} ${POLYFLAGS} -I . -B/usr/share/libhugetlbfs -Wl,--hugetlbfs-align -o $@ $^ -lrbio -lm -lpapi
	$(eval $@_TMP := $(shell mktemp -d))
	$$(nm $@ |grep kernel_spmv|cut -d " " -f 1|awk '{val="0x"$$1;print strtonum(val);}'|sort -n > $@_TMP)
	TEXT_START=$$(head -n 1 < $@_TMP)&& \
	TEXT_END=$$(tail -n 1 < $@_TMP)&& \
	TEXT_SEGMENT_SIZE=$$(bc <<< $${TEXT_END}-$${TEXT_START})&& \
	echo "Recompiling with -DTEXT_SEGMENT_SIZE=$${TEXT_SEGMENT_SIZE}"&& \
	${CC} ${CFLAGS} ${POLYFLAGS} -I . -B/usr/share/libhugetlbfs -Wl,--hugetlbfs-align -DTEXT_SEGMENT_SIZE=$${TEXT_SEGMENT_SIZE} -o $@ $^ -lrbio -lm -lpapi
	$$(rm $@_TMP)

clean:
	rm -f ${BINARIES}
