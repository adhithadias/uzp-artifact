#!/bin/bash

# if PROJECT_ROOT is not set, set it to /uzp-artifact
if [ -z "$PROJECT_ROOT" ]; then
    PROJECT_ROOT="/uzp-artifact"
fi

pushd z_polyhedrator && RUSTUP_HOME="/opt/rustup" CARGO_HOME="/tmp/.cargo_pldi25_artifact" cargo clean && popd

rm -rf /opt/rustup
rm -rf /tmp/.cargo_pldi25_artifact

rm -rf $PROJECT_ROOT/lib/scripts/__pycache__/
rm -rf $PROJECT_ROOT/spmv-executors/__pycache__/

find /tmp -type f -name "*.uzp" -delete

make -C $PROJECT_ROOT/spmv-executors/uzp-genex distclean clean
make -C $PROJECT_ROOT/spmv-executors/spmv-mkl clean
make -C $PROJECT_ROOT/spmv-executors/spmv-csr clean
make -C $PROJECT_ROOT/spmv-executors/spmv-csr5/CSR5_avx2 clean