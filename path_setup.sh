export OPT64=/opt/vegas
export CUDA=/usr/local/cuda
export CUDAC=${CUDA}/bin/nvcc
export SETI_GBT=${PWD%'/src'}
export SETI_GBT=${SETI_GBT%'/bin'}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH//":${SETI_GBT}/lib"/}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SETI_GBT}/lib
export PATH=${PATH//":${CUDA}/bin"/}
export PATH=${PATH}:${CUDA}/bin
export PRESTO=/home/siemion/sw/presto
echo "OPT64="$OPT64
echo "CUDA="$CUDA
echo "CUDAC="$CUDAC
echo "SETI_GBT="$SETI_GBT
echo "LD_LIBRARY_PATH="$LD_LIBRARY_PATH
echo "PRESTO="$PRESTO
echo "PATH="$PATH
