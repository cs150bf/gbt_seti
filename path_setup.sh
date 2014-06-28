export SETI_GBT=${PWD%'/src'}
export SETI_GBT=${SETI_GBT%'/bin'}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH//":${SETI_GBT}/lib"/}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SETI_GBT}/lib
export CUDA_PATH=/usr/local/cuda
export PATH=${PATH//":${CUDA_PATH}/bin"/}
export PATH=${PATH}:${CUDA_PATH}/bin
export PRESTO=/home/siemion/sw/presto
echo "SETI_GBT="$SETI_GBT
echo "LD_LIBRARY_PATH="$LD_LIBRARY_PATH
echo "CUDA_PATH="$CUDA_PATH
echo "PRESTO="$PRESTO
echo "PATH="$PATH
