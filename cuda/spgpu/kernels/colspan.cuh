#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)

#include "cudadebug.h"
#include "cudalang.h"
#include <cuda_runtime.h>
#include "core.h"
#include "cublas_v2.h"
#include "mathbase.cuh"
extern "C"
{
	#include "vector.h"
	//#include "cuda_util.h"
	int getGPUMultiProcessors();
	int getGPUMaxThreadsPerMP();
	cublasHandle_t psb_cudaGetCublasHandle();
}

extern CFIBuffer coeffBuffer;

#undef GEN_CUBLAS_FUNC_NAME
#define GEN_CUBLAS_FUNC_NAME(x, t) CONCAT(CONCAT(cublas, x), t)

#undef GEN_SPGPU_FUNC_NAME
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu, x), colspan)

#undef GEN_SPGPU_KERN_NAME
#define GEN_SPGPU_KERN_NAME(x, d) CONCAT(CONCAT(CONCAT(CONCAT(spgpu, x), colspan), d), _krn)

static const VALUE_TYPE CONCAT(CONCAT(cu, TYPE_SYMBOL), one) = CONCAT(one_, VALUE_TYPE)();
static const VALUE_TYPE CONCAT(CONCAT(cu, TYPE_SYMBOL), zero) = CONCAT(zero_, VALUE_TYPE)();
#define BLOCK_SIZE 512
#define USE_CUBLAS 1

__global__ void GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, 1D)(int n, VALUE_TYPE* x, int pitchX, int countX, const VALUE_TYPE* coeffPtr, VALUE_TYPE* y, bool updFlag)
{
	size_t id = threadIdx.x + BLOCK_SIZE * blockIdx.x;
	const size_t gridSize = blockDim.x * gridDim.x;

	for(; id < n; id += gridSize)
	{
		VALUE_TYPE acc = updFlag ? y[id] : CONCAT(zero_, VALUE_TYPE)();
		for(size_t j = 0; j < countX; ++j)
			acc = CONCAT(VALUE_TYPE, _add)(CONCAT(VALUE_TYPE, _mul)(coeffPtr[j], x[id + (j*pitchX)]), acc);
		y[id] = acc;
	}
}

__global__ void GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, 2D)(int n, VALUE_TYPE* x, int countX, int pitchX, const VALUE_TYPE* coeffPtr, VALUE_TYPE* y, int countY, int pitchY, bool updFlag)
{
	size_t id = threadIdx.x + BLOCK_SIZE * blockIdx.x;
	const size_t gridSize = blockDim.x * gridDim.x;

	for(; id < n; id += gridSize)
		for(size_t k = 0; k < countY; ++k)
		{
			VALUE_TYPE acc = updFlag ? y[id + (k * pitchY)] : CONCAT(zero_, VALUE_TYPE)();
			for(size_t j = 0; j < countX; ++j)
				acc = CONCAT(VALUE_TYPE, _add)(CONCAT(VALUE_TYPE, _mul)(coeffPtr[j + k * countX], x[id + (j * pitchX)]), acc);
			y[id + (k * pitchY)] = acc;
		}
}

void GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(spgpuHandle_t handle, int n, __device VALUE_TYPE* x, int countX, int pitchX, CFI_cdesc_t* coeff, __device VALUE_TYPE* y, int countY, int pitchY, bool updFlag)
{
    CFIBufferLoad(handle, coeff, &coeffBuffer);

	#if USE_CUBLAS
		cublasHandle_t cublasHandle = psb_cudaGetCublasHandle();
		const VALUE_TYPE* beta = updFlag ? &(CONCAT(CONCAT(cu, TYPE_SYMBOL), one)) : &(CONCAT(CONCAT(cu, TYPE_SYMBOL), zero));
		if(coeff->rank == 1) GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL, gemv)(cublasHandle, CUBLAS_OP_N, n, countX, &(CONCAT(CONCAT(cu, TYPE_SYMBOL), one)), x, pitchX, (VALUE_TYPE*) coeffBuffer.devBuffer, 1, beta, y, 1);
		if(coeff->rank == 2) GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL, gemm)(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, n, countY, countX, &(CONCAT(CONCAT(cu, TYPE_SYMBOL), one)), x, pitchX, (VALUE_TYPE*) coeffBuffer.devBuffer, countX, beta, y, pitchY);
	#else
		int num_mp         = getGPUMultiProcessors();
		int max_threads_mp = getGPUMaxThreadsPerMP();
		int num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
		int num_blocks     = num_blocks_mp * num_mp;
		dim3 grid(num_blocks);
		dim3 block(BLOCK_SIZE);
		if(coeff->rank == 1) GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, 1D)<<<grid, block, 0, handle->currentStream>>>(n, x, countX, pitchX, (VALUE_TYPE*) coeffBuffer.devBuffer, y, updFlag);
		if(coeff->rank == 2) GEN_SPGPU_KERN_NAME(TYPE_SYMBOL, 2D)<<<grid, block, 0, handle->currentStream>>>(n, x, countX, pitchX, (VALUE_TYPE*) coeffBuffer.devBuffer, y, countY, pitchY, updFlag);
	#endif
}