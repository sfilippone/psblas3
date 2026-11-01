#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)

#include "stdio.h"
#include "cudalang.h"
#include "cudadebug.h"
#include "core.h"
#include "cublas_v2.h"
#include "mathbase.cuh"
extern "C"
{
	#include "vector.h"
	cublasHandle_t psb_cudaGetCublasHandle();
}

#undef GEN_CUBLAS_FUNC_NAME
#define GEN_CUBLAS_FUNC_NAME(x) CONCAT(CONCAT(cublas, x), nrm2)

#undef GEN_SPGPU_FUNC_NAME
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu, x), nrm2)

#undef GEN_SPGPU_FUNC_NAMES
#define GEN_SPGPU_FUNC_NAMES(x) CONCAT(CONCAT(spgpu, x), msnrm2)

#undef GEN_SPGPU_FUNC_NAMEV
#define GEN_SPGPU_FUNC_NAMEV(x) CONCAT(CONCAT(spgpu, x), mvnrm2)

#undef GEN_SPGPU_KERN_NAME
#define GEN_SPGPU_KERN_NAME(x) CONCAT(CONCAT(CONCAT(spgpu, x), nrm2), _krn)

//#define ASSUME_LOCK_SYNC_PARALLELISM
#define BLOCK_SIZE 512

static __device__ RVALUE_TYPE CONCAT(TYPE_SYMBOL, nrm2ReductionResult)[128];

__global__ void GEN_SPGPU_KERN_NAME(TYPE_SYMBOL)(int n, VALUE_TYPE* x)
{
	__shared__ RVALUE_TYPE sSum[BLOCK_SIZE];

	RVALUE_TYPE res = 0;
	VALUE_TYPE* lastX = x + n;

	x += threadIdx.x + blockIdx.x*BLOCK_SIZE;
	int blockOffset = gridDim.x*BLOCK_SIZE;

	while (x < lastX)
    {
		VALUE_TYPE x1 = x[0];
		res = CONCAT(RVALUE_TYPE, _add)(res, REALOP(CONCAT(VALUE_TYPE, _mul)(x1, CONJOP(x1))));
		
		x += blockOffset;
	}

	if (threadIdx.x >= 32)
		sSum[threadIdx.x] = res;

	__syncthreads();

	// Start reduction!
	if (threadIdx.x < 32) 
	{
		for(int i = 1; i < BLOCK_SIZE/32; ++i)
			res += sSum[i*32 + threadIdx.x];

	//useless (because inter-warp)
#ifndef	ASSUME_LOCK_SYNC_PARALLELISM
	}
	__syncthreads(); 

	if (threadIdx.x < 32) 
	{
#endif	

#ifdef ASSUME_LOCK_SYNC_PARALLELISM
		volatile VALUE_TYPE* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) vsSum[threadIdx.x] += vsSum[threadIdx.x + 16];
		if (threadIdx.x < 8) vsSum[threadIdx.x] += vsSum[threadIdx.x + 8];
		if (threadIdx.x < 4) vsSum[threadIdx.x] += vsSum[threadIdx.x + 4];
		if (threadIdx.x < 2) vsSum[threadIdx.x] += vsSum[threadIdx.x + 2];
		if (threadIdx.x == 0)
			CONCAT(TYPE_SYMBOL, nrm2ReductionResult)[blockIdx.x] = vsSum[0] + vsSum[1];

#else
		RVALUE_TYPE* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) vsSum[threadIdx.x] += vsSum[threadIdx.x + 16];
		__syncthreads();
		if (threadIdx.x < 8) vsSum[threadIdx.x] += vsSum[threadIdx.x + 8];
		__syncthreads();
		if (threadIdx.x < 4) vsSum[threadIdx.x] += vsSum[threadIdx.x + 4];
		__syncthreads();
		if (threadIdx.x < 2) vsSum[threadIdx.x] += vsSum[threadIdx.x + 2];
		__syncthreads();
		if (threadIdx.x == 0)
			CONCAT(TYPE_SYMBOL, nrm2ReductionResult)[blockIdx.x] = vsSum[0] + vsSum[1];
#endif	
	}
}

RVALUE_TYPE GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(spgpuHandle_t handle, int n, VALUE_TYPE* x)
{
	RVALUE_TYPE res = CONCAT(zero_, RVALUE_TYPE)();
	#if 0 	
		int device;
		cudaGetDevice(&device);
		struct cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop,device);	

		int blocks = min(128, min(prop.multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
	#else
		int blocks = min(128, min(handle->multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
	#endif
	RVALUE_TYPE tRes[128];

	GEN_SPGPU_KERN_NAME(TYPE_SYMBOL)<<<blocks, BLOCK_SIZE, 0, handle->currentStream>>>(n, x);
	cudaMemcpyFromSymbol(tRes, CONCAT(TYPE_SYMBOL, nrm2ReductionResult), blocks * sizeof(RVALUE_TYPE));

	for(int i = 0; i < blocks; ++i)
		res += tRes[i];

	cudaCheckError("CUDA error on dnrm2");
	
	return sqrt(res);
}

void GEN_SPGPU_FUNC_NAMES(TYPE_SYMBOL)(spgpuHandle_t handle, RVALUE_TYPE *y, int n, __device VALUE_TYPE *x)
{
	#ifdef USE_CUBLAS
		cublasHandle_t cublasHandle = psb_cudaGetCublasHandle();
		GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL)(cublasHandle, n, x, 1, y);
		cudaDeviceSynchronize();
	#else
		*y = GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, n, x);
	#endif
}

void GEN_SPGPU_FUNC_NAMEV(TYPE_SYMBOL)(spgpuHandle_t handle, CFI_cdesc_t *y, int n, __device VALUE_TYPE *x, int count, int pitch)
{
	for(int i = 0; i < count; ++i)
	{
		#ifdef USE_CUBLAS
			cublasHandle_t cublasHandle = psb_cudaGetCublasHandle();
			RVALUE_TYPE res;
			GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL)(cublasHandle, n, x, 1, &res);
			cudaDeviceSynchronize();
			CFI_AT1(RVALUE_TYPE, y, i) = res;
		#else
			CFI_AT1(RVALUE_TYPE, y, i) = GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, n, x);
		#endif
		x += pitch;
	}
}