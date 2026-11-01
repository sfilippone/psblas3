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

extern CFIBuffer dotBuffer;

#undef GEN_CUBLAS_FUNC_NAME
#define GEN_CUBLAS_FUNC_NAME(x, t) CONCAT(CONCAT(cublas, x), t)

#undef GEN_SPGPU_FUNC_NAME
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu, x), dot)

#undef GEN_SPGPU_FUNC_NAMEM
#define GEN_SPGPU_FUNC_NAMEM(x) CONCAT(CONCAT(spgpu, x), mdot)

#undef GEN_SPGPU_FUNC_NAMEMS
#define GEN_SPGPU_FUNC_NAMEMS(x) CONCAT(CONCAT(spgpu, x), msdot)

#undef GEN_SPGPU_FUNC_NAMEMV
#define GEN_SPGPU_FUNC_NAMEMV(x) CONCAT(CONCAT(spgpu, x), mvdot)

#undef GEN_SPGPU_FUNC_NAMEMM
#define GEN_SPGPU_FUNC_NAMEMM(x) CONCAT(CONCAT(spgpu, x), mmdot)

#undef GEN_SPGPU_KERN_NAME
#define GEN_SPGPU_KERN_NAME(x) CONCAT(CONCAT(CONCAT(spgpu, x), dot), _krn)

#define BLOCK_SIZE 512
// #define ASSUME_LOCK_SYNC_PARALLELISM

static const VALUE_TYPE CONCAT(CONCAT(cu, TYPE_SYMBOL), one) = CONCAT(one_, VALUE_TYPE)();
static const VALUE_TYPE CONCAT(CONCAT(cu, TYPE_SYMBOL), zero) = CONCAT(zero_, VALUE_TYPE)();
static __device__ VALUE_TYPE ddotReductionResult[128];

__global__ void GEN_SPGPU_KERN_NAME(TYPE_SYMBOL)(int n, VALUE_TYPE* x, VALUE_TYPE* y)
{
	__shared__ VALUE_TYPE sSum[BLOCK_SIZE];

	VALUE_TYPE res = CONCAT(zero_, VALUE_TYPE)();
	VALUE_TYPE* lastX = x + n;

	x += threadIdx.x + blockIdx.x * BLOCK_SIZE;
	y += threadIdx.x + blockIdx.x * BLOCK_SIZE;

	int blockOffset = gridDim.x * BLOCK_SIZE;

	while (x < lastX)
    {
		res = CONCAT(VALUE_TYPE, _add)(res, CONCAT(VALUE_TYPE, _mul)(x[0], y[0]));
		
		x += blockOffset;
		y += blockOffset;
	}

	if(threadIdx.x >= 32)
		sSum[threadIdx.x] = res;

	__syncthreads();

	// Start reduction!
	if(threadIdx.x < 32) 
	{
		for(int i = 1; i < BLOCK_SIZE/32; ++i)
			res += sSum[i*32 + threadIdx.x];

	//useless (because inter-warp)
#ifndef	ASSUME_LOCK_SYNC_PARALLELISM
	}
	__syncthreads(); 

	if(threadIdx.x < 32) 
	{
#endif	

#ifdef ASSUME_LOCK_SYNC_PARALLELISM
		volatile VALUE_TYPE* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if(threadIdx.x < 16) vsSum[threadIdx.x] += vsSum[threadIdx.x + 16];
		if(threadIdx.x < 8) vsSum[threadIdx.x] += vsSum[threadIdx.x + 8];
		if(threadIdx.x < 4) vsSum[threadIdx.x] += vsSum[threadIdx.x + 4];
		if(threadIdx.x < 2) vsSum[threadIdx.x] += vsSum[threadIdx.x + 2];
		if(threadIdx.x == 0)
			ddotReductionResult[blockIdx.x] = vsSum[0] + vsSum[1];
#else
		VALUE_TYPE* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if(threadIdx.x < 16) vsSum[threadIdx.x] += vsSum[threadIdx.x + 16];
		__syncthreads();
		if(threadIdx.x < 8) vsSum[threadIdx.x] += vsSum[threadIdx.x + 8];
		__syncthreads();
		if(threadIdx.x < 4) vsSum[threadIdx.x] += vsSum[threadIdx.x + 4];
		__syncthreads();
		if(threadIdx.x < 2) vsSum[threadIdx.x] += vsSum[threadIdx.x + 2];
		__syncthreads();
		if(threadIdx.x == 0)
		ddotReductionResult[blockIdx.x] = vsSum[0] + vsSum[1];
#endif
	}
}

VALUE_TYPE GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(spgpuHandle_t handle, int n, __device VALUE_TYPE* a, __device VALUE_TYPE* b)
{
	VALUE_TYPE res = CONCAT(zero_, VALUE_TYPE)();
	int device;
	cudaGetDevice(&device);
#if 0
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, device);	

	int blocks = min(128, min(prop.multiProcessorCount, (n + BLOCK_SIZE - 1) / BLOCK_SIZE));
#else
	int blocks = min(128, min(handle->multiProcessorCount, (n + BLOCK_SIZE - 1) / BLOCK_SIZE));
#endif
	
	VALUE_TYPE tRes[128];

	GEN_SPGPU_KERN_NAME(TYPE_SYMBOL)<<<blocks, BLOCK_SIZE, 0, handle->currentStream>>>(n, a, b);
	cudaMemcpyFromSymbol(tRes, ddotReductionResult, blocks * sizeof(VALUE_TYPE));

	for(int i = 0; i < blocks; ++i)
		res += tRes[i];

	cudaCheckError("CUDA error on ddot");
	return res;
}

void GEN_SPGPU_FUNC_NAMEM(TYPE_SYMBOL)(spgpuHandle_t handle, VALUE_TYPE* y, int n, __device VALUE_TYPE* a, __device VALUE_TYPE* b, int count, int pitch)
{
	for(int i = 0; i < count; ++i)
	{
		y[i] = GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, n, a, b);
		a += pitch;
		b += pitch;
	}
}

void GEN_SPGPU_FUNC_NAMEMS(TYPE_SYMBOL)(spgpuHandle_t handle, VALUE_TYPE* y, int n, __device VALUE_TYPE* a, __device VALUE_TYPE* b)
{
	#ifdef USE_CUBLAS
		cublasHandle_t cublasHandle = psb_cudaGetCublasHandle();
		GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL, dot)(cublasHandle, n, a, 1, b, 1, y);
	#else
		*y = GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, n, a, b);
	#endif
}

void GEN_SPGPU_FUNC_NAMEMV(TYPE_SYMBOL)(spgpuHandle_t handle, VALUE_TYPE* y, int n, __device VALUE_TYPE* a, __device VALUE_TYPE* b, int countA, int pitchA)
{
	#if 0
		VALUE_TYPE res[countA];
		GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL, gemv)(handle->cublasHandle, CUBLAS_OP_T, n, countA, &c_done, a, pitchA, b, &c_done, &c_dzero, res, 1);
		cudaDeviceSynchronize();
		for(int i = 0; i < countA; ++i)
			y[i] = res[i];
	#else
		//TO DO: optimize with custom gemv kernel?
		for(int i = 0; i < countA; ++i)
			y[i] = GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, n, a + (i * pitchA), b);
	#endif
}

void CONCAT(GEN_SPGPU_FUNC_NAMEMV(TYPE_SYMBOL), _CFI)(spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device VALUE_TYPE* a, __device VALUE_TYPE* b, int countA, int pitchA)
{
	#ifdef USE_CUBLAS
		cublasHandle_t cublasHandle = psb_cudaGetCublasHandle();
		CFIdevBufferAlloc(&dotBuffer, countA * (y->elem_len));
		GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL, gemv)(cublasHandle, CUBLAS_OP_T, n, countA, &CONCAT(CONCAT(cu, TYPE_SYMBOL), one), a, pitchA, b, 1, &CONCAT(CONCAT(cu, TYPE_SYMBOL), zero), (VALUE_TYPE*) dotBuffer.devBuffer, 1);
		cudaDeviceSynchronize();
		CFIBufferStore(handle, y, &dotBuffer, countA * (y->elem_len));
		cudaDeviceSynchronize();
	#else
		//TO DO: optimize with custom gemv kernel?
		for(int i = 0; i < countA; ++i)
			CFI_AT1(VALUE_TYPE, y, i) = GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, n, a + (i * pitchA), b);
	#endif
}

void GEN_SPGPU_FUNC_NAMEMM(TYPE_SYMBOL)(spgpuHandle_t handle, VALUE_TYPE* y, int n, __device VALUE_TYPE* a, __device VALUE_TYPE* b, int countA, int pitchA, int countB, int pitchB)
{
	#if 0
		VALUE_TYPE res[countA*countB];
		GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL, gemm)(handle->cublasHandle, CUBLAS_OP_T, CUBLAS_OP_N, countA, countB, n, &c_done, a, pitchA, b, &c_done, &c_zero, res, countA);
		cudaDeviceSynchronize();
		for(int i = 0; i < countB; ++i)
			for(int j = 0; j < countA; ++j)
				y[i*countA + j] = res[i*countA + j];
	#else
		//TO DO: optimize with custom gemm kernel?
		for(int i = 0; i < countB; ++i)
			for(int j = 0; j < countA; ++j)
				y[i*countA + j] = GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, n, a + (i * pitchA), b + (j * pitchB));
	#endif
}

void CONCAT(GEN_SPGPU_FUNC_NAMEMM(TYPE_SYMBOL), _CFI)(spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device VALUE_TYPE* a, __device VALUE_TYPE* b, int countA, int pitchA, int countB, int pitchB)
{
	#ifdef USE_CUBLAS
		cublasHandle_t cublasHandle = psb_cudaGetCublasHandle();
		CFIdevBufferAlloc(&dotBuffer, countA * countB * (y->elem_len));
		GEN_CUBLAS_FUNC_NAME(TYPE_SYMBOL, gemm)(cublasHandle, CUBLAS_OP_T, CUBLAS_OP_N, countA, countB, n, &CONCAT(CONCAT(cu, TYPE_SYMBOL), one), a, pitchA, b, pitchB, &CONCAT(CONCAT(cu, TYPE_SYMBOL), zero), (VALUE_TYPE*) dotBuffer.devBuffer, countA);
		cudaDeviceSynchronize();
		CFIBufferStore(handle, y, &dotBuffer, countA * countB * (y->elem_len));
		cudaDeviceSynchronize();
	#else
		//TO DO: optimize with custom gemm kernel?
		for(int i = 0; i < countA; ++i)
			for(int j = 0; j < countB; ++j)
				CFI_AT2(VALUE_TYPE, y, i, j) = GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, n, a + (i * pitchA), b + (j * pitchB));
	#endif
}