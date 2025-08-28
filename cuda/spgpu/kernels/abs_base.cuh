/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2015
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */


#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)

#undef GEN_SPGPU_FUNC_NAME
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu,x),abs)

#define BLOCK_SIZE 256

// Define:
//#define RES_VALUE_TYPE
//#define VALUE_TYPE
//#define TYPE_SYMBOL

#include "mathbase.cuh"

__device__ __host__ static inline bool is_one_float(float x) { return (x==1.0f); }
__device__ __host__ static inline bool is_one_cuFloatComplex(cuFloatComplex x) { return ((x.x==1.0f)&&(x.y==0.0f));}


#if (__CUDA_ARCH__ >= 130) || (!__CUDA_ARCH__)
__device__ __host__ static inline bool is_one_double(double x) { return (x==1.0); }
__device__ __host__ static inline bool is_one_cuDoubleComplex(cuDoubleComplex x) { return ((x.x==1.0)&&(x.y==0.0));}
#endif


__global__ void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern_alpha)
	(RES_VALUE_TYPE *y, int n, RES_VALUE_TYPE alpha, VALUE_TYPE* x)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	
	if (id < n)
	{
		// Since y, and x are accessed with the same offset by the same thread,
		// and the write to y follows the read of x, then x could be y.

		y[id] = CONCAT(RES_VALUE_TYPE, _mul)(alpha, CONCAT(VALUE_TYPE, _abs)(x[id]));
	}
}

__global__ void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern)
	(RES_VALUE_TYPE *y, int n, VALUE_TYPE* x)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	
	if (id < n)
	{
		// Since y, and x are accessed with the same offset by the same thread,
		// and the write to y follows the read of x, then x could be y.

		y[id] = CONCAT(VALUE_TYPE, _abs)(x[id]);
	}
}

void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
	(spgpuHandle_t handle, RES_VALUE_TYPE *y, int n, RES_VALUE_TYPE alpha, VALUE_TYPE* x)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;

	dim3 block(BLOCK_SIZE);
	dim3 grid(msize);
	

	if (CONCAT(is_one_,RES_VALUE_TYPE)(alpha))
		CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern)<<<grid, block, 0, handle->currentStream>>>(y, n, x);
	else
		CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern_alpha)<<<grid, block, 0, handle->currentStream>>>(y, n, alpha, x);

}

void 
GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle,
	__device RES_VALUE_TYPE *y,
	int n,
	RES_VALUE_TYPE alpha,
	__device VALUE_TYPE *x)
{
	int maxNForACall = max(handle->maxGridSizeX, BLOCK_SIZE*handle->maxGridSizeX);
	while (n > maxNForACall) //managing large vectors
	{
		CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
			(handle, y, maxNForACall, alpha, x);
		x = x + maxNForACall;
		y = y + maxNForACall;
		n -= maxNForACall;
	}
	
    CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_) (handle, y, n, alpha, x);
	cudaCheckError("CUDA error on abs");
}
