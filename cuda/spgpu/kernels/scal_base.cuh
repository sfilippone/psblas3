/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014
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
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu,x),scal)

#define THREAD_BLOCK 256

// Define:
//#define VALUE_TYPE
//#define TYPE_SYMBOL

#include "mathbase.cuh"

#define BLOCK_SIZE 256

__global__ void
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL), _krn)
(VALUE_TYPE *y, int n, VALUE_TYPE alpha, VALUE_TYPE* x)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;

	// Since x and y are accessed with the same offset by the same thread,
	// and the write to y follows the x read, x and y can share the same base address (in-place computing).	
	if (id < n)
	{
		y[id] = CONCAT(VALUE_TYPE, _mul)(alpha, x[id]);
	}
}

void 
CONCAT(_,GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL))
	(spgpuHandle_t handle, VALUE_TYPE *y, int n, VALUE_TYPE alpha, VALUE_TYPE* x)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;

	dim3 block(BLOCK_SIZE);
	dim3 grid(msize);
	
	CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL), _krn)
		<<<grid, block, 0, handle->currentStream>>>(y, n, alpha, x);
}

void 
GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)
	(spgpuHandle_t handle,
	__device VALUE_TYPE *y,
	int n,
	VALUE_TYPE alpha,
	__device VALUE_TYPE *x)
{

	int maxNForACall = max(handle->maxGridSizeX, THREAD_BLOCK*handle->maxGridSizeX);

	while (n > maxNForACall) //managing large vectors
    	{
		CONCAT(_,GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL))(handle, y, maxNForACall, alpha, x);
		x = x + maxNForACall;
		y = y + maxNForACall;
		n -= maxNForACall;
	}
	
	CONCAT(_,GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL))(handle, y, n, alpha, x);
	
	cudaCheckError("CUDA error on scal");
}
