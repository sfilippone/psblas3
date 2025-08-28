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
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu,x),gath)

#define BLOCK_SIZE 256

// Define:
//#define VALUE_TYPE
//#define TYPE_SYMBOL

#include "mathbase.cuh"

__global__ void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern)
	(VALUE_TYPE* values, int count, const int* indices, int firstIndex, const VALUE_TYPE* vector)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	
	if (id < count)
	{	
		int pos = indices[id]-firstIndex;
		
		if (pos < 0)
			return;
		
		values[id] = vector[pos];
	}
}

void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
	(spgpuHandle_t handle, VALUE_TYPE *xValues, int xNnz,
	const __device int *xIndices, int xBaseIndex, const VALUE_TYPE* y)
{
	int msize = (xNnz+BLOCK_SIZE-1)/BLOCK_SIZE;

	dim3 block(BLOCK_SIZE);
	dim3 grid(msize);

	CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern)<<<grid, block, 0, handle->currentStream>>>(xValues, xNnz, xIndices, xBaseIndex, y);
	
}

void 
GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle,
	__device VALUE_TYPE *xValues,
	int xNnz,
	const __device int *xIndices,
	int xBaseIndex,
	const __device VALUE_TYPE* y)
{
	int maxNForACall = max(handle->maxGridSizeX, BLOCK_SIZE*handle->maxGridSizeX);
	while (xNnz > maxNForACall) //managing large vectors
	{
		CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
			(handle, xValues, maxNForACall, xIndices, xBaseIndex, y);
		xIndices += maxNForACall;
		xValues += maxNForACall;
		xNnz -= maxNForACall;
	}
	
    	CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
			(handle, xValues, xNnz, xIndices, xBaseIndex, y);
	cudaCheckError("CUDA error on gath");
}

