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
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu,x),setscal)

#define BLOCK_SIZE 256

// Define:
//#define VALUE_TYPE
//#define TYPE_SYMBOL

#include "mathbase.cuh"

__global__ void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern)
	(VALUE_TYPE* vector, int count, VALUE_TYPE val)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	
	if (id < count)
	{	
	  vector[id] = val;
	}
}

void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
     (spgpuHandle_t handle,
      int n, VALUE_TYPE val,  __device VALUE_TYPE* y)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;

	dim3 block(BLOCK_SIZE);
	dim3 grid(msize);

	CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern)<<<grid, block, 0, handle->currentStream>>>(y, n, val);
	
}

void 
GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)
     (spgpuHandle_t handle,
      int first,
      int last,
      int baseIndex,
      VALUE_TYPE val,
      __device VALUE_TYPE* y)
{
	int maxNForACall = max(handle->maxGridSizeX, BLOCK_SIZE*handle->maxGridSizeX);
	int n=last-first+1; 
	y += (first-baseIndex);
	while (n > maxNForACall) //managing large vectors
	{
		CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
		  (handle, maxNForACall, val,y);
		y += maxNForACall;
		n -= maxNForACall;
	}
	
    	CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
	  (handle, n, val,y);
	cudaCheckError("CUDA error on scat");
}

