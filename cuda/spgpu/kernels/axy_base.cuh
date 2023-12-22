/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014
 *     Salvatore Filippone - University of Rome Tor Vergata
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
#define GEN_SPGPU_FUNC_NAME(x) CONCAT(CONCAT(spgpu,x),axy)
#define GEN_SPGPU_MFUNC_NAME(x) CONCAT(CONCAT(spgpu,x),maxy)
#define GEN_SPGPU_FUNC_NAME_2(x) CONCAT(CONCAT(spgpu,x),axypbz)
#define GEN_SPGPU_MFUNC_NAME_2(x) CONCAT(CONCAT(spgpu,x),maxypbz)
#define GEN_SPGPU_SCAL_NAME(x) CONCAT(CONCAT(spgpu,x),scal)

#define BLOCK_SIZE 256

// Define:
//#define VALUE_TYPE
//#define TYPE_SYMBOL

#include "mathbase.cuh"

__global__ void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern)
	(VALUE_TYPE *z, int n, VALUE_TYPE alpha, VALUE_TYPE* x, VALUE_TYPE* y)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	
	if (id < n)
	{
		// Since z, x and y are accessed with the same offset by the same thread,
		// and the write to z follows the x and y reads, x, y and z can share the same base address (in-place computing).

		z[id] = CONCAT(VALUE_TYPE, _mul)(alpha, CONCAT(VALUE_TYPE, _mul)(x[id], y[id]));
	}
}

void 
CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
	(spgpuHandle_t handle, VALUE_TYPE *z, int n, VALUE_TYPE alpha, VALUE_TYPE* x, VALUE_TYPE* y)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;

	dim3 block(BLOCK_SIZE);
	dim3 grid(msize);

	CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_kern)<<<grid, block, 0, handle->currentStream>>>(z, n, alpha, x, y);

}

void 
GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle,
	__device VALUE_TYPE *z,
	int n,
	VALUE_TYPE alpha,
	__device VALUE_TYPE *x,
	__device VALUE_TYPE *y)
{
	int maxNForACall = max(handle->maxGridSizeX, BLOCK_SIZE*handle->maxGridSizeX);
	while (n > maxNForACall) //managing large vectors
	{
		CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_)
			(handle, z, maxNForACall, alpha, x, y);
		x = x + maxNForACall;
		y = y + maxNForACall;
		z = z + maxNForACall;
		n -= maxNForACall;
	}
	
    	CONCAT(GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL),_) (handle, z, n, alpha, x, y);
	cudaCheckError("CUDA error on axy");
}

void 
GEN_SPGPU_MFUNC_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle,
	__device VALUE_TYPE *z,
	int n,
	VALUE_TYPE alpha,
	__device VALUE_TYPE* x,
	__device VALUE_TYPE *y,
	int count,
	int pitch)
{
	for (int i=0; i<count; i++)
	{
    		GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)(handle, z, n, alpha, x, y);
		
		x += pitch;
		y += pitch;
		z += pitch;
	}
}

__global__ void 
CONCAT(GEN_SPGPU_FUNC_NAME_2(TYPE_SYMBOL),_kern)
(VALUE_TYPE *w, int n, VALUE_TYPE beta, VALUE_TYPE* z, VALUE_TYPE alpha, VALUE_TYPE* x, VALUE_TYPE* y)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	
	if (id < n)
	{
		// Since w, x and y and z are accessed with the same offset by the same thread,
		// and the write to z follows the x, y and z reads, x, y, z and w can share the same base address (in-place computing).
		w[id] = CONCAT(VALUE_TYPE, _fma) (
				alpha,
				CONCAT(VALUE_TYPE, _mul)(x[id],y[id]), 
				CONCAT(VALUE_TYPE, _mul)(beta,z[id]));
	}
}



void 
CONCAT(GEN_SPGPU_FUNC_NAME_2(TYPE_SYMBOL),_)
	(spgpuHandle_t handle, VALUE_TYPE *w, int n, VALUE_TYPE beta, VALUE_TYPE* z, VALUE_TYPE alpha, VALUE_TYPE* x, VALUE_TYPE* y)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;

	dim3 block(BLOCK_SIZE);
	dim3 grid(msize);
	
	CONCAT(GEN_SPGPU_FUNC_NAME_2(TYPE_SYMBOL),_kern)
		<<<grid, block, 0, handle->currentStream>>>(w, n, beta, z, alpha, x, y);

}

void
GEN_SPGPU_FUNC_NAME_2(TYPE_SYMBOL)
	(spgpuHandle_t handle,
	__device VALUE_TYPE *w,
	int n,
	VALUE_TYPE beta,
	__device VALUE_TYPE *z,
	VALUE_TYPE alpha,
	__device VALUE_TYPE* x,
	__device VALUE_TYPE *y
	)
{

	if (CONCAT(VALUE_TYPE, _isZero(alpha)))
	{
		GEN_SPGPU_SCAL_NAME(TYPE_SYMBOL)
			(handle, w, n, beta, z);
	}
	else if (CONCAT(VALUE_TYPE, _isZero(beta))) {
		GEN_SPGPU_FUNC_NAME(TYPE_SYMBOL)
			(handle, w, n, alpha, x, y);
	} 
	else {
		int maxNForACall = max(handle->maxGridSizeX, BLOCK_SIZE*handle->maxGridSizeX);
		
		while (n > maxNForACall) //managing large vectors
		{
			
			CONCAT(GEN_SPGPU_FUNC_NAME_2(TYPE_SYMBOL),_)
				(handle, w, maxNForACall, beta, z, alpha, x, y);
	
			x = x + maxNForACall;
			y = y + maxNForACall;
			z = z + maxNForACall;
			w = w + maxNForACall;
			n -= maxNForACall;
		}
    
		CONCAT(GEN_SPGPU_FUNC_NAME_2(TYPE_SYMBOL),_)
			(handle, w, n, beta, z, alpha, x, y);
    }	
  
	cudaCheckError("CUDA error on axypbz");
}

void
GEN_SPGPU_MFUNC_NAME_2(TYPE_SYMBOL)
	(spgpuHandle_t handle,
	__device VALUE_TYPE *w,
	int n,
	VALUE_TYPE beta,
	__device VALUE_TYPE *z,
	VALUE_TYPE alpha,
	__device VALUE_TYPE* x,
	__device VALUE_TYPE *y,
	int count,
	int pitch)
{
  for (int i=0; i<count; i++)
    {
	GEN_SPGPU_FUNC_NAME_2(TYPE_SYMBOL)(handle, w, n, beta, z, alpha, x, y);
		
	x += pitch;
	y += pitch;
	z += pitch;
	w += pitch;
    }
}
