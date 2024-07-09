/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
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
#include "cudadebug.h"
#include "cudalang.h"

extern "C"
{
#include "core.h"
#include "vector.h"
  int getGPUMultiProcessors();
  int getGPUMaxThreadsPerMP();
  //#include "cuda_util.h"
}


#include "debug.h"

#define BLOCK_SIZE 512


#if 1
__global__ void spgpuSaxpby_krn(float *z, int n, float beta, float *y, float alpha, float* x)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	if (beta == 0.0f) {
	  for ( ; id < n; id +=gridSize)
	    {
	      // Since z, x and y are accessed with the same offset by the same thread,
	      // and the write to z follows the x and y read, x, y and z can share the same base address (in-place computing).
	      
	      z[id] = PREC_FMUL(alpha,x[id]);
	    }
	} else {
	  for ( ; id < n; id +=gridSize)
	    {
	      // Since z, x and y are accessed with the same offset by the same thread,
	      // and the write to z follows the x and y read, x, y and z can share the same base address (in-place computing).
	      z[id] = PREC_FADD(PREC_FMUL(alpha, x[id]), PREC_FMUL(beta,y[id]));
	    }
	}
}

void spgpuSaxpby(spgpuHandle_t handle,
	__device float *z,
	int n,
	float beta,
	__device float *y,
	float alpha,
	__device float* x)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
	num_mp         = getGPUMultiProcessors();
	max_threads_mp = getGPUMaxThreadsPerMP();
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuSaxpby_krn<<<grid, block, 0, handle->currentStream>>>(z, n, beta, y, alpha, x);
}

#else

__global__ void spgpuSaxpby_krn(float *z, int n, float beta, float *y, float alpha, float* x)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	
	if (id < n)
	{
		// Since z, x and y are accessed with the same offset by the same thread,
		// and the write to z follows the x and y read, x, y and z can share the same base address (in-place computing).

		if (beta == 0.0f)
			z[id] = PREC_FMUL(alpha,x[id]);
		else
			z[id] = PREC_FADD(PREC_FMUL(alpha, x[id]), PREC_FMUL(beta,y[id]));
	}
}



void spgpuSaxpby_(spgpuHandle_t handle,
	__device float *z,
	int n,
	float beta,
	__device float *y,
	float alpha,
	__device float* x)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;

	dim3 block(BLOCK_SIZE);
	dim3 grid(msize);

	spgpuSaxpby_krn<<<grid, block, 0, handle->currentStream>>>(z, n, beta, y, alpha, x);
}

void spgpuSaxpby(spgpuHandle_t handle,
	__device float *z,
	int n,
	float beta,
	__device float *y,
	float alpha,
	__device float* x)
{
	int maxNForACall = max(handle->maxGridSizeX, BLOCK_SIZE*handle->maxGridSizeX);
	while (n > maxNForACall) //managing large vectors
	{
		spgpuSaxpby_(handle, z, maxNForACall, beta, y, alpha, x);
		
		x = x + maxNForACall;
		y = y + maxNForACall;
		z = z + maxNForACall;
		n -= maxNForACall;
	}
	
	spgpuSaxpby_(handle, z, n, beta, y, alpha, x);

	cudaCheckError("CUDA error on saxpby");
}
#endif
void spgpuSmaxpby(spgpuHandle_t handle,
		  __device float *z,
		  int n,
		  float beta,
		  __device float *y,
		  float alpha,
		  __device float* x, 
		  int count, int pitch)
{

  for (int i=0; i<count; i++)
    spgpuSaxpby(handle, z+pitch*i, n, beta, y+pitch*i, alpha, x+pitch*i);
  
}
