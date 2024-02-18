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
#include "cudadebug.h"
#include "cudalang.h"
#include "cuComplex.h"


extern "C"
{
#include "core.h"
#include "vector.h"
}


#include "debug.h"

#define BLOCK_SIZE 512

__global__ void spgpuCaxpby_krn(cuFloatComplex *z, int n, cuFloatComplex beta, cuFloatComplex *y, cuFloatComplex alpha, cuFloatComplex* x)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	if (cuFloatComplex_isZero(beta)) {
	  for ( ; id < n; id +=gridSize)
	    //if (id,n) 
	    {
	      // Since z, x and y are accessed with the same offset by the same thread,
	      // and the write to z follows the x and y read, x, y and z can share the same base address (in-place computing).
	      
	      z[id] = cuCmulf(alpha,x[id]);
	    }
	} else {
	  for ( ; id < n; id +=gridSize)
	    //if (id,n) 
	    {
	      z[id] = cuCfmaf(beta, y[id], cuCmulf(alpha, x[id]));
	    }
	}
}

#if 1
void spgpuCaxpby(spgpuHandle_t handle,
	__device cuFloatComplex *z,
	int n,
	cuFloatComplex beta,
	__device cuFloatComplex *y,
	cuFloatComplex alpha,
	__device cuFloatComplex* x)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, 0);
	num_mp         = deviceProp.multiProcessorCount;
	max_threads_mp = deviceProp.maxThreadsPerMultiProcessor;
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuCaxpby_krn<<<grid, block, 0, handle->currentStream>>>(z, n, beta, y, alpha, x);
}

#else
void spgpuCaxpby_(spgpuHandle_t handle,
	__device cuFloatComplex *z,
	int n,
	cuFloatComplex beta,
	__device cuFloatComplex *y,
	cuFloatComplex alpha,
	__device cuFloatComplex* x)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, 0);
	num_mp         = deviceProp.multiProcessorCount;
	max_threads_mp = deviceProp.maxThreadsPerMultiProcessor;
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuCaxpby_krn<<<grid, block, 0, handle->currentStream>>>(z, n, beta, y, alpha, x);
}

void spgpuCaxpby(spgpuHandle_t handle,
	__device cuFloatComplex *z,
	int n,
	cuFloatComplex beta,
	__device cuFloatComplex *y,
	cuFloatComplex alpha,
	__device cuFloatComplex* x)
{
	int maxNForACall = max(handle->maxGridSizeX, BLOCK_SIZE*handle->maxGridSizeX);

	while (n > maxNForACall) //managing large vectors
	{
		spgpuCaxpby_(handle, z, maxNForACall, beta, y, alpha, x);
		
		x = x + maxNForACall;
		y = y + maxNForACall;
		z = z + maxNForACall;
		n -= maxNForACall;
	}
	
	spgpuCaxpby_(handle, z, n, beta, y, alpha, x);

	cudaCheckError("CUDA error on saxpby");
}
#endif
void spgpuCmaxpby(spgpuHandle_t handle,
		  __device cuFloatComplex *z,
		  int n,
		  cuFloatComplex beta,
		  __device cuFloatComplex *y,
		  cuFloatComplex alpha,
		  __device cuFloatComplex* x, 
		  int count, int pitch)
{

  for (int i=0; i<count; i++)
    spgpuCaxpby(handle, z+pitch*i, n, beta, y+pitch*i, alpha, x+pitch*i);
  
}
