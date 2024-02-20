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
#include <cuda_runtime.h>

extern "C"
{
#include "core.h"
#include "vector.h"
  int getGPUMultiProcessors();
  int getGPUMaxThreadsPerMP();
}


#include "debug.h"

#define BLOCK_SIZE 512

__global__ void spgpuCabgdxyz_krn(int n, cuFloatComplex  alpha, cuFloatComplex  beta,
				  cuFloatComplex  gamma, cuFloatComplex  delta,
				  cuFloatComplex * x, cuFloatComplex  *y, cuFloatComplex  *z)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	cuFloatComplex  t;
	for ( ; id < n; id +=gridSize)
		//if (id,n) 
	{

	  if (cuFloatComplex_isZero(beta)) 
	    t = cuCmulf(alpha,x[id]);
	  else
	    t = cuCfmaf(alpha, x[id], cuCmulf(beta,y[id]));
	  if (cuFloatComplex_isZero(delta))
	    z[id] = cuCmulf(gamma, t);
	  else
	    z[id] = cuCfmaf(gamma, t, cuCmulf(delta,z[id]));
	  y[id] = t;
	}
}


void spgpuCabgdxyz(spgpuHandle_t handle,
		   int n,
		   cuFloatComplex  alpha,
		   cuFloatComplex  beta,
		   cuFloatComplex  gamma,
		   cuFloatComplex  delta,
		   __device cuFloatComplex * x,
		   __device cuFloatComplex * y,
		   __device cuFloatComplex  *z)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
	num_mp         = getGPUMultiProcessors();
	max_threads_mp = getGPUMaxThreadsPerMP();
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuCabgdxyz_krn<<<grid, block, 0, handle->currentStream>>>(n, alpha, beta, gamma, delta,
								   x, y, z);
}

