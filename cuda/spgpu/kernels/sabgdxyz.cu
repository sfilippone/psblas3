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

__global__ void spgpuSabgdxyz_krn(int n, float alpha, float beta, float gamma, float delta,
				  float* x, float *y, float *z)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	float t;
	for ( ; id < n; id +=gridSize)
		//if (id,n) 
	{

	  if (beta == 0.0)
	    t = PREC_FMUL(alpha,x[id]);
	  else
	    t = PREC_FADD(PREC_FMUL(alpha, x[id]), PREC_FMUL(beta,y[id]));
	  if (delta == 0.0)
	    z[id] = gamma * t;
	  else
	    z[id] = PREC_FADD(PREC_FMUL(gamma, t), PREC_FMUL(delta,z[id]));
	  y[id] = t;
	}
}


void spgpuSabgdxyz(spgpuHandle_t handle,
		   int n,
		   float alpha,
		   float beta,
		   float gamma,
		   float delta,
		   __device float* x,
		   __device float* y,
		   __device float *z)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
	num_mp         = getGPUMultiProcessors();
	max_threads_mp = getGPUMaxThreadsPerMP();
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuSabgdxyz_krn<<<grid, block, 0, handle->currentStream>>>(n, alpha, beta, gamma, delta,
								   x, y, z);
}

