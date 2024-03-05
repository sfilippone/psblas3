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

__global__ void spgpuSxyzw_krn(int n, float  a, float  b,
			       float  c, float  d,
			       float  e, float  f,
			       float * x, float  *y,
			       float  *z, float  *w)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	float  ty, tz;
	for ( ; id < n; id +=gridSize)
		//if (id,n) 
	{

	  ty    = PREC_FADD(PREC_FMUL(a, x[id]), PREC_FMUL(b,y[id]));
	  tz    = PREC_FADD(PREC_FMUL(c, ty), PREC_FMUL(d,z[id]));
	  w[id] = PREC_FADD(PREC_FMUL(e, tz), PREC_FMUL(f,w[id]));
	  y[id] = ty;
	  z[id] = tz;
	}
}


void spgpuSxyzw(spgpuHandle_t handle,
		int n,
		float  a, float  b,
		float  c, float  d,
		float  e, float  f,
		__device float * x,
		__device float * y,
		__device float * z,
		__device float *w)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
	num_mp         = getGPUMultiProcessors();
	max_threads_mp = getGPUMaxThreadsPerMP();
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuSxyzw_krn<<<grid, block, 0, handle->currentStream>>>(n, a,b,c,d,e,f,
								  x, y, z,w);
}

