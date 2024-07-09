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

__global__ void spgpuCxyzw_krn(int n, cuFloatComplex  a, cuFloatComplex  b,
			       cuFloatComplex  c, cuFloatComplex  d,
			       cuFloatComplex  e, cuFloatComplex  f,
			       cuFloatComplex * x, cuFloatComplex  *y,
			       cuFloatComplex  *z, cuFloatComplex  *w)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	cuFloatComplex  ty, tz;
	for ( ; id < n; id +=gridSize)
		//if (id,n) 
	{

	  ty    = cuCfmaf(a, x[id], cuCmulf(b,y[id]));
	  tz    = cuCfmaf(c, ty, cuCmulf(d,z[id]));
	  w[id] = cuCfmaf(e, tz, cuCmulf(f,w[id]));
	  y[id] = ty;
	  z[id] = tz;
	}
}


void spgpuCxyzw(spgpuHandle_t handle,
		int n,
		cuFloatComplex  a, cuFloatComplex  b,
		cuFloatComplex  c, cuFloatComplex  d,
		cuFloatComplex  e, cuFloatComplex  f,
		__device cuFloatComplex * x,
		__device cuFloatComplex * y,
		__device cuFloatComplex * z,
		__device cuFloatComplex *w)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
	num_mp         = getGPUMultiProcessors();
	max_threads_mp = getGPUMaxThreadsPerMP();
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuCxyzw_krn<<<grid, block, 0, handle->currentStream>>>(n, a,b,c,d,e,f,
								  x, y, z,w);
}

