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

__global__ void spgpuZxyzw_krn(int n, cuDoubleComplex  a, cuDoubleComplex  b,
			       cuDoubleComplex  c, cuDoubleComplex  d,
			       cuDoubleComplex  e, cuDoubleComplex  f,
			       cuDoubleComplex * x, cuDoubleComplex  *y,
			       cuDoubleComplex  *z, cuDoubleComplex  *w)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	cuDoubleComplex  ty, tz;
	for ( ; id < n; id +=gridSize)
		//if (id,n) 
	{

	  ty    = cuCfma(a, x[id], cuCmul(b,y[id]));
	  tz    = cuCfma(c, ty, cuCmul(d,z[id]));
	  w[id] = cuCfma(e, tz, cuCmul(f,w[id]));
	  y[id] = ty;
	  z[id] = tz;
	}
}


void spgpuZxyzw(spgpuHandle_t handle,
		int n,
		cuDoubleComplex  a, cuDoubleComplex  b,
		cuDoubleComplex  c, cuDoubleComplex  d,
		cuDoubleComplex  e, cuDoubleComplex  f,
		__device cuDoubleComplex * x,
		__device cuDoubleComplex * y,
		__device cuDoubleComplex * z,
		__device cuDoubleComplex *w)
{
	int msize = (n+BLOCK_SIZE-1)/BLOCK_SIZE;
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
	num_mp         = getGPUMultiProcessors();
	max_threads_mp = getGPUMaxThreadsPerMP();
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuZxyzw_krn<<<grid, block, 0, handle->currentStream>>>(n, a,b,c,d,e,f,
								  x, y, z,w);
}

