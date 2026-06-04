/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#include "cudadebug.h"
#include "cudalang.h"
#include <cuda_runtime.h>
#include "core.h"

extern "C"
{
#include "vector.h"
  int getGPUMultiProcessors();
  int getGPUMaxThreadsPerMP();
}


#include "debug.h"

#define BLOCK_SIZE 512

__global__ void spgpuZupd_xyz_krn(int n, cuDoubleComplex  alpha, cuDoubleComplex  beta,
				  cuDoubleComplex  gamma, cuDoubleComplex  delta,
				  cuDoubleComplex * x, cuDoubleComplex  *y, cuDoubleComplex  *z)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	cuDoubleComplex  t;
	for ( ; id < n; id +=gridSize)
		//if (id,n) 
	{

	  if (cuDoubleComplex_isZero(beta)) 
	    t = cuCmul(alpha,x[id]);
	  else
	    t = cuCfma(alpha, x[id], cuCmul(beta,y[id]));
	  if (cuDoubleComplex_isZero(delta))
	    z[id] = cuCmul(gamma, t);
	  else
	    z[id] = cuCfma(gamma, t, cuCmul(delta,z[id]));
	  y[id] = t;
	}
}


void spgpuZupd_xyz(spgpuHandle_t handle,
		   int n,
		   cuDoubleComplex  alpha,
		   cuDoubleComplex  beta,
		   cuDoubleComplex  gamma,
		   cuDoubleComplex  delta,
		   __device cuDoubleComplex * x,
		   __device cuDoubleComplex * y,
		   __device cuDoubleComplex  *z)
{
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
	num_mp         = getGPUMultiProcessors();
	max_threads_mp = getGPUMaxThreadsPerMP();
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuZupd_xyz_krn<<<grid, block, 0, handle->currentStream>>>(n, alpha, beta, gamma, delta,
								   x, y, z);
}

