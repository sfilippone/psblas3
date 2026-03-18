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

__global__ void spgpuDupd_xyz_krn(int n, double alpha, double beta, double gamma, double delta,
				  double* x, double *y, double *z)
{
	int id = threadIdx.x + BLOCK_SIZE*blockIdx.x;
	unsigned int gridSize = blockDim.x * gridDim.x;
	double t;
	for ( ; id < n; id +=gridSize)
		//if (id,n) 
	{

	  if (beta == 0.0)
	    t = PREC_DMUL(alpha,x[id]);
	  else
	    t = PREC_DADD(PREC_DMUL(alpha, x[id]), PREC_DMUL(beta,y[id]));
	  if (delta == 0.0)
	    z[id] = gamma * t;
	  else
	    z[id] = PREC_DADD(PREC_DMUL(gamma, t), PREC_DMUL(delta,z[id]));
	  y[id] = t;
	}
}


void spgpuDupd_xyz(spgpuHandle_t handle,
		   int n,
		   double alpha,
		   double beta,
		   double gamma,
		   double delta,
		   __device double* x,
		   __device double* y,
		   __device double *z)
{
	int num_mp, max_threads_mp, num_blocks_mp, num_blocks;
	dim3 block(BLOCK_SIZE);
	num_mp         = getGPUMultiProcessors();
	max_threads_mp = getGPUMaxThreadsPerMP();
	num_blocks_mp  = max_threads_mp/BLOCK_SIZE;
	num_blocks     = num_blocks_mp*num_mp;
	dim3 grid(num_blocks);

	spgpuDupd_xyz_krn<<<grid, block, 0, handle->currentStream>>>(n, alpha, beta, gamma, delta,
								   x, y, z);
}

