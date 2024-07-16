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

#include "stdio.h"
#include "cudalang.h"
#include "cudadebug.h"
#include "cuComplex.h"



extern "C"
{
#include "core.h"
#include "vector.h"
}


//#define USE_CUBLAS

#define BLOCK_SIZE 512

//#define ASSUME_LOCK_SYNC_PARALLELISM

static __device__ cuDoubleComplex ddotReductionResult[128];

__global__ void spgpuZdot_kern(int n, cuDoubleComplex* x, cuDoubleComplex* y)
{
	__shared__ cuDoubleComplex sSum[BLOCK_SIZE];

	cuDoubleComplex res = make_cuDoubleComplex(0.0, 0.0);

	cuDoubleComplex* lastX = x + n;

	x += threadIdx.x + blockIdx.x*BLOCK_SIZE;
	y += threadIdx.x + blockIdx.x*BLOCK_SIZE;

	int blockOffset = gridDim.x*BLOCK_SIZE;

	while (x < lastX)
    {
		res = cuCfma(x[0], y[0], res);
		
		x += blockOffset;
		y += blockOffset;

	}

	if (threadIdx.x >= 32)
		sSum[threadIdx.x] = res;

	__syncthreads();


	// Start reduction!

	if (threadIdx.x < 32) 
	{
		for (int i=1; i<BLOCK_SIZE/32; ++i)
		{
			res = cuCadd(res, sSum[i*32 + threadIdx.x]);
		}

	//useless (because inter-warp)
#ifndef	ASSUME_LOCK_SYNC_PARALLELISM
	}
	__syncthreads(); 

	if (threadIdx.x < 32) 
	{
#endif	

		cuDoubleComplex* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) vsSum[threadIdx.x] = cuCadd(vsSum[threadIdx.x], vsSum[threadIdx.x + 16]);
		__syncthreads();
		if (threadIdx.x < 8) vsSum[threadIdx.x] = cuCadd(vsSum[threadIdx.x], vsSum[threadIdx.x + 8]);
		__syncthreads();
		if (threadIdx.x < 4) vsSum[threadIdx.x] = cuCadd(vsSum[threadIdx.x], vsSum[threadIdx.x + 4]);
		__syncthreads();
		if (threadIdx.x < 2) vsSum[threadIdx.x] = cuCadd(vsSum[threadIdx.x], vsSum[threadIdx.x + 2]);
		__syncthreads();
		if (threadIdx.x == 0)
		ddotReductionResult[blockIdx.x] = cuCadd(vsSum[0], vsSum[1]);

	}
}

cuDoubleComplex spgpuZdot(spgpuHandle_t handle, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b)
{
#ifdef USE_CUBLAS
	cuDoubleComplex res;
	cublasDdot(n,x,1,y,1,&res);
	cudaDeviceSynchronize();
	
	return res;
#else
	cuDoubleComplex res = make_cuDoubleComplex(0.0, 0.0);

#if 0 	
	int device;
	cudaGetDevice(&device);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop,device);	

	int blocks = min(128, min(prop.multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#else
	int blocks = min(128, min(handle->multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#endif
	
	cuDoubleComplex tRes[128];

	spgpuZdot_kern<<<blocks, BLOCK_SIZE, 0, handle->currentStream>>>(n, a, b);
	cudaMemcpyFromSymbol(tRes, ddotReductionResult,blocks*sizeof(cuDoubleComplex));

	for (int i=0; i<blocks; ++i)
	{
		res = cuCadd(res, tRes[i]);
	}

	cudaCheckError("CUDA error on ddot");
	
	return res;
#endif
}

void spgpuZmdot(spgpuHandle_t handle, cuDoubleComplex* y, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b, int count, int pitch)
{
	for (int i=0; i<count; ++i)
	{
		y[i] = spgpuZdot(handle, n, a, b);
		a += pitch;
		b += pitch;
	}
}
