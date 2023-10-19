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

#define BLOCK_SIZE 320
//#define BLOCK_SIZE 512

//#define ASSUME_LOCK_SYNC_PARALLELISM


static __device__ cuFloatComplex sdotReductionResult[128];

__global__ void spgpuCdot_kern(int n, cuFloatComplex* x, cuFloatComplex* y)
{
	__shared__ cuFloatComplex sSum[BLOCK_SIZE];

	cuFloatComplex res = make_cuFloatComplex(0.0f, 0.0f);

	cuFloatComplex* lastX = x + n;

	x += threadIdx.x + blockIdx.x*BLOCK_SIZE;
	y += threadIdx.x + blockIdx.x*BLOCK_SIZE;

	int blockOffset = gridDim.x*BLOCK_SIZE;

	int numSteps = (lastX - x + blockOffset - 1)/blockOffset;

	// prefetching
	for (int j = 0; j < numSteps / 2; j++)
    {
		cuFloatComplex x1 = x[0]; x += blockOffset;
		cuFloatComplex y1 = y[0]; y += blockOffset;
		cuFloatComplex x2 = x[0]; x += blockOffset;
		cuFloatComplex y2 = y[0]; y += blockOffset;

		res = cuCfmaf(x1, y1, res);
		res = cuCfmaf(x2, y2, res);

	}

	if (numSteps % 2)
	{
		res = cuCfmaf(*x, *y, res);
	}

	if (threadIdx.x >= 32)
		sSum[threadIdx.x] = res;

	__syncthreads();


	// Start reduction!

	if (threadIdx.x < 32) 
	{
		for (int i=1; i<BLOCK_SIZE/32; ++i)
		{
			res = cuCaddf(res, sSum[i*32 + threadIdx.x]);
		}

	//useless (because inter-warp)
#ifndef	ASSUME_LOCK_SYNC_PARALLELISM
	}
	__syncthreads(); 

	if (threadIdx.x < 32) 
	{
#endif	

		cuFloatComplex* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) vsSum[threadIdx.x] = cuCaddf(vsSum[threadIdx.x], vsSum[threadIdx.x + 16]);
		__syncthreads();
		if (threadIdx.x < 8) vsSum[threadIdx.x] = cuCaddf(vsSum[threadIdx.x], vsSum[threadIdx.x + 8]);
		__syncthreads();
		if (threadIdx.x < 4) vsSum[threadIdx.x] = cuCaddf(vsSum[threadIdx.x], vsSum[threadIdx.x + 4]);
		__syncthreads();
		if (threadIdx.x < 2) vsSum[threadIdx.x] = cuCaddf(vsSum[threadIdx.x], vsSum[threadIdx.x + 2]);
		__syncthreads();
		if (threadIdx.x == 0)
			sdotReductionResult[blockIdx.x] = cuCaddf(vsSum[0], vsSum[1]);
	}
}

cuFloatComplex spgpuCdot(spgpuHandle_t handle, int n, __device cuFloatComplex* a, __device cuFloatComplex* b)
{

#ifdef USE_CUBLAS
	cuFloatComplex res;
	cublasSdot(n,x,1,y,1,&res);
	cudaDeviceSynchronize();
	
	return res;
#else
	cuFloatComplex res = make_cuFloatComplex(0.0f, 0.0f);

#if 0 	
	int device;
	cudaGetDevice(&device);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop,device);	

	int blocks = min(128, min(prop.multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#else
	int blocks = min(128, min(handle->multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#endif

	
	cuFloatComplex tRes[128];

	spgpuCdot_kern<<<blocks, (BLOCK_SIZE), 0, handle->currentStream>>>(n, a, b);
	cudaMemcpyFromSymbol(tRes, sdotReductionResult, blocks*sizeof(cuFloatComplex));

	for (int i=0; i<blocks; ++i)
	{
		res = cuCaddf(res, tRes[i]);
	}

	cudaCheckError("CUDA error on sdot (blocks: %i, regs per block: %i)\n", blocks, prop.regsPerBlock);
	
	return res;
#endif
}

void spgpuCmdot(spgpuHandle_t handle, cuFloatComplex* y, int n, __device cuFloatComplex* a, __device cuFloatComplex* b, int count, int pitch)
{
	for (int i=0; i<count; ++i)
	{
		y[i] = spgpuCdot(handle, n, a, b);
		a += pitch;
		b += pitch;
	}
}
