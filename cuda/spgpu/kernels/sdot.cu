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

extern "C"
{
#include "core.h"
#include "vector.h"
}

//#define USE_CUBLAS

#ifdef USE_CUBLAS
#include "cublas.h"
#endif

#define BLOCK_SIZE 320
//#define BLOCK_SIZE 512

//#define ASSUME_LOCK_SYNC_PARALLELISM

#ifndef USE_CUBLAS
static __device__ float sdotReductionResult[128];
#endif

__global__ void spgpuSdot_kern(int n, float* x, float* y)
{
	__shared__ float sSum[BLOCK_SIZE];

	float res = 0;

	float* lastX = x + n;

	x += threadIdx.x + blockIdx.x*BLOCK_SIZE;
	y += threadIdx.x + blockIdx.x*BLOCK_SIZE;

	int blockOffset = gridDim.x*BLOCK_SIZE;

	int numSteps = (lastX - x + blockOffset - 1)/blockOffset;

	// prefetching
	for (int j = 0; j < numSteps / 2; j++)
    {
		float x1 = x[0]; x += blockOffset;
		float y1 = y[0]; y += blockOffset;
		float x2 = x[0]; x += blockOffset;
		float y2 = y[0]; y += blockOffset;

		res = PREC_FADD(res, PREC_FMUL(x1,y1));
		res = PREC_FADD(res, PREC_FMUL(x2,y2));

	}

	if (numSteps % 2)
	{
		res = PREC_FADD(res, PREC_FMUL(*x,*y));
	}

	if (threadIdx.x >= 32)
		sSum[threadIdx.x] = res;

	__syncthreads();


	// Start reduction!

	if (threadIdx.x < 32) 
	{
		for (int i=1; i<BLOCK_SIZE/32; ++i)
		{
			res += sSum[i*32 + threadIdx.x];
		}

	//useless (because inter-warp)
#ifndef	ASSUME_LOCK_SYNC_PARALLELISM
	}
	__syncthreads(); 

	if (threadIdx.x < 32) 
	{
#endif	

#ifdef ASSUME_LOCK_SYNC_PARALLELISM
		volatile float* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) vsSum[threadIdx.x] += vsSum[threadIdx.x + 16];
		if (threadIdx.x < 8) vsSum[threadIdx.x] += vsSum[threadIdx.x + 8];
		if (threadIdx.x < 4) vsSum[threadIdx.x] += vsSum[threadIdx.x + 4];
		if (threadIdx.x < 2) vsSum[threadIdx.x] += vsSum[threadIdx.x + 2];
		if (threadIdx.x == 0)
			sdotReductionResult[blockIdx.x] = vsSum[0] + vsSum[1];

#else
		float* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) vsSum[threadIdx.x] += vsSum[threadIdx.x + 16];
		__syncthreads();
		if (threadIdx.x < 8) vsSum[threadIdx.x] += vsSum[threadIdx.x + 8];
		__syncthreads();
		if (threadIdx.x < 4) vsSum[threadIdx.x] += vsSum[threadIdx.x + 4];
		__syncthreads();
		if (threadIdx.x < 2) vsSum[threadIdx.x] += vsSum[threadIdx.x + 2];
		__syncthreads();
		if (threadIdx.x == 0)
		sdotReductionResult[blockIdx.x] = vsSum[0] + vsSum[1];
#endif
	}
}

float spgpuSdot(spgpuHandle_t handle, int n, __device float* a, __device float* b)
{

#ifdef USE_CUBLAS
	float res;
	cublasSdot(n,a,1,b,1,&res);
	cudaDeviceSynchronize();
	
	return res;
#else
	float res = 0;

#if 0 	
	int device;
	cudaGetDevice(&device);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop,device);	

	int blocks = min(128, min(prop.multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#else
	int blocks = min(128, min(handle->multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#endif
	
	float tRes[128];
	spgpuSdot_kern<<<blocks, (BLOCK_SIZE), 0, handle->currentStream>>>(n, a, b);
	cudaMemcpyFromSymbol(tRes, sdotReductionResult, blocks*sizeof(float));

	for (int i=0; i<blocks; ++i)
	{
		res += tRes[i];
	}

	cudaCheckError("CUDA error on sdot (blocks: %i, regs per block: %i)\n", blocks, prop.regsPerBlock);
	
	return res;
#endif
}

void spgpuSmdot(spgpuHandle_t handle, float* y, int n, __device float* a, __device float* b, int count, int pitch)
{
	for (int i=0; i<count; ++i)
	{
		y[i] = spgpuSdot(handle, n, a, b);
		a += pitch;
		b += pitch;
	}
}
