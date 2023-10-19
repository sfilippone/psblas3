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

#define BLOCK_SIZE 320
//#define BLOCK_SIZE 512

//#define ASSUME_LOCK_SYNC_PARALLELISM

static __device__ float snrm2ReductionResult[128];

__global__ void spgpuSnrm2_kern(int n, float* x)
{
	__shared__ float sSum[BLOCK_SIZE];

	float res = 0;

	float* lastX = x + n;

	x += threadIdx.x + blockIdx.x*BLOCK_SIZE;

	int blockOffset = gridDim.x*BLOCK_SIZE;

	int numSteps = (lastX - x + blockOffset - 1)/blockOffset;

	// prefetching
	for (int j = 0; j < numSteps / 2; j++)
    {
		float x1 = x[0]; x += blockOffset;
		float x2 = x[0]; x += blockOffset;

		res = PREC_FADD(res, PREC_FMUL(x1,x1));
		res = PREC_FADD(res, PREC_FMUL(x2,x2));

	}

	if (numSteps % 2)
	{
		float x1 = x[0];
		res = PREC_FADD(res, PREC_FMUL(x1,x1));
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
			snrm2ReductionResult[blockIdx.x] = vsSum[0] + vsSum[1];

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
		snrm2ReductionResult[blockIdx.x] = vsSum[0] + vsSum[1];
#endif
	}
}

float spgpuSnrm2(spgpuHandle_t handle, int n, __device float* x)
{
#ifdef USE_CUBLAS
	float res;
	cublasSnrm2(n,x,1,&res);
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

	spgpuSnrm2_kern<<<blocks, BLOCK_SIZE, 0, handle->currentStream>>>(n, x);
	cudaMemcpyFromSymbol(tRes, snrm2ReductionResult,blocks*sizeof(float));

	for (int i=0; i<blocks; ++i)
	{
		res += tRes[i];
	}

	cudaCheckError("CUDA error on snrm2");
	
	return sqrtf(res);
#endif
}

void spgpuSmnrm2(spgpuHandle_t handle, float *y, int n, __device float *x, int count, int pitch)
{
	for (int i=0; i < count; ++i)
	{
		y[i] = spgpuSnrm2(handle, n, x);
		x += pitch;
	}
}
