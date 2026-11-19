/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#include "stdio.h"
#include "cudalang.h"
#include "cudadebug.h"

#include "core.h"

extern "C"
{
#include "vector.h"
}


//#define USE_CUBLAS

#define BLOCK_SIZE 512

//#define ASSUME_LOCK_SYNC_PARALLELISM

static __device__ double ddotReductionResult[128];

__global__ void spgpuDdot_kern(int n, double* x, double* y)
{
	__shared__ double sSum[BLOCK_SIZE];

	double res = 0;

	double* lastX = x + n;

	x += threadIdx.x + blockIdx.x*BLOCK_SIZE;
	y += threadIdx.x + blockIdx.x*BLOCK_SIZE;

	int blockOffset = gridDim.x*BLOCK_SIZE;

	while (x < lastX)
    {
		res = PREC_DADD(res, PREC_DMUL(x[0], y[0]));
		
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
		volatile double* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) vsSum[threadIdx.x] += vsSum[threadIdx.x + 16];
		if (threadIdx.x < 8) vsSum[threadIdx.x] += vsSum[threadIdx.x + 8];
		if (threadIdx.x < 4) vsSum[threadIdx.x] += vsSum[threadIdx.x + 4];
		if (threadIdx.x < 2) vsSum[threadIdx.x] += vsSum[threadIdx.x + 2];
		if (threadIdx.x == 0)
			ddotReductionResult[blockIdx.x] = vsSum[0] + vsSum[1];

#else
		double* vsSum = sSum;
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
		ddotReductionResult[blockIdx.x] = vsSum[0] + vsSum[1];
#endif
	}
}

double spgpuDdot(spgpuHandle_t handle, int n, __device double* a, __device double* b)
{
#ifdef USE_CUBLAS
	double res;
	cublasDdot(n,x,1,y,1,&res);
	cudaDeviceSynchronize();
	
	return res;
#else
	double res = 0;

	int device;
	cudaGetDevice(&device);
#if 0 	
	int device;
	cudaGetDevice(&device);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop,device);	

	int blocks = min(128, min(prop.multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#else
	int blocks = min(128, min(handle->multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#endif
	
	double tRes[128];

	spgpuDdot_kern<<<blocks, BLOCK_SIZE, 0, handle->currentStream>>>(n, a, b);
	cudaMemcpyFromSymbol(tRes, ddotReductionResult,blocks*sizeof(double));

	for (int i=0; i<blocks; ++i)
	{
		res += tRes[i];
	}

	cudaCheckError("CUDA error on ddot");
	
	return res;
#endif
}

void spgpuDmdot(spgpuHandle_t handle, double* y, int n, __device double* a, __device double* b, int count, int pitch)
{
	for (int i=0; i<count; ++i)
	{
		y[i] = spgpuDdot(handle, n, a, b);
		a += pitch;
		b += pitch;
	}
}


/*
 * Fused local dot products sharing the first operand:
 *   res[0] = x.y1, res[1] = x.y2, res[2] = x.y3 (y3 may be NULL).
 * One kernel launch and one device->host transfer instead of two or three,
 * as needed by the FCG iteration (PSBLAS psb_dfcg).
 */
#define DOT3_BLOCK_SIZE 256
#define DOT3_MAX_BLOCKS 128

static __device__ double ddot3ReductionResult[3*DOT3_MAX_BLOCKS];

__global__ void spgpuDdot3_kern(int n, const double* __restrict__ x,
				const double* __restrict__ y1,
				const double* __restrict__ y2,
				const double* __restrict__ y3)
{
	__shared__ double s1[DOT3_BLOCK_SIZE];
	__shared__ double s2[DOT3_BLOCK_SIZE];
	__shared__ double s3[DOT3_BLOCK_SIZE];

	double a1 = 0.0, a2 = 0.0, a3 = 0.0;
	int tid = threadIdx.x;
	int stride = gridDim.x*DOT3_BLOCK_SIZE;

	if (y3 != NULL) {
		for (int i = tid + blockIdx.x*DOT3_BLOCK_SIZE; i < n; i += stride) {
			double xi = x[i];
			a1 = PREC_DADD(a1, PREC_DMUL(xi, y1[i]));
			a2 = PREC_DADD(a2, PREC_DMUL(xi, y2[i]));
			a3 = PREC_DADD(a3, PREC_DMUL(xi, y3[i]));
		}
	} else {
		for (int i = tid + blockIdx.x*DOT3_BLOCK_SIZE; i < n; i += stride) {
			double xi = x[i];
			a1 = PREC_DADD(a1, PREC_DMUL(xi, y1[i]));
			a2 = PREC_DADD(a2, PREC_DMUL(xi, y2[i]));
		}
	}
	s1[tid] = a1;
	s2[tid] = a2;
	s3[tid] = a3;
	__syncthreads();

	for (int s = DOT3_BLOCK_SIZE/2; s > 0; s >>= 1) {
		if (tid < s) {
			s1[tid] += s1[tid + s];
			s2[tid] += s2[tid + s];
			s3[tid] += s3[tid + s];
		}
		__syncthreads();
	}

	if (tid == 0) {
		ddot3ReductionResult[3*blockIdx.x]   = s1[0];
		ddot3ReductionResult[3*blockIdx.x+1] = s2[0];
		ddot3ReductionResult[3*blockIdx.x+2] = s3[0];
	}
}

void spgpuDdot3(spgpuHandle_t handle, double* res, int n,
		__device double* x, __device double* y1,
		__device double* y2, __device double* y3)
{
	res[0] = 0.0;
	res[1] = 0.0;
	res[2] = 0.0;
	if (n <= 0) return;

	int blocks = min(DOT3_MAX_BLOCKS,
			 min(2*handle->multiProcessorCount,
			     (n+DOT3_BLOCK_SIZE-1)/DOT3_BLOCK_SIZE));
	double tRes[3*DOT3_MAX_BLOCKS];

	spgpuDdot3_kern<<<blocks, DOT3_BLOCK_SIZE, 0, handle->currentStream>>>(n, x, y1, y2, y3);
	cudaMemcpyFromSymbol(tRes, ddot3ReductionResult, 3*blocks*sizeof(double));

	for (int i = 0; i < blocks; ++i) {
		res[0] += tRes[3*i];
		res[1] += tRes[3*i+1];
		res[2] += tRes[3*i+2];
	}

	cudaCheckError("CUDA error on ddot3");
}
