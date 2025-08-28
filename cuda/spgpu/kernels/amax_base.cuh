/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014
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
 
 
#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)

#undef GEN_SPGPU_AMAX_NAME
#define GEN_SPGPU_AMAX_NAME(x) CONCAT(CONCAT(spgpu,x),amax)

// Define:
//#define VALUE_TYPE
//#define TYPE_SYMBOL

#define BLOCK_SIZE 512

typedef float absType_float;
typedef float absType_cuFloatComplex;
typedef double absType_double;
typedef double absType_cuDoubleComplex;

__device__ __host__ static float zero_float() { return 0.0f; }
__device__ __host__ static cuFloatComplex zero_cuFloatComplex() { return make_cuFloatComplex(0.0, 0.0); }
__device__ __host__ static bool float_isNotZero(float x) { return x != 0.0f; }

__device__ __host__ static float abs_float(float a) { return fabsf(a); }
__device__ __host__ static float abs_cuFloatComplex(cuFloatComplex a) { return cuCabsf(a); }


__device__ static float float_fma(float a, float b, float c) { return PREC_FADD(PREC_FMUL (a, b), c); }
__device__ static float float_add(float a, float b) { return PREC_FADD (a, b); }
__device__ static float float_mul(float a, float b) { return PREC_FMUL (a, b); }

__device__ static cuFloatComplex cuFloatComplex_fma(cuFloatComplex a, cuFloatComplex b, cuFloatComplex c) { return cuCfmaf(a, b, c); } 
__device__ static cuFloatComplex cuFloatComplex_add(cuFloatComplex a, cuFloatComplex b) { return cuCaddf(a, b); }
__device__ static cuFloatComplex cuFloatComplex_mul(cuFloatComplex a, cuFloatComplex b) { return cuCmulf(a, b); }

__device__ static float readValue_float(float fetch) { return fetch; }
__device__ static cuFloatComplex readValue_cuFloatComplex(cuFloatComplex fetch) { return fetch; }

// host or c.c >= 1.3 
#if (__CUDA_ARCH__ >= 130) || (!__CUDA_ARCH__)
__device__ __host__ static double zero_double() { return 0.0; }
__device__ __host__ static cuDoubleComplex zero_cuDoubleComplex() { return make_cuDoubleComplex(0.0, 0.0); }
__device__ __host__ static bool double_isNotZero(double x) { return x != 0.0; }

__device__ __host__ static double abs_double(float a) { return fabs(a); }
__device__ __host__ static double abs_cuDoubleComplex(cuDoubleComplex a) { return cuCabs(a); }

__device__ static double double_fma(double a, double b, double c) { return PREC_DADD(PREC_DMUL (a, b), c); }
__device__ static double double_add(double a, double b) { return PREC_DADD (a, b); }
__device__ static double double_mul(double a, double b) { return PREC_DMUL (a, b); }

__device__ static cuDoubleComplex cuDoubleComplex_fma(cuDoubleComplex a, cuDoubleComplex b, cuDoubleComplex c) { return cuCfma(a, b, c); }
__device__ static cuDoubleComplex cuDoubleComplex_add(cuDoubleComplex a, cuDoubleComplex b) { return cuCadd(a, b); }
__device__ static cuDoubleComplex cuDoubleComplex_mul(cuDoubleComplex a, cuDoubleComplex b) { return cuCmul(a, b); }

__device__ static double readValue_double(int2 fetch) { return __hiloint2double (fetch.y, fetch.x); }
__device__ static cuDoubleComplex readValue_cuDoubleComplex(int4 fetch) 
{
	cuDoubleComplex c;
	c.x = __hiloint2double (fetch.y, fetch.x);
	c.y = __hiloint2double (fetch.w, fetch.z);
	return c;
}
#endif

static __device__ CONCAT(absType_,VALUE_TYPE) CONCAT(TYPE_SYMBOL,amaxReductionResult)[128];

#define MAX(a,b) ((a) > (b) ? (a) : (b))

static __device__ CONCAT(absType_,VALUE_TYPE) amaxvv(VALUE_TYPE a, VALUE_TYPE b)
{
	CONCAT(absType_,VALUE_TYPE) absa = CONCAT(abs_,VALUE_TYPE)(a);
	CONCAT(absType_,VALUE_TYPE) absb = CONCAT(abs_,VALUE_TYPE)(b);
	
	return MAX(absa,absb);
}

static __device__ CONCAT(absType_,VALUE_TYPE) amaxaa(CONCAT(absType_,VALUE_TYPE) a, CONCAT(absType_,VALUE_TYPE) b)
{
	return MAX(a,b);
}

static __device__ CONCAT(absType_,VALUE_TYPE) amaxav(CONCAT(absType_,VALUE_TYPE) a, VALUE_TYPE b)
{
	CONCAT(absType_,VALUE_TYPE) absb = CONCAT(abs_,VALUE_TYPE)(b);
	
	return MAX(a,absb);
}

__global__ void 
CONCAT(GEN_SPGPU_AMAX_NAME(TYPE_SYMBOL),_kern)
(int n, VALUE_TYPE* x)
{
	__shared__ CONCAT(absType_,VALUE_TYPE) sSum[BLOCK_SIZE];

	CONCAT(absType_,VALUE_TYPE) res = 0;

	VALUE_TYPE* lastX = x + n;

	x += threadIdx.x + blockIdx.x*BLOCK_SIZE;

	int blockOffset = gridDim.x*BLOCK_SIZE;

	while (x < lastX)
	{
		VALUE_TYPE x1 = x[0];
		res = amaxav(res, x1);
		
		x += blockOffset;

	}

	if (threadIdx.x >= 32)
		sSum[threadIdx.x] = res;

	__syncthreads();


	// Start reduction!

	if (threadIdx.x < 32) 
	{
		for (int i=1; i<BLOCK_SIZE/32; ++i)
		{
			res = amaxaa(res, sSum[i*32 + threadIdx.x]);
		}

	//useless (because inter-warp)
#ifndef	ASSUME_LOCK_SYNC_PARALLELISM
	}
	__syncthreads(); 

	if (threadIdx.x < 32) 
	{
#endif	

#ifdef ASSUME_LOCK_SYNC_PARALLELISM
		volatile CONCAT(absType_,VALUE_TYPE)* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) amaxaa(vsSum[threadIdx.x], vsSum[threadIdx.x + 16]);
		if (threadIdx.x < 8) amaxaa(vsSum[threadIdx.x], vsSum[threadIdx.x + 8]);
		if (threadIdx.x < 4) amaxaa(vsSum[threadIdx.x], vsSum[threadIdx.x + 4]);
		if (threadIdx.x < 2) amaxaa(vsSum[threadIdx.x], vsSum[threadIdx.x + 2]);
		if (threadIdx.x == 0)
			CONCAT(TYPE_SYMBOL,amaxReductionResult)[blockIdx.x] = amaxaa(vsSum[0], vsSum[1]);

#else
		CONCAT(absType_,VALUE_TYPE)* vsSum = sSum;
		vsSum[threadIdx.x] = res;

		if (threadIdx.x < 16) amaxaa(vsSum[threadIdx.x], vsSum[threadIdx.x + 16]);
		__syncthreads();
		if (threadIdx.x < 8) amaxaa(vsSum[threadIdx.x], vsSum[threadIdx.x + 8]);
		__syncthreads();
		if (threadIdx.x < 4) amaxaa(vsSum[threadIdx.x], vsSum[threadIdx.x + 4]);
		__syncthreads();
		if (threadIdx.x < 2) amaxaa(vsSum[threadIdx.x], vsSum[threadIdx.x + 2]);
		__syncthreads();
		if (threadIdx.x == 0)
			CONCAT(TYPE_SYMBOL,amaxReductionResult)[blockIdx.x] = amaxaa(vsSum[0], vsSum[1]);
#endif	
	}
}

CONCAT(absType_,VALUE_TYPE)
GEN_SPGPU_AMAX_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle, int n, VALUE_TYPE* x)
{
#ifdef USE_CUBLAS
	CONCAT(absType_,VALUE_TYPE) res;
	int id = CONCAT(cublasI,CONCAT(TYPE_SYMBOL,amax))(n,x,1);
	
	cudaError_t err = cudaMemcpy(&res, &(x[id-1]), sizeof(CONCAT(absType_,VALUE_TYPE)), cudaMemcpyDeviceToHost);
	res = CONCAT(abs_,VALUE_TYPE)(res);
	
	return res;

#else
	CONCAT(absType_,VALUE_TYPE) res = 0;

#if 0 	
	int device;
	cudaGetDevice(&device);
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop,device);	

	int blocks = min(128, min(prop.multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#else
	int blocks = min(128, min(handle->multiProcessorCount, (n+BLOCK_SIZE-1)/BLOCK_SIZE));
#endif
	
	CONCAT(absType_,VALUE_TYPE) tRes[128];

	CONCAT(GEN_SPGPU_AMAX_NAME(TYPE_SYMBOL),_kern)<<<blocks, BLOCK_SIZE, 0, handle->currentStream>>>(n, x);;
	cudaMemcpyFromSymbol(tRes, CONCAT(TYPE_SYMBOL,amaxReductionResult), blocks*sizeof(CONCAT(absType_,VALUE_TYPE)));

	for (int i=0; i<blocks; ++i)
	{
		res = MAX(res, tRes[i]);
	}

	cudaCheckError("CUDA error on amax");
	
	return res;
#endif
}

void 
GEN_SPGPU_AMAX_NAME(CONCAT(TYPE_SYMBOL,m))
(spgpuHandle_t handle, CONCAT(absType_,VALUE_TYPE) *y, int n, __device VALUE_TYPE *x, int count, int pitch)
{
	for (int i=0; i < count; ++i)
	{
		y[i] = GEN_SPGPU_AMAX_NAME(TYPE_SYMBOL)(handle, n, x);
		x += pitch;
	}
}
