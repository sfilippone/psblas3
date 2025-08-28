/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2015
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

#undef GEN_SPGPU_HDIA_NAME
#undef X_TEX
#define X_TEX CONCAT(x_tex_, FUNC_SUFFIX)

__device__ __host__ static float zero_float() { return 0.0f; }
__device__ __host__ static cuFloatComplex zero_cuFloatComplex() { return make_cuFloatComplex(0.0, 0.0); }
__device__ __host__ static bool float_isNotZero(float x) { return x != 0.0f; }

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
#if 0
// Texture cache management
texture < TEX_FETCH_TYPE, 1, cudaReadModeElementType > X_TEX;

#define bind_tex_x(x) cudaBindTexture(NULL, X_TEX, x)
#define unbind_tex_x(x) cudaUnbindTexture(X_TEX)

__device__ static VALUE_TYPE 
fetchTex (int pointer)
{
	TEX_FETCH_TYPE fetch = tex1Dfetch (X_TEX, pointer);
	return CONCAT(readValue_,VALUE_TYPE) (fetch);
}
#endif
#define GEN_SPGPU_HDIA_NAME(x) CONCAT(CONCAT(spgpu,x),hdiaspmv_vanilla)
#define GEN_SPGPU_HDIA_NAME_VANILLA(x) CONCAT(CONCAT(spgpu,x),hdiaspmv_vanilla)
#include "hdia_spmv_base_template.cuh"
#if 0 
#undef GEN_SPGPU_HDIA_NAME
#define GEN_SPGPU_HDIA_NAME(x) CONCAT(CONCAT(spgpu,x),hdiaspmv_prefetch)
#define GEN_SPGPU_HDIA_NAME_PREFETCH(x) CONCAT(CONCAT(spgpu,x),hdiaspmv_prefetch)
#undef USE_PREFETCHING
#define USE_PREFETCHING
#include "hdia_spmv_base_template.cuh"
#define ENABLE_CACHE
#undef GEN_SPGPU_HDIA_NAME
#define GEN_SPGPU_HDIA_NAME(x) CONCAT(CONCAT(spgpu,x),hdiaspmv_texcache_prefetch)
#define GEN_SPGPU_HDIA_NAME_TEX_PREFETCH(x) CONCAT(CONCAT(spgpu,x),hdiaspmv_texcache_prefetch)
#include "hdia_spmv_base_template.cuh"
#undef GEN_SPGPU_HDIA_NAME
#undef USE_PREFETCHING
#define GEN_SPGPU_HDIA_NAME(x) CONCAT(CONCAT(spgpu,x),hdiaspmv_texcache)
#define GEN_SPGPU_HDIA_NAME_TEX(x) CONCAT(CONCAT(spgpu,x),hdiaspmv_texcache)
#include "hdia_spmv_base_template.cuh"
#endif
#undef GEN_SPGPU_HDIA_NAME
#define GEN_SPGPU_HDIA_NAME(x) CONCAT(CONCAT(spgpu,x),hdiaspmv)
void
GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle, 
	VALUE_TYPE* z, 
	const VALUE_TYPE *y, 
	VALUE_TYPE alpha, 
	const VALUE_TYPE* dM, 
	const int* offsets, 
	int hackSize, 
	const int* hackOffsets,
	int rows,
	int cols, 
	const VALUE_TYPE *x, 
	VALUE_TYPE beta)
{
	int maxNForACall = max(handle->maxGridSizeX, THREAD_BLOCK*handle->maxGridSizeX);

	// maxNForACall should be a multiple of hackSize
	maxNForACall = (maxNForACall/hackSize)*hackSize;

	while (rows > maxNForACall) //managing large vectors
	{

	CONCAT(_,GEN_SPGPU_HDIA_NAME_VANILLA(TYPE_SYMBOL)) (handle, z, y, alpha, dM, offsets, hackSize, hackOffsets, maxNForACall, cols, x, beta);
	//if (avgDiags < 10 && handle->capabilityMajor > 1)
		//	CONCAT(_,GEN_SPGPU_HDIA_NAME_VANILLA(TYPE_SYMBOL)) (handle, z, y, alpha, dM, offsets, hackSize, hackOffsets, maxNForACall, cols, x, beta);
		//else 
		//if (avgDiags < 20)
		//	CONCAT(_,GEN_SPGPU_HDIA_NAME_TEX(TYPE_SYMBOL)) (handle, z, y, alpha, dM, offsets, hackSize, hackOffsets, maxNForACall, cols, x, beta);
		//else
	//CONCAT(_,GEN_SPGPU_HDIA_NAME_TEX_PREFETCH(TYPE_SYMBOL)) (handle, z, y, alpha, dM, offsets, hackSize, hackOffsets, maxNForACall, cols, x, beta);

		y = y + maxNForACall;
		z = z + maxNForACall;
		hackOffsets += maxNForACall/hackSize;
		
		rows -= maxNForACall;
	}
	CONCAT(_,GEN_SPGPU_HDIA_NAME_VANILLA(TYPE_SYMBOL)) (handle, z, y, alpha, dM, offsets, hackSize, hackOffsets, rows, cols, x, beta);
	//if (avgDiags < 10 && handle->capabilityMajor > 1)
	//	CONCAT(_,GEN_SPGPU_HDIA_NAME_VANILLA(TYPE_SYMBOL)) (handle, z, y, alpha, dM, offsets, dMPitch, rows, cols, diags, x, beta);
	//else 
	//if (avgDiags < 20)
	//	CONCAT(_,GEN_SPGPU_HDIA_NAME_TEX(TYPE_SYMBOL)) (handle, z, y, alpha, dM, offsets, hackSize, hackOffsets, rows, cols, x, beta);
	//else
	//CONCAT(_,GEN_SPGPU_HDIA_NAME_TEX_PREFETCH(TYPE_SYMBOL)) (handle, z, y, alpha, dM, offsets, hackSize, hackOffsets, rows, cols, x, beta);

	cudaCheckError("CUDA error on hdia_spmv");
}

