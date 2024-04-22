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

#undef GEN_SPGPU_ELL_NAME
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
#define GEN_SPGPU_ELL_NAME(x) CONCAT(CONCAT(spgpu,x),ellspmm_vanilla)
#define GEN_SPGPU_ELL_NAME_VANILLA(x) CONCAT(CONCAT(spgpu,x),ellspmm_vanilla)
#include "ell_spmm_base_template.cuh"
#if 0
#undef GEN_SPGPU_ELL_NAME
#define GEN_SPGPU_ELL_NAME(x) CONCAT(CONCAT(spgpu,x),ellspmm_prefetch)
#define GEN_SPGPU_ELL_NAME_PREFETCH(x) CONCAT(CONCAT(spgpu,x),ellspmm_prefetch)
#undef USE_PREFETCHING
#define USE_PREFETCHING
#include "ell_spmm_base_template.cuh"
#define ENABLE_CACHE
#undef GEN_SPGPU_ELL_NAME
#define GEN_SPGPU_ELL_NAME(x) CONCAT(CONCAT(spgpu,x),ellspmm_texcache_prefetch)
#define GEN_SPGPU_ELL_NAME_TEX_PREFETCH(x) CONCAT(CONCAT(spgpu,x),ellspmm_texcache_prefetch)
#include "ell_spmm_base_template.cuh"
#undef GEN_SPGPU_ELL_NAME
#undef USE_PREFETCHING
#define GEN_SPGPU_ELL_NAME(x) CONCAT(CONCAT(spgpu,x),ellspmm_texcache)
#define GEN_SPGPU_ELL_NAME_TEX(x) CONCAT(CONCAT(spgpu,x),ellspmm_texcache)
#include "ell_spmm_base_template.cuh"
#endif
#undef GEN_SPGPU_ELL_NAME
#define GEN_SPGPU_ELL_NAME(x) CONCAT(CONCAT(spgpu,x),ellspmm)

void
GEN_SPGPU_ELL_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle,
	int count,
	VALUE_TYPE* z,
	int zPitch,
	const VALUE_TYPE *y,
	int yPitch,
	VALUE_TYPE alpha, 
	const VALUE_TYPE* cM, 
	const int* rP,
    int cMPitch,
    int rPPitch,
	const __device int* rS,
	const __device int* rIdx, 
    int avgNnzPerRow,
    int maxNnzPerRow,
	int rows, 
	const VALUE_TYPE *x,
	int xPitch,
	VALUE_TYPE beta,
	int baseIndex)
{
// TODO
  VALUE_TYPE *px,*py,*pz;
  int cnt;
  int maxNForACall = max(handle->maxGridSizeX, THREAD_BLOCK*handle->maxGridSizeX);

  int maxShmemSz;
  maxShmemSz=getGPUSharedMemPerBlock();

  while (rows > maxNForACall) {//managing large vectors
    cnt = count;
    px = (VALUE_TYPE *) x;
    py = (VALUE_TYPE *) y;
    pz = (VALUE_TYPE *) z;	  
    while (cnt > MMBSZ) {
      CONCAT(_,GEN_SPGPU_ELL_NAME_VANILLA(TYPE_SYMBOL)) (handle, MMBSZ, pz, zPitch,
							  py, yPitch,
							  alpha, cM, rP,
							  cMPitch, rPPitch,
							  rS, rIdx, avgNnzPerRow,
							  maxNnzPerRow, maxNForACall,
							  px, xPitch, beta, baseIndex);
      px += xPitch*MMBSZ;
      py += yPitch*MMBSZ;
      pz += zPitch*MMBSZ;
      cnt -= MMBSZ;
    }
    if (cnt >0) {
      CONCAT(_,GEN_SPGPU_ELL_NAME_VANILLA(TYPE_SYMBOL)) (handle, cnt, pz, zPitch,
							  py, yPitch,
							  alpha, cM, rP,
							  cMPitch, rPPitch,
							  rS, rIdx, avgNnzPerRow,
							  maxNnzPerRow, maxNForACall,
							  px, xPitch, beta, baseIndex);
    }

    y = y + maxNForACall;
    z = z + maxNForACall;
    cM = cM + maxNForACall;
	rP = rP + maxNForACall;
	rS = rS + maxNForACall;
    rows -= maxNForACall;
  }

  cnt = count;
  px = (VALUE_TYPE *) x;
  py = (VALUE_TYPE *) y;
  pz = (VALUE_TYPE *) z;	  
  while (cnt > MMBSZ) {
    fprintf(stderr,"counts %d %d %d :  pointers: %p %p %p\n",rows,cnt,MMBSZ,px,py,pz);
    CONCAT(_,GEN_SPGPU_ELL_NAME_VANILLA(TYPE_SYMBOL)) (handle, MMBSZ, pz, zPitch,
							  py, yPitch,
							  alpha, cM, rP,
							  cMPitch, rPPitch,
							  rS, rIdx, avgNnzPerRow,
							  maxNnzPerRow, maxNForACall,
							  px, xPitch, beta, baseIndex);
    px += xPitch*MMBSZ;
    py += yPitch*MMBSZ;
    pz += zPitch*MMBSZ;
    cnt -= MMBSZ;
  }
  if (cnt >0) {
    fprintf(stderr,"counts %d %d %d :  pointers: %p %p %p\n",rows,cnt,MMBSZ,px,py,pz);

    CONCAT(_,GEN_SPGPU_ELL_NAME_VANILLA(TYPE_SYMBOL)) (handle, cnt, pz, zPitch,
							  py, yPitch,
							  alpha, cM, rP,
							  cMPitch, rPPitch,
							  rS, rIdx, avgNnzPerRow,
							  maxNnzPerRow, maxNForACall,
							  px, xPitch, beta, baseIndex);
  }
  
  cudaCheckError("CUDA error on hell_spmm");
}