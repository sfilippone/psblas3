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

#undef GEN_SPGPU_HELL_NAME
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
#if __CUDA_ARCH__ < 300
extern __shared__ int dynShrMem[]; 
#endif

#define GEN_SPGPU_HELL_NAME(x) CONCAT(CONCAT(spgpu,x),hellspmm_vanilla)
#define GEN_SPGPU_HELL_NAME_VANILLA(x) CONCAT(CONCAT(spgpu,x),hellspmm_vanilla)
#include "hell_spmm_base_template.cuh"
#undef GEN_SPGPU_HELL_NAME
#if 0
#define GEN_SPGPU_HELL_NAME(x) CONCAT(CONCAT(spgpu,x),hellspmm_prefetch)
#define GEN_SPGPU_HELL_NAME_PREFETCH(x) CONCAT(CONCAT(spgpu,x),hellspmm_prefetch)
#undef USE_PREFETCHING
#define USE_PREFETCHING
#include "hell_spmm_base_template.cuh"
#define ENABLE_CACHE
#undef GEN_SPGPU_HELL_NAME
#define GEN_SPGPU_HELL_NAME(x) CONCAT(CONCAT(spgpu,x),hellspmm_texcache_prefetch)
#define GEN_SPGPU_HELL_NAME_TEX_PREFETCH(x) CONCAT(CONCAT(spgpu,x),hellspmm_texcache_prefetch)
#include "hell_spmm_base_template.cuh"
#undef GEN_SPGPU_HELL_NAME
#undef USE_PREFETCHING
#endif
#define GEN_SPGPU_HELL_NAME(x) CONCAT(CONCAT(spgpu,x),hellspmm_texcache)
#define GEN_SPGPU_HELL_NAME_TEX(x) CONCAT(CONCAT(spgpu,x),hellspmm_texcache)
#include "hell_spmm_base_template.cuh"
#undef GEN_SPGPU_HELL_NAME
#define GEN_SPGPU_HELL_NAME(x) CONCAT(CONCAT(spgpu,x),hellspmm)

#if 0

void
GEN_SPGPU_HELL_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle,
    int count, 
    VALUE_TYPE* z,
    int zPitch,
	const VALUE_TYPE *y,
    int  yPitch,
	VALUE_TYPE alpha, 
	const VALUE_TYPE* cM, 
	const int* rP, 
	int hackSize,
	const __device int* hackOffsets, 
	const __device int* rS,
	const __device int* rIdx, 
	int rows, 
	const VALUE_TYPE *x,
    int xPitch,
	VALUE_TYPE beta, 
    int baseIndex)
{
  VALUE_TYPE *px,*py, *pz;
  int cnt, c1;

  dim3 block (THREAD_BLOCK, 1);
  //  dim3 grid ((rows + THREAD_BLOCK - 1) / THREAD_BLOCK);
  // Should we generalize the code to 1/2/4/8 threads per row?
  // And maybe adjust THREAD_BLOCK size? 
  int shrMemSize,maxShmemSz;
  int numMp=getGPUMultiProcessors();
  int maxThMp=getGPUMaxThreadsPerMP();
  int nmblksMp=maxThMp/THREAD_BLOCK;
  int nmblk=nmblksMp*numMp;
  dim3 grid (nmblk);

  maxShmemSz=getGPUSharedMemPerBlock();
  shrMemSize=MMBSZ*THREAD_BLOCK*sizeof(VALUE_TYPE);
  if (shrMemSize > maxShmemSz) {
    fprintf(stderr,"Fatal error: SHMEM size too large %ld %ld\n",shrMemSize,maxShmemSz);
    return;
  }
  cnt = count;
  px = (VALUE_TYPE *) x;
  py = (VALUE_TYPE *) y;
  pz = (VALUE_TYPE *) z;	  
  while (cnt > 2*MMBSZ) {
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn) 
      <<< grid, block, shrMemSize, handle->currentStream >>> (MMBSZ, pz, zPitch,py, yPitch,
							      alpha, cM, rP, hackSize, hackOffsets,
							      rS, rows, px, xPitch, beta, baseIndex);
    px += xPitch*MMBSZ;
    py += yPitch*MMBSZ;
    pz += zPitch*MMBSZ;
    cnt -= MMBSZ;
  }
  if (cnt > MMBSZ) {
    c1 = cnt/2;
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn) 
      <<< grid, block, shrMemSize, handle->currentStream >>> (c1, pz, zPitch,py, yPitch,
							      alpha, cM, rP, hackSize, hackOffsets,
							      rS, rows, px, xPitch, beta, baseIndex);
    cnt -= c1;
  }
  if (cnt > MMBSZ) {
    fprintf(stderr,"Invalid residual count %d\n",cnt);
  } else if (cnt > 0){
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn) 
      <<< grid, block, shrMemSize, handle->currentStream >>> (cnt, pz, zPitch,py, yPitch,
							      alpha, cM, rP, hackSize, hackOffsets,
							      rS, rows, px, xPitch, beta, baseIndex);
  }
  cudaCheckError("CUDA error on hell_spmm");
}

#endif

void
GEN_SPGPU_HELL_NAME(TYPE_SYMBOL)
(spgpuHandle_t handle,
    int count, 
    VALUE_TYPE* z,
    int zPitch,
	const VALUE_TYPE *y,
    int  yPitch,
	VALUE_TYPE alpha, 
	const VALUE_TYPE* cM, 
	const int* rP, 
	int hackSize,
	const __device int* hackOffsets, 
	const __device int* rS,
	const __device int* rIdx, 
	int rows, 
	const VALUE_TYPE *x,
    int xPitch,
	VALUE_TYPE beta, 
    int baseIndex)
{
  VALUE_TYPE *px,*py, *pz;
  int cnt;
  int maxNForACall = max(handle->maxGridSizeX, THREAD_BLOCK*handle->maxGridSizeX);
  
  // maxNForACall should be a multiple of hackSize
  maxNForACall = (maxNForACall/hackSize)*hackSize;
  int maxShmemSz;
  maxShmemSz=getGPUSharedMemPerBlock();
  //fprintf(stderr,"MaxSHmemSz  %d \n",maxShmemSz);
  while (rows > maxNForACall) {//managing large vectors
    cnt = count;
    px = (VALUE_TYPE *) x;
    py = (VALUE_TYPE *) y;
    pz = (VALUE_TYPE *) z;	  
    while (cnt > MMBSZ) {
      //fprintf(stderr,"counts %d %d %d :  pointers: %p %p %p\n",rows,cnt,MMBSZ,px,py,pz);    
      CONCAT(_,GEN_SPGPU_HELL_NAME_VANILLA(TYPE_SYMBOL)) (handle, MMBSZ, pz, zPitch,
							  py, yPitch,
							  alpha, cM, rP,
							  hackSize, hackOffsets,
							  rS, rIdx,
							  maxNForACall,
							  px, xPitch, beta, baseIndex);
      px += xPitch*MMBSZ;
      py += yPitch*MMBSZ;
      pz += zPitch*MMBSZ;
      cnt -= MMBSZ;
    }
    if (cnt >0) {
      CONCAT(_,GEN_SPGPU_HELL_NAME_VANILLA(TYPE_SYMBOL)) (handle, cnt, pz, zPitch,
							  py, yPitch,
							  alpha, cM, rP,
							  hackSize, hackOffsets,
							  rS, rIdx,
							  maxNForACall,
							  px, xPitch, beta, baseIndex);
    }

    y = y + maxNForACall;
    z = z + maxNForACall;
    hackOffsets = hackOffsets + maxNForACall/hackSize;
    rS = rS + maxNForACall;
    
    rows -= maxNForACall;
  }
  cnt = count;
  px = (VALUE_TYPE *) x;
  py = (VALUE_TYPE *) y;
  pz = (VALUE_TYPE *) z;	  
  while (cnt > MMBSZ) {
    //fprintf(stderr,"counts %d %d %d :  pointers: %p %p %p\n",rows,cnt,MMBSZ,px,py,pz);    
    CONCAT(_,GEN_SPGPU_HELL_NAME_VANILLA(TYPE_SYMBOL)) (handle, MMBSZ, pz, zPitch, py, yPitch,
							alpha, cM, rP, hackSize, hackOffsets,
							rS, rIdx, rows,
							px, xPitch, beta, baseIndex);
    px += xPitch*MMBSZ;
    py += yPitch*MMBSZ;
    pz += zPitch*MMBSZ;
    cnt -= MMBSZ;
  }
  if (cnt >0) {
    CONCAT(_,GEN_SPGPU_HELL_NAME_VANILLA(TYPE_SYMBOL)) (handle, cnt, pz, zPitch,
							py, yPitch,
							alpha, cM, rP,
							hackSize, hackOffsets,
							rS, rIdx,
							rows,
							px, xPitch, beta, baseIndex);
  }
  
  
  cudaCheckError("CUDA error on hell_spmm");
}
