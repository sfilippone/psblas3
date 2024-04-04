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
#define IDX2
#define THREAD_BLOCK 128

__device__ void
CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _ridx_2)
(int i, VALUE_TYPE yVal, int outRow,
 VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int hackSize, const int* hackOffsets, const int* rS, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
  VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();
  
  __shared__ VALUE_TYPE temp[THREAD_BLOCK];
  
  if (i < rows) {
    int hackId = i / hackSize;
    int hackLaneId = i % hackSize;
    
    int hackOffset;
    unsigned int laneId = threadIdx.x % 32;
#if __CUDA_ARCH__ < 300
    // "volatile" used to avoid __syncthreads()
    volatile int* warpHackOffset = dynShrMem;
    
    unsigned int warpId = threadIdx.x / 32;
    
    if (laneId == 0)
      warpHackOffset[warpId] = hackOffsets[hackId];
    
    hackOffset = warpHackOffset[warpId] + hackLaneId;
#elif __CUDA_ARCH__ < 700
    if (laneId == 0)
      hackOffset = hackOffsets[hackId];
    //__syncthreads();
    hackOffset = __shfl(hackOffset, 0) + hackLaneId;
#else
    if (laneId == 0)
      hackOffset = hackOffsets[hackId];
    //__syncthreads();
    hackOffset = __shfl_sync(0xFFFFFFFF,hackOffset, 0) + hackLaneId;		
#endif
    
    rP += hackOffset; 
    cM += hackOffset; 
    
    int rowSize = rS[i]; 
    int rowSizeM = rowSize / 2;
    
    if (threadIdx.y == 0) {
      if (rowSize % 2)
	++rowSizeM;
    } else {
      rP += hackSize; 
      cM += hackSize;
    }
    
    
    for (int j = 0; j < rowSizeM; j++) {
      int pointer;
      VALUE_TYPE value;
      VALUE_TYPE fetch;
      
      pointer = rP[0] - baseIndex;
      rP += hackSize; 
      rP += hackSize;
      
      value = cM[0];
      cM += hackSize;
      cM += hackSize;
      
#ifdef ENABLE_CACHE
      fetch = fetchTex(pointer);
#else
      fetch = x[pointer];
#endif	
      
      // avoid MAD on pre-Fermi
      zProd = CONCAT(VALUE_TYPE, _fma)(value, fetch, zProd);
    }
    
    // Reduction
    if (threadIdx.y == 1)
      temp[threadIdx.x] = zProd;
    
    __syncthreads();
    
    if (threadIdx.y == 0)     {
      zProd = CONCAT(VALUE_TYPE, _add)(zProd, temp[threadIdx.x]);
      // Since z and y are accessed with the same offset by the same thread,
      // and the write to z follows the y read, y and z can share the same base address (in-place computing).
      if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
	z[outRow] = CONCAT(VALUE_TYPE, _fma)(beta, yVal, CONCAT(VALUE_TYPE, _mul) (alpha, zProd));
      else
	z[outRow] = CONCAT(VALUE_TYPE, _mul)(alpha, zProd);
    }
  }
}	

__device__ void
CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _ridx)
  (int i, VALUE_TYPE yVal, int outRow,
   VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP,  int hackSize, const int* hackOffsets, const int* rS, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
  VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();
  
  if (i < rows) {
    
    int hackId = i / hackSize;
    int hackLaneId = i % hackSize;
    
    int hackOffset;
    unsigned int laneId = threadIdx.x % 32;
#if __CUDA_ARCH__ < 300
    // "volatile" used to avoid __syncthreads()
    volatile int* warpHackOffset = dynShrMem;
    
    unsigned int warpId = threadIdx.x / 32;
    
    if (laneId == 0)
      warpHackOffset[warpId] = hackOffsets[hackId];
    
    hackOffset = warpHackOffset[warpId] + hackLaneId;
#elif __CUDA_ARCH__ < 700
    if (laneId == 0)
      hackOffset = hackOffsets[hackId];
    //__syncthreads();    
    hackOffset = __shfl(hackOffset, 0) + hackLaneId;
#else
    if (laneId == 0)
      hackOffset = hackOffsets[hackId];
    //__syncthreads();
    hackOffset = __shfl_sync(0xFFFFFFFF,hackOffset, 0) + hackLaneId;		
#endif
    
    rP += hackOffset; 
    cM += hackOffset; 
    
    int rowSize = rS[i];
    
#ifdef USE_PREFETCHING		
    for (int j = 0; j < rowSize / 2; j++) {
      int pointers1, pointers2;
      VALUE_TYPE values1, values2;
      VALUE_TYPE fetches1, fetches2;
      
      pointers1 = rP[0] - baseIndex;
      rP += hackSize; 
      pointers2 = rP[0] - baseIndex;
      rP += hackSize; 
      
      values1 = cM[0];
      cM += hackSize;
      
      values2 = cM[0];
      cM += hackSize;
      
#ifdef ENABLE_CACHE
      fetches1 = fetchTex(pointers1);
      fetches2 = fetchTex(pointers2);
#else
      fetches1 = x[pointers1];
      fetches2 = x[pointers2];	
#endif
      
      // avoid MAD on pre-Fermi
      zProd = CONCAT(VALUE_TYPE, _fma)(values1, fetches1, zProd);
      zProd = CONCAT(VALUE_TYPE, _fma)(values2, fetches2, zProd);
    }
    
    // odd row size
    if (rowSize % 2) {
      int pointer = rP[0] - baseIndex;
      VALUE_TYPE value = cM[0];
      VALUE_TYPE fetch;
      
#ifdef ENABLE_CACHE
      fetch = fetchTex (pointer);
#else
      fetch = x[pointer];
#endif
      zProd = CONCAT(VALUE_TYPE, _fma)(value, fetch, zProd);
    }
#else
    for (int j = 0; j < rowSize; j++) {
      int pointer;
      VALUE_TYPE value;
      VALUE_TYPE fetch;
      
      pointer = rP[0] - baseIndex;
      rP += hackSize;
      
      value = cM[0];
      cM += hackSize;
      
#ifdef ENABLE_CACHE
      fetch = fetchTex (pointer);
#else
      fetch = x[pointer];
#endif
      zProd = CONCAT(VALUE_TYPE, _fma)(value, fetch, zProd);
    }
#endif	   	
    // Since z and y are accessed with the same offset by the same thread,
    // and the write to z follows the y read, y and z can share the same base address (in-place computing).
    
    if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
      z[outRow] = CONCAT(VALUE_TYPE, _fma)(beta, yVal, CONCAT(VALUE_TYPE, _mul) (alpha, zProd));
    else
      z[outRow] = CONCAT(VALUE_TYPE, _mul)(alpha, zProd);
  }

}

__global__ void
CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn_ridx)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP,  int hackSize, const int* hackOffsets, const int* rS, const int* rIdx, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
  int i = threadIdx.x + blockIdx.x * (THREAD_BLOCK);
  
  VALUE_TYPE yVal = CONCAT(zero_,VALUE_TYPE)();
  int outRow = 0;
  if (i < rows)  {
    
    outRow = rIdx[i];
    if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
      yVal = y[outRow];
  }
#if 1
  if (blockDim.y == 1)
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _ridx)
      (i, yVal, outRow, z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, beta, baseIndex);
  else
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _ridx_2)
      (i, yVal, outRow, z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, beta, baseIndex);
#else
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _ridx)
      (i, yVal, outRow, z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, beta, baseIndex);
#endif  
}


__device__ void
CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _)
  (VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP,  int hackSize, const int* hackOffsets, const int* rS, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
  int i = threadIdx.x + blockIdx.x * (THREAD_BLOCK);
  
  VALUE_TYPE yVal = CONCAT(zero_,VALUE_TYPE)();
  
  if (i < rows)	{
    if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
      yVal = y[i];
    
  }
  
#ifdef IDX2
  if (blockDim.y == 1)
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _ridx)
      (i, yVal, i, z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, beta, baseIndex);
  else
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _ridx_2)
      (i, yVal, i, z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, beta, baseIndex);
#else
      CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _ridx)
      (i, yVal, i, z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, beta, baseIndex);
#endif
  
}

// Force to recompile and optimize with llvm
__global__ void
CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn_b0) 
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP,  int hackSize, const int* hackOffsets, const int* rS, int rows, const VALUE_TYPE *x, int baseIndex)
{
  CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _)
    (z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, CONCAT(zero_,VALUE_TYPE)(), baseIndex);
}

__global__ void
CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn)
  (VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int hackSize, const int* hackOffsets, const int* rS, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
  CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _)
    (z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, beta, baseIndex);
}

void
CONCAT(_,GEN_SPGPU_HELL_NAME(TYPE_SYMBOL))
(spgpuHandle_t handle, VALUE_TYPE* z, const VALUE_TYPE *y, VALUE_TYPE alpha, 
 const VALUE_TYPE* cM, const int* rP, int hackSize, const int* hackOffsets, const int* rS,  
 const __device int* rIdx, int avgNnzPerRow, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
  int avgThreshold;
  
  if (handle->capabilityMajor < 2)
    avgThreshold = 8;
  else if (handle->capabilityMajor < 3)
    avgThreshold = 16;
  else
    avgThreshold = 32;
#ifdef IDX2  
#if defined(HELL_FORCE_THREADS_1)
  dim3 block (THREAD_BLOCK, 1);
#elif defined(HELL_FORCE_THREADS_2)
  dim3 block (THREAD_BLOCK, 2);
#else
  dim3 block (THREAD_BLOCK, avgNnzPerRow >= avgThreshold ? 2 : 1);
#endif
#else
  dim3 block (THREAD_BLOCK, 1);
#endif
  dim3 grid ((rows + THREAD_BLOCK - 1) / THREAD_BLOCK);

  // Should we generalize the code to 1/2/4/8 threads per row?
  // And maybe adjust THREAD_BLOCK size? 
  int shrMemSize;
  shrMemSize=THREAD_BLOCK*sizeof(VALUE_TYPE);
  
#ifdef ENABLE_CACHE
  bind_tex_x ((const TEX_FETCH_TYPE *) x);
#endif
  
  if (rIdx) {
    cudaFuncSetCacheConfig(CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn_ridx), cudaFuncCachePreferL1);
    
    CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn_ridx)
      <<< grid, block, shrMemSize, handle->currentStream >>> (z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rIdx, rows, x, beta, baseIndex);
  } else {
    cudaFuncSetCacheConfig(CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn), cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn_b0), cudaFuncCachePreferL1);
    
    if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
      CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn) 
	<<< grid, block, shrMemSize, handle->currentStream >>> (z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, beta, baseIndex);
    else
      CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn_b0)
	<<< grid, block, shrMemSize, handle->currentStream >>> (z, y, alpha, cM, rP, hackSize, hackOffsets, rS, rows, x, baseIndex);
  }
  
#ifdef ENABLE_CACHE
  unbind_tex_x ((const TEX_FETCH_TYPE *) x);
#endif
  
}
