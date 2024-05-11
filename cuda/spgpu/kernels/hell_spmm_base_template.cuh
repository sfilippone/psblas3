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
#define MMBSZ 8

#if 0

__global__ void
CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn)
  (int count, VALUE_TYPE *z, int zPitch, const VALUE_TYPE *y, int yPitch,
   VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP,
   int hackSize, const int* hackOffsets, const int* rS, int rows,
   const VALUE_TYPE *x, int xPitch, 
   VALUE_TYPE beta, int baseIndex)
{
  VALUE_TYPE *pz,*px,*py;
  VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();
  VALUE_TYPE yVal; 
  __shared__ VALUE_TYPE temp[MMBSZ][THREAD_BLOCK];

  int *rrP;
  VALUE_TYPE *rcM;
  
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int gridSize =  gridDim.x * blockDim.x;
  
  while (i < rows) {
    int j; 
    int hackId = i / hackSize;
    int hackLaneId = i % hackSize;
    
    int hackOffset;
    unsigned int laneId = threadIdx.x % 32;
    if (laneId == 0)
      hackOffset = hackOffsets[hackId];
    //__syncthreads();
    hackOffset = __shfl_sync(0xFFFFFFFF,hackOffset, 0) + hackLaneId;		
    
    rrP = (int *) rP + hackOffset; 
    rcM = (VALUE_TYPE *) cM + hackOffset; 
    
    int rowSize = rS[i];
    for (int k=0; k<count; k++) {
      temp[k][threadIdx.x] = CONCAT(zero_,VALUE_TYPE)();
    }      
    
    for (int j = 0; j < rowSize; j++) {
      int pointer;
      VALUE_TYPE value;
      VALUE_TYPE fetch;
      
      pointer = rrP[0] - baseIndex;
      rrP += hackSize;
      
      value = rcM[0];
      rcM += hackSize;

      px = (VALUE_TYPE *) x;
      for (int k=0; k<count; k++) {
        fetch = px[pointer]; 
        temp[k][threadIdx.x] = CONCAT(VALUE_TYPE, _fma)(value, fetch, temp[k][threadIdx.x]);
	      px = px + xPitch;
      }
    }
    // Since z and y are accessed with the same offset by the same thread,
    // and the write to z follows the y read, y and z can share the same base address (in-place computing).
    py = (VALUE_TYPE *) y;
    pz = z;
    if (CONCAT(VALUE_TYPE, _isNotZero(beta))) {
      for (int k=0; k<count; k++) {
        yVal = py[i];
        pz[i] = CONCAT(VALUE_TYPE, _fma)(beta, yVal, CONCAT(VALUE_TYPE, _mul)(alpha, temp[k][threadIdx.x]));
        py += yPitch;
        pz += zPitch;
      }
    } else {
      for (int k=0; k<count; k++) {
        pz[i] = CONCAT(VALUE_TYPE, _mul)(alpha, temp[k][threadIdx.x]);
        pz += zPitch;
      }
    }
    
    i += gridSize;
  }
}

#endif

__global__ void
CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn)
  (int count, VALUE_TYPE *z, int zPitch, const VALUE_TYPE *y, int yPitch,
   VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP,
   int hackSize, const int* hackOffsets, const int* rS, int rows,
   const VALUE_TYPE *x, int xPitch, 
   VALUE_TYPE beta, int baseIndex)
{
  VALUE_TYPE *pz,*px,*py;
  VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();
  VALUE_TYPE yVal; 
  __shared__ VALUE_TYPE temp[MMBSZ][THREAD_BLOCK];
        
  int i = threadIdx.x + blockIdx.x * (THREAD_BLOCK);
    
  if (i < rows) {
    int j; 
    int hackId = i / hackSize;
    int hackLaneId = i % hackSize;
    
    int hackOffset;
    unsigned int laneId = threadIdx.x % 32;
    if (laneId == 0)
      hackOffset = hackOffsets[hackId];
    //__syncthreads();
    hackOffset = __shfl_sync(0xFFFFFFFF,hackOffset, 0) + hackLaneId;		
    
    rP += hackOffset; 
    cM += hackOffset; 
    
    int rowSize = rS[i];
    for (int k=0; k<count; k++) {
      temp[k][threadIdx.x] = CONCAT(zero_,VALUE_TYPE)();
    }      
    
    for (int j = 0; j < rowSize; j++) {
      int pointer;
      VALUE_TYPE value;
      VALUE_TYPE fetch;
      
      pointer = rP[0] - baseIndex;
      rP += hackSize;
      
      value = cM[0];
      cM += hackSize;

      px = (VALUE_TYPE *) x;
      for (int k=0; k<count; k++) {
        fetch = px[pointer]; 
        temp[k][threadIdx.x] = CONCAT(VALUE_TYPE, _fma)(value, fetch, temp[k][threadIdx.x]);
        px = px + xPitch;
      }
    }
    // Since z and y are accessed with the same offset by the same thread,
    // and the write to z follows the y read, y and z can share the same base address (in-place computing).
    py = (VALUE_TYPE *) y;
    pz = z;
    if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
      for (int k=0; k<count; k++) {
        yVal = py[i];
        pz[i] = CONCAT(VALUE_TYPE, _fma)(beta, yVal, CONCAT(VALUE_TYPE, _mul)(alpha, temp[k][threadIdx.x]));
        py += yPitch;
        pz += zPitch;
      }
    else
      for (int k=0; k<count; k++) {
        pz[i] = CONCAT(VALUE_TYPE, _mul)(alpha, temp[k][threadIdx.x]);
        pz += zPitch;
      }
  }
}


void
CONCAT(_,GEN_SPGPU_HELL_NAME(TYPE_SYMBOL))
  (spgpuHandle_t handle, int count, VALUE_TYPE* z, int zPitch, const VALUE_TYPE *y, int yPitch,
   VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int hackSize, const int* hackOffsets,
   const int* rS,  const __device int* rIdx, int rows,
   const VALUE_TYPE *x, int xPitch, VALUE_TYPE beta, int baseIndex)
{
  dim3 block (THREAD_BLOCK, 1);
  dim3 grid ((rows + THREAD_BLOCK - 1) / THREAD_BLOCK);
  // Should we generalize the code to 1/2/4/8 threads per row?
  // And maybe adjust THREAD_BLOCK size? 
  int shrMemSize,maxShmemSz;
  maxShmemSz=getGPUSharedMemPerBlock();
  shrMemSize=MMBSZ*THREAD_BLOCK*sizeof(VALUE_TYPE);
  CONCAT(GEN_SPGPU_HELL_NAME(TYPE_SYMBOL), _krn) 
    <<< grid, block, shrMemSize, handle->currentStream >>> (count, z, zPitch,y, yPitch,
							    alpha, cM, rP, hackSize, hackOffsets, rS, rows,
							    x, xPitch, beta, baseIndex);
}
