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
 
#define THREAD_BLOCK 128
#define MMBSZ 8

__global__ void
CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _krn)
  (int count, VALUE_TYPE* z, int zPitch, const VALUE_TYPE *y,
	int yPitch, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* hdiaOffsets,
  int hackSize, const int* hackOffsets, int rows, int cols,
  const VALUE_TYPE *x, int xPitch, VALUE_TYPE beta)
{
  VALUE_TYPE *pz,*px,*py;
  VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();
  VALUE_TYPE yVal; 
  __shared__ VALUE_TYPE temp[MMBSZ][THREAD_BLOCK];
        
  int hackCount = (rows + hackSize - 1)/hackSize;

  int i = threadIdx.x + blockIdx.x * (THREAD_BLOCK);

  int hackId = i / hackSize;
  int hackLaneId = i % hackSize;

  // shared between offsetsChunks and warpHackOffsetTemp
	extern __shared__ int dynShrMem[]; 
    
  int hackOffset = 0;
	int nextOffset = 0;

  unsigned int laneId = threadIdx.x % warpSize;
	unsigned int warpId = threadIdx.x / warpSize;

  if (laneId == 0 && i < rows) {
    hackOffset = hackOffsets[hackId];
		nextOffset = hackOffsets[hackId+1];
  }

  hackOffset = __shfl_sync(0xFFFFFFFF,hackOffset, 0);		
	nextOffset = __shfl_sync(0xFFFFFFFF,nextOffset, 0);
    
  if (hackId >= hackCount)
		return;
    
  cM += hackOffset*hackSize + hackLaneId;
	hdiaOffsets += hackOffset;

  for (int k=0; k<count; k++) {
    temp[k][threadIdx.x] = CONCAT(zero_,VALUE_TYPE)();
  }

  // diags for this hack is next hackOffset minus current hackOffset
	int diags = nextOffset - hackOffset;
    
	// Warp oriented
	int rounds = (diags + warpSize - 1)/warpSize;
	
	volatile int *offsetsChunk = dynShrMem + warpId*warpSize;

  for (int r = 0; r < rounds; r++) {
    // in the last round diags will be <= warpSize
		if (laneId < diags)
			offsetsChunk[laneId] = hdiaOffsets[laneId];

    if (i < rows) {
			int dCount = min(diags, warpSize);

      for (int j = 0; j < dCount; ++j) {
        int column = offsetsChunk[j] + i;

				if(column >= 0 && column < cols) {
          px = (VALUE_TYPE *) x;
          for (int k = 0; k < count; k++) {
            VALUE_TYPE xValue = px[column];
            VALUE_TYPE mValue = cM[0];
            temp[k][threadIdx.x] = CONCAT(VALUE_TYPE, _fma)(mValue, xValue, temp[k][threadIdx.x]);
            px = px + xPitch;
          }
        }
        cM += hackSize;
      }
    }
    diags -= warpSize;
    hdiaOffsets += warpSize;
  }

  // Since z and y are accessed with the same offset by the same thread,
	// and the write to z follows the y read, y and z can share the same base address (in-place computing).
	if (i >= rows)
    return;

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

void
CONCAT(_,GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL))
  (spgpuHandle_t handle, int count, VALUE_TYPE* z, int zPitch, const VALUE_TYPE *y,
	int yPitch, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* hdiaOffsets, int hackSize,
	const __device int* hackOffsets, int rows, int cols, const VALUE_TYPE *x,
  int xPitch, VALUE_TYPE beta)
{
  dim3 block (THREAD_BLOCK, 1);
  dim3 grid ((rows + THREAD_BLOCK - 1) / THREAD_BLOCK);
  // Should we generalize the code to 1/2/4/8 threads per row?
  // And maybe adjust THREAD_BLOCK size? 
  int shrMemSize,maxShmemSz;
  maxShmemSz = getGPUSharedMemPerBlock();
  shrMemSize = MMBSZ*THREAD_BLOCK*sizeof(VALUE_TYPE);
  CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _krn) 
    <<< grid, block, shrMemSize, handle->currentStream >>> (count, z, zPitch, y, yPitch,
							    alpha, cM, hdiaOffsets, hackSize, hackOffsets, rows, cols,
							    x, xPitch, beta);
}
