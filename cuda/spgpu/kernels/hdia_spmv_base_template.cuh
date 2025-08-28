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

__device__ void
CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* dM, const int* offsets, int hackSize, const int* hackOffsets, 
	int rows, int cols, const VALUE_TYPE *x, VALUE_TYPE beta, int hackCount)
{
	int i = threadIdx.x + blockIdx.x * (blockDim.x);
	
	VALUE_TYPE yVal = CONCAT(zero_,VALUE_TYPE)();

	if (i < rows && CONCAT(VALUE_TYPE, _isNotZero(beta)))
		yVal = y[i];

	VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();
	
	int hackId = i / hackSize;
	int hackLaneId = i % hackSize;
	
	
	// shared between offsetsChunks and warpHackOffsetTemp
	extern __shared__ int dynShrMem[]; 

	int hackOffset = 0;
	int nextOffset = 0;
	
	unsigned int laneId = threadIdx.x % warpSize;
	unsigned int warpId = threadIdx.x / warpSize;
	
#if __CUDA_ARCH__ < 300	
	int* warpHackOffset = dynShrMem;


	if (laneId == 0 && i < rows)
	{
		warpHackOffset[warpId] = hackOffsets[hackId];
		warpHackOffset[warpId + (blockDim.x / warpSize)] = hackOffsets[hackId+1];
	}
	
	__syncthreads();
	hackOffset = warpHackOffset[warpId];
	nextOffset = warpHackOffset[warpId + blockDim.x / warpSize];
	__syncthreads();
#elif __CUDA_ARCH__ < 700
	if (laneId == 0 && i < rows)
	{
		hackOffset = hackOffsets[hackId];
		nextOffset = hackOffsets[hackId+1];
	}
	
	hackOffset = __shfl(hackOffset, 0);	
	nextOffset = __shfl(nextOffset, 0);
#else
	if (laneId == 0 && i < rows)
	{
		hackOffset = hackOffsets[hackId];
		nextOffset = hackOffsets[hackId+1];
	}
	
	hackOffset = __shfl_sync(0xFFFFFFFF,hackOffset, 0);	
	nextOffset = __shfl_sync(0xFFFFFFFF,nextOffset, 0);
	
#endif
	
	if (hackId >= hackCount)
		return;

	dM += hackOffset*hackSize + hackLaneId;
	offsets += hackOffset;
	
	// diags for this hack is next hackOffset minus current hackOffset
	int diags = nextOffset - hackOffset;
	
	
	// Warp oriented
	int rounds = (diags + warpSize - 1)/warpSize;
	
	volatile int *offsetsChunk = dynShrMem + warpId*warpSize;
	
	for (int r = 0; r < rounds; r++)
	{
		// in the last round diags will be <= warpSize
		if (laneId < diags)
			offsetsChunk[laneId] = offsets[laneId];
	
		if (i < rows)
		{
			int count = min(diags, warpSize);

#ifdef USE_PREFETCHING	
			int j;
			for (j=0; j<=count-2; j += 2)
			{
				// prefetch 2 values
				int column1 = offsetsChunk[j] + i;
				int column2 = offsetsChunk[j+1] + i;			

				VALUE_TYPE xValue1, xValue2;
				VALUE_TYPE mValue1, mValue2;
				
				bool inside1 = column1 >= 0 && column1 < cols;
				bool inside2 = column2 >= 0 && column2 < cols;
				
				if(inside1)
                		{
                			mValue1 = dM[0];
#ifdef ENABLE_CACHE
					xValue1 = fetchTex (column1);
#else
					xValue1 = x[column1];
#endif				
				}
				
				dM += hackSize;
							
				if(inside2)
                		{
                			mValue2 = dM[0];
#ifdef ENABLE_CACHE
					xValue2 = fetchTex (column2);
#else
					xValue2 = x[column2];
#endif					
				}

				dM += hackSize;					
											
				if(inside1)
					zProd = CONCAT(VALUE_TYPE, _fma)(mValue1, xValue1, zProd);
				if(inside2)
					zProd = CONCAT(VALUE_TYPE, _fma)(mValue2, xValue2, zProd);
			}
	
			for (;j<count; ++j)
			{
				int column = offsetsChunk[j] + i;
				
				if(column >= 0 && column < cols)
                		{
					VALUE_TYPE xValue;
#ifdef ENABLE_CACHE
					xValue = fetchTex (column);
#else
					xValue = x[column];
#endif				
                			VALUE_TYPE mValue = dM[0];
					zProd = CONCAT(VALUE_TYPE, _fma)(mValue, xValue, zProd);				
				}
				
				dM += hackSize;
			}
#else
			for (int j=0;j<count; ++j)
			{
				int column = offsetsChunk[j] + i;
				
				if(column >= 0 && column < cols)
				{
					VALUE_TYPE xValue;
#ifdef ENABLE_CACHE
					xValue = fetchTex (column);
#else
					xValue = x[column];
#endif				
					VALUE_TYPE mValue = dM[0];
					zProd = CONCAT(VALUE_TYPE, _fma)(mValue, xValue, zProd);				
				}
				
				dM += hackSize;
			}
#endif
	
		}
	
		diags -= warpSize;
		offsets += warpSize;
	}


	// Since z and y are accessed with the same offset by the same thread,
	// and the write to z follows the y read, y and z can share the same base address (in-place computing).
	
	if (i >= rows)
		return;
	
	if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
		z[i] = CONCAT(VALUE_TYPE, _fma)(beta, yVal, CONCAT(VALUE_TYPE, _mul) (alpha, zProd));
	else
		z[i] = CONCAT(VALUE_TYPE, _mul)(alpha, zProd);
}

// Force to recompile and optimize with llvm
__global__ void
CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _krn_b0)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* dM, const int* offsets, int hackSize, const int* hackOffsets, int rows, int cols, const VALUE_TYPE *x, int hackCount)
{
	CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _)
		(z, y, alpha, dM, offsets, hackSize, hackOffsets, rows, cols, x, CONCAT(zero_,VALUE_TYPE)(), hackCount);
}

__global__ void
CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _krn)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* dM, const int* offsets, int hackSize, const int* hackOffsets, int rows, int cols, const VALUE_TYPE *x, VALUE_TYPE beta, int hackCount)
{
	CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _)
		(z, y, alpha, dM, offsets, hackSize, hackOffsets, rows, cols, x, beta, hackCount);
}

void
CONCAT(_,GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL))
(spgpuHandle_t handle, VALUE_TYPE* z, const VALUE_TYPE *y, VALUE_TYPE alpha, 
	const VALUE_TYPE* dM, const int* offsets, int hackSize, const int* hackOffsets, int rows, int cols,
	const VALUE_TYPE *x, VALUE_TYPE beta)
{
	dim3 block (THREAD_BLOCK);
	dim3 grid ((rows + THREAD_BLOCK - 1) / THREAD_BLOCK);

	int hackCount = (rows + hackSize - 1)/hackSize;
	
#ifdef ENABLE_CACHE
	bind_tex_x ((const TEX_FETCH_TYPE *) x);
#endif

	cudaFuncSetCacheConfig(CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _krn), cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _krn_b0), cudaFuncCachePreferL1);

	if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
		CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _krn) <<< grid, block, block.x*sizeof(int), handle->currentStream >>> (z, y, alpha, dM, offsets, hackSize, hackOffsets, rows, cols, x, beta, hackCount);
	else
		CONCAT(GEN_SPGPU_HDIA_NAME(TYPE_SYMBOL), _krn_b0) <<< grid, block, block.x*sizeof(int), handle->currentStream >>> (z, y, alpha, dM, offsets, hackSize, hackOffsets, rows, cols, x, hackCount);

#ifdef ENABLE_CACHE
  	unbind_tex_x ((const TEX_FETCH_TYPE *) x);
#endif

}

