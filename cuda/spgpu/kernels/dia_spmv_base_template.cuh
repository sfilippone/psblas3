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
CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* dM, const int* offsets, int dMPitch, int rows, int cols, int diags, const VALUE_TYPE *x, VALUE_TYPE beta)
{
	int i = threadIdx.x + blockIdx.x * (blockDim.x);

	VALUE_TYPE yVal = CONCAT(zero_,VALUE_TYPE)();

	if (i < rows && CONCAT(VALUE_TYPE, _isNotZero(beta)))
		yVal = y[i];

	VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();

	dM += i;

	extern __shared__ int offsetsChunk[];

	int rounds = (diags + blockDim.x - 1)/blockDim.x;
	
	for (int r = 0; r < rounds; r++)
	{
		// in the last round diags will be <= blockDim.x
		if (threadIdx.x < diags)
			offsetsChunk[threadIdx.x] = offsets[threadIdx.x];
	
		__syncthreads();
	
		if (i < rows)
		{
			int count = min(diags, blockDim.x );
			

#ifdef USE_PREFETCHING	
			int j;
			for (j=0; j<=count-3; j += 3)
			{
				// Prefetch 3 values
				int column1 = offsetsChunk[j] + i;
				int column2 = offsetsChunk[j+1] + i;
				int column3 = offsetsChunk[j+2] + i;	
				
				bool inside1 = column1 >= 0 && column1 < cols;
				bool inside2 = column2 >= 0 && column2 < cols;
				bool inside3 = column3 >= 0 && column3 < cols;			

				// Anticipate global memory read
				
				VALUE_TYPE xValue1, xValue2, xValue3;
				VALUE_TYPE mValue1, mValue2, mValue3;
				
				if(inside1)
                		{
                			mValue1 = dM[0];
#ifdef ENABLE_CACHE
					xValue1 = fetchTex (column1);
#else
					xValue1 = x[column1];
#endif				
				}
				dM += dMPitch;
							
				if(inside2)
                		{
                			mValue2 = dM[0];
#ifdef ENABLE_CACHE
					xValue2 = fetchTex (column2);
#else
					xValue2 = x[column2];
#endif				
				}
				dM += dMPitch;
							
				if(inside3)
                		{
                			mValue3 = dM[0];
#ifdef ENABLE_CACHE
					xValue3 = fetchTex (column3);
#else
					xValue3 = x[column3];
#endif				
				}
				dM += dMPitch;
				
				if(inside1)
					zProd = CONCAT(VALUE_TYPE, _fma)(mValue1, xValue1, zProd);
				if(inside2)
					zProd = CONCAT(VALUE_TYPE, _fma)(mValue2, xValue2, zProd);
				if(inside3)
					zProd = CONCAT(VALUE_TYPE, _fma)(mValue3, xValue3, zProd);
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
				
				dM += dMPitch;
			}
#else
			for (int j=0; j<count; ++j)
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

				dM += dMPitch;
			}
#endif
		}
		
		diags -= blockDim.x;
		offsets += blockDim.x;
		__syncthreads();
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
CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _krn_b0)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* dM, const int* offsets, int dMPitch, int rows, int cols, int diags, const VALUE_TYPE *x)
{
	CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _)
		(z, y, alpha, dM, offsets, dMPitch, rows, cols, diags, x, CONCAT(zero_,VALUE_TYPE)());
}

__global__ void
CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _krn)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* dM, const int* offsets, int dMPitch, int rows, int cols, int diags, const VALUE_TYPE *x, VALUE_TYPE beta)
{
	CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _)
		(z, y, alpha, dM, offsets, dMPitch, rows, cols, diags, x, beta);
}

void
CONCAT(_,GEN_SPGPU_DIA_NAME(TYPE_SYMBOL))
(spgpuHandle_t handle, VALUE_TYPE* z, const VALUE_TYPE *y, VALUE_TYPE alpha, 
	const VALUE_TYPE* dM, const int* offsets, int dMPitch, int rows, int cols, int diags,
	const VALUE_TYPE *x, VALUE_TYPE beta)
{
	dim3 block (THREAD_BLOCK );
	dim3 grid ((rows + THREAD_BLOCK  - 1) / THREAD_BLOCK );

#ifdef ENABLE_CACHE
	bind_tex_x ((const TEX_FETCH_TYPE *) x);
#endif

	cudaFuncSetCacheConfig(CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _krn), cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _krn_b0), cudaFuncCachePreferL1);

	if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
		CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _krn)
		<<< grid, block, block.x*sizeof(int), handle->currentStream >>> (z, y, alpha, dM, offsets, dMPitch, rows, cols, diags, x, beta);
	else
		CONCAT(GEN_SPGPU_DIA_NAME(TYPE_SYMBOL), _krn_b0) <<< grid, block, block.x*sizeof(int), handle->currentStream >>> (z, y, alpha, dM, offsets, dMPitch, rows, cols, diags, x);

#ifdef ENABLE_CACHE
  	unbind_tex_x ((const TEX_FETCH_TYPE *) x);
#endif

}

