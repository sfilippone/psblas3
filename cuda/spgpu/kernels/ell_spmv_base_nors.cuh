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

__device__ void
CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_4_noRs)
(int i, VALUE_TYPE yVal, int outRow,
	VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, int rows, int maxNnzPerRow, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
	VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();

	__shared__ VALUE_TYPE temp[2][THREAD_BLOCK+1];
	
	if (i < rows)
	{
		rP += i; cM += i;

		int rowSizeM = maxNnzPerRow / 4;
		
				
		if ((maxNnzPerRow % 4) > threadIdx.y)
			++rowSizeM;
		
		rP += threadIdx.y*rPPitch; 
		cM += threadIdx.y*cMPitch;
		
		
		for (int j = 0; j < rowSizeM; j++)
		{
			int pointer;
			VALUE_TYPE value;
			VALUE_TYPE fetch;
		
			pointer = rP[0] - baseIndex;
			rP += 4*rPPitch; 

			value = cM[0];
			cM +=  4*cMPitch;

#ifdef ENABLE_CACHE
			fetch = fetchTex(pointer);
#else
			fetch = x[pointer];
#endif	

			// avoid MAD on pre-Fermi
			zProd = CONCAT(VALUE_TYPE, _fma)(value, fetch, zProd);
		}

		// Reduction
		if (threadIdx.y > 1)
			temp[threadIdx.y - 2][threadIdx.x] = zProd;
	}
	
	__syncthreads();
	
	if (i < rows)
	{
		if (threadIdx.y <= 1)
			zProd = CONCAT(VALUE_TYPE, _add)(zProd, temp[threadIdx.y][threadIdx.x]);
		
		if (threadIdx.y == 1)
			temp[1][threadIdx.x] = zProd;
	}
	
	__syncthreads();
	
	if (i < rows)
	{
		if (threadIdx.y == 0)	
		{
			zProd = CONCAT(VALUE_TYPE, _add)(zProd, temp[1][threadIdx.x]);
		
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
CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_2_noRs)
(int i, VALUE_TYPE yVal, int outRow,
	VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, int maxNnzPerRow, const int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
	VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();

	__shared__ VALUE_TYPE temp[THREAD_BLOCK];
	
	if (i < rows)
	{
		rP += i; cM += i;

		int rowSizeM = maxNnzPerRow / 2;
		
		if (threadIdx.y == 0)
		{
			if (maxNnzPerRow % 2)
				++rowSizeM;
		}
		else
		{
			rP += rPPitch; 
			cM += cMPitch;
		}
		
		
		for (int j = 0; j < rowSizeM; j++)
		{
			int pointer;
			VALUE_TYPE value;
			VALUE_TYPE fetch;
		
			pointer = rP[0] - baseIndex;
			rP += rPPitch; 
			rP += rPPitch;

			value = cM[0];
			cM += cMPitch;
			cM += cMPitch;

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
	}
	
	__syncthreads();
	
	if (i < rows)
	{
		if (threadIdx.y == 0)	
		{
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
CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_noRs)
(int i, VALUE_TYPE yVal, int outRow,
	VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, int maxNnzPerRow, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
	VALUE_TYPE zProd = CONCAT(zero_,VALUE_TYPE)();

	if (i < rows)
	{
		rP += i; cM += i;

#ifdef USE_PREFETCHING		
		for (int j = 0; j < maxNnzPerRow / 2; j++)
		{
			int pointers1, pointers2;
			VALUE_TYPE values1, values2;
			VALUE_TYPE fetches1, fetches2;
		
			pointers1 = rP[0] - baseIndex;
			rP += rPPitch; 
			pointers2 = rP[0] - baseIndex;
			rP += rPPitch; 

			values1 = cM[0];
			cM += cMPitch;
			
			values2 = cM[0];
			cM += cMPitch;

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
		if (maxNnzPerRow % 2)
	    	{
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
		for (int j = 0; j < maxNnzPerRow; j++)
		{
			int pointer;
			VALUE_TYPE value;
			VALUE_TYPE fetch;
		
			pointer = rP[0] - baseIndex;
			rP += rPPitch;

			value = cM[0];
			cM += cMPitch;

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
CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _krn_ridx_noRs)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, const int* rIdx, int maxNnzPerRow, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
	int i = threadIdx.x + blockIdx.x * (THREAD_BLOCK);
	
	VALUE_TYPE yVal = CONCAT(zero_,VALUE_TYPE)();
	int outRow = 0;
	if (i < rows)
	{

		outRow = rIdx[i];
		if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
			yVal = y[outRow];
	}
	
	if (blockDim.y == 1)
		CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_noRs)
			(i, yVal, outRow, z, y, alpha, cM, rP, cMPitch, rPPitch, maxNnzPerRow, rows, x, beta, baseIndex);
	else //if (blockDim.y == 2)
	
		CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_2_noRs)
			(i, yVal, outRow, z, y, alpha, cM, rP, cMPitch, rPPitch, maxNnzPerRow, rows, x, beta, baseIndex);
	/*
	else if (blockDim.y == 4)
	
	 
		CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_4_noRs)
			(i, yVal, outRow, z, y, alpha, cM, rP, cMPitch, rPPitch, maxNnzPerRow, rows, x, beta, baseIndex);
			*/
}


__device__ void
CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _noRs)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, int maxNnzPerRow,  int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
	int i = threadIdx.x + blockIdx.x * (THREAD_BLOCK);
	
	VALUE_TYPE yVal = CONCAT(zero_,VALUE_TYPE)();

	if (i < rows)
	{
		if (CONCAT(VALUE_TYPE, _isNotZero(beta)))
			yVal = y[i];
	
	}
	
	if (blockDim.y == 1)
		CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_noRs)
			(i, yVal, i, z, y, alpha, cM, rP, cMPitch, rPPitch, maxNnzPerRow, rows, x, beta, baseIndex);
	
	else //if (blockDim.y == 2)
	
		CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_2_noRs)
			(i, yVal, i, z, y, alpha, cM, rP, cMPitch, rPPitch, maxNnzPerRow, rows, x, beta, baseIndex);
	/*
	else if (blockDim.y == 4)
	
		CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _ridx_4_noRs)
			(i, yVal, i, z, y, alpha, cM, rP, cMPitch, rPPitch, maxNnzPerRow, rows, x, beta, baseIndex);
			*/
			
}

// Force to recompile and optimize with llvm
__global__ void
CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _krn_b0_noRs) 
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, int maxNnzPerRow, int rows, const VALUE_TYPE *x, int baseIndex)
{
	CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _noRs)
		(z, y, alpha, cM, rP, cMPitch, rPPitch, maxNnzPerRow, rows, x, CONCAT(zero_,VALUE_TYPE)(), baseIndex);
}

__global__ void
CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _krn_noRs)
(VALUE_TYPE *z, const VALUE_TYPE *y, VALUE_TYPE alpha, const VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, int maxNnzPerRow, int rows, const VALUE_TYPE *x, VALUE_TYPE beta, int baseIndex)
{
	CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _noRs)
		(z, y, alpha, cM, rP, cMPitch, rPPitch, maxNnzPerRow, rows, x, beta, baseIndex);
}
