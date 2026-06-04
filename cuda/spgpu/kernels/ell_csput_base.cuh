/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */


#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)

#undef GEN_SPGPU_ELL_NAME
#define GEN_SPGPU_ELL_NAME(x) CONCAT(CONCAT(spgpu,x),ellcsput)

#define THREAD_BLOCK 256

// Define:
//#define VALUE_TYPE
//#define TYPE_SYMBOL

#include "mathbase.cuh"


__global__ void
CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _krn)
(VALUE_TYPE alpha, VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, const int* rS,
	int nnz, int *aI, int *aJ, VALUE_TYPE *aVal, int baseIndex)
{
	int i = threadIdx.x + blockIdx.x * (THREAD_BLOCK);
		
	if (i < nnz)
	{
		int row = aI[i] - baseIndex;
		int column = aJ[i];
		VALUE_TYPE value = aVal[i];
		
		if (row < 0)
			return;
		
		// Select row
		cM += row;
		rP += row;
		rS += row;
		
		// Binary search
		int lower = 0;
		int upper = (*rS - 1);
		
		while(lower <= upper) 
		{
			int medium = (lower + upper) / 2;
			
			int currentColumn = rP[medium*rPPitch];
			
			if(currentColumn == column) 
			{
				cM[medium*cMPitch] = value;
				break;
			}
			else if(currentColumn < column)
			    lower = medium + 1;
			else
			    upper = medium - 1;
		}
	}
}



void
CONCAT(_,GEN_SPGPU_ELL_NAME(TYPE_SYMBOL))
(spgpuHandle_t handle, VALUE_TYPE alpha, VALUE_TYPE* cM, const int* rP, int cMPitch, int rPPitch, const int* rS,
	int nnz, int *aI, int *aJ, VALUE_TYPE *aVal, int baseIndex)
{
	dim3 block (THREAD_BLOCK, 1);
	
	dim3 grid ((nnz + THREAD_BLOCK - 1) / THREAD_BLOCK);

	CONCAT(GEN_SPGPU_ELL_NAME(TYPE_SYMBOL), _krn)
		<<< grid, block, 0, handle->currentStream >>> (alpha, cM, rP, cMPitch, rPPitch, rS, nnz, aI, aJ, aVal, baseIndex);
}

void
GEN_SPGPU_ELL_NAME(TYPE_SYMBOL)
	(spgpuHandle_t handle, 
	VALUE_TYPE alpha, 
	VALUE_TYPE* cM, 
	const int* rP, 
	int cMPitch, 
	int rPPitch, 
	const int* rS,
	int nnz, 
	int *aI, 
	int *aJ, 
	VALUE_TYPE *aVal, 
	int baseIndex)
{
	int maxNForACall = max(handle->maxGridSizeX, THREAD_BLOCK*handle->maxGridSizeX);

	while (nnz > maxNForACall) //managing large vectors
	{
		CONCAT(_,GEN_SPGPU_ELL_NAME(TYPE_SYMBOL))
			(handle, alpha, cM, rP, cMPitch, rPPitch, rS, maxNForACall, aI, aJ, aVal, baseIndex);

		aI = aI + maxNForACall;
		aJ = aJ + maxNForACall;
		aVal = aVal + maxNForACall;
		
		nnz -= maxNForACall;
	}
	
	CONCAT(_,GEN_SPGPU_ELL_NAME(TYPE_SYMBOL))
		(handle, alpha, cM, rP, cMPitch, rPPitch, rS, nnz, aI, aJ, aVal, baseIndex);			
	
	cudaCheckError("CUDA error on ell_csput");
}

