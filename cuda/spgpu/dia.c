#include "dia.h"
#include "dia_conv.h"
#include "stdlib.h"

int computeDiaAllocPitch(int rowsCount)
{
	// returns a pitch good for indices and values
	return ((rowsCount + 31)/32)*32;
}

int computeDiaDiagonalsCount(
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices)
{
	int* diagIds = (int*)malloc((rowsCount + columnsCount - 1)*sizeof(int));
	int diagonalsCount = 0;
	int i;
	
	for (i=0; i<(rowsCount + columnsCount - 1); ++i)
		diagIds[i] = -1;

	for (i=0; i<nonZerosCount; ++i)
	{
		int diagPos = rowsCount - 1 + cooColsIndices[i] - cooRowIndices[i];
		
		if (diagIds[diagPos] < 0)
		{
			diagIds[diagPos] = diagonalsCount++;
		}
	}
	
	free(diagIds);
	
	return diagonalsCount;
}

void coo2dia(
	void* values,
	int* offsets,
	int valuesPitch,	
	int diagonals,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	spgpuType_t valuesType)
{
	int* diagIdsToPos = (int*)malloc((rowsCount + columnsCount - 1)*sizeof(int));
	int i;
	int diagonalsCount = 0;
	
	for (i=0; i<(rowsCount + columnsCount - 1); ++i)
		diagIdsToPos[i] = -1;


	for (i=0; i<nonZerosCount; ++i)
	{
		int rowIdx = cooRowIndices[i];
		int colIdx = cooColsIndices[i];
		int diagId = colIdx - rowIdx;
		int diagPos = rowsCount - 1 + diagId;
		
		if (diagIdsToPos[diagPos] < 0)
		{
			diagIdsToPos[diagPos] = 1;
		}
	}	
	
	// Reorder diags
	for (i=0; i<(rowsCount + columnsCount - 1); ++i)
	{
		if (diagIdsToPos[i] == 1)
		{
			int diagPosInsideOffsets;
			int diagId = i - rowsCount + 1;
			diagIdsToPos[i] = diagPosInsideOffsets = diagonalsCount++;
			offsets[diagPosInsideOffsets] = diagId;
		}
	}
	

	for (i=0; i<nonZerosCount; ++i)
	{
		int rowIdx = cooRowIndices[i];
		int colIdx = cooColsIndices[i];
		int diagId = colIdx - rowIdx;
		
		int diagPosInsideOffsets = diagIdsToPos[rowsCount - 1 + diagId];
		
		size_t elementSize = spgpuSizeOf(valuesType);
		
		void* valAddr = values + elementSize*(rowIdx-cooBaseIndex) + diagPosInsideOffsets*elementSize*valuesPitch;
		
		memcpy(valAddr, (const char*)cooValues + i*elementSize, elementSize);
	}
	
	free(diagIdsToPos);
}

