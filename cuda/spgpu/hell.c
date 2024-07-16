#include "hell.h"
#include "hell_conv.h"

void computeHellAllocSize(
	int* allocationHeight,
	int hackSize,
	int rowsCount,
	const int *ellRowLengths
	)
{
	int totalLen = 0;
	int i;
	int remainings;
	int done;
	int maxLen;

	for (i=0; i<rowsCount/hackSize; ++i)
	{
		int maxLen = 0;
		int j;
		for (j=0; j<hackSize; ++j)
		{
			int row = i*hackSize + j;
			int curLen = ellRowLengths[row];
			if (curLen > maxLen)
				maxLen = curLen;
		}
		totalLen += maxLen;
	}

	remainings = rowsCount % hackSize;
	done = (rowsCount/hackSize)*hackSize;
	maxLen = 0;
	
	for (i=0; i<remainings; ++i)
	{
		int row = done + i;
		int curLen = ellRowLengths[row];
		if (curLen > maxLen)
			maxLen = curLen;
	}
	
	*allocationHeight = totalLen + maxLen;
}

void ellToHell(
	void *hellValues,
	int *hellIndices,
	int* hackOffsets,
	int hackSize,

	const void *ellValues,
	const int *ellIndices,
	int ellValuesPitch,
	int ellIndicesPitch,
	int *ellRowLengths,
	int rowsCount,
	spgpuType_t valuesType
	)
{

	size_t elementSize = spgpuSizeOf(valuesType);
	
	int hacks = (rowsCount + hackSize - 1)/hackSize;
	
	char* currValPos = (char*)hellValues;
	int* currIndPos = hellIndices;

	int hackOffset = 0;
	int i;
	for (i=0; i<hacks; ++i)
	{
		int maxLen = 0;
		int j;
		hackOffsets[i] = hackOffset;

		for (j=0; j<hackSize; ++j)
		{
			int row = i*hackSize + j;
			int rowLen;
			int k;

			if (row >= rowsCount)
				break;

			rowLen = ellRowLengths[row];

			if (rowLen > maxLen)
				maxLen = rowLen;

			for (k=0; k<rowLen; ++k)
			{
				memcpy(currValPos + (j + k*hackSize)*elementSize,
				 (((char*)ellValues + k*ellValuesPitch*elementSize) + row*elementSize),
				 elementSize);
				currIndPos[j + k*hackSize] = *(ellIndices + k*ellIndicesPitch + row);
			}
		}

		hackOffset += hackSize*maxLen;
		currValPos = currValPos + hackSize*maxLen*elementSize;
		currIndPos += hackSize*maxLen;
	}
}
