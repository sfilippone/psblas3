#include "hdia_conv.h"
#include "stdlib.h"
#include "string.h"

#include <vector>
#include <map>

int getHdiaHacksCount(int hackSize, int rowsCount)
{
	return (rowsCount + hackSize - 1)/hackSize;
}

void computeHdiaHackOffsets(
	int *allocationHeight,
	int *hackOffsets,
	int hackSize,
	const void* diaValues,
	int diaValuesPitch,	
	int diagonals,
	int rowsCount,
	spgpuType_t valuesType
	)
{
	int i,r,s, hack;
	int hackCount = getHdiaHacksCount(hackSize, rowsCount);
	
	size_t elementSize = spgpuSizeOf(valuesType);
	
	int hackHeight = 0;
		
	hackOffsets[0] = 0;
	for (hack=0; hack<hackCount; ++hack)
	{
		for (i=0; i<diagonals; ++i)
		{
			int found = 0;
			for (r=0; r<hackSize; ++r)
			{
				int row = hack*hackSize + r;
				
				if (row >= rowsCount)
					break;
			
				const char* val = (char*)diaValues + elementSize*(row + i*diaValuesPitch);
				
				for (s=0; s<elementSize; ++s)
					if (*(val+s) != 0)
					{
						found = 1;
						goto hackTest1;
					}
			}
hackTest1:			
			if (found != 0)
				++hackHeight;
		}	
		hackOffsets[hack+1] = hackHeight;
	}
	
	*allocationHeight = hackOffsets[hackCount];	
}






void diaToHdia(
	void *hdiaValues,
	int *hdiaOffsets,
	const int *hackOffsets,
	int hackSize,
	const void* diaValues,
	const int* diaOffsets,
	int diaValuesPitch,	
	int diagonals,
	int rowsCount,
	spgpuType_t valuesType
	)
{
	int i,r,s;
	int hack;
	int hackCount = getHdiaHacksCount(hackSize, rowsCount);
	
	// Compute offsets
	int hackOffsetsSize = hackCount + 1;
	
	size_t elementSize = spgpuSizeOf(valuesType);
	
	for (hack=0; hack<hackCount; ++hack)
	{
		int posOffset = hackOffsets[hack];
		
		int hackHeight = 0;
		for (i=0; i<diagonals; ++i)
		{
			int found = 0;
			for (r=0; r<hackSize; ++r)
			{
				int row = hack*hackSize + r;
				
				if (row >= rowsCount)
					break;
			
				const char* val = (const char*)diaValues + elementSize*(row + i*diaValuesPitch);
				
				for (s=0; s<elementSize; ++s)
					if (*(val+s) != 0)
					{
						found = 1;
						goto hackTest2;
					}
			}
hackTest2:			
			if (found != 0)
			{
				// use hdiaOffsets to temporarely store i, instead of diaOffsets[i]
				hdiaOffsets[posOffset + hackHeight++] = i;

			}
		}
	}
	
	// Copy values
	for (hack=0; hack<hackCount; ++hack)
	{
		// get diagonal offset
		int posOffset = hackOffsets[hack];
		int hackDiags = hackOffsets[hack+1] - posOffset;
		
		for (i=0; i<hackDiags; ++i)
		{
			int diagPosInsideDia = hdiaOffsets[posOffset + i];
			int diagOffset;
			
			// reupdate hdiOffsets with the correct value
			hdiaOffsets[posOffset + i] = diagOffset = diaOffsets[diagPosInsideDia];
			
			for (r=0; r<hackSize; ++r)
			{
				int row = hack*hackSize + r;
				
				if (row >= rowsCount)
					break;
			
				char* dest = (char*)hdiaValues + elementSize*((posOffset + i)*hackSize + r);
				const char* src = (const char*)diaValues + elementSize*(row + diagPosInsideDia*diaValuesPitch);
			
				memcpy(dest, src, elementSize);
			}
		}
	}
}







void computeHdiaHackOffsetsFromCoo(
	int *allocationHeight,
	int *hackOffsets,
	int hackSize,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	int cooBaseIndex
	)
{	

	int i,j,h;
	
	int hackCount = getHdiaHacksCount(hackSize, rowsCount);
	
	
	// Find rows per hack
	std::vector<int> *rowsPerHack = new std::vector<int> [hackCount];
	
	for (i=0; i<nonZerosCount; ++i)
	{
		int rowIdx = cooRowIndices[i];
			
		// It is inside current hack
		int h = (rowIdx-cooBaseIndex)/hackSize;
		rowsPerHack[h].push_back(i);
	}
			

	// Use hackOffsets to deduce hack's heights
	std::map<int, int> diagIdsToPos;

	hackOffsets[0] = 0;
	for (h=0; h<hackCount; ++h)
	{
		diagIdsToPos.clear();
			
		int diagonalsCount = 0;
	
		std::vector<int> *hackRows = &rowsPerHack[h];
		int hackRowsSize = hackRows->size();
		
		for (j=0; j<hackRowsSize; ++j)
		{
			i = hackRows->at(j);
			int rowIdx = cooRowIndices[i];
			int colIdx = cooColsIndices[i];
			int diagId = (colIdx-cooBaseIndex) - ((rowIdx-cooBaseIndex) % hackSize);
			int diagPos = hackSize - 1 + diagId;
		
			std::map<int,int>::iterator it = diagIdsToPos.find(diagPos);
			
			if(it == diagIdsToPos.end())
			{
				diagIdsToPos[diagPos] = 1;
				++diagonalsCount;
			}
		}		
				
		hackOffsets[h+1] = hackOffsets[h] + diagonalsCount;
	} 
	
	*allocationHeight = hackOffsets[hackCount];
	
	delete[] rowsPerHack;
}

void cooToHdia_size(
	void *hdiaValues,
	int *hdiaOffsets,
	const int *hackOffsets,
	int hackSize,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	size_t elementSize
	)
{	
	int i,j,h;
	
	int hackCount = getHdiaHacksCount(hackSize, rowsCount);

	// Find rows per hack
	std::vector<int> *rowsPerHack = new std::vector<int> [hackCount];
	
	for (i=0; i<nonZerosCount; ++i)
	{
		int rowIdx = cooRowIndices[i];
			
		// It is inside current hack
		int h = (rowIdx-cooBaseIndex)/hackSize;
		rowsPerHack[h].push_back(i);
	}
		

	std::map<int, int> hackDiagIdsToPos;
	
	for (h=0; h<hackCount; ++h)
	{
		int diagonalsCount = 0;
	
		hackDiagIdsToPos.clear();

		std::vector<int> *hackRows = &rowsPerHack[h];
		int hackRowsSize = hackRows->size();
		
		for (j=0; j<hackRowsSize; ++j)
		{
			i = hackRows->at(j);
			
			int rowIdx = cooRowIndices[i];
			int colIdx = cooColsIndices[i];
			int globalDiagId = colIdx - rowIdx;
			int diagId = (colIdx - cooBaseIndex) - ((rowIdx - cooBaseIndex) % hackSize);
			int diagPos = hackSize - 1 + diagId;
		
			std::map<int,int>::iterator it = hackDiagIdsToPos.find(diagPos);
			
			if(it == hackDiagIdsToPos.end())
			{
				hackDiagIdsToPos[diagPos] = globalDiagId;
			}
		}	
	
		// Reorder diags
		for (std::map<int, int>::iterator it = hackDiagIdsToPos.begin(); it != hackDiagIdsToPos.end(); ++it)
		{
			int i = it->first;
			
			int globalDiagId = it->second;
			int diagPosInsideOffsets;
			int diagId = i - hackSize + 1;
			hackDiagIdsToPos[i] = diagPosInsideOffsets = diagonalsCount++;
			hdiaOffsets[diagPosInsideOffsets] = globalDiagId;
		}

	
		hdiaOffsets += diagonalsCount;
		
		for (j=0; j<hackRowsSize; ++j)
		{
			i = hackRows->at(j);
			int rowIdx = cooRowIndices[i];
			int colIdx = cooColsIndices[i];
			int diagId = (colIdx - cooBaseIndex) - ((rowIdx - cooBaseIndex) % hackSize);
		
			int diagPosInsideOffsets = hackDiagIdsToPos[hackSize - 1 + diagId];
		
			char* valAddr = (char*)hdiaValues + 
				elementSize*(((rowIdx - cooBaseIndex) % hackSize) 
					+ hackSize* (hackOffsets[h] + diagPosInsideOffsets));
		
			memcpy(valAddr, (const char*)cooValues + i*elementSize, elementSize);
		}
	}
	
	delete[] rowsPerHack;
}



void cooToHdia(
	void *hdiaValues,
	int *hdiaOffsets,
	const int *hackOffsets,
	int hackSize,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	spgpuType_t valuesType
	)
{
	size_t elementSize = spgpuSizeOf(valuesType);
	
	cooToHdia_size(hdiaValues, hdiaOffsets,
		hackOffsets, hackSize, rowsCount,
		columnsCount, nonZerosCount, 
		cooRowIndices, cooColsIndices, cooValues, cooBaseIndex, elementSize);
}

void bcooToBhdia(
	void *hdiaValues,
	int *hdiaOffsets,
	const int *hackOffsets,
	int hackSize,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	spgpuType_t valuesType,
	int blockSize
	)
{
	size_t elementSize = blockSize*spgpuSizeOf(valuesType);
	
	cooToHdia_size(hdiaValues, hdiaOffsets,
		hackOffsets, hackSize, rowsCount,
		columnsCount, nonZerosCount, 
		cooRowIndices, cooColsIndices, cooValues, cooBaseIndex, elementSize);
}

