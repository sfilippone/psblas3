#include "coo_conv.h"
#include "core.h"

#include <map>

// returns the number of non-zero blocks
int computeBcooSize(int blockRows, int blockCols, const int* rows, const int* cols, int nonZeros)
{
	// use a map to count al non zero blocks
	std::map<unsigned long long int, int> blocksPositions;
	
	int blockCount = 0;
	
	int i;
	for (i=0; i<nonZeros; ++i)
	{
		int rowId = rows[i];
		int colId = cols[i];
		
		unsigned long long int blockRowId = rowId/blockRows;
		unsigned long long int blockColId = colId/blockCols;
		
		unsigned long long int blockId = blockRowId | (blockColId << 32);
		
		std::map<unsigned long long ,int>::iterator it = blocksPositions.find(blockId);
		
		// not found	
		if(it == blocksPositions.end())
		{
			blocksPositions[blockId] = blockCount;
			++blockCount;
		}
	}
	
	return blockCount;
}


void cooToBcoo(int* bRows, int* bCols, void* blockValues, /*int isBlockColumnMajor,*/ int blockRows, int blockCols, 
	const int* rows, const int* cols, const void* values, int nonZeros, spgpuType_t valuesType)
{
	// use a map to count al non zero blocks
	std::map<unsigned long long int, int> blocksPositions;
	
	int blockCount = 0;
	
	size_t elementSize = spgpuSizeOf(valuesType);
	size_t blockElementSize = elementSize*blockRows*blockCols;
	int i;
	for (i=0; i<nonZeros; ++i)
	{
		int rowId = rows[i];
		int colId = cols[i];
		
		unsigned long long int blockRowId = rowId/blockRows;
		unsigned long long int blockColId = colId/blockCols;
		
		unsigned long long int blockId = blockRowId | (blockColId << 32);
		
		std::map<unsigned long long int,int>::iterator it = blocksPositions.find(blockId);
		
		int blockPos;
		
		// not found	
		if(it == blocksPositions.end())
		{
			blocksPositions[blockId] = blockCount;
			blockPos = blockCount;
			
			bRows[blockCount] = blockRowId;
			bCols[blockCount] = blockColId;
			
			memset((char*)blockValues + blockCount*blockElementSize, 0, blockElementSize);
			
			++blockCount;
		}
		else
			blockPos = it->second;
		
		int blockRowOffset = rowId % blockRows;
		int blockColOffset = colId % blockCols;
		
		int blockOffset;
		
		//if (isBlockColumnMajor)
			blockOffset = blockRowOffset + blockColOffset*blockRows;
		/*else
			blockOffset = blockRowOffset*blockCols + blockColOffset;
		*/
		
		char* dest = (char*)blockValues + blockPos*blockElementSize + blockOffset*elementSize;
		const char* src = (const char*)values + elementSize*i;
			
		memcpy(dest, src, elementSize);
	}
}


