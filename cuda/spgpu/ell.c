#include "ell.h"
#include "ell_conv.h"
#include "stdlib.h"

void computeEllRowLenghts(
	int *ellRowLengths,
	int *ellMaxRowSize,
	int rowsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	int cooBaseIndex
	)
{
	// find the max number of non zero per row
	int maxRowSize = 0;
	int i;
	for (i=0; i<rowsCount; i++) 
		ellRowLengths[i] = 0;

	for (i=0; i<nonZerosCount; i++)
		++ellRowLengths[cooRowIndices[i] - cooBaseIndex];

	for (i=0; i<rowsCount; i++)
	{
		int currCount = ellRowLengths[i];
		if (currCount > maxRowSize) 
			maxRowSize = currCount;
	}

	*ellMaxRowSize = maxRowSize;
}

int computeEllAllocPitch(int rowsCount)
{
	// returns a pitch good for indices and values
	return ((rowsCount + 31)/32)*32;
}

void cooToEll(
	void *ellValues,
	int *ellIndices,
	int ellValuesPitch,
	int ellIndicesPitch,
	int ellMaxRowSize,
	int ellBaseIndex,
	int rowsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	spgpuType_t valuesType
	)
{	

	size_t elementSize = spgpuSizeOf(valuesType);
		
	// fill values and indices
	int* currentPos = (int*)malloc(rowsCount*sizeof(int));
	int i;
	
	for (i=0; i<rowsCount; i++)
		currentPos[i] = 0;

	for (i=0; i<nonZerosCount; i++)
	{
		int argRow = cooRowIndices[i] - cooBaseIndex;

		void* currentCm = ((char*)ellValues + argRow*elementSize) + currentPos[argRow]*ellValuesPitch*elementSize;
		void* currentRp = ((char*)&ellIndices[argRow]) + currentPos[argRow]*ellIndicesPitch*sizeof(int);

		*((int*)currentRp) = cooColsIndices[i] - cooBaseIndex + ellBaseIndex;
		
		memcpy(currentCm, (char*)cooValues + i*elementSize, elementSize);

		currentPos[argRow]++;

	}
	free(currentPos);
}




void merge(int *app_dstRs, int *app_rIdx, int *dstRs, int *rIdx, int start, int center, int end, int size) {
	int i, j, k; 

	i = start;
	j = center+1;
	k = 0;
 
	while ((i<=center) && (j<=end)) {
	
		if(dstRs[i] > dstRs[j]) {
			app_dstRs[k] = dstRs[i];
			app_rIdx[k] = rIdx[i];
			
			++k; ++i;
		} else {
			app_dstRs[k] = dstRs[j];
			app_rIdx[k] = rIdx[j];
			
			++k; ++j;
		}
	}
 
	while (i<=center) 
	{
		app_dstRs[k] = dstRs[i];
		app_rIdx[k] = rIdx[i];
			
		++k; ++i;
	}
 
	while (j<=end) 
	{
		app_dstRs[k] = dstRs[j];
		app_rIdx[k] = rIdx[j];
			
		++k; ++j;
	}
 
	for (k=start; k<=end; k++)
	{
		dstRs[k] = app_dstRs[k-start];
		rIdx[k] = app_rIdx[k-start];
	}
}
 
void mergesort(int *dstRs, int *rIdx, int size) {
	int* app_dstRs = (int*)malloc(size*sizeof(int));
	int* app_rIdx = (int*)malloc(size*sizeof(int));

	int sizetomerge=size-1;
	size--;
	int i;
	int n=2;
 
	while (n<sizetomerge*2) {
		for (i=0; (i+n-1)<=sizetomerge; i+=n) {
			merge(app_dstRs, app_rIdx, dstRs, rIdx, i,(i+i+n-1)/2,i+(n-1),sizetomerge); 
		}
 
		i--;
		if ((sizetomerge+1)%n!=0) {
			if (size>sizetomerge)
				merge (app_dstRs, app_rIdx, dstRs, rIdx, sizetomerge -((sizetomerge)%n),sizetomerge,size,size);
			sizetomerge=sizetomerge-((sizetomerge+1)%n);}
		n=n*2;
	}
 
	if (size>sizetomerge) 
		merge (app_dstRs, app_rIdx, dstRs,rIdx,0,size-(size-sizetomerge),size,size);
		
	free(app_dstRs);
	free(app_rIdx);
}



void ellToOell(
	int *rIdx,
	void *dstEllValues,
	int *dstEllIndices,
	int *dstRs,
	const void *srcEllValues,
	const int *srcEllIndices,
	const int *srcRs,
	int ellValuesPitch,
	int ellIndicesPitch,
	int rowsCount,
	spgpuType_t valuesType
	)
{
	size_t elementSize = spgpuSizeOf(valuesType);

	int i,j;
	for (i=0; i<rowsCount; ++i)
	{
		rIdx[i] = i;
		dstRs[i] = srcRs[i];
	}
	// sort..
	mergesort(dstRs, rIdx, rowsCount);
	
	for(i=0;i<rowsCount;i++) 
	{ 
		//Copy a row
		int srcId = rIdx[i];
		int srcLen = srcRs[srcId];

		int k;
		for (k=0; k<srcLen; ++k)
		{
			void* dstCm = ((char*)dstEllValues + i*elementSize) + k*ellValuesPitch*elementSize;
			void* srcCm = ((char*)srcEllValues + srcId*elementSize) + k*ellValuesPitch*elementSize;
			memcpy(dstCm, srcCm, elementSize);

			dstEllIndices[i + k*ellIndicesPitch] = srcEllIndices[srcId + k*ellIndicesPitch];
		}
	}
}
