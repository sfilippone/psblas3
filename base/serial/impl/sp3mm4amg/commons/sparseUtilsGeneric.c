/*
 *              Sp3MM_for_AlgebraicMultiGrid
 *    (C) Copyright 2021-2022
 *        Andrea Di Iorio      
 * 
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *    2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions, and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *    3. The name of the Sp3MM_for_AlgebraicMultiGrid or the names of its contributors may
 *       not be used to endorse or promote products derived from this
 *       software without specific written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE Sp3MM_for_AlgebraicMultiGrid GROUP OR ITS CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */ 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include "macros.h"
#include "sparseMatrix.h"
#include "utils.h"

///////////////   no offset --> SINGLE implementation functions
#ifndef SPARSE_UTILS_C
#define SPARSE_UTILS_C
///DENSE accumulator utils
/*
 * alloc threads' aux arrays once and split them in threads' structures
 * so free them once from the first thread struct, with the original pointers returned from the alloc
 */
ACC_DENSE* _initAccVectors_monoalloc(idx_t num,idx_t size){ //TODO PERF WITH NEXT
	ACC_DENSE* out	 = NULL;
	double* vAll	 = NULL;
	idx_t* vAllNzIdx = NULL;
	if (!(out = calloc(num,sizeof(*out)))){
		ERRPRINT("_initAccVectors aux struct alloc failed\n");
		return NULL;
	}
	if (!(vAll = calloc(num*size,sizeof(*vAll)))) {
		ERRPRINT("_initAccVectors aux dense vectors alloc failed\n");
		goto err;
	}
	if (!(vAllNzIdx = calloc(num*size,sizeof(*vAllNzIdx)))) {
		ERRPRINT("_initAccVectors aux dense vectors' idx alloc failed\n");
		goto err;
	}
	for (idx_t i=0; i<num; i++){
		out[i].vLen	= size; //TODO USELESS INFO?
		out[i].v	= vAll	  + i*size;  
		out[i].nnzIdx	= vAllNzIdx + i*size;;
		out[i].nnzIdxMap.len  = 0;
	}
	return out;
	
	err:
	free(out);
	free(vAll);
	free(vAllNzIdx);
	return NULL;
}
int allocAccDense(ACC_DENSE* v,idx_t size){
		v->vLen = size; 
		if (!(v->v = calloc(size,sizeof(*(v->v))))) {
			ERRPRINT("_initAccVectors aux dense vector alloc failed\n");
			return EXIT_FAILURE;
		}
		if (!(v->nnzIdx = calloc(size,sizeof(*(v->nnzIdx))))) {
			ERRPRINT("_initAccVectors aux dense vector' idx alloc failed\n");
			return EXIT_FAILURE;
		}
		if (initSpVectIdxDenseAcc(size, &v->nnzIdxMap))	return EXIT_FAILURE;

		return EXIT_SUCCESS;
}

ACC_DENSE* _initAccVectors(idx_t num,idx_t size){
	ACC_DENSE* out	= NULL;
	if (!(out = calloc(num,sizeof(*out)))){
		ERRPRINT("_initAccVectors aux struct alloc failed\n");
		return NULL;
	}
	for (idx_t i=0; i<num; i++){
		if (allocAccDense(out+i,size))   goto _err;
	}
	return out;
	
	_err:
	for (idx_t i=0; i<num; i++){
		if (out[i].v)	   free(out[i].v);
		if (out[i].nnzIdx)  free(out[i].nnzIdx);
	}
	free(out);
	return NULL;
}

void freeAccsDense(ACC_DENSE* vectors,idx_t num){
	for (idx_t i=0; i<num; i++){
		free(vectors[i].v);
		free(vectors[i].nnzIdx);
		free(vectors[i].nnzIdxMap.idxsMap);
	}
	free(vectors);
}
void _freeAccsDenseChecks(ACC_DENSE* vectors,idx_t num){ 
	if (!vectors)   return;
	for (idx_t i=0; i<num; i++){
		if(vectors[i].v)	free(vectors[i].v);
		if(vectors[i].nnzIdx)   free(vectors[i].nnzIdx);
	}
	free(vectors);
}

int initSpVectIdxDenseAcc(idx_t idxMax,SPVECT_IDX_DENSE_MAP* vectIdxsMap){
	vectIdxsMap->len = 0;
	nnz_idxs_flags_t* idxMaps = &vectIdxsMap->idxsMap;
	#if SPVECT_IDX_BITWISE == TRUE //nnz presence falgs as bitflags in limbs
	vectIdxsMap->idxsMapN = INT_DIV_CEIL(idxMax, LIMB_SIZE_BIT);
	#else	//nnz presence falgs in a single array
	vectIdxsMap->idxsMapN = idxMax;
	#endif //SPVECT_IDX_BITWISE == TRUE
	if (!(*idxMaps = calloc(vectIdxsMap->idxsMapN, sizeof(**idxMaps)))) {
		ERRPRINT("initSpVectIdxDenseAcc\tidxMaps SPVECT_IDX_BITWISE callc err\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

void checkOverallocPercent(idx_t* forecastedSizes,spmat* AB){
	for (idx_t r=0,rSize,forecastedSize; r < AB->M; r++){
		forecastedSize = forecastedSizes[r];
		#ifdef ROWLENS
		rSize = AB->RL[r];
		#else
		rSize = AB->IRP[r+1] - AB->IRP[r];
		#endif
		DEBUGCHECKS{
			if ( forecastedSize < rSize ){
				ERRPRINT("BAD FORECASTING\n");
				assert(forecastedSize >= rSize );
			}
		}
		DEBUGPRINT
			printf("extra forecastedSize of row: %d\t=\t%lf %% \n",
			  r,100*(forecastedSize-rSize) / (double) forecastedSize);
	}
	idx_t extraMatrix = forecastedSizes[AB->M] - AB->NZ;
	printf("extra forecastedSize of the matrix: \t%d\t = %lf %% \n",
	  extraMatrix, 100*extraMatrix /(double) forecastedSizes[AB->M]);
}
int spmatDiff(spmat* A, spmat* B){
	if (A->NZ != B->NZ){
		ERRPRINT("NZ differ\n");
		return EXIT_FAILURE;
	}
	if (memcmp(A->IRP,B->IRP,A->M)){
		ERRPRINT("IRP differ\n");
		return EXIT_FAILURE;
	}
	if (doubleVectorsDiff(A->AS,B->AS,A->NZ,NULL)){
		ERRPRINT("AS DIFFER\n");
		return EXIT_FAILURE;
	}
	if (memcmp(A->JA,B->JA,A->NZ)){
		ERRPRINT("JA differ\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

double* CSRToDense(spmat* sparseMat){
	double* denseMat;
	idx_t i,j,idxNZ,denseSize;
	if (__builtin_smul_overflow(sparseMat->M,sparseMat->N,&denseSize)){
		ERRPRINT("overflow in dense allocation\n");
		return NULL;
	}
	if (!(denseMat = calloc(denseSize, sizeof(*denseMat)))){
		ERRPRINT("dense matrix alloc failed\n");
		return  NULL;
	}
	for (i=0;i<sparseMat->M;i++){
		for (idxNZ=sparseMat->IRP[i]; idxNZ<sparseMat->IRP[i+1]; ++idxNZ){
			 j = sparseMat->JA[idxNZ];
			 //converting sparse item into dense entry
			 denseMat[(idx_t) IDX2D(i,j,sparseMat->N)] = sparseMat->AS[idxNZ]; 
		}
	}
	return denseMat;
}
void printSparseMatrix(spmat* spMatrix,char justNZMarkers){
	double* denseMat = CSRToDense(spMatrix);
	if (!denseMat)  return;
	printMatrix(denseMat,spMatrix->M,spMatrix->N,justNZMarkers);
	free(denseMat);
}



static inline int _colsPartitioningUnifRanges_init(spmat* A,int gridCols,
  spmat** colParts,idx_t** colPartsLens){

	spmat* colPart;
	idx_t _colBlock = A->N/gridCols, _colBlockRem = A->N%gridCols, *tmpJA;
	///alloc/init partitions structures
	if (!(*colParts = calloc(gridCols, sizeof(**colParts)))){
		ERRPRINT("colsPartitioningUnifRanges\tcolumns partitions of A calloc fail\n");
		return EXIT_FAILURE;
	}
	for (idx_t i=0,colBlock; i<gridCols; i++){
		colBlock = UNIF_REMINDER_DISTRI(i,_colBlock,_colBlockRem);
		colPart  = *colParts + i;
		if (allocSpMatrixInternal(A->M,colBlock,colPart)){
			ERRPRINT("colsPartitioningUnifRanges\tallocSpMatrixInternal partition err\n");
			return EXIT_FAILURE;
		}
		//TODO TODO overalloc A cols partitions NZ arrays, then realloc
		if (!(colPart->AS = malloc(A->NZ * sizeof(*A->AS)))){
			ERRPRINT("colPart of A overalloc of AS errd\n");
			return EXIT_FAILURE;
		}
		if (!(colPart->JA = malloc(A->NZ * sizeof(*A->JA)))){
			ERRPRINT("colPart of A overalloc of JA errd\n");
			return EXIT_FAILURE;
		}
	}
	//for each A col partition -> last copied nz index = nnz copied ammount
	if (! (*colPartsLens = calloc(gridCols, sizeof(**colPartsLens))) ) {
		ERRPRINT("colsPartitioningUnifRanges: colPartsLens calloc errd\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

static inline int _colsPartitioningUnifRanges_finalRealloc(spmat* A,int gridCols,
  spmat* colParts,idx_t* colPartsLens){

	spmat* colPart;
	double* tmpAS; idx_t* tmpJA;
	//realloc overallcd A parts NZ arrays (TODO -> downsizing -> nofails?)
	for (idx_t i=0,partLen; i<gridCols; i++){
		colPart = colParts + i;
		partLen = colPartsLens[i];
		if (!(tmpAS = realloc(colPart->AS,partLen*sizeof(*(colPart->AS))))){
			ERRPRINT("realloc overallocated cols partition AS array\n");
			return EXIT_FAILURE;
		}
		colPart->AS = tmpAS;
		if (!(tmpJA = realloc(colPart->JA,partLen*sizeof(*(colPart->JA))))){
			ERRPRINT("realloc overallocated cols partition JA array\n");
			return EXIT_FAILURE;
		}
		colPart->JA		 = tmpJA;
		colPart->NZ		 = partLen;
		colPart->IRP[A->M]  = partLen;
	}
	return EXIT_SUCCESS;
}
#endif

 
#ifndef OFF_F
	#error generic implementation requires OFF_F defined
#endif
////////////////////////  CSR SPECIFIC -- TODO RENAME //////////////////
///SPARSE MATRIX PARTITIONING
idx_t* CAT(colsOffsetsPartitioningUnifRanges_,OFF_F)(spmat* A,int gridCols){
	idx_t subRowsN = A->M * gridCols;
	idx_t _colBlock = A->N/gridCols, _colBlockRem = A->N%gridCols;
	idx_t* offsets = malloc( (subRowsN+1) * sizeof(*offsets) );
	if (!offsets)  {
		ERRPRINT("colsOffsetsPartitioningUnifRanges:\toffsets malloc errd\n");
		return NULL;
	}
	///OFFSETS COMPUTE FOR COL GROUPS -> O( A.NZ )
	for (idx_t r=0, j=0;	 r<A->M;	 j=A->IRP[++r]-OFF_F){
		offsets[ IDX2D(r,0,gridCols) ] = j;  //row's first gc start is costrained
		//navigate column groups inside current row
		for (idx_t gc=1,gcStartCol;  gc<gridCols;  gc++){
			gcStartCol = UNIF_REMINDER_DISTRI_STARTIDX(gc,_colBlock,_colBlockRem);
			//goto GroupCols start entry,keeping A's nnz entries navigation (idx j)
			//for (idx_t c=A->JA[j]-OFF_F; c<gcStartCol && j < A->IRP[r+1]-OFF_F; c=A->JA[++j]-OFF_F);
			while (j < A->IRP[r+1]-OFF_F &&  A->JA[j]-OFF_F < gcStartCol)  
				j++;
			offsets[ IDX2D(r,gc,gridCols) ] = j;  //row's gc group startIdx
		}
	}
	offsets[subRowsN] = A->NZ;  //last row's partition end costrained
	return offsets;
}

spmat* CAT(colsPartitioningUnifRangesOffsetsAux_,OFF_F)(spmat* A,int gridCols,
	   idx_t** colPartsOffsets)
{
	spmat *colParts = NULL, *colPart;
	idx_t _colBlock = A->N/gridCols, _colBlockRem = A->N%gridCols;
	idx_t *colPartsLens=NULL, *tmpJA;
	double* tmpAS;

	///alloc/init partitions structures
	idx_t* colOffsets  = NULL;
	if (!(colOffsets = CAT(colsOffsetsPartitioningUnifRanges_,OFF_F)(A,gridCols))) 
		goto _err;
	if (_colsPartitioningUnifRanges_init(A,gridCols,&colParts,&colPartsLens))
		goto _err;
	//OFFSET BASED COPY OF A.COL_GROUPS -> O( A.NZ )
	for (idx_t r=0,gcId=0;	 r<A->M;	r++){
		for (idx_t gc=0,gcStartIdx=0,gLen=0;  gc<gridCols; gc++,gcId++){
			//gcId		= IDX2D(r,gc,gridCols); //seqent. scan of parts
			colPart		= colParts + gc;
			gcStartIdx 	= colOffsets[ gcId   ];
			gLen		= colOffsets[ gcId+1 ] - gcStartIdx;
			colPart->IRP[r] = colPartsLens[gc];	//new line for the col partition
			//actual copy of nnz entries to colPartitions
			memcpy(colPart->AS+colPart->IRP[r], A->AS+gcStartIdx, gLen*sizeof(*A->AS));
			memcpy(colPart->JA+colPart->IRP[r], A->JA+gcStartIdx, gLen*sizeof(*A->JA));
			colPartsLens[gc] += gLen;
			#ifdef ROWLENS
			colPart->RL[r] = i;
			#endif
		}
	}

	//realloc overallcd A parts NZ arrays
	if(_colsPartitioningUnifRanges_finalRealloc(A,gridCols,colParts,colPartsLens))
		goto  _err;

	free(colPartsLens);
	if (colPartsOffsets)	*colPartsOffsets = colOffsets;	//save for the caller
	else			free(colOffsets);

	return colParts;
	_err:
	free(*colPartsOffsets);
	for (idx_t i=0; i<gridCols; i++)   	freeSpmatInternal(colParts+i);
	free(colOffsets);
	free(colParts);
	free(colPartsLens);
	return NULL;
}
spmat* CAT(colsPartitioningUnifRanges_,OFF_F)(spmat* A,int gridCols){
	spmat *colParts, *colPart;
	idx_t _colBlock = A->N/gridCols, _colBlockRem = A->N%gridCols, *colPartsLens=NULL, *tmpJA;
	double* tmpAS;
	///alloc/init partitions structures
	if (_colsPartitioningUnifRanges_init(A,gridCols,&colParts,&colPartsLens))
		goto _err;
	/* TODO
	 * Parallelize: 2for collapse OMP, gcEndCol -> startIdxize, ...
	 * oppure wrappare cio in static inline 
	 */
	for (idx_t r=0, j=0;	 r<A->M;	 j=A->IRP[++r]-OFF_F){
		//navigate column groups inside current row
		for (idx_t gc=0,gcEndCol=0,i;  gc<gridCols ;  gc++,j+=i){
			i = 0;  //@i=len current subpartition of row @r to copy
			colPart = colParts + gc;
			//NB: not shifting IRP because is handled as internal implementation component
 			//But JA idx memcopied -> kept as they were originally, handled with shift in functions
			colPart->IRP[r] = colPartsLens[gc];	
			gcEndCol += UNIF_REMINDER_DISTRI(gc,_colBlock,_colBlockRem);
			//goto next GroupCols,keeping A's nnz entries navigation ( index j+i )
			//for (idx_t c=A->JA[j+i]-OFF_F; c<gcEndCol && j+i  < A->IRP[r+1]-OFF_F; c=A->JA[j+ ++i]-OFF_F);
			while ( j+i < A->IRP[r+1]-OFF_F && A->JA[j+i]-OFF_F < gcEndCol ) i++;
			memcpy(colPart->AS+colPart->IRP[r], A->AS+j, i*sizeof(*A->AS));
			memcpy(colPart->JA+colPart->IRP[r], A->JA+j, i*sizeof(*A->JA));
			
			colPartsLens[gc] += i;
			#ifdef ROWLENS
			colPart->RL[r] = i;
			#endif
		}
	}
	//realloc overallcd A parts NZ arrays
	if(_colsPartitioningUnifRanges_finalRealloc(A,gridCols,colParts,colPartsLens))
		goto  _err;
	free(colPartsLens);
	return colParts;
	_err:
	for (idx_t i=0; i<gridCols; i++)   freeSpmatInternal(colParts+i);
	if(colParts)		free(colParts);
	if(colPartsLens)	free(colPartsLens);
	return NULL;
}


void CAT(checkOverallocRowPartsPercent_,OFF_F)(idx_t* forecastedSizes,spmat* AB,
  				   idx_t gridCols,idx_t* bColOffsets){
	idx_t* abColOffsets = CAT(colsOffsetsPartitioningUnifRanges_,OFF_F)(AB, gridCols);
	assert(abColOffsets);	//partitioning error
	for (idx_t rLen,forecast,partId=0; partId<AB->M*gridCols; partId++){
		forecast = forecastedSizes[partId];
		rLen = abColOffsets[partId+1] - abColOffsets[partId];
		DEBUGCHECKS	assert(forecast >= rLen);
	}
	idx_t extraMatrix = forecastedSizes[AB->M] - AB->NZ;
	printf("extra forecastedSize of the matrix: \t%d\t = %lf %% \n",
		extraMatrix, 100*extraMatrix /(double) forecastedSizes[AB->M]);

	free(abColOffsets);
	
}


#ifdef SPARSEUTILS_MAIN_TEST	///unit test embbeded

///inline export here 
//SPMV_CHUNKS_DISTR spmvChunksFair; 
spmat* allocSpMatrix(idx_t rows, idx_t cols);
int allocSpMatrixInternal(idx_t rows, idx_t cols, spmat* mat);
void freeSpmatInternal(spmat* mat);
void freeSpmat(spmat* mat);

////INTERNAL TEST FUNCTIONS
//test that each row's partition from colsOffsetsPartitioningUnifRanges is in the correct index range
#include <alloca.h>
int testColsOffsetsPartitioningUnifRanges(spmat* mat,idx_t gridCols,idx_t* partsOffs){
	idx_t _colBlock = mat->N/gridCols, _colBlockRem = mat->N%gridCols;
	idx_t j=0;	//CSR scanning nnz idx
	idx_t* colPartsPopulations = alloca(gridCols * sizeof(*colPartsPopulations));
	memset(colPartsPopulations,0,gridCols * sizeof(*colPartsPopulations));
	for (idx_t r=0,pId=0; r<mat->M; r++){
		for (idx_t gc=0,pStartIdx,pEndIdx; gc<gridCols; gc++,pId++){
			pStartIdx = UNIF_REMINDER_DISTRI_STARTIDX(gc,_colBlock,_colBlockRem);
			pEndIdx   = UNIF_REMINDER_DISTRI_STARTIDX(gc+1,_colBlock,_colBlockRem)-1; 
			//pId=IDX2D(r,gc,gridCols);
			for (idx_t idx=partsOffs[pId],c; idx<partsOffs[pId+1]; idx++,j++){
				c = mat->JA[idx];
				assert(j == idx); //consecutive index in partitioning
				assert(pStartIdx <= c && c <= pEndIdx);	//colRange
				assert(mat->IRP[r] <= idx && idx <= mat->IRP[r+1] ); //rowRange
			}
			colPartsPopulations[gc] += partsOffs[pId+1] - partsOffs[pId]; 
		}
	}
	assert(j == mat->NZ);
	idx_t s=0;
	for (idx_t gc=0,partSize; gc < gridCols; gc++,s+=partSize){
		partSize = colPartsPopulations[gc];
		double partShare=partSize/(double)mat->NZ,partsAvg=1/(double)gridCols;
		double partShareAvgDiff = partShare - partsAvg;
		printf("colPartition %d has:\t%d = %lf of NNZ\t\t .. %lf\tAVG diff\n",
		  gc,partSize,partShare,partShareAvgDiff);
	}
	assert(s == mat->NZ); //TODO DUPLICATED
	return EXIT_SUCCESS;
}

CONFIG Conf = {
	.gridRows = 8,
	.gridCols = 8,
};

#include "parser.h"
int main(int argc, char** argv){
	int out=EXIT_FAILURE;
	if (init_urndfd())  return out;
	if (argc < 2 )  {ERRPRINT("COO MATRIX FOR TEST"); return out;}
	////parse sparse matrix and dense vector
	spmat* mat;
	char* trgtMatrix = TMP_EXTRACTED_MARTIX;
	if (extractInTmpFS(argv[1],TMP_EXTRACTED_MARTIX) < 0)
		trgtMatrix = argv[1];
	if (!(mat = MMtoCSR(trgtMatrix))){
		ERRPRINT("err during conversion MM -> CSR\n");
		return out;
	}
	////partitioning test
	idx_t* colsPartitions = colsOffsetsPartitioningUnifRanges_0(mat,Conf.gridCols);
	if (!colsPartitions)	goto _free;
	if (testColsOffsetsPartitioningUnifRanges(mat,Conf.gridCols,colsPartitions))
		goto _free;

	out=EXIT_SUCCESS;
	printf("testColsOffsetsPartitioningUnifRanges passed with "
		   "mat: %dx%d-%dNNZ\tgrid: %dx%d\n",
			mat->M,mat->N,mat->NZ,Conf.gridRows,Conf.gridCols);
	_free:
	if (colsPartitions) free(colsPartitions);

	return out;
}
#endif //SPARSEUTILS_MAIN_TEST
