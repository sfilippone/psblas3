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
#include "sparseMatrix.h"
#include "config.h"

#ifndef OFF_F
	#pragma error OFF_F required
#endif
//////////////////////////////// CSR SPECIFIC /////////////////////////////////
///SPARSE MATRIX PARTITIONING
/*
 * partition CSR sparse matrix @A in @gridCols columns partitions 
 * returning an offsets matrix out[i][j] = start of jth colPartition of row i
 * subdivide @A columns in uniform cols ranges in the output 
 */
idx_t* CAT(colsOffsetsPartitioningUnifRanges_,OFF_F)(spmat* A,uint gridCols);

/*
 * partition CSR sparse matrix @A in @gridCols columns partitions as 
 * indipended and allocated sparse matrixes and return them
 * subdivide @A columns in uniform cols ranges in the output 
 */
spmat* CAT(colsPartitioningUnifRanges_,OFF_F)(spmat* A,uint gridCols);
//same as above but with (aux) use of offsets partitoning (also returned if colOffsets!=NULL
spmat* CAT(colsPartitioningUnifRangesOffsetsAux_,OFF_F)(spmat* A,uint gridCols,idx_t** colPartsOffsets);

//same as checkOverallocPercent but with 2D partitioning - CSR col partitioning
void CAT(checkOverallocRowPartsPercent_,OFF_F)(ulong* forecastedSizes,spmat* AB,
				  	      idx_t gridCols,idx_t* bColOffsets);
///////////////////////////////////////////////////////////////////////////////

///Single implementations headers
#ifndef SPARSEUTILS_H_COMMON_IDX_IMPLS
#define SPARSEUTILS_H_COMMON_IDX_IMPLS

//shift every index about the sparse data for use the matrix in a fortran app
inline void C_FortranShiftIdxs(spmat* m){
	for(ulong r=0; r<m->M+1; m -> IRP[r]++,r++);
	for(ulong i=0; i<m->NZ;  m -> JA[i]++, i++);
}
//shift every index about the sparse data for use the matric in a C app
inline void Fortran_C_ShiftIdxs(spmat* m){	//TODO DBG ONLY and compleatness
	for(ulong r=0; r<m->M+1; m -> IRP[r]--,r++);
	for(ulong i=0; i<m->NZ;  m -> JA[i]--, i++);
}

/*
 * check SpMM resulting matrix @AB = A * B nnz distribution in rows
 * with the preallocated, forecasted size in @forecastedSizes 
 * in @forecastedSizes there's for each row -> forecasted size 
 * and in the last entry the cumulative of the whole matrix
 */
void checkOverallocPercent(ulong* forecastedSizes,spmat* AB);
/*  
	check if sparse matrixes A<->B differ up to 
	DOUBLE_DIFF_THREASH per element
*/
int spmatDiff(spmat* A, spmat* B);
////dyn alloc of spMM output matrix
/*
///size prediction of AB = @A * @B
inline ulong SpMMPreAlloc(spmat* A,spmat* B){
	//TODO BETTER PREALLOC HEURISTICS HERE 
	return MAX(A->NZ,B->NZ);
}
//init a sparse matrix AB=@A * @B with a initial allocated space by an euristic
inline spmat* initSpMatrixSpMM(spmat* A, spmat* B){
	spmat* out;
	if (!(out = allocSpMatrix(A->M,B->N)))  return NULL;
	out -> NZ = SpMMPreAlloc(A,B);
	if (!(out->AS = malloc(out->NZ*sizeof(*(out->AS))))){
		ERRPRINT("initSpMatrix: out->AS malloc errd\n");
		free(out);
		return NULL;
	}
	if (!(out->JA = malloc(out->NZ*sizeof(*(out->JA))))){
		ERRPRINT("initSpMatrix: out->JA malloc errd\n");
		freeSpmat(out);
		return NULL;
	}
	return out;
}

#define REALLOC_FACTOR  1.5
//realloc sparse matrix NZ arrays
inline int reallocSpMatrix(spmat* mat,ulong newSize){
	mat->NZ *= newSize;
	void* tmp;
	if (!(tmp = realloc(mat->AS,mat->NZ * sizeof(*(mat->AS))))){
		ERRPRINT("reallocSpMatrix:  realloc AS errd\n");
		return EXIT_FAILURE;
	}
	mat->AS = tmp;
	if (!(tmp = realloc(mat->JA,mat->NZ * sizeof(*(mat->JA))))){
		ERRPRINT("reallocSpMatrix:  realloc JA errd\n");
		return EXIT_FAILURE;
	}
	mat->JA = tmp;
	return EXIT_SUCCESS;
}
*/
////MISC
//print useful information about 3SPMM about to compute
void print3SPMMCore(spmat* R,spmat* AC,spmat* P,CONFIG* conf);
void printSparseMatrix(spmat* sparseMat,char justNZMarkers);
/*convert @sparseMat sparse matrix in dense matrix returned*/
double* CSRToDense(spmat* sparseMat);

void freeAccsDense(ACC_DENSE* vectors,ulong num);

#endif //SPARSEUTILS_H_COMMON_IDX_IMPLS 
