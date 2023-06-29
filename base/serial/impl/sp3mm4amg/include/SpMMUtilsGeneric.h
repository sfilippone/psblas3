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
///scalar-vector multiply
//TODO che guadagno si ha ad utilizzare solo la versione generica delle successive 2 funzioni
/*
 * Sparse vector part <->scalar multiplication in a dense output
 * sparse vector part will hold nnz values in @vectVals 		 (AS subvector)
 * with corresponding indexes in @vectIdxs in range [0,@vectLen] (JA subvector)
 * resulting vector accumulated in a dense array in @aux->v, along with nnzIdx
 * both of accumulator's dense array and nnzIdx in @aux and has to be big @vectLen
 */
inline void CAT(scSparseVectMul_,OFF_F)
  (double scalar,double* vectVals,ulong* vectIdxs,ulong vectLen, ACC_DENSE* aux){
	for (ulong i=0,j; i<vectLen; i++){
		j = vectIdxs[i]-OFF_F;
		DEBUGCHECKS{
			if (j>=aux->vLen){
				fprintf(stderr,"index %lu outside vLen %lu\n",j,aux->vLen);
				assert(j < aux->vLen);
			}
		}
		aux->v[j] += vectVals[i] * scalar;  //accumulate
		//append new nonzero index to auxVNNZeroIdxs for quick sparsify
		//if (!(aux->v[j]))	aux->nnzIdx[ aux->nnzIdxMap.len++ ] = j; //TODO numerical zero may then cause readd
		if(!spVect_idx_in(j,&aux->nnzIdxMap))
			aux->nnzIdx[ aux->nnzIdxMap.len-1 ] = j;
	}
}

/* Same as above scSparseVectMul_ , but consider an initial offset to remove from each idx
 * Sparse vector part <->scalar multiplication in a dense output
 * @vectVals:	sparse vector part values	( from spmat.AS )
 * with [at least] @vectLen corresponding 
 * target idxs in @vectIdxs (from spmat.JA ) starting from @startIdx 
 *
 * Resulting vector accumulated in a dense array in @aux->v, along with nnzIdx
 * all nnz values indexes will be shifted back of @startIdx in @aux
 * both of accumulator's dense array and nnzIdx in @aux and has to be big @vectLen
 */
inline void CAT(scSparseVectMulPart_,OFF_F)(double scalar,double* vectVals,
  ulong* vectIdxs,ulong vectLen,ulong startIdx,ACC_DENSE* aux){
	for (ulong i=0,j; i<vectLen; i++){
		j = vectIdxs[i]-OFF_F - startIdx;
		DEBUGCHECKS{
			if (j>=aux->vLen){
				fprintf(stderr,"index %lu outside vLen %lu\n",j,aux->vLen);
				assert(j < aux->vLen);
			}
		}
		aux->v[j] += vectVals[i] * scalar;  //accumulate
		//append new nonzero index to auxVNNZeroIdxs for quick sparsify
		//if (!(aux->v[j]))	aux->nnzIdx[ aux->nnzIdxMap.len++ ] = j; 
		////TODO numerical zero may then cause readd
		if(!spVect_idx_in(j,&aux->nnzIdxMap))
			aux->nnzIdx[ aux->nnzIdxMap.len-1 ] = j;
	}
}


///MultiImplementations functions
#ifndef OFF_F
	#error generic implementations requires OFF_F defined
#endif

///TODO IMPLEMENT SCALAR<->ROW MUL AS GENERIC SPARSE VECTOR<->SCALAR MUL
//inline void scSparseRowMul(double scalar,spmat* mat,ulong trgtR, ACC_DENSE* aux){
inline void CAT(scSparseRowMul_,OFF_F)(double scalar,spmat* mat,ulong trgtR, ACC_DENSE* aux){
	ulong  rowStartIdx = mat->IRP[trgtR]-OFF_F,rowLen;
	#ifdef ROWLENS
	rowLen = mat->RL[trgtR];
	#else
	rowLen = mat->IRP[trgtR+1] - mat->IRP[trgtR];
	#endif
	CAT(scSparseVectMul_,OFF_F)(scalar,mat->AS+rowStartIdx,mat->JA+rowStartIdx,rowLen,aux);
	//TODO check impact of generic version use
}


///OUTPUT SIZE PREDICTION	

/*			O(A.NZ [+B.M] )
 * return array of upper bounded row sizes of the product @A * @B
 * also appended at the end for the cumulative total size of the matrix AB = A*B
 */
inline idx_t* CAT(spMMSizeUpperbound_,OFF_F)(spmat* A,spmat* B){
	AUDIT_INTERNAL_TIMES	Start = omp_get_wtime();
	idx_t* rowSizes = calloc((A->M+1),  sizeof(*rowSizes));
	if (!rowSizes){
		ERRPRINT("spMMSizeUpperbound: rowSizes calloc errd\n");
		return NULL;
	}
	idx_t fullMatBound = 0;
	#pragma omp parallel for schedule(static) reduction(+:fullMatBound)
	for (idx_t r=0;  r<A->M; r++){
		for (idx_t jj=A->IRP[r] - OFF_F,j,rlen;  jj<A->IRP[r+1] - OFF_F; jj++){
			j = A->JA[jj] - OFF_F;
			#ifdef ROWLENS
			rlen = B->RL[j];
			#else
			rlen = B->IRP[j+1] - B->IRP[j]; //OFF_F: delta is offset independent
			#endif 
			rowSizes[r]	 += rlen;
			fullMatBound	+= rlen;
			//rowSizes[A->M]  += rlen;	//just below omp reduction sum
		}
	}
	rowSizes[A->M]  = fullMatBound;
	AUDIT_INTERNAL_TIMES	End= omp_get_wtime();
	VERBOSE 
		printf("spMMSizeUpperbound:%lu\t%le s\n",rowSizes[A->M],End-Start);
	return rowSizes;
}

/*	O(A.NZ + B.NZ)
 * return matrix @A.M x gridCols of upper bound 
 * for each of the @gridCols col partitions of the output matrix AB = @A * @B 
 * also appended at the end for the cumulative total size of the matrix AB 
 */
inline idx_t* CAT(spMMSizeUpperboundColParts_,OFF_F)
  (spmat* A,spmat* B,ushort gridCols,idx_t* bColPartOffsets){
	AUDIT_INTERNAL_TIMES	Start = omp_get_wtime();
	
	idx_t* rowPartsSizes = calloc((A->M*gridCols +1),  sizeof(*rowPartsSizes));
	if (!rowPartsSizes){
		ERRPRINT("spMMSizeUpperbound: rowPartsSizes calloc errd\n");
		return NULL;
	}

	idx_t fullMatBound = 0;
	#pragma omp parallel for schedule(static) reduction(+:fullMatBound)
	for (idx_t r=0;  r<A->M; r++){
		//for each A.row -> sum B.colParts lens		
		for (idx_t jj=A->IRP[r]-OFF_F,j,rlen;  jj<A->IRP[r+1]-OFF_F; jj++){
			j = A->JA[jj] - OFF_F;
			for (idx_t gc=0,bPartID=IDX2D(j,gc,gridCols);  gc < gridCols; gc++,bPartID++){
				rlen = bColPartOffsets[bPartID+1] - bColPartOffsets[bPartID];
				rowPartsSizes[ IDX2D(r,gc,gridCols) ]	+= rlen;
				fullMatBound	+= rlen;
			}
		}
	}
	rowPartsSizes[ A->M*gridCols ]  = fullMatBound;
	AUDIT_INTERNAL_TIMES	End= omp_get_wtime();
	VERBOSE 
	 printf("spMMSizeUpperboundColParts_:%lu\t%le s\n",rowPartsSizes[A->M],End-Start);
	return rowPartsSizes;
}


////Single implementations funcs
#ifndef SPMMUTILS_H_SINGLE_IMPLEMENTATION
#define SPMMUTILS_H_SINGLE_IMPLEMENTATION

#include "utils.h"
///Allocs - Free
//SpMM holder of accumulats 
inline SPMM_ACC* initSpMMAcc(ulong entriesNum, ulong accumulatorsNum){
	SPMM_ACC* out = calloc(1,sizeof(*out));
	if (!out){
		ERRPRINT("initSpMMAcc:	out calloc errd\n");
		return NULL;
	}
	out->size = entriesNum;
	if (!(out->JA = malloc(entriesNum * sizeof(*(out->JA))))){
		ERRPRINT("initSpMMAcc:	JA malloc errd\n");
		goto _err;
	}
	if (!(out->AS = malloc(entriesNum * sizeof(*(out->AS))))){
		ERRPRINT("initSpMMAcc:	AS malloc errd\n");
		goto _err;
	}
	if (!(out->accs = malloc(accumulatorsNum * sizeof(*(out->accs))))){
		ERRPRINT("initSpMMAcc:	accs malloc errd\n");
		goto _err;
	}
	return out;

	_err:
	if (out->JA)	free(out->JA);
	if (out->AS)	free(out->AS);
	if (out->accs)  free(out->accs);
	if (out)		free(out);
	return NULL;
}
inline void freeSpMMAcc(SPMM_ACC* acc){
	free(acc->JA);
	free(acc->AS);
	free(acc->accs);
	free(acc);
}
///// dense acc sparsify functions
/*
 * sparsify dense accumulated vector @accV (with shifted of @startColAcc) 
 * into sparse accumulator @accSparse that'll use space for nnz entries from @acc
*/
///DIRECT OUT MATRIX SPARSIFY
/*
 * sparsify @accV directly inside row @row of matrix @m
 * considering, if given not NULL, 2D partitioning with 
 * @gridCols cols groups and colGroups offsets per row matrix @colPartsOffsets
 * :Returns	inplace modify of @m
 */
static inline void sparsifyDirect(ACC_DENSE* accV,spmat* m,idx_t row){
	idx_t nnz = accV->nnzIdxMap.len, sparsifyStart = m->IRP[row], sparsifyEnd = m->IRP[row+1];
	sort_idx_t(accV->nnzIdx,nnz); //sort nnz idx for ordered write
	DEBUGCHECKS	assertArrNoRepetitions(accV->nnzIdx,nnz);
	DEBUGCHECKS	assert(nnz <= (sparsifyEnd - sparsifyStart));

	for (idx_t i=0,j;	i < nnz;   i++){
		j =  accV->nnzIdx[i];
		m -> JA[sparsifyStart + i] = j; //+ startColAcc;
		m -> AS[sparsifyStart + i] = accV->v[j];
	}
}
//TODO 2D partitioning - colParts version of above...best to use multi impl trick of symbStep
static inline void sparsifyDirectColParts(ACC_DENSE* accV,spmat* m,idx_t row,
  ushort colGroupId,ushort gridCols, idx_t* colPartsOffsets,idx_t startCol){
	idx_t nnz = accV->nnzIdxMap.len;
	ushort partID = IDX2D(row,colGroupId,gridCols);
	idx_t sparsifyStart = colPartsOffsets[partID], sparsifyEnd = colPartsOffsets[partID+1];
	sort_idx_t(accV->nnzIdx,nnz); //sort nnz idx for ordered write
	DEBUGCHECKS		assertArrNoRepetitions(accV->nnzIdx,nnz);
	DEBUGCHECKS		assert( nnz == sparsifyEnd - sparsifyStart );
	sort_idx_t(accV->nnzIdx,nnz); //sort nnz idx for ordered write

	for (idx_t i=0,j;	i < nnz;   i++){ 
		j =  accV->nnzIdx[i];
		m -> JA[sparsifyStart + i] = j + startCol;
		m -> AS[sparsifyStart + i] = accV->v[j];
	}
	
}	

///UB SPACE SPARSIFY
//internal sparsivy dense acc inside (prepared) sparse acc struct
static inline void _sparsifyUB(ACC_DENSE* accV,SPACC* accSparse,idx_t startColAcc){
	idx_t nnz = accV->nnzIdxMap.len;
	sort_idx_t(accV->nnzIdx,nnz); //sort nnz idx for ordered write
	DEBUGCHECKS		assertArrNoRepetitions(accV->nnzIdx,nnz);
	for (idx_t i=0,j;	i < nnz;   i++){ 
		j =  accV->nnzIdx[i];
		accSparse -> JA[i] = j + startColAcc;
		accSparse -> AS[i] = accV->v[j];
	}
	accSparse -> len = nnz;
}

//row[Part] sparsified in a thread safe (exactly long) reserved area using atomics
static inline void sparsifyUBNoPartsBounds
  (SPMM_ACC* acc,ACC_DENSE* accV,SPACC* accSparse, ulong startColAcc){
	//sort nnz indexes of dense accumulator
	idx_t nnz = accV->nnzIdxMap.len;
	idx_t sparsifyStartV;		 //start index(inside @accSparse) of @accV to sparsify
	//sparsifyStartV = __atomic_fetch_add(&(acc->lastAssigned),nnz,__ATOMIC_ACQ_REL); 
	#pragma omp atomic capture
	{   //fetch and add like .... 
		sparsifyStartV = acc->lastAssigned;
		acc->lastAssigned += nnz;
	}
	DEBUGCHECKS{
		if (acc->lastAssigned >= acc->size){
			ERRPRINT("OMP ATOMIC OR SG ERRD IN SPACE ASSIGNMENTS...\n");
			assert(acc->lastAssigned < acc->size);
		}
	}
	//
	accSparse -> AS = acc->AS + sparsifyStartV; 
	accSparse -> JA = acc->JA + sparsifyStartV; 
	_sparsifyUB(accV,accSparse,startColAcc);
}
////output-gather functons
/*
 * merge @conf->gridCols*@mat->M sparse rows partitions into @mat
 * EXPECTED rowsParts @rowsParts to be sorted in accord to the 
 * 2D rowMajor computing grid given in @conf
 * allocd arrays to hold non zero values and indexes into @mat
 */
inline int mergeRowsPartitions(SPACC* rowsParts,spmat* mat,
  CONFIG* conf){
	ulong nzNum=0,j,rLen,idx,partsNum = mat->M * conf->gridCols;
	//TODO PARALLEL MERGE - SAVE STARTING offset OF EACH PARTITION IN THE OUT MATRIX
	ulong* rowsPartsOffsets=alloca(partsNum*sizeof(*rowsPartsOffsets));
	///count nnz entries and alloc arrays for them
	for (ulong r=0; r<mat->M; r++){
		//for each partition ->get len -> outMat.IRP and aux offsets  
		for (j=0,rLen=0; j<conf->gridCols; j++){
			idx = IDX2D(r,j,conf->gridCols);
			rowsPartsOffsets[idx]=nzNum+rLen;//part start=prev accumulated end
			rLen += rowsParts[idx].len;
		}
		nzNum += rLen;
		mat->IRP[r+1] = nzNum;
		#ifdef ROWLENS
		mat->RL[r] = rLen;
		#endif
	}
	mat->NZ = nzNum;
	if (!(mat->AS = malloc(nzNum * sizeof(*(mat->AS))))){
		ERRPRINT("merged sparse matrix AS alloc errd\n");
		return EXIT_FAILURE;
	}  
	if (!(mat->JA = malloc(nzNum * sizeof(*(mat->JA))))){
		ERRPRINT("merged sparse matrix JA alloc errd\n");
		return EXIT_FAILURE;
	}
	///popolate with rows nnz values and indexes
	ulong pLen; //omp for aux vars
	#pragma omp parallel for schedule(static) private(pLen)
	for (ulong i=0;  i<partsNum; i++){
		pLen = rowsParts[i].len;
		memcpy(mat->AS + rowsPartsOffsets[i],rowsParts[i].AS,pLen*sizeof(*(mat->AS)));
		memcpy(mat->JA + rowsPartsOffsets[i],rowsParts[i].JA,pLen*sizeof(*(mat->JA)));
	}
	CONSISTENCY_CHECKS{ //TODO REMOVE written nnz check manually
		for (ulong i=0,w=0; i<mat->M; i++){
			if (mat->IRP[i] != w) 
				{ERRPRINT("MERGE ROW ERR IRP\n");return -1;}
			for (j=0; j<conf->gridCols; j++){
				SPACC r = rowsParts[IDX2D(i,j,conf->gridCols)];
				for (ulong jj=0; jj<r.len; jj++,w++){
					if (mat->AS[w]!= r.AS[jj]){
						ERRPRINT("MERGE ROW ERR AS\n"); return -1;}
					if (mat->JA[w]!= r.JA[jj]){
						ERRPRINT("MERGE ROW ERR JA\n"); return -1;}
				}
			}
		}
	}
	return EXIT_SUCCESS;
}

/*
 * merge @mat->M sparse rows @rows in sparse matrix @mat
 * EXPECTED @rows to be sorted in accord to trgt matrix row index @r
 * allocd arrays to hold non zero values and indexes into @mat
 */
inline int mergeRows(SPACC* rows,spmat* mat){
	ulong nzNum=0;
	//count nnz entries and alloc arrays for them
	for (ulong r=0;   r<mat->M;   ++r){
		nzNum += rows[r].len;
		mat->IRP[r+1] = nzNum;
		#ifdef ROWLENS
		mat->RL[r]  = rows[r].len
		#endif
	}
	mat->NZ = nzNum;
	if (!(mat->AS = malloc(nzNum * sizeof(*(mat->AS))))){
		ERRPRINT("merged sparse matrix AS alloc errd\n");
		return EXIT_FAILURE;
	}  
	if (!(mat->JA = malloc(nzNum * sizeof(*(mat->JA))))){
		ERRPRINT("merged sparse matrix JA alloc errd\n");
		return EXIT_FAILURE;
	}
	///POPOLATE WITH ROWS NNZ VALUES AND INDEXES
	//TODO PARALLEL COPY
	#pragma omp parallel for schedule(static)
	for (ulong r=0; r<mat->M; r++){
		memcpy(mat->AS+mat->IRP[r], rows[r].AS, rows[r].len*sizeof(*(mat->AS)));
		memcpy(mat->JA+mat->IRP[r], rows[r].JA, rows[r].len*sizeof(*(mat->JA)));
	}
	CONSISTENCY_CHECKS{ //TODO REMOVE written nnz check manually
		for (ulong r=0,i=0; r<mat->M; r++){
			if (i != mat->IRP[r])
				{ERRPRINT("MERGE ROW ERR IRP\n");return -1;}
			for (ulong j=0; j<rows[r].len; j++,i++){
				if (mat->AS[i]!= rows[r].AS[j]){
					ERRPRINT("MERGE ROW ERR AS\n"); return -1;}
				if (mat->JA[i]   != rows[r].JA[j]){
					ERRPRINT("MERGE ROW ERR JA\n"); return -1;}
			}
		}
	}
	return EXIT_SUCCESS;
} 

#endif	//SPMMUTILS_H_SINGLE_IMPLEMENTATION

