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

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sparseMatrix.h"
#include "mmio.h"
#include "parser.h"
#include "macros.h"
#include "utils.h"

////COO PARSE
int MMCheck(MM_typecode mcode) {
	if (!mm_is_matrix(mcode)){  //consistency checks among flags in @mcode
		ERRPRINT("invalid matrix: not a matrix\n");
		return EXIT_FAILURE;
	}
	if (mm_is_dense(mcode) ){   //|| mm_is_array(mcode) ){
		ERRPRINT("invalid matrix: not a supported sparse matrix\tDENSE MAT\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

entry* MMtoCOO(ulong* NZ, FILE *fp, MM_typecode mcode,ulong* rowLens){
	int scanndRet=0;
	ulong nzTrgt=*NZ,nzIdx=0;	//expanded num of nz (in case of sym matrix)
	ulong diagEntries=0, row = 0, col = 0;//current entry's row,col from MM -> 1 based
	double val = 0;
	entry* entries = NULL;		//COO parsed entries
	///init
	if (mm_is_symmetric(mcode)){
		nzTrgt = 2* (*NZ); //upscale max num of nz in the matrix
		VERBOSE	 printf("MMtoCOO:\tparsing a simmetric matrix\n");
	}
	if (!(entries	 = malloc(nzTrgt * sizeof(*entries)))){
		ERRPRINT("MMtoCOO:  entries malloc errd\n");
		return NULL;
	}
	///parse MM fp lines into COOordinate entries
	while (1) { // Reading the fp until EOF
		if (mm_is_pattern(mcode)){
			scanndRet = fscanf(fp, "%lu %lu\n", &row, &col);
			val = 1.0;
		} else if (mm_is_real(mcode) || (mm_is_integer(mcode))){
			scanndRet = fscanf(fp, "%lu %lu %lf\n", &row, &col, &val);
		}
		

		if (scanndRet == EOF){ //TODO more strict check with type&ret?
			if (ferror(fp)){
				perror("fscanf EOF");
				goto _err; 
			} else  break;
		}
		CONSISTENCY_CHECKS{ 
			//TODO USELESS ? ? ?
			if ((mm_is_pattern(mcode) && scanndRet != 2) || 
				(!mm_is_pattern(mcode) && scanndRet != 3)){
				ERRPRINT("invalid matrix: not consistent entry scannable\n");
				goto _err; 
			}
		}
		////ADD THE CURRENT MATRIX ENTRY
		rowLens[row-1]++;
		entries[nzIdx++]=(entry) { .row=row-1, .col=col-1, .val=val };
		//also mirrored entry if sym.matrix with reflected idx inside matrix limits
		if (mm_is_symmetric(mcode) && row != col ){
			//TODO COSTRAINED FORMAT ?&& row <= mat->N && col <= mat->M ){
			SWAP(row,col);
			rowLens[row-1]++;
			entries[nzIdx++]=(entry) { .row=row-1, .col=col-1, .val=val };
		}
		else	diagEntries++;	//for CONSISTENCY_CHECKS only
	}
 
	//CONSISTENCY_CHECKS 
	nzTrgt = *NZ;
	if(mm_is_symmetric(mcode))	nzTrgt = 2*(*NZ) - diagEntries;
	assert( nzIdx == nzTrgt );
	
	//update NZ 
	*NZ = nzIdx;
	return entries;
	
	_err:
	free(entries);  return NULL;
} 

void freeMatrixMarket(MatrixMarket* mm){
	if (!mm)	return;
	free(mm->entries);
	free(mm->rowLens);
	free(mm);
}
MatrixMarket* MMRead(char* matPath){
	FILE* fp = fopen(matPath, "r");
	if (!fp){
		perror("fopen");
		return NULL;
	}
	MatrixMarket* out = calloc(1,sizeof(*out));
	if (!out){
		ERRPRINT("MMRead out malloc errd\n");
		goto err;
	}
	//banner -> parse  matrix specs
	if (mm_read_banner(fp, &out->mcode) != 0) {
		fprintf(stderr,"mm_read_banner err at:%s\n",matPath);
		goto err;
	}
	//assert matrix is compatible with this app scope
	if (MMCheck(out->mcode))	 goto err;

	//parse sizes
	//TODO OVERCOME uint limitation?
	if(mm_read_mtx_crd_size(fp,(uint*) &out->M, (uint*) &out->N, (uint*) &out->NZ)){
		fprintf(stderr,"mm_read_mtx_crd_size err at %s:\n",matPath);
		goto err;
	}
	if (!(out->rowLens = calloc(out->M,sizeof(*(out->rowLens))))){
		ERRPRINT("MMRead:\trowLens calloc errd\n");
		goto err;
	}
	if (!(out->entries = MMtoCOO(&out->NZ, fp, out->mcode,out->rowLens))){
		ERRPRINTS("MAT PARSE TO CSR ERR at:%s\n",matPath);
		goto err;
	}
	goto _end;

	err:
	freeMatrixMarket(out);
	out = NULL;
	_end:
	fclose(fp);
	return out;
}


////COO -> ANYTHING ELSE CONVERSION
int COOtoCSR(entry* entries, spmat* mat,ulong* rowLens){ 
	int out = EXIT_FAILURE;
	ulong idx;
	long* _rowsLastCol = NULL;	//for each row -> last added entry's columnIdx 
	ulong* rowsNextIdx = NULL;	//for each row -> next entry progressive idx
	if (!(rowsNextIdx = calloc(mat->M,sizeof(*rowsNextIdx)))){
		ERRPRINT("MMtoCOO:  rowsNextIdx calloc errd\n");
		goto _end;
	}
	CONSISTENCY_CHECKS{ //alloc and init aux arr for entries sort check
		if (!(_rowsLastCol = malloc(mat->M*sizeof(*_rowsLastCol)))){
			ERRPRINT("MMtoCOO:  _rowsLastCol malloc errd\n");
			goto _end;
		}
		memset(_rowsLastCol,-1,mat->M*sizeof(*_rowsLastCol));
	}
	/*TODO OLD
	 * //get rowLens->IRP (partial), TODO moved MMtoCOO to avoid FULL rescan entries
	 * for (ulong i=0; i<mat->NZ; i++)	 mat->IRP[entries[i].row+1]++;
	 * memcpy(mat->RL,mat->IRP + 1,sizeof(*mat->IRP) * mat->M); //TODO in next ifdef
	 * for (ulong i=2; i<mat->M+1; i++)	mat->IRP[i] += mat->IRP[i-1];
	 * OLD2: rowLens memcpy ... no just moved the pointer
	 * #ifdef ROWLENS
	 * memcpy(mat->RL,rowLens,sizeof(*rowLens) * mat->M); //TODO in next ifdef
	 * #endif
	 */
	//IRP: trasform rows lens as increments to build row index "pointer"
	//0th -> 0 mandatory; 1th = 0th row len, ...., M+1th = end of Mth row
	memcpy(mat->IRP+1,rowLens,sizeof(*rowLens) * mat->M);//init IRP with rows lens
	for (ulong i=2; i<mat->M+1; i++)	mat->IRP[i] += mat->IRP[i-1];
	CONSISTENCY_CHECKS  assert(mat->IRP[mat->M] == mat->NZ);
	///FILL 
	//TODO EXPECTED entries with .col entries -> CONSISTENCY_CHECKS
	//sorted for each row (even nn sequential in @entries)
	//entries write in CSR format
	entry* e;
	for (ulong i=0; i<mat->NZ; i++) {
		e = entries+i;
		CONSISTENCY_CHECKS{ //TODO CHECK IF COO ENTRIES ARE SORTED
			/*#pragma message("COO sorting check enabled")*/
			if (_rowsLastCol[e->row] >= (long) e->col){
				ERRPRINTS("not sorted entry:%ld,%ld,%lf",e->row,e->col,e->val);
				goto _end;
			}
			_rowsLastCol[e->row] = e->col;
		}
		idx = mat -> IRP[e->row] + rowsNextIdx[e->row]++;
		mat -> AS[idx] = e->val;
		mat -> JA[idx] = e->col;
	}
	
	out = EXIT_SUCCESS;

	_end:
	if(rowsNextIdx)	 free(rowsNextIdx);
	if(_rowsLastCol)	free(_rowsLastCol);

	return out;
}

int COOtoELL(entry* entries, spmat* mat, ulong* rowLens){
	int out=EXIT_FAILURE;
	ulong maxRow = 0, col, _ellEntriesTot, *rowsNextCol;
	long* _rowsLastCol=NULL;
	entry* e;
	for (ulong i=0; i<mat->M; i++)  maxRow = MAX(maxRow,rowLens[i]); 
	_ellEntriesTot = 2*mat->M*maxRow;
	#ifdef LIMIT_ELL_SIZE
	if ( _ellEntriesTot > ELL_MAX_ENTRIES ){
		ERRPRINTS("Required entries %lu -> %lu uMB for the matrix exceed the "
		  "designated threashold of: %lu  -> %lu MB for ellpack\n",
		  _ellEntriesTot,(sizeof(double)*_ellEntriesTot) >> 20,
		  ELL_MAX_ENTRIES,(sizeof(double)*ELL_MAX_ENTRIES) >> 20);
		return EXIT_FAILURE;
	}
	#endif
	//malloc aux vects
	if (!(rowsNextCol = calloc(mat->M,sizeof(*rowsNextCol)))){
		ERRPRINT("MMtoELL:\trowsNextCol calloc errd\n");
		goto _end;
	}
	CONSISTENCY_CHECKS{ //alloc and init aux arr for entries SORT CHECK
		if (!(_rowsLastCol = malloc(mat->M*sizeof(*_rowsLastCol)))){
			ERRPRINT("MMtoELL:\trowsLastCol malloc errd\n");
			goto _end;
		}
		memset(_rowsLastCol,-1,mat->M*sizeof(*_rowsLastCol));
	}
	///malloc dependant to MAX ROW LEN, err free in the caller
	if (!(mat->AS = calloc(mat->M * maxRow,  sizeof(*(mat->AS))))){
		ERRPRINT("MMtoELL:\tELL->AS calloc errd\n");
		goto _end;
	} //zero init for auto rows residual fill with 0
	if (!(mat->JA = calloc(mat->M * maxRow,  sizeof(*(mat->JA))))){
		ERRPRINT("MMtoELL:\tELL->JA calloc errd\n");
		goto _end;
	}

	mat->MAX_ROW_NZ = maxRow;
	/*#ifdef ROWLENS
	 *memcpy(mat->RL,rowLens,sizeof(*rowLens) * mat->M); //TODO in next ifdef
	 *#endif 
	 */
	///FILL NZ
	//TODO EXPECTED entries with .col entries -> CONSISTENCY_CHECKS
	//sorted for each row (even nn sequential in @entries)
	for (ulong i=0; i<mat->NZ; i++){
		e = entries + i;
		CONSISTENCY_CHECKS{ //TODO CHECK IF COO ENTRIES ARE COLS SORTED for righe successive
			/*#pragma message("COO sorting check enabled")*/
			if (_rowsLastCol[e->row] >= (long) e->col){
				ERRPRINTS("not sorted entry:%ld,%ld,%lf",
					  e->row,e->col,e->val);
				goto _end;
			}
			_rowsLastCol[e->row] = e->col;
		}
		col = rowsNextCol[e->row]++;	//place entry in its row's sequent spot
		mat->AS[ IDX2D(e->row,col,maxRow) ] = e->val;
		mat->JA[ IDX2D(e->row,col,maxRow) ] = e->col;
		 
	}
	///FILL PAD
	ulong padded = 0,paddedEntries = mat->M*mat->MAX_ROW_NZ;
	for (ulong r=0; r<mat->M; r++){
		for (ulong c=rowLens[r],j=IDX2D(r,c,maxRow); c<maxRow; c++,j++,padded++){
			//mat->AS[j] = ELL_AS_FILLER; //TODO ALREADY DONE IN CALLOC
			//mat->JA[j] = mat->JA[rowLens[r]-1]; //ELL_JA_FILLER; //TODO calloc CUDA benefit?
		}
	}
	VERBOSE{ 
	  printf("padded %lu entries = %lf%% of NZ\n",padded,100*padded/(double) mat->NZ);
	  printf("ELL matrix of: %lu paddedEntries -> %lu MB of JA+AS\n",
		paddedEntries,(paddedEntries*sizeof(*(mat->JA))+paddedEntries*sizeof(*(mat->AS))) >> 20);
	}
	out = EXIT_SUCCESS;
	_end:
	if(rowsNextCol)  free(rowsNextCol);
	if(_rowsLastCol) free(_rowsLastCol);
	return out;
}
////wrapper MM -> specialized target
spmat* MMtoCSR(char* matPath){
	spmat* mat	= NULL;
	MatrixMarket* mm = MMRead(matPath);
	if (!mm){
		ERRPRINT("MMtoCSR parse err\n");
		return NULL;
	}
	if (!(mat = calloc(1,sizeof(*mat)))){
		ERRPRINT("MMtoCSR: mat struct alloc errd");
		goto err;
	}
	mat -> M = mm->M;
	mat -> N = mm->N;
	mat -> NZ= mm->NZ;
	//alloc sparse matrix components
	if (!(mat->IRP = calloc(mat->M+1,sizeof(*(mat->IRP))))){
		ERRPRINT("MMtoCSR: IRP calloc err\n");
		goto err;
	}
	////alloc core struct of CSR
	if(!(mat->JA = malloc(mat->NZ*sizeof(*(mat->JA))))){
		ERRPRINT("MMtoCSR: JA malloc err\n");
		goto err; 
	}
	if(!(mat->AS = malloc(mat->NZ*sizeof(*(mat->AS))))){
		ERRPRINT("MMtoCSR: AS malloc err\n");
		goto err;  
	}
	if (COOtoCSR(mm->entries,mat,mm->rowLens))  goto err;
	#ifdef ROWLENS
	mat->RL = mm->rowLens;
	mm->rowLens = NULL; //avoid free in @freeMatrixMarket
	#endif

	VERBOSE
	printf("MMtoCSR: %lu NZ entries-> %lu MB of AS+JA+IRP\n",mat->NZ,
	(mat->NZ*sizeof(*(mat->AS))+mat->NZ*sizeof(*(mat->JA))+(1+mat->M*sizeof(*(mat->IRP))))>>20);
	goto _free;


	err:
	if (mat)	freeSpmat(mat);
	mat = NULL;
	_free:
	freeMatrixMarket(mm);
	return mat;
}


spmat* MMtoELL(char* matPath){
	spmat* mat	= NULL;
	MatrixMarket* mm = MMRead(matPath);
	if (!mm){
		ERRPRINT("MMtoELL: parse err\n");
		return NULL;
	}
	if (!(mat = calloc(1,sizeof(*mat)))){
		ERRPRINT("MMtoELL:  mat struct alloc errd");
		goto err;
	}
	////alloc core struct of CSR
	mat -> M = mm->M;
	mat -> N = mm->N;
	mat -> NZ= mm->NZ;
	if (COOtoELL(mm->entries,mat,mm->rowLens))  goto err;
	#ifdef ROWLENS
	mat->RL = mm->rowLens;
	mm->rowLens = NULL; //avoid free in @freeMatrixMarket
	#endif

	goto _free;

	err:
	if(mat) freeSpmat(mat);
	mat = NULL;
	_free:
	freeMatrixMarket(mm);
	return mat;
}
