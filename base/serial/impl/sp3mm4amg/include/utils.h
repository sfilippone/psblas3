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
#ifndef UTILS
#define UTILS

//#include <stddef.h> 
#include "macros.h"

#include "linuxK_rbtree_minimalized.h"
extern int urndFd;	//file pointer to DRNG_DEVFILE O_RDONLY opened
int init_urndfd(); // wrap init urndFd

///IO
//UNBUFFERED IO
/*
 * urndFd usage template to populate random timeout
	if(_read_wrap(urndFd,(char*)&timeout,sizeof(timeout))<0){
		fprintf(stderr,"rnd read for thread's timeout failed");
		ret=EXIT_FAILURE;
		goto end;
	}
 */
//wrap read cycle over @fd
int read_wrap(int fd,void* dst,size_t count);
int fread_wrap(FILE* fp,void* dst,size_t count);
//dual of read_wrap
int write_wrap(int fd,void* src,size_t count);
//create or open file at @outFpath for write
int createNewFile(char* const outFpath);
///STRUCTURED DATA IO
#define DOUBLE_STR_FORMAT	"%25le\n"
//write double vector @v as row sequence of double at @fpath
//e.g. read with od -tf8 -w8 fpath : OCTALOFFSET:   DOUBLE FULL DIGITS
int writeDoubleVector(char* fpath,double* v,ulong size);
/*
 * read vector of double [Str] of arbitrary size from @fpath, true lenght in *size
 * if size point to a nnz value, the initial allocation will be of *size
 * eventual successive reallocation done multipling *size with VECTOR_STEP_REALLOC
 */
double* readDoubleVector(char* fpath,ulong* size);
double* readDoubleVectorStr(char* fpath,ulong* size);

///STRUCTURED DATA IO -- BUFFERED: FSCANF - FPRINTF
//dual of readDoubleVectorVector
int writeDoubleVectorAsStr(char* fpath,double* v,ulong size);
int MPI_Dims_create(int nnodes, int ndims, int dims[]);	//commons/ompi_dims_create/ompi_dims_create.c

#include "config.h"
///config from ENV
#define GRID_ROWS   "GRID_ROWS"
#define GRID_COLS   "GRID_COLS"
//parse configuration from env
int getConfig(CONFIG* conf);

//append only list implemented with a reallocated array
typedef struct{
	ulong* a;
	ulong  size;
	ulong  lastIdx;
} APPENDARRAY;
//append @val to @list, reallocating if reached end
//TODO inline int appendArr(ulong val,APPENDARRAY* list);

void sortuint(uint* arr, uint len);     //sort uint array @arr of @len elements
void sort_idx_t(idx_t* arr, idx_t len); 
void sortulong(ulong* arr, ulong len);
void sortRbNode(rbNode* arr,idx_t len);

///ranges functions 
/*
 * considering @rangesN uniform ranges from [0,rangesEnd)
 * return the rangeId that element @i is in
 */
inline ushort matchingUnifRangeIdxLinear(idx_t i,idx_t rangesEnd,ushort rangesN){
	double rangeW = rangesEnd / (double) rangesN;
	for(idx_t rEnd = rangeW,j=0; j < rangesN; rEnd = rangeW * (j++)){
		if( i < rEnd )	return j;
	}
	assert( FALSE );	//i should be in a range!
	return EXIT_FAILURE;
}
/*
 * find which range @idx match among 
 * a uniform range divion of @size element in @rangesN ranges with carried reminder
 * return 0based idx of matched range
 */
inline ushort matchingUnifRangeIdx(idx_t idx, idx_t size, ushort rangesN){
	idx_t rangeW = size/rangesN, rangeRem = size%rangesN;
	idx_t searchStart,searchEnd;
	//!IN_RANGE(idx,unifRemShareBlock(idx,rangeW,rangeRem)) && IN_RANGE(r,0,rangesN);){
	for(ushort r=rangesN/2-1, hMoveWidth=rangesN/2; hMoveWidth>0;hMoveWidth/=2){
		searchStart = unifRemShareStart(idx,rangeW,rangeRem);
		searchEnd	= unifRemShareEnd(idx,rangeW,rangeRem);
		if 		(IN_RANGE(idx,searchStart,searchEnd))	return r;
		else if	(idx < searchStart)						r = AVG(0,r);
		else											r = AVG(r,rangesN-1);
	}
	assert(FALSE);
}

///reductionUtils
inline idx_t reductionSumSeq(idx_t* arr,idx_t arrLen){
	idx_t i,out;
	for(i=0,out=0; i<arrLen; out += arr[i++] );
	return out;
}
inline idx_t reductionMaxSeq(idx_t* arr,idx_t arrLen){
	idx_t i,out;
	for( i=0,out=0; i<arrLen; i++,out=out<arr[i]?arr[i]:out );
	return out;
}
inline idx_t reductionSumOmp(idx_t* arr,idx_t arrLen){
	idx_t out;
	#pragma omp parallel for reduction(+:out)
	for(idx_t i=0; i<arrLen; i++){
		out += arr[i];
	}
	return out;
}
inline idx_t reductionMaxOmp(idx_t* arr,idx_t arrLen){
	idx_t out;
	#pragma omp parallel for reduction(max:out)
	for(idx_t i=0; i<arrLen; i++){	out=out<arr[i]?arr[i]:out;	}
	return out;	
}

/*
 * return 0 if vectors a and b has elements that differ at most of DOUBLE_DIFF_THREASH 
 * if diffMax!=NULL save there the max difference value  
 *  of the 2 entries of the input vectors, signed as a[i] - b[i] (as well as dump prints)
 * CONVENTION:	@a = true result, @b vector to check with the true result 
 */
int doubleVectorsDiff(double* a, double* b, ulong n,double* diffMax);
//fill a random vector in @v long @size doubles
int fillRndVector(ulong size, double* v);
//read vector as a sequence of space separated double from file at @fpath 
#define VECTOR_STEP_MALLOC 100

/* 
 * decompress file at @path into @tmpFsDecompressPath, 
 * decompression command obtanined first looking at the extension
 * then matching it with a list of avaible decompression cmd
 * that can be make as shell cmd adding @path > @tmpFsDecompressPath
 * e.g. decompress xz -d -c @path > @tmpFsDecompressPath
 * Returns: -1 if decompression wasn't possible otherwise decompress command exti status
 */
int extractInTmpFS(char* path, char* tmpFsDecompressPath);
//compute E[@values] in @out[0] and VAR[@values] in @out[1] of @numVals values
void statsAvgVar(double* values,uint numVals, double* out);
void printMatrix(double* mat,ulong m,ulong n,char justNZMarkers);
void printVector(double* v,ulong size);

void assertArrNoRepetitions(idx_t* arrSorted, idx_t arrLen);
#endif
