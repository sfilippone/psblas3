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

#ifndef PARSER
#define PARSER

#include "mmio.h" 
#include "sparseMatrix.h" 

typedef struct{
	ulong row;
	ulong col;
	double val;
} entry;	 //MatrixMarket COO entry

typedef struct{
	MM_typecode mcode;
	entry* entries;
	ulong* rowLens;
	ulong M,N,NZ;   //spmat sizes
} MatrixMarket;

////COO PARSE
//parse and check MatrixMarket mat in @matPath file
MatrixMarket* MMRead(char* matPath);
void freeMatrixMarket(MatrixMarket* mm);
//basic check for sparse matrix compliance to the app, return posix bool
int MMCheck(MM_typecode typecode);
/* 
 * parse MatrixMarket matrix entries in @fp, of type @mcode
 * into COOrdinate list of entries
 *  -> expand simmetric matrixes into a normal matrix with both parts
 *	  so NZ will be inplace doubled
 * return allocated and filled COO entries with the NNZ number into
 * NO SORT CHECKING HERE
 */
entry* MMtoCOO(ulong* NZ, FILE *fp, MM_typecode mcode,ulong* rowLens);

////COO -> ANYTHING ELSE CONVERSION
/*
 * write COO entries in @entries inside sparse matrix @mat in ELL format
 * EXPECTED: CSR arrays allocated, @entries col sorted in (not madatory consecut) rows
 * [simmetrical parts explicitly rappresented --> not important here]
 */
int COOtoCSR(entry* entries, spmat* mat,ulong* rowLens);
/*
 * write COO entries in @entries inside sparse matrix @mat in ELL format
 * EXPECTED: @entries col sorted in (not madatory consecut) rows
 * ELL internal array allocated in this function, not freed in case of error
 */
int COOtoELL(entry* entries, spmat* mat, ulong* rowLens);
////wrapper MM -> specialized target
/*
 * Parse MatrixMarket matrix stored in file at @matPath
 * IMPLEMENTED WRAPPING: MMtoCOO -> COOtoCSR
 * Returns: allocated spmat sparse matrix with all field allocated
 * symetric matrix are expanded in a full matrix
 */
spmat* MMtoCSR(char* matPath);
spmat* MMtoELL(char* matPath);


#endif
