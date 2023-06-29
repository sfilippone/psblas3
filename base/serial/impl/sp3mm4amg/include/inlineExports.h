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

//inline export - single implmentation functions
void cleanRbNodes(rbRoot* root,rbNode* nodes,idx_t nodesNum);
int rbInsertNewKey(rbRoot *root,rbNode *node, idx_t key);
//void C_FortranShiftIdxs(spmat* outMat);
//void Fortran_C_ShiftIdxs(spmat* m);
ACC_DENSE* _initAccVectors(ulong num,ulong size);
ACC_DENSE* _initAccVectors_monoalloc(ulong num,ulong size); //TODO PERF WITH NEXT
SPMM_ACC* initSpMMAcc(ulong entriesNum, ulong accumulatorsNum);
idx_t reductionMaxSeq(idx_t* arr,idx_t arrLen);
int _allocAccDense(ACC_DENSE* v,ulong size);
int mergeRows(SPACC* rows,spmat* mat);
int mergeRowsPartitions(SPACC* rowsParts,spmat* mat,CONFIG* conf);
void _freeAccsDenseChecks(ACC_DENSE* vectors,ulong num); 
void _resetAccVect(ACC_DENSE* acc);
void _resetIdxMap(SPVECT_IDX_DENSE_MAP* acc);
void assertArrNoRepetitions(idx_t* arrSorted, idx_t arrLen);
void freeAccsDense(ACC_DENSE* vectors,ulong num);
void freeSpMMAcc(SPMM_ACC* acc);
void sparsifyDenseVect(SPMM_ACC* acc,ACC_DENSE* accV,SPACC* accSparse, ulong startColAcc);
