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
//CORE IMPLEMENTATIONS HEADER
#ifndef SP3MM_CSR_OMP_MULTI_H
#define SP3MM_CSR_OMP_MULTI_H

///commons single implementation stuff
#include "macros.h"
#include "sparseMatrix.h"

///aux structures
//hold SPMM result over a unpartitionated space among threads-row[s' blocks]
typedef struct{
	//space to hold SPMM output
	ulong*  JA;
	double* AS;
	ulong   size;			//num of entries allocated -> only dbg checks
	ulong   lastAssigned;	//last JA&AS assigned index to an accumulator(atom)
	SPACC*  accs;			//SPARSIFIED ACC POINTERS
	uint	accsNum;	
} SPMM_ACC; //accumulator for SPMM
///compute function interface and its pointer definitions
typedef spmat* ( SPMM)		(spmat*,spmat*,CONFIG*);
typedef spmat* (*SPMM_INTERF)	(spmat*,spmat*,CONFIG*);
typedef spmat* ( SP3MM)		(spmat*,spmat*,spmat*,CONFIG*,SPMM_INTERF);
typedef spmat* (*SP3MM_INTERF)  (spmat*,spmat*,spmat*,CONFIG*,SPMM_INTERF);

typedef enum {
	_1D_DIRECT,
	_1D_BLOCKS,
	_2D_OFFSET,
	_2D_ALLOCD 
} SPMM_IMPL_TYPE;
///-- commons single implementation stuff

///includes
#include "linuxK_rbtree_minimalized.h"
#include "Sp3MM_CSR_OMP_SymbStep_Multi.h"

//extern char TRGT_IMPL_START_IDX; //multi implementation switch
#include "sparseUtilsMulti.h"
#ifdef OFF_F	//save "includer" OFF_F value before overwriting it
	#pragma push_macro("OFF_F")
	#define _OFF_F_OLD
	#undef  OFF_F
#endif


#define OFF_F 0
#include "Sp3MM_CSR_OMP_UB_Generic.h"
#include "Sp3MM_CSR_OMP_Num_Generic.h"
#undef OFF_F

#define OFF_F 1
#include "Sp3MM_CSR_OMP_UB_Generic.h"
#include "Sp3MM_CSR_OMP_Num_Generic.h"
#undef OFF_F



#ifdef _OFF_F_OLD
	#pragma pop_macro("OFF_F")
	#undef  _OFF_F_OLD
#endif


#endif	//SP3MM_CSR_OMP_MULTI_H 
