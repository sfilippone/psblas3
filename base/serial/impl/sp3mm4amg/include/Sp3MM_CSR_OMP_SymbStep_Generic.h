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
/*
 * CSR Sp[3]MM Symbolic step implementations
 * target: 	compute the output matrix size and the row lens for preallocation
 * 			direct write out partial results
 * See interfaces in respective header
 *
 * MultiImplementations functions with all parameters
 * several functions (lower level) has config based on 
 * 	OUT_IDXS for returing or using also symbMul Output idxs somehow
 * 	COL_PARTS for returing or using also distribution of symb out idxs in col partition
 */

///MultiImplementations
///setup aux macros for different signatures implementation via #if arith expr
#ifdef OUT_IDXS 
	#define _OUT_IDXS  	TRUE
#else   
	#define _OUT_IDXS 	FALSE
	#define OUT_IDXS	_UNDEF
#endif
#ifdef COL_PARTS
	#define _COL_PARTS	TRUE
#else   
	#define _COL_PARTS	FALSE
	#define COL_PARTS	_UNDEF
#endif
///
#if !defined OFF_F 
	#error generic implementations requires OFF_F 
#endif




/*
 * Compute symbolic product of (nnz indexes of) row @aRowJA and matrix @b
 * insert nnz indexes of the mul. result row as nodes in a rbtree rooted at @root
 * with nodes in @nodes which have to be enough for the mul result row (use an UB)
 * Retuns: multiplication result row NNZ number,se CONFIG_MACROS below for more 
 *
 * CONFIG_MACROS:
 * if _OUT_IDXS == TRUE	return mul.result row nnz idxs in @outIdxs
 *	  ifdef: OUT_IDXS_RBTREE_NODES: nnz indexes returned inplace sorting rbtree
 * 								as nnz indexes(JA) of the mul result row
 * else:							stop at returning the mul. result row lenght
 * if _COL_PARTS == TRUE	return the number of nonzero elements in
 *		in each of the @gridCols column partitions inside @rowColPartsLens
 * OFF_F:	offset back indexes from fortran
 * 		TODO also output indexes are shifted (see c_b )

idx_t CAT4(SpMM_Row_Symb_Rbtree,OUT_IDXS,COL_PARTS,OFF_F) (idx_t* aRowJA,idx_t aRowLen,
  spmat* b,rbRoot* root,rbNode* nodes
	#if _OUT_IDXS  == TRUE
	#ifndef OUT_IDXS_RBTREE_NODES 
	,idx_t* outIdxs
	#endif
	#endif
	#if _COL_PARTS == TRUE
	,ushort gridCols,idx_t* rowColPartsOffsets
	#endif
  );
 */

/*
 * SPVECT_IDX_DENSE_MAP based, as SpMM_Row_Symb_Rbtree but with idxMap aux idx keeping
 * CONFIG_MACROS (new)
 * IDX_RMUL_SYMB_RBTREE && ( _OUT_IDXS == T || _COL_PARTS == T ):
 * 	 (tmp) symb mult out indexes will be kept via a rbtree
 * 	 otherwise directly in the out array appending them and then sorting them
 * 	 		(potentially same n log n)

idx_t CAT4(SpMM_Row_Symb_IdxMap,OUT_IDXS,COL_PARTS,OFF_F)  
  (
   idx_t* aRowJA, idx_t aRowLen, spmat* b, SPVECT_IDX_DENSE_MAP* idxsMapAcc
   #if _OUT_IDXS  == TRUE 
   ,idx_t* outIdxs
	#if IDX_RMUL_SYMB_RBTREE == TRUE
	,rbRoot* root, rbNode* nodes
	#endif
   #endif	// _OUT_IDXS == TRUE
   #if _COL_PARTS == TRUE
   ,ushort gridCols,idx_t* rowColPartsLens
   #endif
  );
 */

/*
 * Compute symbolic product of sparse matrixes @a * @b
 * Alloc aux structures based on a upper bounded allocation
 * Returns array of exact row lens of the result mul matrix c=@a*@b
 * 		plus an extra entry for the cumulative num of matrix nnz 
 * 		(interoperability with upper bound implementations)
 *
 * CONFIG_MACROS:
 * 	if _OUT_IDXS == TRUE:	
 * 		return in *@outIdxs pointers to overallocd JAs nnz indexes
 * 		*outIdxs[0] --> start of first row (its size in 0th of returned array)
 * 						also is the malloc returned address for later free
 * 							[NB: rows are not contiguos in the allocation]
 * 	if _COL_PARTS == TRUE: 	return in *@rowColPartsLens
 * 		a matrix @a->M * @gridCols of offsets 
 * 		for each column partition in each row
 */
idx_t* CAT4(SpMM_Symb_,OUT_IDXS,COL_PARTS,OFF_F) 
  (
   ROW_MMSYM_IMPL_MODE symbRowImplID, spmat* a, spmat* b
   #if _OUT_IDXS  == TRUE
   ,idx_t*** outIdxs
   #endif
   #if _COL_PARTS == TRUE
   ,ushort gridCols, idx_t** rowColPartsLens
   #endif
  );

//////Sp3MM - rowByrowByrow
#define SP3MM_UB_HEURISTIC	2
#if !defined COL_PARTS && defined OUT_IDXS
///MultiImplementations functions without COL_PARTS
//NB required OUT_IDXS for initial row-by-row step...
/*
 * as earlier but meant for work with 
 * triple product as rob-by-row-by-row forwarding:
 * @abRowJATmp is a tmp storage for first row-by-row sym product
 * 	has to be big enough to store all the possible nonzero elements
 * @nodes has to be big enough to store all possible nnz of ab and abc row
 * 	(not only ab row as earlier version)
 *
 */
idx_t CAT3(Sp3MM_Row_Symb_,OUT_IDXS,OFF_F) 
  (
	ROW_MMSYM_IMPL_MODE symbMMRowImplID, idx_t* aRowJA,idx_t aRowLen,
	spmat* b,spmat* c, rbRoot* root,rbNode* nodes, idx_t* abRowJATmp
	#if _OUT_IDXS  == TRUE
	#ifndef OUT_IDXS_RBTREE_NODES 
	,idx_t* outIdxs
	#endif
	#endif
  );

/*
 * triple product @a*@b*@c as rob-by-row-by-row forwarding:
 * Returns: resulting matrix rows' sizes exactly in an array
 * 			plus an extra element for the res.matrix's total size
 * CONFIG_MACROS:
 * 	if _OUT_IDXS == TRUE:	
 * 		return in *@outIdxs pointers to overallocd JAs nnz indexes
 * 		*outIdxs[0] --> start of first row (its size in 0th of returned array)
 * 						also is the malloc returned address for later free
 * 							[NB: rows are not contiguos in the allocation]
 */ 
//TODO 	HEURISTICS IN _OUT_IDXS to avoid serializing after first sym.product
//		see HEURISTICS_UB
idx_t* CAT3(Sp3MM_Symb_,OUT_IDXS,OFF_F) 
  (
	ROW_MMSYM_IMPL_MODE symbMMRowImplID, spmat* a, spmat* b, spmat* c
	#if _OUT_IDXS  == TRUE
	,idx_t*** outIdxs
	#endif
  );

#endif
	
#undef OUT_IDXS
#undef _OUT_IDXS
#undef COL_PARTS
#undef _COL_PARTS
