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
 */


/*#pragma message( "compiling Sp3MM_CSR_OMP_Symb_Generic.c with config as:" \
	STR(OFF_F) " - " STR(OUT_IDXS) " - " STR(COL_PARTS) )*/
#ifndef OFF_F
    #error generic implementation requires OFF_F defined
#endif

///setup aux macros for different signatures implementation via #if arith expr
#pragma push_macro("OUT_IDXS")
#pragma push_macro("_OUT_IDXS")
#pragma push_macro("COL_PARTS")
#pragma push_macro("_COL_PARTS")

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

//////SpMM - rowByrow
///1row->matrix->outRow
//RBTREE based
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
 */
static inline idx_t CAT4(SpMM_Row_Symb_Rbtree,OUT_IDXS,COL_PARTS,OFF_F)  
  (
   idx_t* aRowJA, idx_t aRowLen, spmat* b,rbRoot* root, rbNode* nodes
   #if _OUT_IDXS  == TRUE && !defined OUT_IDXS_RBTREE_NODES 
   ,idx_t* outIdxs
   #endif
   #if _COL_PARTS == TRUE
   ,ushort gridCols,idx_t* rowColPartsLens
   #endif
  )
{
	//Compute resulting ab's row non zero indexes and total lenght
	idx_t abRowLen = 0;		//mul.result row len,	return value
	for ( idx_t i=0,c_a,inserted; i < aRowLen; i++ ){	//for each entry in a's row
		c_a = aRowJA[i]-OFF_F;
		//gather diffrent nnz indexes in corresponding b's `c_a`th row 
		for ( idx_t j = b->IRP[c_a]-OFF_F,c_b; j < b->IRP[c_a+1]-OFF_F; j++ ){ 
			c_b = b->JA[j]-OFF_F;
			//check if c_b is nonzero index for mul.result row
			inserted 	 =  rbInsertNewKey (root, nodes+abRowLen, c_b);
			abRowLen	+=	inserted; //inserted needed just after this
			/*	LESS EFFICIENT THEN BELOW (here no memory of last colPart)
			#if _COL_PARTS == TRUE //keep track of which col partition is c_b in
			if (inserted)
			 rowColPartsLens[ matchingUnifRangeIdx(c_b, b->N, gridCols) ]++;
			#endif */
		}
	}
	#if 	_OUT_IDXS == T  && defined OUT_IDXS_RBTREE_NODES
	/* return the non zero indexes of the mul.result row
 	 * sorting inplace the nodes inserted in the rbtree */
	sortRbNode(nodes,abRowLen);
	#elif	_OUT_IDXS == T || _COL_PARTS == T
	int i=0;
	idx_t k;
	#if _COL_PARTS == T
	//colParts aux vars
	idx_t _colBlock = abRowLen / gridCols, _colBlockRem = abRowLen % gridCols;
	ushort gc=0;
	idx_t gcStartCol = unifRemShareStart(gc,_colBlock,_colBlockRem);
	idx_t gcEndCol = unifRemShareEnd(gc,_colBlock,_colBlockRem);
	#endif	//_COL_PARTS == T
  	for (struct rb_node* n = rb_first(&root->rb_root); n; n = rb_next(n)){
		k = rb_entry(n,rbNode,rb)->key;
		#if _OUT_IDXS == T
		//return the mul.result nnz index inside the rbNodes
		outIdxs[ i++ ] = k;
		#endif
		#if _COL_PARTS == T
		while (k >= gcEndCol ){	//see if the idx is in another col partition
								//	TODO also = since gcEndCol as k is 0based
			gcEndCol = 	unifRemShareEnd(gc ,_colBlock, _colBlockRem);
			gc++;
			DEBUGCHECKS{ assert( gc < gridCols ); }
		}
		rowColPartsLens[gc]++;
		#endif //_COL_PARTS == T
	}
	#endif	//_OUT_IDXS == T ... _COL_PARTS == T
	/*DEBUGCHECKS{	//TODO PRINT NNZ INDEXES FOR MANUAL (PAINFUL CHECK)
		idx_t k;
  		for (struct rb_node* n = rb_first(&root->rb_root); n; n = rb_next(n)){
			k = rb_entry(n,rbNode,rb)->key;
			printf("%d, ",k);
		}
		printf("\n");
	}*/

	return abRowLen;
}

//SPVECT_IDX_DENSE_MAP based	TODO double implementation for trick syntax folding here...
/*
 * SPVECT_IDX_DENSE_MAP based, as SpMM_Row_Symb_Rbtree but with idxMap aux idx keeping
 * CONFIG_MACROS (new)
 * IDX_RMUL_SYMB_RBTREE && ( _OUT_IDXS == T || _COL_PARTS == T ):
 * 	 (tmp) symb mult out indexes will be kept via a rbtree
 * 	 otherwise directly in the out array appending them and then sorting them
 * 	 		(potentially same n log n)
 */
static inline idx_t CAT4(SpMM_Row_Symb_IdxMap,OUT_IDXS,COL_PARTS,OFF_F)  
  (
   idx_t* aRowJA, idx_t aRowLen, spmat* b, SPVECT_IDX_DENSE_MAP* idxsMapAcc
   #if _OUT_IDXS  == TRUE 
   ,idx_t* outIdxs
   #endif
   #if ( _OUT_IDXS == TRUE && IDX_RMUL_SYMB_RBTREE == T ) || _COL_PARTS == T
   ,rbRoot* root, rbNode* nodes
   #endif	// _OUT_IDXS == TRUE
   #if _COL_PARTS == TRUE
   ,ushort gridCols,idx_t* rowColPartsLens
   #endif
  )
{
	//Compute resulting ab's row non zero indexes and total lenght
	idx_t abRowLen = 0;		//mul.result row len,	return value
	for ( idx_t i=0,c_a,inserted; i < aRowLen; i++ ){	//for each entry in a's row
		c_a = aRowJA[i]-OFF_F;
		//gather diffrent nnz indexes in corresponding b's `c_a`th row 
		for ( idx_t j = b->IRP[c_a]-OFF_F,c_b; j < b->IRP[c_a+1]-OFF_F; j++ ){ 
			c_b = b->JA[j]-OFF_F;
			//check if c_b is nonzero index for mul.result row
			inserted  = spVect_idx_in(c_b,idxsMapAcc);
			#if	_OUT_IDXS == T || _COL_PARTS == T //idxs HAS TO be accumulated
			if (inserted)
				#if IDX_RMUL_SYMB_RBTREE == T || _OUT_IDXS == F 	//add it in a RBTREE struct
				rbInsertNewKey (root, nodes+idxsMapAcc->len, c_b);
				#else						//append it, then sort
				outIdxs[idxsMapAcc->len] = c_b;
				#endif	//IDX_RMUL_SYMB_RBTREE == T //how accumulated key c_b
			#endif	//#if	_OUT_IDXS == T || _COL_PARTS == T
		}
	}
	abRowLen = idxsMapAcc->len;
	//gather idxs or their sparsity struct in output row
	#if	_OUT_IDXS == T || _COL_PARTS == T 
	idx_t j = 0,k;
	#if _COL_PARTS == T 
	//colParts aux vars
	idx_t _colBlock = abRowLen / gridCols, _colBlockRem = abRowLen % gridCols;
	ushort gc = 0;
	idx_t  gcStartCol = unifRemShareStart(gc,_colBlock,_colBlockRem);
	idx_t  gcEndCol = unifRemShareEnd(gc,_colBlock,_colBlockRem);
	#endif
	#if IDX_RMUL_SYMB_RBTREE == T || _OUT_IDXS == F	///idxs recorded in a aux rbtree
  	for (struct rb_node* n = rb_first(&root->rb_root); n; n = rb_next(n)){
		k = rb_entry(n,rbNode,rb)->key;
		#if _OUT_IDXS == T
		outIdxs[ j++ ] = k;	//record ordered key sotred from aux rbtree
		#endif
	#else				///idxs recorded in aux append array
	sort_idx_t(outIdxs,abRowLen);
  	for (; j < abRowLen; j++){
		k = outIdxs[j];		//(OSS) already ordered in outIndexes arr
	#endif	//IDX_RMUL_SYMB_RBTREE == T
		#if _COL_PARTS == T
		while (k >= gcEndCol ){	//see if the idx is in another col partition
				//	TODO also = since gcEndCol as k is 0based
			gcEndCol = unifRemShareEnd(gc ,_colBlock, _colBlockRem);
			gc++;
			DEBUGCHECKS{ assert( gc < gridCols ); }
		}
		rowColPartsLens[gc]++;
		#endif //_COL_PARTS == T
	}
	#endif	//_OUT_IDXS == T ... _COL_PARTS == T
	return abRowLen;
}

//switch among 2 row_symb_XXX implemenetation aux
/*
 * SpMM single row symbolic computation
 * select one implementation via @implID 
 * among SpMM_Row_Symb_Rbtree or SpMM_Row_Symb_IdxMap
 * args will be forwared accordingly
 */
static inline idx_t CAT4(SpMM_Row_Symb_,OUT_IDXS,COL_PARTS,OFF_F)  
  (
   ROW_MMSYM_IMPL_MODE implID, idx_t* aRowJA, idx_t aRowLen, spmat* b,
   rbRoot* root, rbNode* nodes, SPVECT_IDX_DENSE_MAP* idxsMapAcc
   #if _OUT_IDXS  == TRUE 
   ,idx_t* outIdxs
   #endif
   #if _COL_PARTS == TRUE
   ,ushort gridCols,idx_t* rowColPartsLens
   #endif
  )
{
	if (implID == RBTREE)	{
		return CAT4(SpMM_Row_Symb_Rbtree,OUT_IDXS,COL_PARTS,OFF_F)  
			(
				aRowJA,aRowLen,b,root,nodes
   				#if _OUT_IDXS  == TRUE && !defined OUT_IDXS_RBTREE_NODES 
				,outIdxs
				#endif
				#if _COL_PARTS == TRUE
				,gridCols,rowColPartsLens
				#endif
			);
	}
	else { //IDXMAP
		return CAT4(SpMM_Row_Symb_IdxMap,OUT_IDXS,COL_PARTS,OFF_F)  
			(
				aRowJA,aRowLen,b,idxsMapAcc
   				#if		 _OUT_IDXS  == TRUE 
				,outIdxs
				#endif
   				#if ( _OUT_IDXS == TRUE && IDX_RMUL_SYMB_RBTREE == T ) || _COL_PARTS == T
   				,root, nodes
				#endif
				#if _COL_PARTS == TRUE
				,gridCols, rowColPartsLens
				#endif
			);
	}
}
///SpMM row-by-row
idx_t* CAT4(SpMM_Symb_,OUT_IDXS,COL_PARTS,OFF_F) 
  (
   ROW_MMSYM_IMPL_MODE symbRowImplID, spmat* a, spmat* b
   #if _OUT_IDXS  == TRUE
   ,idx_t*** outIdxs
   #endif
   #if _COL_PARTS == TRUE
   ,ushort gridCols, idx_t** rowColPartsLens
   #endif
  )
{

	///initial allocations
	rbRoot* rbRoots = NULL;	
	rbNode* rbNodes	= NULL; 
	SPVECT_IDX_DENSE_MAP* idxsMapAccs = NULL;
	idx_t  *rowLens=NULL,*upperBoundedRowsLens=NULL,*upperBoundedSymMat=NULL;
	idx_t  maxRowLen=0;
	int rbTreeUsed = (symbRowImplID == RBTREE || 
			  (IDX_RMUL_SYMB_RBTREE && (_COL_PARTS || _OUT_IDXS)) );
	
	if ( !(rowLens = malloc(sizeof(*rowLens) * (a->M+1))) ){ 
		ERRPRINT("SpMM_Symb_ rowLens malloc errd\n");
		goto _err;
	}
	if (_OUT_IDXS == TRUE || rbTreeUsed ){
		if (!(upperBoundedRowsLens = CAT(spMMSizeUpperbound_,OFF_F)(a,b)))
			goto _err;
	}
	#if _OUT_IDXS  == TRUE
	if (!(*outIdxs = malloc(sizeof(**outIdxs) * a->M))){
		ERRPRINT("SpMM_Symb_ outIdxs malloc errd\n");
		goto _err;
	}
	if (!(upperBoundedSymMat   = malloc(
	  sizeof(*upperBoundedSymMat)*upperBoundedRowsLens[a->M]))){
		ERRPRINT("SpMM_Symb_ upperBoundedSymMat malloc errd\n");
		goto _err;
	}
	//write rows' start pointer from full matrix JA allocated
	for (idx_t i=0,cumul=0; i<a->M; cumul += upperBoundedRowsLens[i++])
		*outIdxs[i] = upperBoundedSymMat + cumul; 
	#endif	//#if _OUT_IDXS  == TRUE
	#if _COL_PARTS == TRUE
	if (!(*rowColPartsLens = malloc(a->M * gridCols * sizeof(**rowColPartsLens)))){
		ERRPRINT("SpMM_Symb_ rowColPartsLens malloc errd\n");
		goto _err;
	}
	#endif //_COL_PARTS
	int maxThreads;
	maxThreads = omp_get_max_threads();	//TODO FROM CFG
	//index keeping aux struct
	//rbtree implementation or idxMap with aux of symbTree for tmp outIdx keeping
	if ( rbTreeUsed ){
		maxRowLen = reductionMaxSeq(upperBoundedRowsLens, a->M);
		//rbTrees for index keeping
		rbRoots = malloc(maxThreads * sizeof(*rbRoots));
		rbNodes	= calloc(maxThreads * maxRowLen, sizeof(*rbNodes));
		if( !rbRoots || !rbNodes ){
			ERRPRINT("SpMM_Symb_ threads' aux rbTree mallocs errd\n");
			goto _err;
		}
		//init roots
		for (int i=0; i<maxThreads; i++)
			rbRoots[i] = RB_ROOT_CACHED;
	} 
	if (symbRowImplID == IDXMAP){  //idxMap implementation && not outIdxs via rbtree
		idxsMapAccs 	= calloc(maxThreads, sizeof(*idxsMapAccs));
		if (!idxsMapAccs){
			ERRPRINT("SpMM_Symb_ idxsMapAccs calloc errd\n");
			goto _err;
		}
		//init idxs maps
		for (int i=0; i<maxThreads; i++){
			if(initSpVectIdxDenseAcc(b->N,idxsMapAccs+i))	goto _err;
		}
	}
	///rows parallel compute
	idx_t* aRow;
	idx_t  aRowLen,rLen,abCumulLen=0;
	int tid;
	rbRoot* tRoot;	rbNode* tNodes; 
	SPVECT_IDX_DENSE_MAP* tIdxsMapAcc = NULL;
	#pragma omp parallel for schedule(static) \
	  private(aRow,aRowLen,rLen, tRoot,tNodes,tid) reduction(+:abCumulLen)
	for(idx_t r=0; r<a->M; r++){
		aRow		= a->JA + a->IRP[r]-OFF_F;
		aRowLen		= a->IRP[r+1] - a->IRP[r];
		tid 		= omp_get_thread_num();
		//TODO low overhead pointer airth can be avoided with if (symbRowImplID .. && )
		tIdxsMapAcc	= idxsMapAccs + tid;
		tRoot 		= rbRoots + tid;
		tNodes  	= rbNodes + tid * maxRowLen;

		rLen = CAT4(SpMM_Row_Symb_,OUT_IDXS,COL_PARTS,OFF_F)  
		(
			symbRowImplID, aRow, aRowLen, b, tRoot,tNodes, tIdxsMapAcc
   			#if _OUT_IDXS  == TRUE
			,*outIdxs[r]
			#endif
			#if _COL_PARTS == TRUE
			,gridCols, (*rowColPartsLens) + IDX2D(r,0,gridCols)
			#endif
		);
		rowLens[r]  = rLen;
		abCumulLen += rLen;
		///reset symb idxs keeping aux structs
		if (symbRowImplID==RBTREE || (IDX_RMUL_SYMB_RBTREE && (_COL_PARTS || _OUT_IDXS))){
			*tRoot = RB_ROOT_CACHED;
			memset(tNodes,0,rLen * sizeof(*tNodes));
		}
		if (symbRowImplID == IDXMAP)	_resetIdxMap(tIdxsMapAcc);
		
	}
	rowLens[a->M]	= abCumulLen;
	goto _free;

	_err:
	free(rowLens);
	#if _OUT_IDXS  == T
	if  (outIdxs)			free(*outIdxs);
	#endif
	#if _COL_PARTS == T
	if  (rowColPartsLens)	free(*rowColPartsLens);
	#endif
	rowLens = NULL;
	_free:
	free(upperBoundedRowsLens);
	free(upperBoundedSymMat);
	free(rbRoots);
	free(rbNodes);
	if (idxsMapAccs){
		for (int i=0; i<maxThreads; i++)	free(idxsMapAccs[i].idxsMap);
		free(idxsMapAccs);
	}

	return rowLens;
}


///Sp3MM - rowByrowByrow	////////////////////////////////////////////////////////
//TODO COL_PARTS VERSION?
#if _OUT_IDXS == TRUE && _COL_PARTS == FALSE
idx_t CAT3(Sp3MM_Row_Symb_,OUT_IDXS,OFF_F) 
  (
	ROW_MMSYM_IMPL_MODE symbMMRowImplID,idx_t* aRowJA,idx_t aRowLen,spmat* b,spmat* c, 
	rbRoot* root,rbNode* nodes, SPVECT_IDX_DENSE_MAP* idxsMapAcc,idx_t* abRowJATmp
	#if _OUT_IDXS  == TRUE && !defined  OUT_IDXS_RBTREE_NODES 
	,idx_t* outIdxs
	#endif
  )
{
	struct rb_node* n;
	idx_t abRowLen	= CAT4(SpMM_Row_Symb_,OUT_IDXS,COL_PARTS,OFF_F) 
	  (symbMMRowImplID,aRowJA,aRowLen,b,root,nodes,idxsMapAcc,abRowJATmp);
	cleanRbNodes(root, nodes, abRowLen);
	idx_t abcRowLen	= CAT4(SpMM_Row_Symb_,OUT_IDXS,COL_PARTS,OFF_F)
	  (symbMMRowImplID,abRowJATmp,abRowLen,c,root,nodes,idxsMapAcc,abRowJATmp);

	#if 	_OUT_IDXS == TRUE
	#ifndef OUT_IDXS_RBTREE_NODES	
	//return the mul.result nnz index inside the rbNodes
	int i=0;
	//rbNodeOrderedVisit(n,root)	outIdxs[ i++ ] = rb_entry(n,rbNode,rb)->key;
  	for (struct rb_node* n = rb_first(&root->rb_root); n; n = rb_next(n)){
		outIdxs[ i++ ] = rb_entry(n,rbNode,rb)->key;
	}
	#else
	/* return the non zero indexes of the mul.result row
 	 * sorting inplace the nodes inserted in the rbtree */
	sortRbNode(nodes,abcRowLen);
	#endif
	#endif
	return abcRowLen;
}

idx_t* CAT3(Sp3MM_Symb_,OUT_IDXS,OFF_F) 
  (
	ROW_MMSYM_IMPL_MODE symbMMRowImplID, spmat* a, spmat* b, spmat* c
	#if _OUT_IDXS  == TRUE
	,idx_t*** outIdxs
	#endif
  )
{
	//idxs keeping aux buffs
	idx_t*  abRowsJATmp	= NULL;
	rbRoot* rbRoots 	= NULL; rbNode* rbNodes	= NULL; 
	SPVECT_IDX_DENSE_MAP* idxsMapAccs = NULL;
	
	idx_t *abUpperBoundedRowsLens = NULL, *upperBoundedSymMat = NULL;
	///initial allocations
	idx_t* rowLens = malloc(sizeof(*rowLens) * (a->M +1) ); //to return
	if (!rowLens){
		ERRPRINT("SpMM_Symb_ rowLens malloc errd\n");
		goto _err;
	}
	
	#if _OUT_IDXS  == TRUE
	if (!(*outIdxs = malloc(sizeof(**outIdxs) * a->M))){
		ERRPRINT("SpMM_Symb_ outIdxs malloc errd\n");
		goto _err;
	}
	if (!(abUpperBoundedRowsLens = CAT(spMMSizeUpperbound_,OFF_F)(a,b)))	
		goto _err;
  	/*TODO TODO instead of doing one sym product first to have a correct UB
	 *	use an heuristics here to get output matrix size
	 */
	idx_t abcUBSize = abUpperBoundedRowsLens[a->M] * SP3MM_UB_HEURISTIC;
	if (!(upperBoundedSymMat=malloc(sizeof(*upperBoundedSymMat)*abcUBSize))){
		ERRPRINT("SpMM_Symb_ upperBoundedSymMat malloc errd\n");
		goto _err;
	}
	//TODO heuristic TO UB rows bounds ... require compacting copy
	for (idx_t i=0,cumul=0; i<a->M; 
	  cumul += SP3MM_UB_HEURISTIC * abUpperBoundedRowsLens[i++])
		*outIdxs[i] = upperBoundedSymMat + cumul; 
  	#endif	//#if _OUT_IDXS  == TRUE
	//rbTrees for index keeping
	int maxThreads 	= omp_get_max_threads(); //TODO FROM CFG
	idx_t abMaxRowLen 	= reductionMaxSeq(abUpperBoundedRowsLens, a->M);
	#ifdef 	  HEURISTICS_UB
	idx_t maxRowLenUB 	= abMaxRowLen * SP3MM_UB_HEURISTIC; //TODO UB HEURISTC
	#else
	idx_t maxRowLenUB 	= c->N;
	#endif	//HEURISTICS_UB

	if (!(abRowsJATmp	= malloc(maxThreads*maxRowLenUB*sizeof(*abRowsJATmp)) ) ){
		ERRPRINT("Sp3MM_Symb_ abRowsJATmp malloc errd\n");
		goto _err;
	}
	if (symbMMRowImplID == RBTREE || IDX_RMUL_SYMB_RBTREE ){
		rbNodes	= malloc(maxThreads * maxRowLenUB * sizeof(*rbNodes));
		rbRoots = malloc(maxThreads * sizeof(*rbRoots));
		if (!rbRoots || !rbNodes ){
			ERRPRINT("Sp3MM_Symb_ rbRoots || rbNodes malloc errd\n");
			goto _err;
		}
		//init roots
		for (int i=0; i<maxThreads; i++)	rbRoots[i] = RB_ROOT_CACHED;
		if (!( rbNodes	= malloc(maxThreads * maxRowLenUB * sizeof(*rbNodes)) )){
			ERRPRINT("Sp3MM_Symb_ rbNodes malloc errdn\n");
			goto _err;
		}
	}
	if (symbMMRowImplID == IDXMAP){
		if (!( idxsMapAccs = calloc(maxThreads,sizeof(*idxsMapAccs)) )){
			ERRPRINT("Sp3MM_Symb_ idxsMapAccs calloc errd\n");
			goto _err;
		}	
		//init idxs maps
		for (int i=0; i<maxThreads; i++){
			if(initSpVectIdxDenseAcc(c->N,idxsMapAccs+i))	goto _err;
		}
	}
	///rows parallel compute
	idx_t* aRow;
	idx_t  aRowLen,rLen,outCumulLen=0;
	//threads local pointers
	int tid;
	rbRoot* tRoot;
	rbNode* tNodes;
	SPVECT_IDX_DENSE_MAP* tIdxsMapAcc;
	idx_t* tABRowJATmp;

	#pragma omp parallel for schedule(static) \
	private(aRow,aRowLen,rLen, tRoot,tNodes,tid) reduction(+:outCumulLen)
	for(idx_t r=0; r<a->M; r++){
		aRow		= a->JA + a->IRP[r]-OFF_F;
		aRowLen		= a->IRP[r+1] - a->IRP[r];
		tid 		= omp_get_thread_num();
		tRoot 		= rbRoots + tid;
		tNodes  	= rbNodes + tid * maxRowLenUB;
		tIdxsMapAcc = NULL;	//TODO
		tABRowJATmp	= abRowsJATmp + tid * maxRowLenUB;
		rLen = CAT3(Sp3MM_Row_Symb_,OUT_IDXS,OFF_F)
		  (symbMMRowImplID, aRow,aRowLen,b,c,tRoot,tNodes,tIdxsMapAcc,tABRowJATmp,
		  	#if _OUT_IDXS == TRUE
			*outIdxs[r]
			#endif
		  );
		outCumulLen += rLen;
		rowLens[r] = rLen;
		
	}
	goto _free;
	_err:
	free(rowLens);
	#if _OUT_IDXS == T
	if (*outIdxs)	free(*outIdxs);
	#endif
	rowLens = NULL;
	_free:
	free(abUpperBoundedRowsLens);
	free(upperBoundedSymMat);
	free(rbRoots);
	free(rbNodes);
	free(abRowsJATmp);

	return rowLens;
}

#endif	//#if !defined COL_PARTS && defined OUT_IDXS


///restore aux macros entry state
//#undef _OUT_ID
//#undef _COL_PARTS
#pragma pop_macro("OUT_IDXS")
#pragma pop_macro("_OUT_IDXS")
#pragma pop_macro("COL_PARTS")
#pragma pop_macro("_COL_PARTS")




