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

#ifndef OFF_F
    /*#pragma message("generic implementation requires OFF_F defined")*/
    #error generic implementation requires OFF_F defined
#endif

/////SYMBOLIC - NUMERIC IMPLEMENTATIONS
///SP3MM FUNCTIONS
/*
 *  triple matrix multiplication among @R * @AC * @P using gustavson parallel implementation
 *  implmented as a pair of subsequent spmm operations
 *  if @conf->spmm != NULL, it will be used as spmm function, otherwise euristics will be 
 *  used to decide wich implementation to use
 */
SP3MM CAT(sp3mmRowByRowPair_SymbNum_,OFF_F);

/*
 * row-by-row-by-row implementation: forwarding @R*@AC rth row to P for row-by-row
 * accumulation in preallocated space, TODO exactly determined
 * basic parallelization: 1thread per @R's rows that will also forward the result to P
 */
SP3MM CAT(sp3mmRowByRowMerged_SymbNum_,OFF_F);

///SUB FUNCTIONS
///SPMM FUNCTIONS
/*
 * sparse parallel implementation of @A * @B parallelizing Gustavson row-by-row
 * formulation using an aux dense vector @_auxDense
 * return resulting product matrix
 */
SPMM CAT(spmmRowByRow_SymbNum_,OFF_F);
/*
 * sparse parallel implementation of @A * @B parallelizing Gustavson 
 * with partitioning of @A in @conf->gridRows blocks of rows  
 * return resulting product matrix
 */
SPMM CAT(spmmRowByRow1DBlocks_SymbNum_,OFF_F);

/* 
 * sparse parallel implementation of @A * @B as Gustavson parallelizzed in 2D
 * with partitioning of
 * @A into rows groups, uniform rows division
 * @B into cols groups, uniform cols division, accessed by aux offsets
 */
SPMM CAT(spmmRowByRow2DBlocks_SymbNum_,OFF_F);

/* 
 * sparse parallel implementation of @A * @B as Gustavson parallelizzed in 2D
 * with partitioning of
 * @A into rows groups, uniform rows division
 * @B into cols groups, uniform cols division, ALLOCATED as CSR submatrixes
 */
SPMM CAT(spmmRowByRow2DBlocksAllocated_SymbNum_,OFF_F);

///implementation wrappers as static array of function pointers 
//sp3mm as pair of spmm
static SPMM_INTERF  CAT(Spmm_SymbNum_Funcs_,OFF_F)[] = {
	& CAT(spmmRowByRow_SymbNum_,OFF_F),
	& CAT(spmmRowByRow1DBlocks_SymbNum_,OFF_F),
	//& CAT(spmmRowByRow2DBlocks_SymbNum_,OFF_F),
	//& CAT(spmmRowByRow2DBlocksAllocated_SymbNum_,OFF_F)
};
//sp3mm as pair of spmm
static SP3MM_INTERF CAT(Sp3mm_SymbNum_Funcs_,OFF_F)[] = {
	//& CAT(sp3mmRowByRowMerged_SymbNum_,OFF_F)
};
