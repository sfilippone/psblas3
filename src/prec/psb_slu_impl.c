/*
 *                    MD2P4
 *    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
 *                              for 
 *                    Parallel Sparse BLAS  v2.0
 *
 *   (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
 *                      Alfredo Buttari
 *                      Daniela Di Serafino    II University of Naples
 *                      Pasqua D'Ambra         ICAR-CNR
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the MD2P4 group or the names of its contributors may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
/* This file is an interface to the SuperLU routines for sparse
   factorization.  It was obtaned by modifying the
   c_fortran_dgssv.c file from the SuperLU source distribution;
   original copyright terms reproduced below. 
   
   PSBLAS v 2.0  */ 


/*		=====================

Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer. 
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution. 
(3) Neither the name of Lawrence Berkeley National Laboratory, U.S. Dept. of
Energy nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
  
*/

/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

#ifdef Have_SLU_
#include "dsp_defs.h"

#define HANDLE_SIZE  8
/* kind of integer to hold a pointer.  Use int.
   This might need to be changed on 64-bit systems. */
#ifdef LargeFptr
typedef long long fptr;  /* 32-bit by default */
#else
typedef int fptr;  /* 32-bit by default */
#endif

typedef struct {
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_c;
    int *perm_r;
} factors_t;


#else

#include <stdio.h>

#endif


#ifdef Add_
#define psb_dslu_factor_ psb_dslu_factor_
#define psb_dslu_solve_  psb_dslu_solve_
#define psb_dslu_free_   psb_dslu_free_
#endif
#ifdef AddDouble_
#define psb_dslu_factor_ psb_dslu_factor__
#define psb_dslu_solve_  psb_dslu_solve__
#define psb_dslu_free_   psb_dslu_free__
#endif
#ifdef NoChange
#define psb_dslu_factor_ psb_dslu_factor
#define psb_dslu_solve_  psb_dslu_solve
#define psb_dslu_free_   psb_dslu_free
#endif




void
psb_dslu_factor_(int *n, int *nnz,
                 double *values, int *rowptr, int *colind,
#ifdef Have_SLU_		 
		 fptr *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
#else 
		 void *f_factors,
#endif
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *  performs LU decomposition.
 *
 * f_factors (input/output) fptr* 
 *      On  output contains the pointer pointing to
 *       the structure of the factored matrices.
 *
 */
 
#ifdef Have_SLU_
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, panel_size, permc_spec, relax;
    trans_t  trans;
    double   drop_tol = 0.0;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors;

    trans = NOTRANS;

    
    /* Set the default input options. */
    set_default_options(&options);
    
    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    /* Adjust to 0-based indexing */
    for (i = 0; i < *nnz; ++i) --colind[i];
    for (i = 0; i <= *n; ++i) --rowptr[i];
    
    dCreate_CompRow_Matrix(&A, *n, *n, *nnz, values, colind, rowptr,
			   SLU_NR, SLU_D, SLU_GE);
    L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
    if ( !(perm_r = intMalloc(*n)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(*n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(etree = intMalloc(*n)) ) ABORT("Malloc fails for etree[].");
    
    /*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = 0: natural ordering 
	 *   permc_spec = 1: minimum degree on structure of A'*A
	 *   permc_spec = 2: minimum degree on structure of A'+A
	 *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	 */    	
    options.ColPerm=2;
    permc_spec = options.ColPerm;
    get_perm_c(permc_spec, &A, perm_c);
    
    sp_preorder(&options, &A, perm_c, etree, &AC);
    
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    
    dgstrf(&options, &AC, drop_tol, relax, panel_size, 
	   etree, NULL, 0, perm_c, perm_r, L, U, &stat, info);
    
    if ( *info == 0 ) {
      Lstore = (SCformat *) L->Store;
      Ustore = (NCformat *) U->Store;
      dQuerySpace(L, U, &mem_usage);
#if 0
      printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
      printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
      printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
      printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	     mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	     mem_usage.expansions);
#endif
    } else {
      printf("dgstrf() error returns INFO= %d\n", *info);
      if ( *info <= *n ) { /* factorization completes */
	dQuerySpace(L, U, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	       mem_usage.expansions);
      }
    }
    
    /* Restore to 1-based indexing */
    for (i = 0; i < *nnz; ++i) ++colind[i];
    for (i = 0; i <= *n; ++i) ++rowptr[i];
    
    /* Save the LU factors in the factors handle */
    LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
    LUfactors->L = L;
    LUfactors->U = U;
    LUfactors->perm_c = perm_c;
    LUfactors->perm_r = perm_r;
    *f_factors = (fptr) LUfactors;
    
    /* Free un-wanted storage */
    SUPERLU_FREE(etree);
    Destroy_SuperMatrix_Store(&A);
    Destroy_CompCol_Permuted(&AC);
    StatFree(&stat);
#else
    fprintf(stderr," SLU Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif
}


void
psb_dslu_solve_(int *itrans, int *n, int *nrhs, 
                 double *b, int *ldb,
#ifdef Have_SLU_		 
		 fptr *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
#else 
		 void *f_factors,
#endif
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *      performs triangular solve
 *
 */
#ifdef Have_SLU_ 
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, panel_size, permc_spec, relax;
    trans_t  trans;
    double   drop_tol = 0.0;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors;

    if (*itrans == 0) {
      trans = NOTRANS;
    } else if (*itrans ==1) {
      trans = TRANS;
    } else if (*itrans ==2) {
      trans = CONJ;
    } else {
      trans = NOTRANS;
    }
    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    /* Extract the LU factors in the factors handle */
    LUfactors = (factors_t*) *f_factors;
    L = LUfactors->L;
    U = LUfactors->U;
    perm_c = LUfactors->perm_c;
    perm_r = LUfactors->perm_r;
    
    dCreate_Dense_Matrix(&B, *n, *nrhs, b, *ldb, SLU_DN, SLU_D, SLU_GE);
    /* Solve the system A*X=B, overwriting B with X. */
    dgstrs (trans, L, U, perm_c, perm_r, &B, &stat, info);

    Destroy_SuperMatrix_Store(&B);
    StatFree(&stat);
#else
    fprintf(stderr," SLU Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif
    
}


void
psb_dslu_free_(
#ifdef Have_SLU_		 
 fptr *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
#else 
		 void *f_factors,
#endif
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *
 *      free all storage in the end
 *
 */
#ifdef Have_SLU_ 
    SuperMatrix A, AC, B;
    SuperMatrix *L, *U;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int *etree;  /* column elimination tree */
    SCformat *Lstore;
    NCformat *Ustore;
    int      i, panel_size, permc_spec, relax;
    trans_t  trans;
    double   drop_tol = 0.0;
    mem_usage_t   mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;
    factors_t *LUfactors; 

    trans = NOTRANS;
    /* Free the LU factors in the factors handle */
    LUfactors = (factors_t*) *f_factors;
    SUPERLU_FREE (LUfactors->perm_r);
    SUPERLU_FREE (LUfactors->perm_c);
    Destroy_SuperNode_Matrix(LUfactors->L);
    Destroy_CompCol_Matrix(LUfactors->U);
    SUPERLU_FREE (LUfactors->L);
    SUPERLU_FREE (LUfactors->U);
    SUPERLU_FREE (LUfactors);
    *info = 0;
#else
    fprintf(stderr," SLU Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif
}


