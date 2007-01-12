/*
 *                    MD2P4
 *    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
 *                              for 
 *                    Parallel Sparse BLAS  v2.0
 *
 *   (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
 *                      Alfredo Buttari
 *                      Daniela di Serafino    Second University of Naples
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
/* This file is an interface to the UMFPACK  routines for 
   factorization. It was obtained by adapting  umfpack_zi_demo
   under the original copyright terms reproduced below.   
   PSBLAS v 2.0  */ 


/*		=====================
UMFPACK Version 4.4 (Jan. 28, 2005), Copyright (c) 2005 by Timothy A.
Davis.  All Rights Reserved.

UMFPACK License:

    Your use or distribution of UMFPACK or any modified version of
    UMFPACK implies that you agree to this License.

    THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
    EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

    Permission is hereby granted to use or copy this program, provided
    that the Copyright, this License, and the Availability of the original
    version is retained on all copies.  User documentation of any code that
    uses UMFPACK or any modified version of UMFPACK code must cite the
    Copyright, this License, the Availability note, and "Used by permission."
    Permission to modify the code and to distribute modified code is granted,
    provided the Copyright, this License, and the Availability note are
    retained, and a notice that the code was modified is included.  This
    software was developed with support from the National Science Foundation,
    and is provided to you free of charge.

Availability:

    http://www.cise.ufl.edu/research/sparse/umfpack

*/



#ifdef Add_
#define psb_zumf_factor_ psb_zumf_factor_
#define psb_zumf_solve_  psb_zumf_solve_
#define psb_zumf_free_   psb_zumf_free_
#endif
#ifdef AddDouble_
#define psb_zumf_factor_ psb_zumf_factor__
#define psb_zumf_solve_  psb_zumf_solve__
#define psb_zumf_free_   psb_zumf_free__
#endif
#ifdef NoChange
#define psb_zumf_factor_ psb_zumf_factor
#define psb_zumf_solve_  psb_zumf_solve
#define psb_zumf_free_   psb_zumf_free
#endif


#include <stdio.h>
#ifdef Have_UMF_		 
#include "umfpack.h"
#endif

#ifdef LargeFptr
typedef long long fptr;  /* 64-bit*/
#else
typedef int fptr;  /* 32-bit by default */
#endif

void
psb_zumf_factor_(int *n, int *nnz,
                 double *values, int *rowind, int *colptr,
#ifdef Have_UMF_		 
		 fptr *symptr, 
		 fptr *numptr, 
		 
#else 
		 void *symptr,
		 void *numptr,
#endif
		 int *info)

{
 
#ifdef Have_UMF_
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
  void *Symbolic, *Numeric ;
  int i;
  
  
  umfpack_zi_defaults(Control);
  
  for (i = 0; i <= *n;  ++i) --colptr[i];
  for (i = 0; i < *nnz; ++i) --rowind[i];
  *info = umfpack_zi_symbolic (*n, *n, colptr, rowind, values, NULL, &Symbolic,
				Control, Info);
  
    
  if ( *info == UMFPACK_OK ) {
    *info = 0;
  } else {
    printf("umfpack_zi_symbolic() error returns INFO= %d\n", *info);
    *info = -11;
    *numptr = (fptr) NULL; 
    return;
  }
    
  *symptr = (fptr) Symbolic; 
  
  *info = umfpack_zi_numeric (colptr, rowind, values, NULL, Symbolic, &Numeric,
				Control, Info) ;
  
    
  if ( *info == UMFPACK_OK ) {
    *info = 0;
    *numptr = (fptr) Numeric; 
  } else {
    printf("umfpack_zi_numeric() error returns INFO= %d\n", *info);
    *info = -12;
    *numptr = (fptr) NULL; 
  }
    
  for (i = 0; i <= *n;  ++i) ++colptr[i];
  for (i = 0; i < *nnz; ++i) ++rowind[i];
#else
    fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif    
}


void
psb_zumf_solve_(int *itrans, int *n,  
                 double *x,  double *b, int *ldb,
#ifdef Have_UMF_		 
		 fptr *numptr, 
		 
#else 
		 void *numptr,
#endif
		 int *info)

{
#ifdef Have_UMF_ 
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
  void *Symbolic, *Numeric ;
  int i,trans;
  
  
  umfpack_di_defaults(Control);
  Control[UMFPACK_IRSTEP]=0;


  if (*itrans == 0) {
    trans = UMFPACK_A;
  } else if (*itrans ==1) {
    trans = UMFPACK_At;
  } else {
    trans = UMFPACK_A;
  }

  *info = umfpack_zi_solve(trans,NULL,NULL,NULL,NULL,
			   x,NULL,b,NULL,(void *) *numptr,Control,Info);
  
#else
    fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif
    
}


void
psb_zumf_free_(
#ifdef Have_UMF_		 
		 fptr *symptr, 
		 fptr *numptr, 
		 
#else 
		 void *symptr,
		 void *numptr,
#endif
		 int *info)

{
#ifdef Have_UMF_ 
  void *Symbolic, *Numeric ;
  Symbolic = (void *) *symptr;
  Numeric  = (void *) *numptr;
  
  umfpack_zi_free_numeric(&Numeric);
  umfpack_zi_free_symbolic(&Symbolic);
  *info=0;
#else
    fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif
}


