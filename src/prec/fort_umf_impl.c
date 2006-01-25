/* This file is an interface to the UMFPACK  routines for sparse
   factorization.   
   
   PSBLAS v 2.0, rc1,  May 03, 2005  */ 


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
#define fort_umf_factor_ fort_umf_factor_
#define fort_umf_solve_  fort_umf_solve_
#define fort_umf_free_   fort_umf_free_
#endif
#ifdef AddDouble_
#define fort_umf_factor_ fort_umf_factor__
#define fort_umf_solve_  fort_umf_solve__
#define fort_umf_free_   fort_umf_free__
#endif
#ifdef NoChange
#define fort_umf_factor_ fort_umf_factor
#define fort_umf_solve_  fort_umf_solve
#define fort_umf_free_   fort_umf_free
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
fort_umf_factor_(int *n, int *nnz,
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
/* 
 * This routine can be called from Fortran.
 *  performs LU decomposition.
 *
 * f_factors (input/output) fptr* 
 *      On  output contains the pointer pointing to
 *       the structure of the factored matrices.
 *
 */
 
#ifdef Have_UMF_
  double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
  void *Symbolic, *Numeric ;
  int i;
  
  
  umfpack_di_defaults(Control);
  
  for (i = 0; i <= *n;  ++i) --colptr[i];
  for (i = 0; i < *nnz; ++i) --rowind[i];
  *info = umfpack_di_symbolic (*n, *n, colptr, rowind, values, &Symbolic,
				Control, Info);
  
    
  if ( *info == UMFPACK_OK ) {
    *info = 0;
  } else {
    printf("umfpack_di_symbolic() error returns INFO= %d\n", *info);
  }
    
  *symptr = (fptr) Symbolic; 
  
  *info = umfpack_di_numeric (colptr, rowind, values, Symbolic, &Numeric,
				Control, Info) ;
  
    
  if ( *info == UMFPACK_OK ) {
    *info = 0;
  } else {
    printf("umfpack_di_numeric() error returns INFO= %d\n", *info);
  }
    
  *numptr = (fptr) Numeric; 
  for (i = 0; i <= *n;  ++i) ++colptr[i];
  for (i = 0; i < *nnz; ++i) ++rowind[i];
#else
    fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif    
}


void
fort_umf_solve_(int *itrans, int *n,  
                 double *x,  double *b, int *ldb,
#ifdef Have_UMF_		 
		 fptr *numptr, 
		 
#else 
		 void *numptr,
#endif
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *      performs triangular solve
 *
 */
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

  *info = umfpack_di_solve(trans,NULL,NULL,NULL,
			   x,b,(void *) *numptr,Control,Info);
  
#else
    fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif
    
}


void
fort_umf_free_(
#ifdef Have_UMF_		 
		 fptr *symptr, 
		 fptr *numptr, 
		 
#else 
		 void *symptr,
		 void *numptr,
#endif
		 int *info)

{
/* 
 * This routine can be called from Fortran.
 *
 *      free all storage in the end
 *
 */
#ifdef Have_UMF_ 
  void *Symbolic, *Numeric ;
  Symbolic = (void *) *symptr;
  Numeric  = (void *) *numptr;
  
  umfpack_di_free_numeric(&Numeric);
  umfpack_di_free_symbolic(&Symbolic);
  *info=0;
#else
    fprintf(stderr," UMF Not Configured, fix make.inc and recompile\n");
    *info=-1;
#endif
}


