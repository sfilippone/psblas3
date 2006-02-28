/*
 *             Parallel Sparse BLAS  v2.0
 *   (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
 *                      Alfredo Buttari        University of Rome Tor Vergata
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the PSBLAS group or the names of its contributors may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

/*
* This file includes the standard C libraries, as well as system
* dependent include files.  All PSBLAS routines include this file.
*/
#include <string.h>

#ifndef PSBLASH
#define PSBLASH
/*
* ========================================================================
* Machine Specific PBLAS macros
* ========================================================================
*/
/* This is a debugging option. 
   #define PS_CONTROL_LEVEL   */

#define _HAL_           0
#define _T3D_           1

#ifdef T3D
#define _MACH_          _T3D_
#endif

#ifndef _MACH_
#define _MACH_          _HAL_
#endif

/*
* ========================================================================
* Include files
* ========================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if( _MACH_ == _T3D_ )
#include <fortran.h>
#endif

#ifdef USE_FBLACS
#ifndef CTOF_BLACS
#include "ctof_blacs.h"
#endif
#endif



/*
* ========================================================================
* FORTRAN <-> C interface
* ========================================================================
*
* These macros define how the PBLAS will be called. _F2C_ADD_ assumes
* that they will be called by FORTRAN, which expects C routines to have
* an underscore postfixed to the name (Suns, and Intel machines expect
* this). _F2C_NOCHANGE indicates that FORTRAN will be calling, and that
* it expects the name called by FORTRAN to be identical to that compiled
* by the C (RS6K's do this).  _F2C_UPCASE says it expects C routines
* called by FORTRAN to be in all upcase (CRAY wants this).
*/

#define _F2C_ADD_       0
#define _F2C_NOCHANGE   1
#define _F2C_UPCASE     2

#ifdef UpCase
#define _F2C_CALL_      _F2C_UPCASE
#endif

#ifdef NoChange
#define _F2C_CALL_      _F2C_NOCHANGE
#endif

#ifdef Add_
#define _F2C_CALL_      _F2C_ADD_
#endif

#ifndef _F2C_CALL_
#define _F2C_CALL_      _F2C_ADD_
#endif

/*
* ========================================================================
* TYPE DEFINITIONS AND CONVERSION UTILITIES
* ========================================================================
*/

typedef struct { float  re, im; } complex;
typedef struct { double re, im; } complex16;

#if( _MACH_ == _T3D_ )
                       /* Type of character argument in a FORTRAN call */
#define F_CHAR          _fcd
                                     /* Character conversion utilities */
#define F2C_CHAR(a)     ( _fcdtocp( (a) ) )
#define C2F_CHAR(a)     ( _cptofcd( (a), 1 ) )
                                          /* Type of FORTRAN functions */
#define F_VOID_FCT      void   fortran                   /* Subroutine */
#define F_INTG_FCT      int    fortran             /* INTEGER function */
#define F_DBLE_FCT      double fortran    /* DOUBLE PRECISION function */

#else
                       /* Type of character argument in a FORTRAN call */
typedef char *          F_CHAR;
                                     /* Character conversion utilities */
#define F2C_CHAR(a)     (a)
#define C2F_CHAR(a)     (a)
                                          /* Type of FORTRAN functions */
#define F_VOID_FCT      void                             /* Subroutine */
#define F_INTG_FCT      int                        /* INTEGER function */
#define F_DBLE_FCT      double            /* DOUBLE PRECISION function */

#endif

/*
* ======================================================================
* FUNCTIONS PROTOTYPES
* ======================================================================
*/
void DVSct(int n, int k,int idx[],int flag, double X[], int lx,
	   double beta, double Y[], int ly);
void DVGth(int n, int k,int idx[],int flag, double X[], int lx,double Y[], int ly);
void IVSct(int n, int k,int idx[],int flag, int X[], int lx,
	   int beta, int Y[], int ly);
void IVGth(int n, int k,int idx[],int flag, int X[], int lx,int Y[], int ly);

void PSI_dSwapData(int iflag, int n, double beta, double Y[], int ly,
		   int desc_data[], int desc_halo[], 
		   double *work, int *lwork, int *ierror);

void PSI_dSwapTran(int flag, int n, double beta, double Y[], int ly,
		   int desc_data[], int desc_halo[],
		   double *work, int *lwork, int *ierror);

void PSI_zSwapData(int n, double Y[], int ly, int desc_data[], int desc_halo[], 
		   double *work, int *lwork, int *ierror);

void PSI_zSwapOverlap(double Y[], double Sum_Ovrlap[], int desc_data[], 
		      int desc_ovrlap[], double work[], int *lwork, int *ierror);
void PSI_iSwapData(int iflag, int n, int beta, int Y[], int ly,
		   int desc_data[], int desc_halo[], 
		   int *work, int *lwork, int *ierror);

void PSI_iSwapTran(int flag, int n, int beta, int Y[], int ly,
		   int desc_data[], int desc_halo[], 
		   int *work, int *lwork, int *ierror);

/*
* ========================================================================
* #DEFINE MACRO CONSTANTS
* ========================================================================
*/
/* MACRO max */
#define max(x,y) ((x)>(y)?(x):(y))

/*MACRO for ovrlap update*/
#define NOHALO_         0
#define HALO_           4
#define NONE_           0
#define SUM_            1
#define AVG_            2
#define SQUARE_ROOT_    3 

/*  Bit fields to control swapdata/ovrlap behaviour.
    BEWARE: check consistency with tools_const.f.
    Should it be automated? */
#define SWAP_SEND       1
#define SWAP_RECV       2
#define SWAP_SYNC       4
#define SWAP_MPI        8


/* Macro for MATRIX_DATA array */
#define DEC_TYPE_	0		     /* The type of decomposition of global
   					        matrix A. */
#define M_	        1		     /* Number of equations */
#define N_              2		     /* Number of variables */
#define N_ROW_		3		     /* The number of row of local matrix. */
#define N_COL_		4		     /* The number of columns of local
   					        matrix. */
#define CTXT_		5		     /* The BLACS context handle, indicating
   					        the global context of the operation
   					        on the matrix.
   					        The context itself is global. */
#define LOC_TO_GLOB_    6                    /* The pointer to the array 
						loc_to_glob */
#define MPI_C_		8		     /* The MPI Fortran handle */
/* values for DEC_TYPE_  */
#define DESC_ASB     3099
#define DESC_BLD     (DESC_ASB+1)

/* Macro for HALO array */
#define PROC_ID_	0		     /* The identifier of domain. */	
#define N_ELEM_RECV_	1		     /* The number of elements to receive*/
#define ELEM_RECV_	2		     /* The first index of local elements */
#define N_ELEM_SEND_    2		     /* The number of elements to send */
#define ELEM_SEND_	3		     /* The first index of local elements */

/* Macro for OVERLAP array */
#define N_OVRLP_ELEM_	1		     /* The number of overlap elements to recv/send */
#define OVRLP_ELEM_TO_	2		     /* The first index of local elements */

/* Macro for OVR_ELEM_D array */
#define OVRLP_ELEM_	0
#define N_DOM_OVR_	1

#define    BROADCAST    "B"              /* Blacs operation definitions */
#define    COMBINE      "C"

#define    ALL          "A"                        /* Scope definitions */
#define    COLUMN       "C"
#define    ROW          "R"

#define    TOPDEF       " " /* Default BLACS topology, PB-BLAS routines */
#define    CTOPDEF      ' '
#define    TOPGET       "!"

#define    YES          "Y"
#define    NO           "N"

#define    MULLENFAC    2

#define    ONE          1.0
#define    ZERO         0.0

/* Integer values for error checking */
#define    no_err       0
#define    act_ret      0
#define    act_abort    1


/*
* ========================================================================
* PREPROCESSOR MACRO FUNCTIONS USED FOR OPTIMIZATION & CONVENIENCE
* ========================================================================
*/

#define ABS(a)   ((a > 0) ? (a) : (-a))

#define MIN(a,b) ((a < b) ? (a) : (b))

#define MAX(a,b) ((a > b) ? (a) : (b))

#define CEIL(a,b) ( (a+b-1) / (b) )

#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )

#define Mupcase(C) ( ((C) > 96 && (C) < 123) ? (C) & 0xDF : (C) )

#define INDXG2L( iglob, nb, iproc, isrcproc, nprocs )\
    ( (nb) * ( ( (iglob)-1) / ( (nb) * (nprocs) ) ) +\
      ( ( (iglob) - 1 ) % (nb) ) + 1 )

#define INDXL2G( iloc, nb, iproc, isrcproc, nprocs )\
    ( (nprocs) * (nb) * ( ( (iloc) - 1 ) / (nb) ) +\
      ( ( (iloc) - 1 ) % (nb) ) +\
      ( ( (nprocs) + (iproc) - (isrcproc) ) % (nprocs) ) * (nb) + 1 )

#define INDXG2P( iglob, nb, iproc, isrcproc, nprocs ) \
    ( ( (isrcproc) + ( (iglob) - 1 ) / (nb) ) % (nprocs) )

#define MYROC0( nblocks, n, nb, nprocs )\
  ( ( (nblocks) % (nprocs) ) ? ( ( (nblocks) / (nprocs) ) * (nb) + (nb) )\
                   : ( ( (nblocks) / (nprocs) )* (nb) + ( (n) % (nb) ) ) )

#if( _F2C_CALL_ == _F2C_ADD_ )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in).
* No redefinition necessary to have following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void pdgemm_(...)
*
* This is the default.
*/
#define pbchkvectf pbchkvectf_
#define fcpsb_errcomm            fcpsb_errcomm_
#define fcpsb_erractionsave      fcpsb_erractionsave_
#define fcpsb_erractionrestore   fcpsb_erractionrestore_
#define fcpsb_perror             fcpsb_perror_
#define fcpsb_serror             fcpsb_serror_
#define fcpsb_errpush            fcpsb_errpush_
#endif

