/* ---------------------------------------------------------------------
*
*  -- PSBLAS routine (version 1.0) --
*
*  ---------------------------------------------------------------------
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
#define dcsmm   dcsmm_
#define dcssm   dcssm_
#define dcsnmi  dcsnmi_
#define idamax  idamax_
#define izamax  izamax_
#define ddot    ddot_
#define dasum   dasum_
#define daxpby  daxpby_
#define dscal   dscal_
#define zcsmm   zcsmm_
#define zcssm   zcssm_
#define zcsnmi  zcsnmi_
#define zdot    zdot_
#define dzasum   dzasum_
#define zaxpby  zaxpby_
#define zscal   zscal_
#define pbchkvectf pbchkvectf_
#define fcpsb_errcomm            fcpsb_errcomm_
#define fcpsb_erractionsave      fcpsb_erractionsave_
#define fcpsb_erractionrestore   fcpsb_erractionrestore_
#define fcpsb_perror             fcpsb_perror_
#define fcpsb_serror             fcpsb_serror_
#define fcpsb_errpush            fcpsb_errpush_
#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in)
* following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void PDGEMM(...)
*/
#define pbchkvectf        PBCHKVECTF                          /* PSBLAS */
#define psddot_           PSDDOT
#define psdmdot_           PSDMDOT
#define psddot_sub_       PSDDOT_SUB
#define psdaxpby_         PSDAXPBY
#define psdamax_          PSDAMAX
#define psdmamax_          PSDMAMAX
#define psdasum_          PSDASUM
#define psdnrm2_          PSDNRM2
#define psdnrmi_          PSDNRMI
#define psdnrmisym_       PSDNRMISYM
#define psdhalo_          PSDHALO
#define psihalo_          PSIHALO
#define psdhred_          PSDHRED
#define psdovrl_          PSDOVRL
#define psdspmm_          PSDSPMM
#define psdswaptran_      PSDSWAPTRAN
#define psdspmmsym_       PSDSPMMSYM
#define psdspsm_          PSDSPSM
#define psderror_         PSDERROR
#define psdverify_        PSDVERIFY
#define psdscatterm_      PSDSCATTERM
#define psdgatherm        PSDGATHERM
                                                        /* PSBLAS */
#define pszdotc_          PSZDOTC
#define pszdotu_          PSZDOTU
#define pszmdot_           PSZMDOT
#define pszaxpby_         PSZAXPBY
#define pszamax_          PSZAMAX
#define pszmamax_         PSZMAMAX
#define pszasum_          PSZASUM
#define psznrm2_          PSZNRM2
#define psznrmi_          PSZNRMI
#define psznrmisym_       PSZNRMISYM
#define pszhalo_          PSZHALO
#define pszovrl_          PSZOVRL
#define pszspmm_          PSZSPMM
#define pszspmmsym_       PSZSPMMSYM
#define pszspsm_          PSZSPSM
#define pszerror_         PSZERROR
#define pszverify_        PSZVERIFY
#define pszscatterm_      PSZSCATTERM
#define pszgatherm_       PSZGATHERM
                                                   /* BLACS */
#define blacs_abort_      BLACS_ABORT
#define blacs_gridinfo_   BLACS_GRIDINFO

#define igesd2d_          IGESD2D
#define igebs2d_          IGEBS2D
#define itrsd2d_          ITRSD2D
#define itrbs2d_          ITRBS2D
#define igerv2d_          IGERV2D
#define igebr2d_          IGEBR2D
#define itrrv2d_          ITRRV2D
#define itrbr2d_          ITRBR2D
#define igamx2d_          IGAMX2D
#define igamn2d_          IGAMN2D
#define igsum2d_          IGSUM2D

#define sgesd2d_          SGESD2D
#define sgebs2d_          SGEBS2D
#define strsd2d_          STRSD2D
#define strbs2d_          STRBS2D
#define sgerv2d_          SGERV2D
#define sgebr2d_          SGEBR2D
#define strrv2d_          STRRV2D
#define strbr2d_          STRBR2D
#define sgamx2d_          SGAMX2D
#define sgamn2d_          SGAMN2D
#define sgsum2d_          SGSUM2D

#define dgesd2d_          DGESD2D
#define dgebs2d_          DGEBS2D
#define dtrsd2d_          DTRSD2D
#define dtrbs2d_          DTRBS2D
#define dgerv2d_          DGERV2D
#define dgebr2d_          DGEBR2D
#define dtrrv2d_          DTRRV2D
#define dtrbr2d_          DTRBR2D
#define dgamx2d_          DGAMX2D
#define dgamn2d_          DGAMN2D
#define dgsum2d_          DGSUM2D

#define cgesd2d_          CGESD2D
#define cgebs2d_          CGEBS2D
#define ctrsd2d_          CTRSD2D
#define ctrbs2d_          CTRBS2D
#define cgerv2d_          CGERV2D
#define cgebr2d_          CGEBR2D
#define ctrrv2d_          CTRRV2D
#define ctrbr2d_          CTRBR2D
#define cgamx2d_          CGAMX2D
#define cgamn2d_          CGAMN2D
#define cgsum2d_          CGSUM2D

#define zgesd2d_          ZGESD2D
#define zgebs2d_          ZGEBS2D
#define ztrsd2d_          ZTRSD2D
#define ztrbs2d_          ZTRBS2D
#define zgerv2d_          ZGERV2D
#define zgebr2d_          ZGEBR2D
#define ztrrv2d_          ZTRRV2D
#define ztrbr2d_          ZTRBR2D
#define zgamx2d_          ZGAMX2D
#define zgamn2d_          ZGAMN2D
#define zgsum2d_          ZGSUM2D
                                                     /* Level-1 BLAS */
#define srotg_            SROTG
#define srotmg_           SROTMG
#define srot_             SROT
#define srotm_            SROTM
#define sswap_            SSWAP
#define sscal_            SSCAL
#define scopy_            SCOPY
#define saxpy_            SAXPY
#define ssdot_            SSDOT
#define isamax_           ISAMAX

#define drotg_            DROTG
#define drotmg_           DROTMG
#define drot_             DROT
#define drotm_            DROTM
#define dswap_            DSWAP
#define dscal_            DSCAL
#define dcopy_            DCOPY
#define daxpy_            DAXPY
#define dddot_            DDDOT
#define dnrm2_            DNRM2
#define dsnrm2_           DSNRM2
#define dasum_            DASUM
#define dsasum_           DSASUM
#define idamax_           IDAMAX
#define daxpby_           DAXPBY

#define zaxpby_           ZAXPBY      /* to match added internal function */

#define cswap_            CSWAP
#define cscal_            CSCAL
#define csscal_           CSSCAL
#define ccopy_            CCOPY
#define caxpy_            CAXPY
#define ccdotu_           CCDOTU
#define ccdotc_           CCDOTC
#define icamax_           ICAMAX

#define zswap_            ZSWAP
#define zscal_            ZSCAL
#define zdscal_           ZDSCAL
#define zcopy_            ZCOPY
#define zaxpy_            ZAXPY
#define zzdotu_           ZZDOTU
#define zzdotc_           ZZDOTC
#define dscnrm2_          DSCNRM2
#define dznrm2_           DZNRM2
#define dscasum_          DSCASUM
#define dzasum_           DZASUM
#define izamax_           IZAMAX
                                                     /* Level-2 BLAS */
#define sgemv_            SGEMV
#define ssymv_            SSYMV
#define strmv_            STRMV
#define strsv_            STRSV
#define sger_             SGER
#define ssyr_             SSYR
#define ssyr2_            SSYR2

#define dgemv_            DGEMV
#define dsymv_            DSYMV
#define dtrmv_            DTRMV
#define dtrsv_            DTRSV
#define dger_             DGER
#define dsyr_             DSYR
#define dsyr2_            DSYR2

#define cgemv_            CGEMV
#define chemv_            CHEMV
#define ctrmv_            CTRMV
#define ctrsv_            CTRSV
#define cgeru_            CGERU
#define cgerc_            CGERC
#define cher_             CHER
#define cher2_            CHER2

#define zgemv_            ZGEMV
#define zhemv_            ZHEMV
#define ztrmv_            ZTRMV
#define ztrsv_            ZTRSV
#define zgeru_            ZGERU
#define zgerc_            ZGERC
#define zher_             ZHER
#define zher2_            ZHER2
                                                     /* Level-3 BLAS */
#define sgemm_            SGEMM
#define ssymm_            SSYMM
#define ssyrk_            SSYRK
#define ssyr2k_           SSYR2K
#define strmm_            STRMM
#define strsm_            STRSM

#define dgemm_            DGEMM
#define dsymm_            DSYMM
#define dsyrk_            DSYRK
#define dsyr2k_           DSYR2K
#define dtrmm_            DTRMM
#define dtrsm_            DTRSM

#define cgemm_            CGEMM
#define chemm_            CHEMM
#define csymm_            CSYMM
#define csyrk_            CSYRK
#define cherk_            CHERK
#define csyr2k_           CSYR2K
#define cher2k_           CHER2K
#define ctrmm_            CTRMM
#define ctrsm_            CTRSM

#define zgemm_            ZGEMM
#define zhemm_            ZHEMM
#define zsymm_            ZSYMM
#define zsyrk_            ZSYRK
#define zherk_            ZHERK
#define zsyr2k_           ZSYR2K
#define zher2k_           ZHER2K
#define ztrmm_            ZTRMM
#define ztrsm_            ZTRSM
                                                 /* Auxilliary PBLAS */
#define pberror_          PBERROR
#define pbfreebuf_        PBFREEBUF

#define dcsmm   DCSMM
#define dcssm   DCSSM
#define dcsnmi  DCSNMI
#define zcsnmi  ZCSNMI

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in)
* for following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void pdgemm(...)
*/


                                                        /* PSBLAS */
#define psddot_           psddot
#define psdmdot_          psdmdot
#define psdaxpby_         psdaxpby
#define psdamax_          psdamax
#define psdmamax_         psdmamax
#define psdasum_          psdasum
#define psdnrm2_          psdnrm2
#define psdnrmi_          psdnrmi
#define psdnrmisym_       psdnrmisym
#define psdhalo_          psdhalo
#define psihalo_          psihalo
#define psdhred_          psdhred
#define psdovrl_          psdovrl
#define psdspmm_          psdspmm
#define psdswaptran_      psdswaptran
#define psdspmmsym_       psdspmmsym
#define psdspsm_          psdspsm
#define psderror_         psderror
#define psdverify_        psdverify
#define psdscatterm_      psdscatterm
#define psdgatherm_       psdgatherm

#define pszmdot_          pszmdot
#define pszdotc_          pszdotc
#define pszdotu_          pszdotu
#define pszaxpby_         pszaxpby
#define pszamax_          pszamax
#define pszmamax_         pszmamax
#define pszasum_          pszasum
#define psznrm2_          psznrm2
#define psznrmi_          psznrmi
#define psznrmisym_       psznrmisym
#define pszhalo_          pszhalo
#define pszovrl_          pszovrl
#define pszspmm_          pszspmm
#define pszspmmsym_       pszspmmsym
#define pszspsm_          pszspsm
#define pszerror_         pszerror
#define pszverify_        pszverify
#define pszscatterm_      pszscatterm
#define pszgatherm_       pszgatherm


                                                            /* BLACS */
#define blacs_abort_      blacs_abort
#define blacs_gridinfo_   blacs_gridinfo

#define igesd2d_          igesd2d
#define igebs2d_          igebs2d
#define itrsd2d_          itrsd2d
#define itrbs2d_          itrbs2d
#define igerv2d_          igerv2d
#define igebr2d_          igebr2d
#define itrrv2d_          itrrv2d
#define itrbr2d_          itrbr2d
#define igamx2d_          igamx2d
#define igamn2d_          igamn2d
#define igsum2d_          igsum2d

#define sgesd2d_          sgesd2d
#define sgebs2d_          sgebs2d
#define strsd2d_          strsd2d
#define strbs2d_          strbs2d
#define sgerv2d_          sgerv2d
#define sgebr2d_          sgebr2d
#define strrv2d_          strrv2d
#define strbr2d_          strbr2d
#define sgamx2d_          sgamx2d
#define sgamn2d_          sgamn2d
#define sgsum2d_          sgsum2d

#define dgesd2d_          dgesd2d
#define dgebs2d_          dgebs2d
#define dtrsd2d_          dtrsd2d
#define dtrbs2d_          dtrbs2d
#define dgerv2d_          dgerv2d
#define dgebr2d_          dgebr2d
#define dtrrv2d_          dtrrv2d
#define dtrbr2d_          dtrbr2d
#define dgamx2d_          dgamx2d
#define dgamn2d_          dgamn2d
#define dgsum2d_          dgsum2d

#define cgesd2d_          cgesd2d
#define cgebs2d_          cgebs2d
#define ctrsd2d_          ctrsd2d
#define ctrbs2d_          ctrbs2d
#define cgerv2d_          cgerv2d
#define cgebr2d_          cgebr2d
#define ctrrv2d_          ctrrv2d
#define ctrbr2d_          ctrbr2d
#define cgamx2d_          cgamx2d
#define cgamn2d_          cgamn2d
#define cgsum2d_          cgsum2d

#define zgesd2d_          zgesd2d
#define zgebs2d_          zgebs2d
#define ztrsd2d_          ztrsd2d
#define ztrbs2d_          ztrbs2d
#define zgerv2d_          zgerv2d
#define zgebr2d_          zgebr2d
#define ztrrv2d_          ztrrv2d
#define ztrbr2d_          ztrbr2d
#define zgamx2d_          zgamx2d
#define zgamn2d_          zgamn2d
#define zgsum2d_          zgsum2d
                                                     /* Level-1 BLAS */
#define srotg_            srotg
#define srotmg_           srotmg
#define srot_             srot
#define srotm_            srotm
#define sswap_            sswap
#define sscal_            sscal
#define scopy_            scopy
#define saxpy_            saxpy
#define ssdot_            ssdot
#define isamax_           isamax

#define drotg_            drotg
#define drotmg_           drotmg
#define drot_             drot
#define drotm_            drotm
#define dswap_            dswap
#define dscal_            dscal
#define dcopy_            dcopy
#define daxpy_            daxpy
#define dddot_            dddot
#define dnrm2_            dnrm2
#define dsnrm2_           dsnrm2
#define dasum_            dasum
#define dsasum_           dsasum
#define idamax_           idamax
#define daxpby_           daxpby

#define zaxpby_           zaxpby

#define cswap_            cswap
#define cscal_            cscal
#define csscal_           csscal
#define ccopy_            ccopy
#define caxpy_            caxpy
#define ccdotu_           ccdotu
#define ccdotc_           ccdotc
#define icamax_           icamax

#define zswap_            zswap
#define zscal_            zscal
#define zdscal_           zdscal
#define zcopy_            zcopy
#define zaxpy_            zaxpy
#define zzdotu_           zzdotu
#define zzdotc_           zzdotc
#define dscnrm2_          dscnrm2
#define dznrm2_           dznrm2
#define dscasum_          dscasum
#define dzasum_           dzasum
#define izamax_           izamax
                                                     /* Level-2 BLAS */
#define sgemv_            sgemv
#define ssymv_            ssymv
#define strmv_            strmv
#define strsv_            strsv
#define sger_             sger
#define ssyr_             ssyr
#define ssyr2_            ssyr2

#define dgemv_            dgemv
#define dsymv_            dsymv
#define dtrmv_            dtrmv
#define dtrsv_            dtrsv
#define dger_             dger
#define dsyr_             dsyr
#define dsyr2_            dsyr2

#define cgemv_            cgemv
#define chemv_            chemv
#define ctrmv_            ctrmv
#define ctrsv_            ctrsv
#define cgeru_            cgeru
#define cgerc_            cgerc
#define cher_             cher
#define cher2_            cher2

#define zgemv_            zgemv
#define zhemv_            zhemv
#define ztrmv_            ztrmv
#define ztrsv_            ztrsv
#define zgeru_            zgeru
#define zgerc_            zgerc
#define zher_             zher
#define zher2_            zher2
                                                     /* Level-3 BLAS */
#define sgemm_            sgemm
#define ssymm_            ssymm
#define ssyrk_            ssyrk
#define ssyr2k_           ssyr2k
#define strmm_            strmm
#define strsm_            strsm

#define dgemm_            dgemm
#define dsymm_            dsymm
#define dsyrk_            dsyrk
#define dsyr2k_           dsyr2k
#define dtrmm_            dtrmm
#define dtrsm_            dtrsm

#define cgemm_            cgemm
#define chemm_            chemm
#define csymm_            csymm
#define csyrk_            csyrk
#define cherk_            cherk
#define csyr2k_           csyr2k
#define cher2k_           cher2k
#define ctrmm_            ctrmm
#define ctrsm_            ctrsm

#define zgemm_            zgemm
#define zhemm_            zhemm
#define zsymm_            zsymm
#define zsyrk_            zsyrk
#define zherk_            zherk
#define zsyr2k_           zsyr2k
#define zher2k_           zher2k
#define ztrmm_            ztrmm
#define ztrsm_            ztrsm
                                                 /* Auxilliary PBLAS */
#define pberror_          pberror
#define pbfreebuf_        pbfreebuf

#endif
#endif





void pbchkvect( int, int, int, int, int, int, int *, int, int, int *, int *, 
	       int *) ;

void pbchkmat( int, int, int, int, int, int, int *, int, int, int *, int *, int *);
















