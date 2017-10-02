#ifndef PSB_KRYL_CBIND_
#define PSB_KRYL_CBIND_

#include "psb_base_cbind.h"
#include "psb_prec_cbind.h"


#ifdef __cplusplus
extern "C" {
#endif

/* Object handle related routines */
/* No new handles for Krylov methods. */
/* Here's a choice: define a struct to hold the options */
/* Drawback: we end up defining defaults in two places  */
/* Note: must be interoperable */ 
typedef struct psb_c_solveroptions {
  int iter;       /* On exit how many iterations were performed */
  int itmax;      /* On entry maximum number of iterations      */
  int itrace;     /* On entry print an info message every itrace iterations */
  int irst;       /* Restart depth for RGMRES or BiCGSTAB(L)    */
  int istop;      /* Stopping criterion: 1:backward error 2: ||r||_2/||b||_2 */
  double eps;     /* Stopping tolerance */ 
  double err;     /* Convergence indicator on exit */
} psb_c_SolverOptions; 

int psb_c_DefaultSolverOptions(psb_c_SolverOptions *opt);

int psb_c_skrylov(const char *method, psb_c_sspmat *ah, psb_c_sprec *ph, 
		  psb_c_svector *bh, psb_c_svector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);

int psb_c_dkrylov(const char *method, psb_c_dspmat *ah, psb_c_dprec *ph, 
		  psb_c_dvector *bh, psb_c_dvector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);

int psb_c_ckrylov(const char *method, psb_c_cspmat *ah, psb_c_cprec *ph, 
		  psb_c_cvector *bh, psb_c_cvector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);

int psb_c_zkrylov(const char *method, psb_c_zspmat *ah, psb_c_zprec *ph, 
		  psb_c_zvector *bh, psb_c_zvector *xh,
		  psb_c_descriptor *cdh, psb_c_SolverOptions *opt);

#define PSB_VALID_KRYLOV_METHODS_STRINGS "CG","CGS","BICG","BICGSTAB","RGMRES","BICGSTABL","FCG","GCR"
#define PSB_VALID_KRYLOV_METHODS_STRING  "CG CGS BICG BICGSTAB RGMRES BICGSTABL FCG GCR"

#ifdef __cplusplus
}
#endif  /* __cplusplus */
#endif
