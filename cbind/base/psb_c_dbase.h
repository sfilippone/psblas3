#ifndef PSB_C_DBASE_
#define PSB_C_DBASE_
#include "psb_c_base.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PSB_C_DVECTOR {
  void *dvector;
} psb_c_dvector; 

typedef struct PSB_C_DSPMAT {
  void *dspmat;
} psb_c_dspmat; 


/* dense vectors */
psb_c_dvector* psb_c_new_dvector();
int    psb_c_dvect_get_nrows(psb_c_dvector *xh);
double *psb_c_dvect_get_cpy( psb_c_dvector *xh);
int    psb_c_dvect_f_get_cpy(double *v, psb_c_dvector *xh);
int    psb_c_dvect_zero(psb_c_dvector *xh);

int    psb_c_dgeall(psb_c_dvector *xh, psb_c_descriptor *cdh);
int    psb_c_dgeins(int nz, const int *irw, const double *val,
		    psb_c_dvector *xh, psb_c_descriptor *cdh);
int    psb_c_dgeins_add(int nz, const int *irw, const double *val,
			psb_c_dvector *xh, psb_c_descriptor *cdh);
int    psb_c_dgeasb(psb_c_dvector *xh, psb_c_descriptor *cdh);
int    psb_c_dgefree(psb_c_dvector *xh, psb_c_descriptor *cdh);

/* sparse matrices*/
psb_c_dspmat* psb_c_new_dspmat();
int    psb_c_dspall(psb_c_dspmat *mh, psb_c_descriptor *cdh);
int    psb_c_dspasb(psb_c_dspmat *mh, psb_c_descriptor *cdh);
int    psb_c_dspfree(psb_c_dspmat *mh, psb_c_descriptor *cdh);
int    psb_c_dspins(int nz, const int *irw, const int *icl, const double *val, 
		    psb_c_dspmat *mh, psb_c_descriptor *cdh);
int    psb_c_dmat_get_nrows(psb_c_dspmat *mh);
int    psb_c_dmat_get_ncols(psb_c_dspmat *mh);

/* int    psb_c_dspasb_opt(psb_c_dspmat *mh, psb_c_descriptor *cdh,  */
/* 			const char *afmt, int upd, int dupl); */
int    psb_c_dsprn(psb_c_dspmat *mh, psb_c_descriptor *cdh, _Bool clear);
/* int    psb_c_dspprint(psb_c_dspmat *mh); */

/* psblas computational routines */
double psb_c_dgedot(psb_c_dvector *xh, psb_c_dvector *yh, psb_c_descriptor *cdh);
double psb_c_dgenrm2(psb_c_dvector *xh, psb_c_descriptor *cdh);
double psb_c_dgeamax(psb_c_dvector *xh, psb_c_descriptor *cdh);
double psb_c_dgeasum(psb_c_dvector *xh, psb_c_descriptor *cdh);
double psb_c_dspnrmi(psb_c_dvector *xh, psb_c_descriptor *cdh);
int    psb_c_dgeaxpby(double alpha, psb_c_dvector *xh, 
		      double beta, psb_c_dvector *yh, psb_c_descriptor *cdh);
int    psb_c_dspmm(double alpha, psb_c_dspmat *ah, psb_c_dvector *xh, 
		   double beta, psb_c_dvector *yh, psb_c_descriptor *cdh);
int    psb_c_dspsm(double alpha, psb_c_dspmat *th, psb_c_dvector *xh, 
		   double beta, psb_c_dvector *yh, psb_c_descriptor *cdh);
#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
