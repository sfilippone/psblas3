#ifndef PSB_C_DBASE_
#define PSB_C_DBASE_
#include "psb_c_base.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PSB_C_CVECTOR {
  void *cvector;
} psb_c_cvector; 

typedef struct PSB_C_CSPMAT {
  void *cspmat;
} psb_c_cspmat; 


/* dense vectors */
psb_c_cvector* psb_c_new_cvector();
psb_i_t    psb_c_cvect_get_nrows(psb_c_cvector *xh);
psb_c_t *psb_c_cvect_get_cpy( psb_c_cvector *xh);
psb_i_t    psb_c_cvect_f_get_cpy(psb_c_t *v, psb_c_cvector *xh);
psb_i_t    psb_c_cvect_zero(psb_c_cvector *xh);

psb_i_t    psb_c_dgeall(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeins(psb_i_t nz, const psb_i_t *irw, const psb_c_t *val,
		    psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeins_add(psb_i_t nz, const psb_i_t *irw, const psb_c_t *val,
			psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeasb(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgefree(psb_c_cvector *xh, psb_c_descriptor *cdh);

/* sparse matrices*/
psb_c_cspmat* psb_c_new_cspmat();
psb_i_t    psb_c_dspall(psb_c_cspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspasb(psb_c_cspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspfree(psb_c_cspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspins(psb_i_t nz, const psb_i_t *irw, const psb_i_t *icl, const psb_c_t *val, 
		    psb_c_cspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dmat_get_nrows(psb_c_cspmat *mh);
psb_i_t    psb_c_dmat_get_ncols(psb_c_cspmat *mh);

/* psb_i_t    psb_c_dspasb_opt(psb_c_cspmat *mh, psb_c_descriptor *cdh,  */
/* 			const char *afmt, psb_i_t upd, psb_i_t dupl); */
psb_i_t    psb_c_dsprn(psb_c_cspmat *mh, psb_c_descriptor *cdh, _Bool clear);
/* psb_i_t    psb_c_dspprint(psb_c_cspmat *mh); */

/* psblas computational routines */
psb_c_t psb_c_dgedot(psb_c_cvector *xh, psb_c_cvector *yh, psb_c_descriptor *cdh);
psb_s_t psb_c_dgenrm2(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_dgeamax(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_dgeasum(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_dspnrmi(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t psb_c_dgeaxpby(psb_c_t alpha, psb_c_cvector *xh, 
		       psb_c_t beta, psb_c_cvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_cspmm(psb_c_t alpha, psb_c_cspmat *ah, psb_c_cvector *xh, 
		    psb_c_t beta, psb_c_cvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_dspsm(psb_c_t alpha, psb_c_cspmat *th, psb_c_cvector *xh, 
		      psb_c_t beta, psb_c_cvector *yh, psb_c_descriptor *cdh);
#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
