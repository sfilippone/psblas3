#ifndef PSB_C_DBASE_
#define PSB_C_DBASE_
#include "psb_c_base.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PSB_C_SVECTOR {
  void *svector;
} psb_c_svector; 

typedef struct PSB_C_SSPMAT {
  void *sspmat;
} psb_c_sspmat; 


/* dense vectors */
psb_c_svector* psb_c_new_svector();
psb_i_t    psb_c_svect_get_nrows(psb_c_svector *xh);
psb_s_t *psb_c_svect_get_cpy( psb_c_svector *xh);
psb_i_t    psb_c_svect_f_get_cpy(psb_s_t *v, psb_c_svector *xh);
psb_i_t    psb_c_svect_zero(psb_c_svector *xh);

psb_i_t    psb_c_dgeall(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeins(psb_i_t nz, const psb_i_t *irw, const psb_s_t *val,
		    psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeins_add(psb_i_t nz, const psb_i_t *irw, const psb_s_t *val,
			psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeasb(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgefree(psb_c_svector *xh, psb_c_descriptor *cdh);

/* sparse matrices*/
psb_c_sspmat* psb_c_new_sspmat();
psb_i_t    psb_c_dspall(psb_c_sspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspasb(psb_c_sspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspfree(psb_c_sspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspins(psb_i_t nz, const psb_i_t *irw, const psb_i_t *icl, const psb_s_t *val, 
		    psb_c_sspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dmat_get_nrows(psb_c_sspmat *mh);
psb_i_t    psb_c_dmat_get_ncols(psb_c_sspmat *mh);

/* psb_i_t    psb_c_dspasb_opt(psb_c_sspmat *mh, psb_c_descriptor *cdh,  */
/* 			const char *afmt, psb_i_t upd, psb_i_t dupl); */
psb_i_t    psb_c_dsprn(psb_c_sspmat *mh, psb_c_descriptor *cdh, _Bool clear);
/* psb_i_t    psb_c_dspprint(psb_c_sspmat *mh); */

/* psblas computational routines */
psb_s_t psb_c_dgedot(psb_c_svector *xh, psb_c_svector *yh, psb_c_descriptor *cdh);
psb_s_t psb_c_dgenrm2(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_dgeamax(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_dgeasum(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_dspnrmi(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t psb_c_dgeaxpby(psb_s_t alpha, psb_c_svector *xh, 
		       psb_s_t beta, psb_c_svector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_sspmm(psb_s_t alpha, psb_c_sspmat *ah, psb_c_svector *xh, 
		    psb_s_t beta, psb_c_svector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_dspsm(psb_s_t alpha, psb_c_sspmat *th, psb_c_svector *xh, 
		      psb_s_t beta, psb_c_svector *yh, psb_c_descriptor *cdh);
#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
