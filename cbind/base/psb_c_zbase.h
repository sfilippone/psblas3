#ifndef PSB_C_DBASE_
#define PSB_C_DBASE_
#include "psb_c_base.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct PSB_C_ZVECTOR {
  void *zvector;
} psb_c_zvector; 

typedef struct PSB_C_ZSPMAT {
  void *zspmat;
} psb_c_zspmat; 


/* dense vectors */
psb_c_zvector* psb_c_new_zvector();
psb_i_t    psb_c_zvect_get_nrows(psb_c_zvector *xh);
psb_z_t *psb_c_zvect_get_cpy( psb_c_zvector *xh);
psb_i_t    psb_c_zvect_f_get_cpy(psb_z_t *v, psb_c_zvector *xh);
psb_i_t    psb_c_zvect_zero(psb_c_zvector *xh);

psb_i_t    psb_c_dgeall(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeins(psb_i_t nz, const psb_i_t *irw, const psb_z_t *val,
		    psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeins_add(psb_i_t nz, const psb_i_t *irw, const psb_z_t *val,
			psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeasb(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgefree(psb_c_zvector *xh, psb_c_descriptor *cdh);

/* sparse matrices*/
psb_c_zspmat* psb_c_new_zspmat();
psb_i_t    psb_c_dspall(psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspasb(psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspfree(psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspins(psb_i_t nz, const psb_i_t *irw, const psb_i_t *icl, const psb_z_t *val, 
		    psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dmat_get_nrows(psb_c_zspmat *mh);
psb_i_t    psb_c_dmat_get_ncols(psb_c_zspmat *mh);

/* psb_i_t    psb_c_dspasb_opt(psb_c_zspmat *mh, psb_c_descriptor *cdh,  */
/* 			const char *afmt, psb_i_t upd, psb_i_t dupl); */
psb_i_t    psb_c_dsprn(psb_c_zspmat *mh, psb_c_descriptor *cdh, _Bool clear);
/* psb_i_t    psb_c_dspprint(psb_c_zspmat *mh); */

/* psblas computational routines */
psb_z_t psb_c_dgedot(psb_c_zvector *xh, psb_c_zvector *yh, psb_c_descriptor *cdh);
psb_d_t psb_c_dgenrm2(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_dgeamax(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_dgeasum(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_dspnrmi(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t psb_c_dgeaxpby(psb_z_t alpha, psb_c_zvector *xh, 
		       psb_z_t beta, psb_c_zvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_zspmm(psb_z_t alpha, psb_c_zspmat *ah, psb_c_zvector *xh, 
		    psb_z_t beta, psb_c_zvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_dspsm(psb_z_t alpha, psb_c_zspmat *th, psb_c_zvector *xh, 
		      psb_z_t beta, psb_c_zvector *yh, psb_c_descriptor *cdh);
#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
