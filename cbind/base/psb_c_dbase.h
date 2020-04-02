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
psb_i_t    psb_c_dvect_get_nrows(psb_c_dvector *xh);
psb_d_t   *psb_c_dvect_get_cpy( psb_c_dvector *xh);
psb_i_t    psb_c_dvect_f_get_cpy(psb_d_t *v, psb_c_dvector *xh);
psb_i_t    psb_c_dvect_zero(psb_c_dvector *xh);
psb_d_t	  *psb_c_dvect_f_get_pnt( psb_c_dvector *xh);

psb_i_t    psb_c_dgeall(psb_c_dvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeins(psb_i_t nz, const psb_l_t *irw, const psb_d_t *val,
		    psb_c_dvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeins_add(psb_i_t nz, const psb_l_t *irw, const psb_d_t *val,
			psb_c_dvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgeasb(psb_c_dvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dgefree(psb_c_dvector *xh, psb_c_descriptor *cdh);

/* sparse matrices*/
psb_c_dspmat* psb_c_new_dspmat();
psb_i_t    psb_c_dspall(psb_c_dspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspasb(psb_c_dspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspfree(psb_c_dspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dspins(psb_i_t nz, const psb_l_t *irw, const psb_l_t *icl,
			const psb_d_t *val, psb_c_dspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_dmat_get_nrows(psb_c_dspmat *mh);
psb_i_t    psb_c_dmat_get_ncols(psb_c_dspmat *mh);

/* psb_i_t    psb_c_dspasb_opt(psb_c_dspmat *mh, psb_c_descriptor *cdh,  */
/* 			const char *afmt, psb_i_t upd, psb_i_t dupl); */
psb_i_t    psb_c_dsprn(psb_c_dspmat *mh, psb_c_descriptor *cdh, _Bool clear);
psb_i_t    psb_c_dmat_name_print(psb_c_dspmat *mh, char *name);

/* psblas computational routines */
psb_d_t psb_c_dgedot(psb_c_dvector *xh, psb_c_dvector *yh, psb_c_descriptor *cdh);
psb_d_t psb_c_dgenrm2(psb_c_dvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_dgeamax(psb_c_dvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_dgeasum(psb_c_dvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_dgenrmi(psb_c_dvector *xh, psb_c_descriptor *cdh);
psb_i_t psb_c_dgeaxpby(psb_d_t alpha, psb_c_dvector *xh,
		       psb_d_t beta, psb_c_dvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_dgeaxpbyz(psb_d_t alpha, psb_c_dvector *xh,
					 		       psb_d_t beta, psb_c_dvector *yh, psb_c_dvector *zh, psb_c_descriptor *cdh);
psb_i_t psb_c_dspmm(psb_d_t alpha, psb_c_dspmat *ah, psb_c_dvector *xh,
		    psb_d_t beta, psb_c_dvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_dspmm_opt(psb_d_t alpha, psb_c_dspmat *ah, psb_c_dvector *xh,
			psb_d_t beta, psb_c_dvector *yh, psb_c_descriptor *cdh,
			char *trans, bool doswap);
psb_i_t psb_c_dspsm(psb_d_t alpha, psb_c_dspmat *th, psb_c_dvector *xh,
		      psb_d_t beta, psb_c_dvector *yh, psb_c_descriptor *cdh);
/* Additional computational routines */
psb_i_t psb_c_dgemlt(psb_c_dvector *xh,psb_c_dvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_dgemlt2(psb_d_t alpha, psb_c_dvector *xh, psb_c_dvector *yh, psb_d_t beta, psb_c_dvector *zh, psb_c_descriptor *cdh);
psb_i_t psb_c_dgediv(psb_c_dvector *xh,psb_c_dvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_dgediv_check(psb_c_dvector *xh,psb_c_dvector *yh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_dgediv2(psb_c_dvector *xh,psb_c_dvector *yh,psb_c_dvector *zh,psb_c_descriptor *cdh);
psb_i_t psb_c_dgediv2_check(psb_c_dvector *xh,psb_c_dvector *yh,psb_c_dvector *zh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_dgeinv(psb_c_dvector *xh,psb_c_dvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_dgeinv_check(psb_c_dvector *xh,psb_c_dvector *yh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_dgeabs(psb_c_dvector *xh,psb_c_dvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_dgecmp(psb_c_dvector *xh,psb_d_t ch,psb_c_dvector *zh,psb_c_descriptor *cdh);
psb_i_t psb_c_dgeaddconst(psb_c_dvector *xh,psb_d_t bh,psb_c_dvector *zh,psb_c_descriptor *cdh);
psb_d_t psb_c_dgenrm2_weight(psb_c_dvector *xh,psb_c_dvector *wh,psb_c_descriptor *cdh);
psb_d_t psb_c_dgenrm2_weightmask(psb_c_dvector *xh,psb_c_dvector *wh,psb_c_dvector *idvh,psb_c_descriptor *cdh);
psb_i_t psb_c_dmask(psb_c_dvector *ch,psb_c_dvector *xh,psb_c_dvector *mh, void *t, psb_c_descriptor *cdh);
psb_d_t psb_c_dgemin(psb_c_dvector *xh,psb_c_descriptor *cdh);
psb_d_t psb_c_dminquotient(psb_c_dvector *xh,psb_c_dvector *yh, psb_c_descriptor *cdh);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
