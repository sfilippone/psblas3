#ifndef PSB_C_SBASE_
#define PSB_C_SBASE_
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
psb_s_t   *psb_c_svect_get_cpy( psb_c_svector *xh);
psb_i_t    psb_c_svect_f_get_cpy(psb_s_t *v, psb_c_svector *xh);
psb_i_t    psb_c_svect_zero(psb_c_svector *xh);
psb_s_t	  *psb_c_svect_f_get_pnt( psb_c_svector *xh);

psb_i_t    psb_c_sgeall(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_sgeins(psb_i_t nz, const psb_l_t *irw, const psb_s_t *val,
		    psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_sgeins_add(psb_i_t nz, const psb_l_t *irw, const psb_s_t *val,
			psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_sgeasb(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_sgefree(psb_c_svector *xh, psb_c_descriptor *cdh);

/* sparse matrices*/
psb_c_sspmat* psb_c_new_sspmat();
psb_i_t    psb_c_sspall(psb_c_sspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_sspasb(psb_c_sspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_sspfree(psb_c_sspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_sspins(psb_i_t nz, const psb_l_t *irw, const psb_l_t *icl,
			const psb_s_t *val, psb_c_sspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_smat_get_nrows(psb_c_sspmat *mh);
psb_i_t    psb_c_smat_get_ncols(psb_c_sspmat *mh);
psb_l_t    psb_c_snnz(psb_c_sspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_sis_matupd(psb_c_sspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_sis_matasb(psb_c_sspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_sis_matbld(psb_c_sspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_sset_matupd(psb_c_sspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_sset_matasb(psb_c_sspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_sset_matbld(psb_c_sspmat *mh,psb_c_descriptor *cdh);
psb_i_t		 psb_c_scopy_mat(psb_c_sspmat *ah,psb_c_sspmat *bh,psb_c_descriptor *cdh);

/* psb_i_t    psb_c_sspasb_opt(psb_c_sspmat *mh, psb_c_descriptor *cdh,  */
/* 			const char *afmt, psb_i_t upd, psb_i_t dupl); */
psb_i_t    psb_c_ssprn(psb_c_sspmat *mh, psb_c_descriptor *cdh, _Bool clear);
psb_i_t    psb_c_smat_name_print(psb_c_sspmat *mh, char *name);

/* psblas computational routines */
psb_s_t psb_c_sgedot(psb_c_svector *xh, psb_c_svector *yh, psb_c_descriptor *cdh);
psb_s_t psb_c_sgenrm2(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_sgeamax(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_sgeasum(psb_c_svector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_sgenrmi(psb_c_sspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_sgeaxpby(psb_s_t alpha, psb_c_svector *xh,
		       psb_s_t beta, psb_c_svector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_sgeaxpbyz(psb_s_t alpha, psb_c_svector *xh,
		       psb_s_t beta, psb_c_svector *yh, psb_c_svector *zh, psb_c_descriptor *cdh);
psb_i_t psb_c_sspmm(psb_s_t alpha, psb_c_sspmat *ah, psb_c_svector *xh,
		    psb_s_t beta, psb_c_svector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_sspmm_opt(psb_s_t alpha, psb_c_sspmat *ah, psb_c_svector *xh,
			psb_s_t beta, psb_c_svector *yh, psb_c_descriptor *cdh,
			char *trans, bool doswap);
psb_i_t psb_c_sspsm(psb_s_t alpha, psb_c_sspmat *th, psb_c_svector *xh,
		      psb_s_t beta, psb_c_svector *yh, psb_c_descriptor *cdh);
/* Additional computational routines */
psb_i_t psb_c_sgemlt(psb_c_svector *xh,psb_c_svector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_sgemlt2(psb_s_t alpha, psb_c_svector *xh, psb_c_svector *yh, psb_s_t beta, psb_c_svector *zh, psb_c_descriptor *cdh);
psb_i_t psb_c_sgediv(psb_c_svector *xh,psb_c_svector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_sgediv_check(psb_c_svector *xh,psb_c_svector *yh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_sgediv2(psb_c_svector *xh,psb_c_svector *yh,psb_c_svector *zh,psb_c_descriptor *cdh);
psb_i_t psb_c_sgediv2_check(psb_c_svector *xh,psb_c_svector *yh,psb_c_svector *zh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_sgeinv(psb_c_svector *xh,psb_c_svector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_sgeinv_check(psb_c_svector *xh,psb_c_svector *yh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_sgeabs(psb_c_svector *xh,psb_c_svector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_sgecmp(psb_c_svector *xh,psb_s_t ch,psb_c_svector *zh,psb_c_descriptor *cdh);
bool psb_c_sgecmpmat(psb_c_sspmat *ah,psb_c_sspmat *bh,psb_s_t tol,psb_c_descriptor *cdh);
bool psb_c_sgecmpmat_val(psb_c_sspmat *ah,psb_s_t val,psb_s_t tol,psb_c_descriptor *cdh);
psb_i_t psb_c_sgeaddconst(psb_c_svector *xh,psb_s_t bh,psb_c_svector *zh,psb_c_descriptor *cdh);
psb_s_t psb_c_sgenrm2_weight(psb_c_svector *xh,psb_c_svector *wh,psb_c_descriptor *cdh);
psb_s_t psb_c_sgenrm2_weightmask(psb_c_svector *xh,psb_c_svector *wh,psb_c_svector *idvh,psb_c_descriptor *cdh);
psb_i_t psb_c_smask(psb_c_svector *ch,psb_c_svector *xh,psb_c_svector *mh, void *t, psb_c_descriptor *cdh);
psb_s_t psb_c_sgemin(psb_c_svector *xh,psb_c_descriptor *cdh);
psb_i_t psb_c_sspscal(psb_s_t alpha, psb_c_sspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_sspscalpid(psb_s_t alpha, psb_c_sspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_sspaxpby(psb_s_t alpha, psb_c_sspmat *ah, psb_s_t beta, psb_c_sspmat *bh, psb_c_descriptor *cdh);
#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
