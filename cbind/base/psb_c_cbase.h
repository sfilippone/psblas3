#ifndef PSB_C_CBASE_
#define PSB_C_CBASE_
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
psb_c_t   *psb_c_cvect_get_cpy( psb_c_cvector *xh);
psb_i_t    psb_c_cvect_f_get_cpy(psb_c_t *v, psb_c_cvector *xh);
psb_i_t    psb_c_cvect_zero(psb_c_cvector *xh);
psb_i_t   *psb_c_cvect_f_get_pnt(psb_c_cvector *xh);

psb_i_t    psb_c_cgeall(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_cgeins(psb_i_t nz, const psb_l_t *irw, const psb_c_t *val,
		    psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_cgeins_add(psb_i_t nz, const psb_l_t *irw, const psb_c_t *val,
			psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_cgeasb(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_cgefree(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_c_t    psb_c_cgetelem(psb_c_cvector *xh,psb_l_t index,psb_c_descriptor *cd);

/* sparse matrices*/
psb_c_cspmat* psb_c_new_cspmat();
psb_i_t    psb_c_cspall(psb_c_cspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_cspasb(psb_c_cspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_cspfree(psb_c_cspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_cspins(psb_i_t nz, const psb_l_t *irw, const psb_l_t *icl,
			const psb_c_t *val, psb_c_cspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_cmat_get_nrows(psb_c_cspmat *mh);
psb_i_t    psb_c_cmat_get_ncols(psb_c_cspmat *mh);
psb_l_t    psb_c_cnnz(psb_c_cspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_cis_matupd(psb_c_cspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_cis_matasb(psb_c_cspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_cis_matbld(psb_c_cspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_cset_matupd(psb_c_cspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_cset_matasb(psb_c_cspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_cset_matbld(psb_c_cspmat *mh,psb_c_descriptor *cdh);
psb_i_t		 psb_c_ccopy_mat(psb_c_cspmat *ah,psb_c_cspmat *bh,psb_c_descriptor *cdh);

/* psb_i_t    psb_c_cspasb_opt(psb_c_cspmat *mh, psb_c_descriptor *cdh,  */
/* 			const char *afmt, psb_i_t upd, psb_i_t dupl); */
psb_i_t    psb_c_csprn(psb_c_cspmat *mh, psb_c_descriptor *cdh, _Bool clear);
psb_i_t    psb_c_cmat_name_print(psb_c_cspmat *mh, char *name);

/* psblas computational routines */
psb_c_t psb_c_cgedot(psb_c_cvector *xh, psb_c_cvector *yh, psb_c_descriptor *cdh);
psb_s_t psb_c_cgenrm2(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_cgeamax(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_cgeasum(psb_c_cvector *xh, psb_c_descriptor *cdh);
psb_s_t psb_c_cgenrmi(psb_c_cspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_cgeaxpby(psb_c_t alpha, psb_c_cvector *xh,
		       psb_c_t beta, psb_c_cvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_cgeaxpbyz(psb_c_t alpha, psb_c_cvector *xh,
 		       psb_c_t beta, psb_c_cvector *yh, psb_c_cvector *zh, psb_c_descriptor *cdh);
psb_i_t psb_c_cspmm(psb_c_t alpha, psb_c_cspmat *ah, psb_c_cvector *xh,
		    psb_c_t beta, psb_c_cvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_cspmm_opt(psb_c_t alpha, psb_c_cspmat *ah, psb_c_cvector *xh,
			psb_c_t beta, psb_c_cvector *yh, psb_c_descriptor *cdh,
			char *trans, bool doswap);
psb_i_t psb_c_cspsm(psb_c_t alpha, psb_c_cspmat *th, psb_c_cvector *xh,
		      psb_c_t beta, psb_c_cvector *yh, psb_c_descriptor *cdh);
/* Additional computational routines */
psb_i_t psb_c_cgemlt(psb_c_cvector *xh,psb_c_cvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_cgemlt2(psb_c_t alpha, psb_c_cvector *xh, psb_c_cvector *yh, psb_c_t beta, psb_c_cvector *zh, psb_c_descriptor *cdh);
psb_i_t psb_c_cgediv(psb_c_cvector *xh,psb_c_cvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_cgediv_check(psb_c_cvector *xh,psb_c_cvector *yh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_cgediv2(psb_c_cvector *xh,psb_c_cvector *yh,psb_c_cvector *zh,psb_c_descriptor *cdh);
psb_i_t psb_c_cgediv2_check(psb_c_cvector *xh,psb_c_cvector *yh,psb_c_cvector *zh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_cgeinv(psb_c_cvector *xh,psb_c_cvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_cgeinv_check(psb_c_cvector *xh,psb_c_cvector *yh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_cgeabs(psb_c_cvector *xh,psb_c_cvector *yh,psb_c_cvector *cdh);
psb_i_t psb_c_cgecmp(psb_c_cvector *xh,psb_s_t ch,psb_c_cvector *zh,psb_c_descriptor *cdh);
bool psb_c_cgecmpmat(psb_c_cspmat *ah,psb_c_cspmat *bh,psb_s_t tol,psb_c_descriptor *cdh);
bool psb_c_cgecmpmat_val(psb_c_cspmat *ah,psb_c_t val,psb_s_t tol,psb_c_descriptor *cdh);
psb_i_t psb_c_cgeaddconst(psb_c_cvector *xh,psb_c_t bh,psb_c_cvector *zh,psb_c_descriptor *cdh);
psb_s_t psb_c_cgenrm2_weight(psb_c_cvector *xh,psb_c_cvector *wh,psb_c_descriptor *cdh);
psb_s_t psb_c_cgenrm2_weightmask(psb_c_cvector *xh,psb_c_cvector *wh,psb_c_cvector *idvh,psb_c_descriptor *cdh);
psb_i_t psb_c_cspscal(psb_c_t alpha, psb_c_cspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_cspscalpid(psb_c_t alpha, psb_c_cspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_cspaxpby(psb_c_t alpha, psb_c_cspmat *ah, psb_c_t beta, psb_c_cspmat *bh, psb_c_descriptor *cdh);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
