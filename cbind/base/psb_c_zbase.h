#ifndef PSB_C_ZBASE_
#define PSB_C_ZBASE_
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
psb_z_t   *psb_c_zvect_get_cpy( psb_c_zvector *xh);
psb_i_t    psb_c_zvect_f_get_cpy(psb_z_t *v, psb_c_zvector *xh);
psb_i_t    psb_c_zvect_zero(psb_c_zvector *xh);
psb_z_t	  *psb_c_zvect_f_get_pnt( psb_c_zvector *xh);

psb_i_t    psb_c_zgeall(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zgeall_remote(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zgeins(psb_i_t nz, const psb_l_t *irw, const psb_z_t *val,
		    psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zgeins_add(psb_i_t nz, const psb_l_t *irw, const psb_z_t *val,
			psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zgeasb(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zgefree(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_z_t    psb_c_zgetelem(psb_c_zvector *xh,psb_l_t index,psb_c_descriptor *cd);

/* sparse matrices*/
psb_c_zspmat* psb_c_new_zspmat();
psb_i_t    psb_c_zspall(psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zspall_remote(psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zspasb(psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zspfree(psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zspins(psb_i_t nz, const psb_l_t *irw, const psb_l_t *icl,
			const psb_z_t *val, psb_c_zspmat *mh, psb_c_descriptor *cdh);
psb_i_t    psb_c_zmat_get_nrows(psb_c_zspmat *mh);
psb_i_t    psb_c_zmat_get_ncols(psb_c_zspmat *mh);
psb_l_t    psb_c_znnz(psb_c_zspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_zis_matupd(psb_c_zspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_zis_matasb(psb_c_zspmat *mh,psb_c_descriptor *cdh);
bool    	 psb_c_zis_matbld(psb_c_zspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_zset_matupd(psb_c_zspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_zset_matasb(psb_c_zspmat *mh,psb_c_descriptor *cdh);
psb_i_t    psb_c_zset_matbld(psb_c_zspmat *mh,psb_c_descriptor *cdh);
psb_i_t		 psb_c_zcopy_mat(psb_c_zspmat *ah,psb_c_zspmat *bh,psb_c_descriptor *cdh);


/* psb_i_t    psb_c_zspasb_opt(psb_c_zspmat *mh, psb_c_descriptor *cdh,  */
/* 			const char *afmt, psb_i_t upd, psb_i_t dupl); */
psb_i_t    psb_c_zsprn(psb_c_zspmat *mh, psb_c_descriptor *cdh, _Bool clear);
psb_i_t    psb_c_zmat_name_print(psb_c_zspmat *mh, char *name);
psb_i_t		 psb_c_zvect_set_scal(psb_c_zvector *xh, psb_z_t val);
psb_i_t		 psb_c_zvect_set_vect(psb_c_zvector *xh, psb_z_t *val, psb_i_t n);

/* psblas computational routines */
psb_z_t psb_c_zgedot(psb_c_zvector *xh, psb_c_zvector *yh, psb_c_descriptor *cdh);
psb_d_t psb_c_zgenrm2(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_zgeamax(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_zgeasum(psb_c_zvector *xh, psb_c_descriptor *cdh);
psb_d_t psb_c_zgenrmi(psb_c_zspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_zgeaxpby(psb_z_t alpha, psb_c_zvector *xh,
		       psb_z_t beta, psb_c_zvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_zgeaxpbyz(psb_z_t alpha, psb_c_zvector *xh,
		       psb_z_t beta, psb_c_zvector *yh, psb_c_zvector *zh, psb_c_descriptor *cdh);
psb_i_t psb_c_zspmm(psb_z_t alpha, psb_c_zspmat *ah, psb_c_zvector *xh,
		    psb_z_t beta, psb_c_zvector *yh, psb_c_descriptor *cdh);
psb_i_t psb_c_zspmm_opt(psb_z_t alpha, psb_c_zspmat *ah, psb_c_zvector *xh,
			psb_z_t beta, psb_c_zvector *yh, psb_c_descriptor *cdh,
			char *trans, bool doswap);
psb_i_t psb_c_zspsm(psb_z_t alpha, psb_c_zspmat *th, psb_c_zvector *xh,
		      psb_z_t beta, psb_c_zvector *yh, psb_c_descriptor *cdh);
/* Additional computational routines */
psb_i_t psb_c_zgemlt(psb_c_zvector *xh,psb_c_zvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_zgemlt2(psb_z_t alpha, psb_c_zvector *xh, psb_c_zvector *yh, psb_z_t beta, psb_c_zvector *zh, psb_c_descriptor *cdh);
psb_i_t psb_c_zgediv(psb_c_zvector *xh,psb_c_zvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_zgediv_check(psb_c_zvector *xh,psb_c_zvector *yh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_zgediv2(psb_c_zvector *xh,psb_c_zvector *yh,psb_c_zvector *zh,psb_c_descriptor *cdh);
psb_i_t psb_c_zgediv2_check(psb_c_zvector *xh,psb_c_zvector *yh,psb_c_zvector *zh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_zgeinv(psb_c_zvector *xh,psb_c_zvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_zgeinv_check(psb_c_zvector *xh,psb_c_zvector *yh,psb_c_descriptor *cdh, bool flag);
psb_i_t psb_c_zgeabs(psb_c_zvector *xh,psb_c_zvector *yh,psb_c_descriptor *cdh);
psb_i_t psb_c_zgecmp(psb_c_zvector *xh,psb_d_t ch,psb_c_zvector *zh,psb_c_descriptor *cdh);
bool psb_c_zgecmpmat(psb_c_zspmat *ah,psb_c_zspmat *bh,psb_d_t tol,psb_c_descriptor *cdh);
bool psb_c_zgecmpmat_val(psb_c_zspmat *ah,psb_z_t val,psb_d_t tol,psb_c_descriptor *cdh);
psb_i_t psb_c_zgeaddconst(psb_c_zvector *xh,psb_z_t bh,psb_c_zvector *zh,psb_c_descriptor *cdh);
psb_d_t psb_c_zgenrm2_weight(psb_c_zvector *xh,psb_c_zvector *wh,psb_c_descriptor *cdh);
psb_d_t psb_c_zgenrm2_weightmask(psb_c_zvector *xh,psb_c_zvector *wh,psb_c_zvector *idvh,psb_c_descriptor *cdh);
psb_i_t psb_c_zspscal(psb_z_t alpha, psb_c_zspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_zspscalpid(psb_z_t alpha, psb_c_zspmat *ah, psb_c_descriptor *cdh);
psb_i_t psb_c_zspaxpby(psb_z_t alpha, psb_c_zspmat *ah, psb_z_t beta, psb_c_zspmat *bh, psb_c_descriptor *cdh);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
