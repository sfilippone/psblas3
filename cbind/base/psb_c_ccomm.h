#ifndef PSB_C_CCOMM_
#define PSB_C_CCOMM_
#include "psb_c_cbase.h"

#ifdef __cplusplus
extern "C" {
#endif

  psb_i_t       psb_c_chalo(psb_c_cvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_chalo_opt(psb_c_cvector *xh, psb_c_descriptor *cdh,
				char *trans, psb_i_t mode);
  psb_i_t       psb_c_covrl(psb_c_cvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_covrl_opt(psb_c_cvector *xh, psb_c_descriptor *cdh,
				psb_i_t update, psb_i_t mode);
  psb_i_t       psb_c_cvscatter(psb_l_t ng, psb_c_t *gx, psb_c_cvector *xh, psb_c_descriptor *cdh);

  psb_c_t*      psb_c_cvgather(psb_c_cvector *xh, psb_c_descriptor *cdh);
  psb_c_cspmat* psb_c_cspgather(psb_c_cspmat *ah, psb_c_descriptor *cdh);
  
  psb_i_t       psb_c_cvgather_f(psb_c_t* gv, psb_c_cvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_cspgather_f(psb_c_cspmat* ga, psb_c_cspmat *ah, psb_c_descriptor *cdh);


#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
