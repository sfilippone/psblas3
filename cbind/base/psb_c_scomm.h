#ifndef PSB_C_SCOMM_
#define PSB_C_SCOMM_
#include "psb_c_sbase.h"

#ifdef __cplusplus
extern "C" {
#endif

  psb_i_t       psb_c_shalo(psb_c_svector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_shalo_opt(psb_c_svector *xh, psb_c_descriptor *cdh,
				char *trans, psb_i_t mode);
  psb_i_t       psb_c_sovrl(psb_c_svector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_sovrl_opt(psb_c_svector *xh, psb_c_descriptor *cdh,
				psb_i_t update, psb_i_t mode);
  psb_i_t       psb_c_svscatter(psb_i_t ng, psb_s_t *gx, psb_c_svector *xh, psb_c_descriptor *cdh);

  psb_s_t*      psb_c_svgather(psb_c_svector *xh, psb_c_descriptor *cdh);
  psb_c_sspmat* psb_c_sspgather(psb_c_sspmat *ah, psb_c_descriptor *cdh);
  
  psb_i_t       psb_c_svgather_f(psb_s_t* gv, psb_c_svector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_sspgather_f(psb_c_sspmat* ga, psb_c_sspmat *ah, psb_c_descriptor *cdh);


#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
