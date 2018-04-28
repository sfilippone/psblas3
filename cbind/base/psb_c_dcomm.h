#ifndef PSB_C_DCOMM_
#define PSB_C_DCOMM_
#include "psb_c_dbase.h"

#ifdef __cplusplus
extern "C" {
#endif

  psb_i_t       psb_c_dhalo(psb_c_dvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_dhalo_opt(psb_c_dvector *xh, psb_c_descriptor *cdh,
				char *trans, psb_i_t mode);
  psb_i_t       psb_c_dovrl(psb_c_dvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_dovrl_opt(psb_c_dvector *xh, psb_c_descriptor *cdh,
				psb_i_t update, psb_i_t mode);
  psb_i_t       psb_c_dvscatter(psb_l_t ng, psb_d_t *gx, psb_c_dvector *xh, psb_c_descriptor *cdh);

  psb_d_t*      psb_c_dvgather(psb_c_dvector *xh, psb_c_descriptor *cdh);
  psb_c_dspmat* psb_c_dspgather(psb_c_dspmat *ah, psb_c_descriptor *cdh);
  
  psb_i_t       psb_c_dvgather_f(psb_d_t* gv, psb_c_dvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_dspgather_f(psb_c_dspmat* ga, psb_c_dspmat *ah, psb_c_descriptor *cdh);


#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
