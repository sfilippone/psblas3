#ifndef PSB_C_ZCOMM_
#define PSB_C_ZCOMM_
#include "psb_c_zbase.h"

#ifdef __cplusplus
extern "C" {
#endif

  psb_i_t       psb_c_zhalo(psb_c_zvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_zhalo_opt(psb_c_zvector *xh, psb_c_descriptor *cdh,
				char *trans, psb_i_t mode);
  psb_i_t       psb_c_zovrl(psb_c_zvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_zovrl_opt(psb_c_zvector *xh, psb_c_descriptor *cdh,
				psb_i_t update, psb_i_t mode);
  psb_i_t       psb_c_zvscatter(psb_i_t ng, psb_z_t *gx, psb_c_zvector *xh, psb_c_descriptor *cdh);

  psb_z_t*      psb_c_zvgather(psb_c_zvector *xh, psb_c_descriptor *cdh);
  psb_c_zspmat* psb_c_zspgather(psb_c_zspmat *ah, psb_c_descriptor *cdh);
  
  psb_i_t       psb_c_zvgather_f(psb_z_t* gv, psb_c_zvector *xh, psb_c_descriptor *cdh);
  psb_i_t       psb_c_zspgather_f(psb_c_zspmat* ga, psb_c_zspmat *ah, psb_c_descriptor *cdh);


#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
