#ifndef PSB_C_DCOMM_
#define PSB_C_DCOMM_
#include "psb_c_dbase.h"

#ifdef __cplusplus
extern "C" {
#endif

  psb_i_t       psb_c_dvscatter(psb_c_dvector *xh, psb_c_descriptor *cdh);
  psb_d_t*      psb_c_dvgather(psb_c_dvector *xh, psb_c_descriptor *cdh);
  psb_c_dspmat* psb_c_dspgather(psb_c_dspmat *ah, psb_c_descriptor *cdh);


#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
