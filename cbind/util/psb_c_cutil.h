#ifndef PSB_C_CUTIL_
#define PSB_C_CUTIL_
#include "psb_base_cbind.h"

#ifdef __cplusplus
extern "C" {
#endif

/* I/O Routine */
psb_i_t psb_c_cmm_mat_write(psb_c_cspmat *ah, char *matrixtitle, char *filename);
psb_i_t psb_c_cglobal_mat_write(psb_c_cspmat *ah,psb_c_descriptor *cdh);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
