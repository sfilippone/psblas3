#ifndef PSB_C_SUTIL_
#define PSB_C_SUTIL_
#include "psb_base_cbind.h"

#ifdef __cplusplus
extern "C" {
#endif

/* I/O Routine */
psb_i_t psb_c_smm_mat_write(psb_c_sspmat *ah, char *matrixtitle, char *filename);
psb_i_t psb_c_sglobal_mat_write(psb_c_sspmat *ah,psb_c_descriptor *cdh);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
