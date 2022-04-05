#ifndef PSB_C_ZUTIL_
#define PSB_C_ZUTIL_
#include "psb_base_cbind.h"

#ifdef __cplusplus
extern "C" {
#endif

/* I/O Routine */
psb_i_t psb_c_zmm_mat_write(psb_c_zspmat *ah, char *matrixtitle, char *filename);
psb_i_t psb_c_zglobal_mat_write(psb_c_zspmat *ah,psb_c_descriptor *cdh);
psb_i_t psb_c_zglobal_vec_write(psb_c_zvector *vh,psb_c_descriptor *cdh);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
