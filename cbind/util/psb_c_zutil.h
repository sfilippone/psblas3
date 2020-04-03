#ifndef PSB_C_ZUTIL_
#define PSB_C_ZUTIL_
#include "psb_base_cbind.h"

#ifdef __cplusplus
extern "C" {
#endif

/* I/O Routine */
psb_i_t psb_c_zmm_mat_write(psb_c_zspmat *ah, char *matrixtitle, char *filename);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
