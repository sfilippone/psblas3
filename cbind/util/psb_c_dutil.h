#ifndef PSB_C_DUTIL_
#define PSB_C_DUTIL_
#include "psb_base_cbind.h"

#ifdef __cplusplus
extern "C" {
#endif

/* I/O Routine */
psb_i_t psb_c_dmm_mat_write(psb_c_dspmat *ah, char *matrixtitle, char *filename);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
