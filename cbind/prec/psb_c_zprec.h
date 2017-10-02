#ifndef PSB_C_ZPREC_
#define PSB_C_ZPREC_
#include "psb_base_cbind.h"
/* Object handle related routines */
/* Note:  psb_get_XXX_handle returns:  <= 0  unsuccessful */
/*                                     >0    valid handle */
#ifdef __cplusplus
extern "C" {
#endif

typedef struct PSB_C_ZPREC {
  void *zprec;
} psb_c_zprec; 

psb_c_zprec* psb_c_new_zprec();

psb_i_t  psb_c_zprecinit(psb_c_zprec *ph, const char *ptype);
psb_i_t  psb_c_zprecbld(psb_c_zspmat *ah, psb_c_descriptor *cdh, psb_c_zprec *ph);
psb_i_t  psb_c_zprecfree(psb_c_zprec *ph);
#ifdef __cplusplus
}
#endif

#endif
