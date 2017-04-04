#ifndef PSB_C_DPREC_
#define PSB_C_DPREC_
#include "psb_base_cbind.h"
/* Object handle related routines */
/* Note:  psb_get_XXX_handle returns:  <= 0  unsuccessful */
/*                                     >0    valid handle */
#ifdef __cplusplus
extern "C" {
#endif

typedef struct PSB_C_DPREC {
  void *dprec;
} psb_c_dprec; 

psb_c_dprec* psb_c_new_dprec();

int  psb_c_dprecinit(psb_c_dprec *ph, const char *ptype);
int  psb_c_dprecbld(psb_c_dspmat *ah, psb_c_descriptor *cdh, psb_c_dprec *ph);
int  psb_c_dprecfree(psb_c_dprec *ph);
#ifdef __cplusplus
}
#endif

#endif
