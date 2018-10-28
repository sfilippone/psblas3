#ifndef PSB_C_SPREC_
#define PSB_C_SPREC_
#include "psb_base_cbind.h"
/* Object handle related routines */
/* Note:  psb_get_XXX_handle returns:  <= 0  unsuccessful */
/*                                     >0    valid handle */
#ifdef __cplusplus
extern "C" {
#endif

  typedef struct PSB_C_SPREC {
    void *sprec;
  } psb_c_sprec; 
  
  psb_c_sprec* psb_c_new_sprec();
  
  psb_i_t  psb_c_sprecinit(psb_i_t ictxt, psb_c_sprec *ph, const char *ptype);
  psb_i_t  psb_c_sprecbld(psb_c_sspmat *ah, psb_c_descriptor *cdh, psb_c_sprec *ph);
  psb_i_t  psb_c_sprecfree(psb_c_sprec *ph);
#ifdef __cplusplus
}
#endif

#endif
