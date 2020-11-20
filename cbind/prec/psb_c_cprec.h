#ifndef PSB_C_CPREC_
#define PSB_C_CPREC_
#include "psb_base_cbind.h"
/* Object handle related routines */
/* Note:  psb_get_XXX_handle returns:  <= 0  unsuccessful */
/*                                     >0    valid handle */
#ifdef __cplusplus
extern "C" {
#endif

  typedef struct PSB_C_CPREC {
    void *cprec;
  } psb_c_cprec; 
  
  psb_c_cprec* psb_c_new_cprec(); 
  
  psb_i_t  psb_c_cprecinit(psb_c_ctxt cctxt,psb_c_cprec *ph, const char *ptype);
  psb_i_t  psb_c_cprecbld(psb_c_cspmat *ah, psb_c_descriptor *cdh, psb_c_cprec *ph);
  psb_i_t  psb_c_cprecfree(psb_c_cprec *ph);
#ifdef __cplusplus
}
#endif

#endif
