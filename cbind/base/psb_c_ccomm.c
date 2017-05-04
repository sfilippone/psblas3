#include <stdlib.h>
#include "psb_c_ccomm.h"
#include "psb_c_cbase.h"


psb_c_t* psb_c_cvgather(psb_c_cvector *xh, psb_c_descriptor *cdh)
{ 
  psb_c_t *temp=NULL;
  psb_i_t vsize=0; 
  
  if ((vsize=psb_c_cd_get_global_rows(cdh))<0) 
    return(temp);
  
  if (vsize==0) 
    vsize=1;
  
  if ((temp=(psb_c_t *)malloc(vsize*sizeof(psb_c_t)))!=NULL)
    psb_c_cvgather_f(temp,xh,cdh);
 
 return(temp);

}


psb_c_cspmat* psb_c_cspgather(psb_c_cspmat *ah, psb_c_descriptor *cdh)
{
  psb_c_cspmat* temp=psb_c_new_cspmat();

  if (temp != NULL) 
    psb_c_cspgather_f(temp, ah, cdh);
  return(temp);
}


