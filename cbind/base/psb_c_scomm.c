#include <stdlib.h>
#include "psb_c_scomm.h"
#include "psb_c_sbase.h"


psb_s_t* psb_c_svgather(psb_c_svector *xh, psb_c_descriptor *cdh)
{ 
  psb_s_t *temp=NULL;
  psb_i_t vsize=0; 
  
  if ((vsize=psb_c_cd_get_global_rows(cdh))<0) 
    return(temp);
  
  if (vsize==0) 
    vsize=1;
  
  if ((temp=(psb_s_t *)malloc(vsize*sizeof(psb_s_t)))!=NULL)
    psb_c_svgather_f(temp,xh,cdh);
 
 return(temp);

}


psb_c_sspmat* psb_c_sspgather(psb_c_sspmat *ah, psb_c_descriptor *cdh)
{
  psb_c_sspmat* temp=psb_c_new_sspmat();

  if (temp != NULL) 
    psb_c_sspgather_f(temp, ah, cdh);
  return(temp);
}


