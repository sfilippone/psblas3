#include <stdlib.h>
#include "psb_c_dcomm.h"

#include <stdlib.h>
#include "psb_c_dbase.h"


psb_d_t* psb_c_dvgather(psb_c_dvector *xh, psb_c_descriptor *cdh)
{ 
  psb_d_t *temp=NULL;
  psb_i_t vsize=0; 
  
  if ((vsize=psb_c_cd_get_global_rows(cdh))<0) 
    return(temp);
  
  if (vsize==0) 
    vsize=1;
  
  if ((temp=(psb_d_t *)malloc(vsize*sizeof(psb_d_t)))!=NULL)
    psb_c_dvgather_f(temp,xh,cdh);
 
 return(temp);

}


psb_c_dspmat* psb_c_dspgather(psb_c_dspmat *ah, psb_c_descriptor *cdh)
{
  psb_c_dspmat* temp=psb_c_new_dspmat();

  if (temp != NULL) 
    psb_c_dspgather_f(temp, ah, cdh);
  return(temp);
}


