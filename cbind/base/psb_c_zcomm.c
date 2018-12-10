#include <stdlib.h>
#include "psb_c_zcomm.h"
#include "psb_c_zbase.h"


psb_z_t* psb_c_zvgather(psb_c_zvector *xh, psb_c_descriptor *cdh)
{ 
  psb_z_t *temp=NULL;
  psb_l_t vsize=0; 
  
  if ((vsize=psb_c_cd_get_global_rows(cdh))<0) 
    return(temp);
  
  if (vsize==0) 
    vsize=1;
  
  if ((temp=(psb_z_t *)malloc(vsize*sizeof(psb_z_t)))!=NULL)
    psb_c_zvgather_f(temp,xh,cdh);
 
 return(temp);

}


psb_c_zspmat* psb_c_zspgather(psb_c_zspmat *ah, psb_c_descriptor *cdh)
{
  psb_c_zspmat* temp=psb_c_new_zspmat();

  if (temp != NULL) 
    psb_c_zspgather_f(temp, ah, cdh);
  return(temp);
}


