#include <stdlib.h>
#include "psb_c_sbase.h"

psb_c_svector* psb_c_new_svector()
{
  psb_c_svector* temp;
  
  temp=(psb_c_svector *) malloc(sizeof(psb_c_svector));
  temp->svector=NULL;
  return(temp);
}

psb_s_t* psb_c_svect_get_cpy(psb_c_svector *xh)
{ 
  psb_s_t *temp=NULL;
  psb_i_t vsize=0; 
  
  if ((vsize=psb_c_svect_get_nrows(xh))<0) 
    return(temp);
  
  if (vsize==0) 
    vsize=1;
  
  if ((temp=(psb_s_t *)malloc(vsize*sizeof(psb_s_t)))!=NULL)
    psb_c_svect_f_get_cpy(temp,xh);
 
 return(temp);

}


psb_c_sspmat* psb_c_new_sspmat()
{
  psb_c_sspmat* temp;
  
  temp=(psb_c_sspmat *) malloc(sizeof(psb_c_sspmat));
  temp->sspmat=NULL;
  return(temp);
}
