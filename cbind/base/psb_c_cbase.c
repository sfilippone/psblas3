#include <stdlib.h>
#include "psb_c_cbase.h"

psb_c_cvector* psb_c_new_cvector()
{
  psb_c_cvector* temp;
  
  temp=(psb_c_cvector *) malloc(sizeof(psb_c_cvector));
  temp->cvector=NULL;
  return(temp);
}

psb_c_t* psb_c_cvect_get_cpy(psb_c_cvector *xh)
{ 
  psb_c_t *temp=NULL;
  psb_i_t vsize=0; 
  
  if ((vsize=psb_c_cvect_get_nrows(xh))<0) 
    return(temp);
  
  if (vsize==0) 
    vsize=1;
  
  if ((temp=(psb_c_t *)malloc(vsize*sizeof(psb_c_t)))!=NULL)
    psb_c_cvect_f_get_cpy(temp,xh);
 
 return(temp);

}


psb_c_cspmat* psb_c_new_cspmat()
{
  psb_c_cspmat* temp;
  
  temp=(psb_c_cspmat *) malloc(sizeof(psb_c_cspmat));
  temp->cspmat=NULL;
  return(temp);
}
