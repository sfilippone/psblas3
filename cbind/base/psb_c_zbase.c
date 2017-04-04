#include <stdlib.h>
#include "psb_c_zbase.h"

psb_c_zvector* psb_c_new_zvector()
{
  psb_c_zvector* temp;
  
  temp=(psb_c_zvector *) malloc(sizeof(psb_c_zvector));
  temp->zvector=NULL;
  return(temp);
}

psb_z_t* psb_c_zvect_get_cpy(psb_c_zvector *xh)
{ 
  psb_z_t *temp=NULL;
  psb_i_t vsize=0; 
  
  if ((vsize=psb_c_zvect_get_nrows(xh))<0) 
    return(temp);
  
  if (vsize==0) 
    vsize=1;
  
  if ((temp=(psb_z_t *)malloc(vsize*sizeof(psb_z_t)))!=NULL)
    psb_c_zvect_f_get_cpy(temp,xh);
 
 return(temp);

}


psb_c_zspmat* psb_c_new_zspmat()
{
  psb_c_zspmat* temp;
  
  temp=(psb_c_zspmat *) malloc(sizeof(psb_c_zspmat));
  temp->zspmat=NULL;
  return(temp);
}
