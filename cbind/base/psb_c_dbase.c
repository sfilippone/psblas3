#include <stdlib.h>
#include "psb_c_dbase.h"

psb_c_dvector* psb_c_new_dvector()
{
  psb_c_dvector* temp;
  
  temp=(psb_c_dvector *) malloc(sizeof(psb_c_dvector));
  temp->dvector=NULL;
  return(temp);
}

double* psb_c_dvect_get_cpy(psb_c_dvector *xh)
{ 
  double *temp=NULL;
  int vsize=0; 
  
  if ((vsize=psb_c_dvect_get_nrows(xh))<0) 
    return(temp);
  
  if (vsize==0) 
    vsize=1;
  
  if ((temp=(double *)malloc(vsize*sizeof(double)))!=NULL)
    psb_c_dvect_f_get_cpy(temp,xh);
 
 return(temp);

}


psb_c_dspmat* psb_c_new_dspmat()
{
  psb_c_dspmat* temp;
  
  temp=(psb_c_dspmat *) malloc(sizeof(psb_c_dspmat));
  temp->dspmat=NULL;
  return(temp);
}
