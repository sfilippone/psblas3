#include <stdlib.h>
#include "psb_c_dbase.h"

psb_c_dvector* psb_c_new_dvector()
{
  psb_c_dvector* temp;

  temp=(psb_c_dvector *) malloc(sizeof(psb_c_dvector));
  temp->dvector=NULL;
  return(temp);
}

psb_d_t* psb_c_dvect_get_cpy(psb_c_dvector *xh)
{
  psb_d_t *temp=NULL;
  psb_i_t vsize=0;

  if ((vsize=psb_c_dvect_get_nrows(xh))<0)
    return(temp);

  if (vsize==0)
    vsize=1;

  if ((temp=(psb_d_t *)malloc(vsize*sizeof(psb_d_t)))!=NULL)
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
