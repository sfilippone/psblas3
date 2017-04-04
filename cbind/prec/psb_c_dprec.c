#include <stdlib.h>
#include "psb_c_dprec.h"

psb_c_dprec* psb_c_new_dprec()
{
  psb_c_dprec* temp;
  
  temp=(psb_c_dprec *) malloc(sizeof(psb_c_dprec));
  temp->dprec=NULL;
  return(temp);
}

