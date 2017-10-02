#include <stdlib.h>
#include "psb_c_cprec.h"

psb_c_cprec* psb_c_new_cprec()
{
  psb_c_cprec* temp;
  
  temp=(psb_c_cprec *) malloc(sizeof(psb_c_cprec));
  temp->cprec=NULL;
  return(temp);
}

