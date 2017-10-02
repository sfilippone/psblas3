#include <stdlib.h>
#include "psb_c_sprec.h"

psb_c_sprec* psb_c_new_sprec()
{
  psb_c_sprec* temp;
  
  temp=(psb_c_sprec *) malloc(sizeof(psb_c_sprec));
  temp->sprec=NULL;
  return(temp);
}

