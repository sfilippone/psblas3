#include <stdlib.h>
#include "psb_c_zprec.h"

psb_c_zprec* psb_c_new_zprec()
{
  psb_c_zprec* temp;
  
  temp=(psb_c_zprec *) malloc(sizeof(psb_c_zprec));
  temp->zprec=NULL;
  return(temp);
}

