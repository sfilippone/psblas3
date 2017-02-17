#include <stdlib.h>
#include <stdio.h>
#include "psb_c_base.h"

psb_c_descriptor* psb_c_new_descriptor()
{
  psb_c_descriptor* temp;
  
  temp=(psb_c_descriptor *) malloc(sizeof(psb_c_descriptor));
  temp->descriptor=NULL;
  return(temp);
}


void psb_c_print_errmsg()
{ 
  char *mesg; 
  
  for (mesg = psb_c_pop_errmsg(); mesg != NULL; mesg = psb_c_pop_errmsg()) {
    fprintf(stderr,"%s\n",mesg); 
    free(mesg);
  }

}


#define PSB_MAX_ERRLINE_LEN 132
#define PSB_MAX_ERR_LINES   4
static int maxlen=PSB_MAX_ERR_LINES*(PSB_MAX_ERRLINE_LEN+2);
char *psb_c_pop_errmsg()
{ 
  char *tmp;
  tmp = (char*) malloc(maxlen*sizeof(char));
  if (psb_c_f2c_errmsg(tmp,maxlen)<=0) {
    free(tmp); tmp = NULL;
  }
  return(tmp);
}
    
