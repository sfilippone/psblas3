#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>

double wtime() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec)*1.0e3 + ((double)tt.tv_usec)*1.0e-3;
  return(temp);
}

double timef_() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec)*1.0e3 + ((double)tt.tv_usec)*1.0e-3;
  return(temp);
}

double timef() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec)*1.0e3 + ((double)tt.tv_usec)*1.0e-3;
  return(temp);
}

double etime() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec) + ((double)tt.tv_usec)*1.0e-6;
  return(temp);
}

double etime_() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec) + ((double)tt.tv_usec)*1.0e-6;
  return(temp);
}

double etimef() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec) + ((double)tt.tv_usec)*1.0e-6;
  return(temp);
}

double etimef_() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec) + ((double)tt.tv_usec)*1.0e-6;
  return(temp);
}



