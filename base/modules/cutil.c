#include <sys/time.h>
#include <stdio.h>
#include <string.h>

void psi_c_diffadd(void *p1, void *p2, int *ret)
{
  *ret = (int)((char *)p2-(char *)p1);
  return;
}
