#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include "psb_internals.h"

#ifdef LowerUnderscore
#define psi_c_diffadd   psi_c_diffadd_
#endif
#ifdef LowerDoubleUnderscore
#define psi_c_diffadd   psi_c_diffadd__
#endif
#ifdef LowerCase
#define psi_c_diffadd   psi_c_diffadd
#endif
#ifdef UpperUnderscore
#define psi_c_diffadd   PSI_C_DIFFADD_
#endif
#ifdef UpperDoubleUnderscore 
#define psi_c_diffadd   PSI_C_DIFFADD__
#endif
#ifdef UpperCase
#define psi_c_diffadd   PSI_C_DIFFADD
#endif

void psi_c_diffadd(char *p1, char *p2, int *ret)
{
  *ret = (int)(p2-p1);
  return;
}
