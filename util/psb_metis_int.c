#include <stdio.h>
#if defined(HAVE_METIS_)
#include "metis.h"

typedef int32_t psb_m_t;

#if defined(IPK4) &&  defined(LPK4)
typedef int32_t psb_i_t;
typedef int32_t psb_l_t;
#elif defined(IPK4) &&  defined(LPK8)
typedef int32_t psb_i_t;
typedef int64_t psb_l_t;
#elif defined(IPK8) &&  defined(LPK8)
typedef int64_t psb_i_t;
typedef int64_t psb_l_t;
#else
#endif
typedef int64_t psb_e_t;

typedef float  psb_s_t;
typedef double psb_d_t;
typedef float  complex psb_c_t;
typedef double complex psb_z_t;


int metis_PartGraphKway_C(int *n, int *ixadj, int *iadj, int *ivwg, 
				int *iajw, int *nparts, float *weights, 
				int *graphpart)
{
  int res = -1;
  idx_t objval = 0;
  idx_t options[METIS_NOPTIONS];
  //printf("Inside Metis/C interface\n");
  idx_t ncon=1;
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 1;
  //printf("n:%p ncon:%p ixadj:%p iadj:%p npart:%p weights:%p options:%p objval:%p graphpart: %p\n",n,&ncon,ixadj,iadj,nparts,NULL,options,&objval,graphpart);
  /* fprintf(stderr,"From metis_int: %f\n",weights[0]); */
  if (weights[0] == -1.0) {
    res = METIS_PartGraphKway((idx_t*)n,(idx_t *)&ncon,(idx_t *)ixadj,(idx_t *)iadj,
				   NULL,NULL,NULL,(idx_t *)nparts,NULL,NULL,options,
				   &objval,(idx_t *)graphpart);
  } else {
    /* res = METIS_PartGraphKway((idx_t*)n,(idx_t *)&ncon,(idx_t *)ixadj,(idx_t *)iadj, */
    /* 				   NULL,NULL,NULL,(idx_t *)nparts,NULL,NULL,NULL, */
    /* 				   &objval,(idx_t *)graphpart); */
    res = METIS_PartGraphKway((idx_t*)n,(idx_t *)&ncon,(idx_t *)ixadj,(idx_t *)iadj,
    				   NULL,NULL,NULL,(idx_t *)nparts,weights,NULL,options,
    				   &objval,(idx_t *)graphpart);
  }
  if (res == METIS_OK) {
    return(0);
  }  else {
    return res;
  }
}


#else 

int metis_PartGraphKway_C(int *n, int *ixadj, int *iadj, int *ivwg, 
				int *iajw, int *nparts, float *weights, 
				int *graphpart)
{
  return(-1);
}
#endif
