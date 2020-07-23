#include <stdio.h>
#if defined(HAVE_METIS_)
#include "psb_metis_int.h"

int metis_PartGraphKway_C(idx_t *n, idx_t *ixadj, idx_t *iadj, idx_t *ivwg, 
				idx_t *iajw, idx_t *nparts, float *weights, 
				idx_t *graphpart)
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

int metis_PartGraphKway_C(idx_t *n, idx_t *ixadj, idx_t *iadj, idx_t *ivwg, 
				idx_t *iajw, idx_t *nparts, float *weights, 
				idx_t *graphpart)
{
  return(-1);
}
#endif
