#include <stdio.h>
#if defined(HAVE_METIS_)
#include "metis.h"

/* extern int METIS_PartGraphKway(int *, int *, int *, int *, int *, int *, int *, int *, float *, float, int *, int *, int *); */


int metis_PartGraphKway_C(int *n, int *ixadj, int *iadj, int *ivwg, 
				int *iajw, int *nparts, float *weights, 
				int *graphpart)
{
  int res = -1;
#if defined(METIS_5) 
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
#elif defined(METIS_4)
  idxtype objval = 0;
  int options[8];
  //printf("Inside Metis/C interface\n");
  idxtype ncon=1;
  int wflag=0;
  int numflag=1;
  int ecut;
  options[0]=0;
  METIS_PartGraphKway((int *)n,(idxtype *)ixadj,(idxtype *)iadj,
		      NULL,NULL,&wflag,&numflag,nparts,options,
		      &ecut,(idxtype *)graphpart);
  return(0); 
    
#else
  choke on me!
#endif
}


#else 

int metis_PartGraphKway_C(int *n, int *ixadj, int *iadj, int *ivwg, 
				int *iajw, int *nparts, float *weights, 
				int *graphpart)
{
  return(-1);
}
#endif
