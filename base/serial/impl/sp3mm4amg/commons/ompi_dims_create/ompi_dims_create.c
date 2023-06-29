/*
 * Copyright (c) 2004-2007 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2020 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2014 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2012      Los Alamos National Security, LLC.  All rights
 *                         reserved.
 * Copyright (c) 2014      Intel, Inc. All rights reserved
 * Copyright (c) 2015 Cisco Systems, Inc.  All rights reserved.
 * Copyright (c) 2015-2016 Research Organization for Information Science
 *                         and Technology (RIST). All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 *
 * extraction from openMPI, mocking not needed dependencies for MPI_Dims_create function export: 
 * Copyright (c) 2022 		Andrea Di Iorio
 */

//#include "ompi_config.h"
#include "ompi_config_minimal.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <math.h>

//#include "ompi/mpi/c/bindings.h"
//#include "ompi/runtime/params.h"
//#include "ompi/communicator/communicator.h"
//#include "ompi/errhandler/errhandler.h"


static const char FUNC_NAME[] = "MPI_Dims_create";

/* static functions */
static int assignnodes(int ndim, int nfactor, int *pfacts,int **pdims);
static int getfactors(int num, int *nfators, int **factors);


/*
 * This is a utility function, no need to have anything in the lower
 * layer for this at all
 * ****
 *	original manpage ... https://www.open-mpi.org/doc/v3.1/man3/MPI_Dims_create.3.php
 *  For Cartesian topologies, the function MPI_Dims_create helps the user select a balanced distribution of processes per coordinate direction, 
 *  depending on the number of processes in the group to be balanced and optional constraints that can be specified by the user. 
 *  One use is to partition all the processes (the size of MPI_COMM_WORLD’s group) into an n-dimensional topology.
 *	The entries in the array dims are set to describe a Cartesian grid with ndims dimensions and a total of nnodes nodes. 
 *	The dimensions are set to be as close to each other as possible, using an appropriate divisibility algorithm. 
 *
 *	The caller may further constrain the operation of this routine by specifying elements of array dims. 
 *	If dims[i] is set to a positive number, the routine will not modify the number of nodes in dimension i; only those entries where dims[i] = 0 are modified by the call.
 *	Negative input values of dims[i] are erroneous. An error will occur if nnodes is not a multiple of ((pi) over (i, dims[i] != 0)) dims[i].
 *	For dims[i] set by the call, dims[i] will be ordered in nonincreasing order. Array dims is suitable for use as input to routine MPI_Cart_create. MPI_Dims_create is local.
 *	Example:
 *	
 *	dims
 *	before                    dims
 *	call        function call        on return
 *	-----------------------------------------------------
 *	(0,0)    MPI_Dims_create(6, 2, dims)    (3,2)
 *	(0,0)    MPI_Dims_create(7, 2, dims)     (7,1)
 *	(0,3,0)    MPI_Dims_create(6, 3, dims)    (2,3,1)
 *	(0,3,0)    MPI_Dims_create(7, 3, dims)    erroneous call
 *	------------------------------------------------------
 */
int MPI_Dims_create(int nnodes, int ndims, int dims[])
{
    int i;
    int freeprocs;
    int freedims;
    int nfactors;
    int *factors;
    int *procs;
    int *p;
    int err;

    /*if (MPI_PARAM_CHECK) {
        OMPI_ERR_INIT_FINALIZE(FUNC_NAME);

        if (0 > ndims) {
            return OMPI_ERRHANDLER_INVOKE (MPI_COMM_WORLD,
                                           MPI_ERR_DIMS, FUNC_NAME);
        }

        if ((0 != ndims) && (NULL == dims)) {
            return OMPI_ERRHANDLER_INVOKE (MPI_COMM_WORLD,
                                           MPI_ERR_ARG, FUNC_NAME);
        }

        if (1 > nnodes) {
            return OMPI_ERRHANDLER_INVOKE (MPI_COMM_WORLD,
                                           MPI_ERR_DIMS, FUNC_NAME);
        }
    }*/
	assert( ndims > 0 && nnodes > 1 && NULL != dims ); 
    /* Get # of free-to-be-assigned processes and # of free dimensions */
    freeprocs = nnodes;
    freedims = 0;
    for (i = 0, p = dims; i < ndims; ++i,++p) {
        if (*p == 0) {
            ++freedims;
        } else if ((*p < 0) || ((nnodes % *p) != 0)) {
            return OMPI_ERRHANDLER_INVOKE (MPI_COMM_WORLD, MPI_ERR_DIMS,
                                           FUNC_NAME);
        } else {
            freeprocs /= *p;
        }
    }

    if (freedims == 0) {
       if (freeprocs == 1) {
          return MPI_SUCCESS;
       }
       return OMPI_ERRHANDLER_NOHANDLE_INVOKE(MPI_ERR_DIMS,
                                     FUNC_NAME);
    }

    if (freeprocs == 1) {
        for (i = 0; i < ndims; ++i, ++dims) {
            if (*dims == 0) {
               *dims = 1;
            }
        }
        return MPI_SUCCESS;
    }

    /* Factor the number of free processes */
    if (MPI_SUCCESS != (err = getfactors(freeprocs, &nfactors, &factors))) {
       return OMPI_ERRHANDLER_NOHANDLE_INVOKE(err,
                                     FUNC_NAME);
    }

    /* Assign free processes to free dimensions */
    if (MPI_SUCCESS != (err = assignnodes(freedims, nfactors, factors, &procs))) {
       free(factors);
       return OMPI_ERRHANDLER_NOHANDLE_INVOKE(err,
                                     FUNC_NAME);
    }

    /* Return assignment results */
    p = procs;
    for (i = 0; i < ndims; ++i, ++dims) {
        if (*dims == 0) {
           *dims = *p++;
        }
    }

    free((char *) factors);
    free((char *) procs);

    /* all done */
    return MPI_SUCCESS;
}

/*
 *  assignnodes
 *
 *  Function:   - assign processes to dimensions
 *          - get "best-balanced" grid
 *          - greedy bin-packing algorithm used
 *          - sort dimensions in decreasing order
 *          - dimensions array dynamically allocated
 *  Accepts:    - # of dimensions
 *          - # of prime factors
 *          - array of prime factors
 *          - ptr to array of dimensions (returned value)
 *  Returns:    - 0 or ERROR
 */
static int
assignnodes(int ndim, int nfactor, int *pfacts, int **pdims)
{
    int *bins;
    int i, j;
    int n;
    int f;
    int *p;
    int *pmin;

    if (0 >= ndim) {
       return MPI_ERR_DIMS;
    }

    /* Allocate and initialize the bins */
    bins = (int *) malloc((unsigned) ndim * sizeof(int));
    if (NULL == bins) {
       return MPI_ERR_NO_MEM;
    }
    *pdims = bins;

    for (i = 0, p = bins; i < ndim; ++i, ++p) {
        *p = 1;
     }

    /* Loop assigning factors from the highest to the lowest */
    for (j = nfactor - 1; j >= 0; --j) {
        f = pfacts[j];
        /* Assign a factor to the smallest bin */
        pmin = bins;
        for (i = 1, p = pmin + 1; i < ndim; ++i, ++p) {
            if (*p < *pmin) {
                pmin = p;
            }
        }
        *pmin *= f;
     }

     /* Sort dimensions in decreasing order (O(n^2) for now) */
     for (i = 0, pmin = bins; i < ndim - 1; ++i, ++pmin) {
         for (j = i + 1, p = pmin + 1; j < ndim; ++j, ++p) {
             if (*p > *pmin) {
                n = *p;
                *p = *pmin;
                *pmin = n;
             }
         }
     }

     return MPI_SUCCESS;
}

/*
 *  getfactors
 *
 *  Function:   - factorize a number
 *  Accepts:    - number
 *          - # prime factors
 *          - array of prime factors
 *  Returns:    - MPI_SUCCESS or ERROR
 */
static int
getfactors(int num, int *nfactors, int **factors) {
    int size;
    int d;
    int i;
    int sqrtnum;

    if(num  < 2) {
        (*nfactors) = 0;
        (*factors) = NULL;
        return MPI_SUCCESS;
    }
    /* Allocate the array of prime factors which cannot exceed log_2(num) entries */
    sqrtnum = ceil(sqrt(num));
    size = ceil(log(num) / log(2));
    *factors = (int *) malloc((unsigned) size * sizeof(int));

    i = 0;
    /* determine all occurences of factor 2 */
    while((num % 2) == 0) {
        num /= 2;
        (*factors)[i++] = 2;
    }
    /* determine all occurences of uneven prime numbers up to sqrt(num) */
    d = 3;
    for(d = 3; (num > 1) && (d <= sqrtnum); d += 2) {
        while((num % d) == 0) {
            num /= d;
            (*factors)[i++] = d;
        }
    }
    /* as we looped only up to sqrt(num) one factor > sqrt(num) may be left over */
    if(num != 1) {
        (*factors)[i++] = num;
    }
    (*nfactors) = i;
    return MPI_SUCCESS;
}


#ifdef TEST_MAIN
int main(int argc,char** argv){
	if (argc != 2)	{fprintf(stderr,"usage: <N> ->  to fit in a 2D grid");return EXIT_FAILURE;}
	int dims2[2] = {0,0};
	int n = atoi(argv[1]);
	if (MPI_Dims_create(n,2,dims2)){return EXIT_FAILURE;} printf("%d 2D division:\t%d-%d\n",n,dims2[0],dims2[1]);
}
#endif //TEST_MAIN
