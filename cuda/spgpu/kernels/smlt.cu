/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

extern "C"
{
	//#include "cuda_util.h"
	int getGPUMultiProcessors();
	int getGPUMaxThreadsPerMP();
}


#define VALUE_TYPE float
#define TYPE_SYMBOL S
#include "mlt.cuh"

