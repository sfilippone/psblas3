/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#include "cudadebug.h"
#include "cudalang.h"
#include <cuda_runtime.h>
#include "core.h"

extern "C"
{
	#include "vector.h"
	int getGPUMultiProcessors();
	int getGPUMaxThreadsPerMP();
	//#include "cuda_util.h"
}

#include "debug.h"

#define VALUE_TYPE cuFloatComplex
#define TYPE_SYMBOL C
#include "axpby.cuh"