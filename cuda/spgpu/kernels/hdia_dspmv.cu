/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#include "cudadebug.h"
#include "cudalang.h"
#include "core.h"

extern "C"
{
#include "hdia.h"
}

#include "debug.h"

//#define ENABLE_CACHE
#define VALUE_TYPE double
#define TYPE_SYMBOL D
//#define TEX_FETCH_TYPE int2
#include "hdia_spmv_base.cuh"

