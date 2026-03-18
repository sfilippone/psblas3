/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2015
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */
 
#include "cudadebug.h"
#include "cudalang.h"
#include "cuComplex.h"
#include "core.h"

extern "C"
{
#include "hdia.h"
}

#include "debug.h"

#define VALUE_TYPE cuFloatComplex
#define TYPE_SYMBOL C
#define TEX_FETCH_TYPE cuFloatComplex
#include "hdia_spmv_base.cuh"
