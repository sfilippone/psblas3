/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2015
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#include "cudadebug.h"
#include "cudalang.h"
#include "core.h"

extern "C"
{
#include "dia.h"
}

#include "debug.h"

#define VALUE_TYPE double
#define TYPE_SYMBOL D
#define TEX_FETCH_TYPE int2
#include "dia_spmv_base.cuh"

