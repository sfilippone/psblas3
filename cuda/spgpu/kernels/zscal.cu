/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */
 
#include "stdio.h"
#include "cudadebug.h"
#include "cudalang.h"
#include "core.h"

extern "C"
{
#include "vector.h"
}

#include "debug.h"

#define VALUE_TYPE cuDoubleComplex
#define TYPE_SYMBOL Z
#include "scal_base.cuh"
