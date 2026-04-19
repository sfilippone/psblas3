/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2015
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

#define VALUE_TYPE cuFloatComplex
#define RES_VALUE_TYPE cuFloatComplex
#define TYPE_SYMBOL C
#include "abs_base.cuh"

