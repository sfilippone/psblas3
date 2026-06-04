/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */
 
#include "cudadebug.h"
#include "cudalang.h"
#include "cuComplex.h"
#include "core.h"

extern "C"
{
#include "ell.h"
}

#include "debug.h"

#define VALUE_TYPE cuDoubleComplex
#define TYPE_SYMBOL Z
#include "ell_csput_base.cuh"

