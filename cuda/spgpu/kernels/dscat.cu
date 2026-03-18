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

#define VALUE_TYPE double
#define TYPE_SYMBOL D
#include "scat_base.cuh"
