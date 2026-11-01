/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#define VALUE_TYPE cuDoubleComplex
#define TYPE_SYMBOL Z
#undef USE_CUBLAS
#include "dot.cuh"