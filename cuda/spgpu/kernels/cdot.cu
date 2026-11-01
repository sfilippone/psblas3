/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#define VALUE_TYPE cuFloatComplex
#define TYPE_SYMBOL C
#undef USE_CUBLAS
#include "dot.cuh"