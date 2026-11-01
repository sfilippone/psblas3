/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#define VALUE_TYPE cuDoubleComplex
#define RVALUE_TYPE double
#define TYPE_SYMBOL Z
#define REALOP(x) cuCreal(x)
#define CONJOP(x) cuConj(x)

#undef USE_CUBLAS
#include "nrm2.cuh"