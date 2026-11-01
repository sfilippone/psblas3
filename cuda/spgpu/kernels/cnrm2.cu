/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#define VALUE_TYPE cuFloatComplex
#define RVALUE_TYPE float
#define TYPE_SYMBOL C
#define REALOP(x) cuCrealf(x)
#define CONJOP(x) cuConjf(x)

#undef USE_CUBLAS
#include "nrm2.cuh"