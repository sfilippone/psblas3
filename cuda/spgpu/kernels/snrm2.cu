/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */

#define VALUE_TYPE float
#define RVALUE_TYPE float
#define TYPE_SYMBOL S
#define REALOP(x) x
#define CONJOP(x) x

#define USE_CUBLAS
#include "nrm2.cuh"