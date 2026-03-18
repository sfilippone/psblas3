#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */
 
#include "cuda_runtime.h"
#include "stdio.h"

#ifdef DEBUG
#define cudaCheckError(...) \
	{ \
		cudaThreadSynchronize(); \
		cudaError_t lastError = cudaGetLastError(); \
		if (lastError != 0)	 \
		{  \
			printf(__VA_ARGS__); \
			printf("Error code: %i (%s)\n",lastError,cudaGetErrorString(lastError)); exit(0); \
		} \
	}
#else
#define cudaCheckError(...)
#endif
