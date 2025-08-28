#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
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
