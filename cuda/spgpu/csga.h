#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2013 
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

#include "core.h"
#include "cuComplex.h"

/** \addtogroup diaFun DIA/HDIA Format
 *  @{
 */
 
#ifdef __cplusplus
extern "C" {
#endif

int dCSGAMV(spgpuHandle_t handle, 
	    double beta,
	    double* y, 
	    double alpha, 
	    const double* as, 
	    const int* ja,
	    const int* irp,
	    int m,
	    int n,
	    int ncol,
	    int  numBlocks,
	    const int* rowBlocks,
	    const double *x,
	    int baseIndex,
	    int *rb);
	
/** @}*/
#ifdef __cplusplus
}
#endif
