#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014 
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
 
#include <string.h>
#include "core.h"

/** \addtogroup conversionRoutines Conversion Routines
 *  @{
 */

#ifdef __cplusplus
extern "C" {
#endif

int computeBcooSize(int blockRows, int blockCols, const int* rows, const int* cols, int nonZeros);

// column-major format for blocks
void cooToBcoo(int* bRows, int* bCols, void* blockValues, /*int isBlockColumnMajor,*/ int blockRows, int blockCols, 
	const int* rows, const int* cols, const void* values, int nonZeros, spgpuType_t valuesType);

#ifdef __cplusplus
}
#endif

/** @}*/
