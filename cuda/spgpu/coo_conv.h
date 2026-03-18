#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2014 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
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
