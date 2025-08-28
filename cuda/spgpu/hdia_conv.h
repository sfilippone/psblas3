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


/** \addtogroup conversionRoutines Conversion Routines
 *  @{
 */
  
#ifdef __cplusplus
extern "C" {
#endif

int getHdiaHacksCount(int hackSize, int rowsCount);

void computeHdiaHackOffsets(
	int *allocationHeight,
	int *hackOffsets,
	int hackSize,
	const void* diaValues,
	int diaValuesPitch,	
	int diagonals,
	int rowsCount,
	spgpuType_t valuesType);
	
void diaToHdia(
	void *hdiaValues,
	int *hdiaOffsets,
	const int *hackOffsets,
	int hackSize,
	const void* diaValues,
	const int* diaOffsets,
	int diaValuesPitch,	
	int diagonals,
	int rowsCount,
	spgpuType_t valuesType
	);

void computeHdiaHackOffsetsFromCoo(
	int *allocationHeight,
	int *hackOffsets,
	int hackSize,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	int cooBaseIndex
	);

void cooToHdia(
	void *hdiaValues,
	int *hdiaOffsets,
	const int *hackOffsets,
	int hackSize,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	spgpuType_t valuesType
	);

void bcooToBhdia(
	void *hdiaValues,
	int *hdiaOffsets,
	const int *hackOffsets,
	int hackSize,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	spgpuType_t valuesType,
	int blockSize);
	
#ifdef __cplusplus
}
#endif


/** @}*/
