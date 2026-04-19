#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2013
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */
 
#include "dia.h"
#include <string.h>

/** \addtogroup conversionRoutines Conversion Routines
 *  @{
 */

#ifdef __cplusplus
extern "C" {
#endif

int computeDiaDiagonalsCount(
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices);


void coo2dia(
	void* values,
	int* offsets,
	int valuesPitch,	
	int diagonals,
	int rowsCount,
	int columnsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	spgpuType_t valuesType);
	
	
/** 
* \fn int computeDiaAllocPitch(int rowsCount)
 * This function returns a pitch (in number of elements) that can be used to allocate the values array for DIA matrix format.
 * \param rowsCount the rows count
 * \return the pitch for an DIA matrix of rowsCount rows.
*/
int computeDiaAllocPitch(int rowsCount);

	
#ifdef __cplusplus
}
#endif

/** @}*/
