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

#include "hell.h"
#include <string.h>

/** \addtogroup conversionRoutines Conversion Routines
 *  @{
 */
  
#ifdef __cplusplus
extern "C" {
#endif

/** 
* \fn void computeHellAllocSize(int* allocationHeight, int hackSize, int rowsCount, const int *ellRowLengths)
 * Compute the HELL format allocation's height for the Hell allocation 
 * (the resulting size should be allocationHeight*hackSize*sizeof(elementsType)).
 * \param allocationHeight outputs the total allocation's height
 * \param hackSize the hack size for this matrix (32 or 64 are good choices)
 * \param rowsCount the rows count
 * \param ellRowLengths the row lengths array from the ell matrix to convert
*/
void computeHellAllocSize(
	int* allocationHeight,
	int hackSize,
	int rowsCount,
	const int *ellRowLengths
	);

/** 
* \fn void ellToHell(void *hellValues, int *hellIndices, int* hackOffsets, int hackSize, const void *ellValues, const int *ellIndices, int ellValuesPitch, int ellIndicesPitch, int *ellRowLengths, int rowsCount, spgpuType_t valuesType)
 * Convert a matrix from the ELL format to the HELL format.
 * \param hellValues pointer to the area that will be filled by the non zero coefficients
 * \param hellIndices pointer to the area that will be filled by the non zero indices
 * \param hackOffsets 
 * \param hackSize the hack size used to allocate hellValues and hellIndices (32 or 64 are good choices)
 * \param ellValues the input matrix coefficients
 * \param ellIndices the input matrix indices
 * \param ellValuesPitch the input values allocation pitch (in number of elements)
 * \param ellIndicesPitch the input indices allocation pitch (in number of elements)
 * \param ellRowLengths the row lengths array of the input matrix
 * \param rowsCount the rows count
 * \param valuesType the type of hellValues and ellValues elements (i.e. SPGPU_TYPE_FLOAT or SPGPU_TYPE_DOUBLE)
*/
void ellToHell(
	void *hellValues,
	int *hellIndices,
	int* hackOffsets,
	int hackSize,
	const void *ellValues,
	const int *ellIndices,
	int ellValuesPitch,
	int ellIndicesPitch,
	int *ellRowLengths,
	int rowsCount,
	spgpuType_t valuesType
	);

#ifdef __cplusplus
}
#endif


/** @}*/
