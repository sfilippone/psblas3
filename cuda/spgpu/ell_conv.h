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
 
#include "ell.h"
#include <string.h>

/** \addtogroup conversionRoutines Conversion Routines
 *  @{
 */

#ifdef __cplusplus
extern "C" {
#endif

/** 
* \fn void computeEllRowLenghts(int *ellRowLengths, int *ellMaxRowSize, int rowsCount, int nonZerosCount, const int* cooRowIndices, int cooBaseIndex)
 * Compute the Ell row lengths array (and the greatest row size) from the COO matrix format.
 * \param ellRowLengths Array of length rowsCount to be filled by the non zeros count for every matrix row
 * \param ellMaxRowSize outputs the greatest row size (in non zeros)
 * \param rowsCount the number of rows of the coo matrix to convert
 * \param nonZerosCount the non zeros count of the coo matrix to convert
 * \param cooRowIndices the row indices array for the coo matrix to convert
 * \param cooBaseIndex the input base index (e.g. 0 for C, 1 for Fortran)
 */
void computeEllRowLenghts(
	int *ellRowLengths,
	int *ellMaxRowSize,
	int rowsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	int cooBaseIndex
	);

/** 
* \fn int computeEllAllocPitch(int rowsCount)
 * This function returns a pitch (in number of elements) that can be used to allocate both indices and values arrays for ELL matrix format.
 * \param rowsCount the rows count
 * \return the pitch for an ELL matrix of rowsCount rows.
*/
int computeEllAllocPitch(int rowsCount);


/** 
* \fn void cooToEll(void *ellValues,int *ellIndices,int ellValuesPitch,int ellIndicesPitch,int ellMaxRowSize,int ellBaseIndex,int rowsCount,int nonZerosCount,const int* cooRowIndices,const int* cooColsIndices,const void* cooValues,int cooBaseIndex, spgpuType_t valuesType)
 * Convert a matrix in COO format to a matrix in ELL format.
 * The matrix is stored in column-major format.  The ellValues and ellIndices sizes are ellMaxRowSize * pitch (pitch is in bytes).
 * \param ellValues pointer to the area that will be filled by the non zero coefficients
 * \param ellIndices pointer to the area that will be filled by the non zero indices
 * \param ellValuesPitch the column-major allocation's pitch of ellValues (in number of elements)
 * \param ellIndicesPitch the column-major allocation's pitch of ellIndices (in number of elements)
 * \param ellMaxRowSize the greatest row size
 * \param ellBaseIndex the desired base index for the ELL matrix (e.g. 0 for C, 1 for Fortran)
 * \param rowsCount input matrix rows count
 * \param nonZerosCount input matrix non zeros count
 * \param cooRowIndices input matrix row indices pointer 
 * \param cooColsIndices input matrix column indices pointer 
 * \param cooValues input matrix non zeros values pointer
 * \param cooBaseIndex input matrix base index
 * \param valuesType the type for elements in ellValues and cooValues (i.e. SPGPU_TYPE_FLOAT or SPGPU_TYPE_DOUBLE)
 */
void cooToEll(
	void *ellValues,
	int *ellIndices,
	int ellValuesPitch,
	int ellIndicesPitch,
	int ellMaxRowSize,
	int ellBaseIndex,
	int rowsCount,
	int nonZerosCount,
	const int* cooRowIndices,
	const int* cooColsIndices,
	const void* cooValues,
	int cooBaseIndex,
	spgpuType_t valuesType
	);

void ellToOell(
	int *rIdx,
	void *dstEllValues,
	int *dstEllIndices,
	int *dstRs,
	const void *srcEllValues,
	const int *srcEllIndices,
	const int *srcRs,
	int ellValuesPitch,
	int ellIndicesPitch,
	int rowsCount,
	spgpuType_t valuesType
	);

#ifdef __cplusplus
}
#endif

/** @}*/
