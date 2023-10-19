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


/** \addtogroup diaFun DIA/HDIA Format
 *  @{
 */
 
#ifdef __cplusplus
extern "C" {
#endif


/** 
* \fn spgpuShdiaspmv (spgpuHandle_t handle, float* z, const float *y, float alpha, const float* dM, const int* offsets, int hackSize, const int* hackOffsets, int rows, int cols, const float *x, float beta)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in Hacked Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param dM The stacked HDIA non zero values allocation pointer
 * \param offsets The stacked HDIA diagonals offsets vector
 * \param hackSize The constant size of every hack (must be a multiple of 32) 
 * \param hackOffsets the array of base index offset for every hack of HDIA offsets vector, plus a last value equal to the size of the offsets vector
 * \param rows the rows count
 * \param cols the columns count
 * \param x the x vector
 * \param beta the beta scalar
 */
void 
spgpuShdiaspmv (spgpuHandle_t handle, 
	float* z, 
	const float *y, 
	float alpha, 
	const float* dM, 
	const int* offsets, 
	int hackSize, 
	const int* hackOffsets,
	int rows,
	int cols, 
	const float *x, 
	float beta);


/** 
* \fn spgpuDhdiaspmv (spgpuHandle_t handle, double* z, const double *y, double alpha, const double* dM, const int* offsets, int hackSize, const int* hackOffsets, int rows, int cols, const double *x, double beta)
 * Computes double precision z = alpha*A*x + beta*y, with A stored in Hacked Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param dM The stacked HDIA non zero values allocation pointer
 * \param offsets The stacked HDIA diagonals offsets vector
 * \param hackSize The constant size of every hack (must be a multiple of 32) 
 * \param hackOffsets the array of base index offset for every hack of HDIA offsets vector, plus a last value equal to the size of the offsets vector
 * \param rows the rows count
 * \param cols the columns count
 * \param x the x vector
 * \param beta the beta scalar
 */
void 
spgpuDhdiaspmv (spgpuHandle_t handle, 
	double* z, 
	const double *y, 
	double alpha, 
	const double* dM, 
	const int* offsets, 
	int hackSize, 
	const int* hackOffsets,
	int rows,
	int cols, 
	const double *x, 
	double beta);


/** 
* \fn spgpuChdiaspmv (spgpuHandle_t handle, cuFloatComplex* z, const cuFloatComplex *y, cuFloatComplex alpha, const cuFloatComplex* dM, const int* offsets, int hackSize, const int* hackOffsets, int rows, int cols, const cuFloatComplex *x, cuFloatComplex beta)
 * Computes single precision complex z = alpha*A*x + beta*y, with A stored in Hacked Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param dM The stacked HDIA non zero values allocation pointer
 * \param offsets The stacked HDIA diagonals offsets vector
 * \param hackSize The constant size of every hack (must be a multiple of 32) 
 * \param hackOffsets the array of base index offset for every hack of HDIA offsets vector, plus a last value equal to the size of the offsets vector
 * \param rows the rows count
 * \param cols the columns count
 * \param x the x vector
 * \param beta the beta scalar
 */
void 
spgpuChdiaspmv (spgpuHandle_t handle, 
	cuFloatComplex* z, 
	const cuFloatComplex *y, 
	cuFloatComplex alpha, 
	const cuFloatComplex* dM, 
	const int* offsets, 
	int hackSize, 
	const int* hackOffsets,
	int rows,
	int cols, 
	const cuFloatComplex *x, 
	cuFloatComplex beta);


/** 
* \fn spgpuZhdiaspmv (spgpuHandle_t handle, cuDoubleComplex* z, const cuDoubleComplex *y, cuDoubleComplex alpha, const cuDoubleComplex* dM, const int* offsets, int hackSize, const int* hackOffsets, int rows, int cols, const cuDoubleComplex *x, cuDoubleComplex beta)
 * Computes double precision complex z = alpha*A*x + beta*y, with A stored in Hacked Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param dM The stacked HDIA non zero values allocation pointer
 * \param offsets The stacked HDIA diagonals offsets vector
 * \param hackSize The constant size of every hack (must be a multiple of 32) 
 * \param hackOffsets the array of base index offset for every hack of HDIA offsets vector, plus a last value equal to the size of the offsets vector
 * \param rows the rows count
 * \param cols the columns count
 * \param x the x vector
 * \param beta the beta scalar
 */
void 
spgpuZhdiaspmv (spgpuHandle_t handle, 
	cuDoubleComplex* z, 
	const cuDoubleComplex *y, 
	cuDoubleComplex alpha, 
	const cuDoubleComplex* dM, 
	const int* offsets, 
	int hackSize, 
	const int* hackOffsets,
	int rows,
	int cols, 
	const cuDoubleComplex *x, 
	cuDoubleComplex beta);
		
/** @}*/

#ifdef __cplusplus
}
#endif

