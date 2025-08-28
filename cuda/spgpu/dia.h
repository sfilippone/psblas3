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


// DIA/HDIA Compressed Matrix Format routines

/// This is the pitch alignment that must be fullfilled by the coefficients and the row pointers allocations.
#define DIA_PITCH_ALIGN_BYTE 128

/** 
* \fn void spgpuSdiaspmv (spgpuHandle_t handle,__device float *z,const __device float *y, float alpha, const __device float* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device float *x, float beta)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param dM The DIA non zero values allocation pointer (stored column-major)
 * \param offsets The DIA diagonals offsets vector
 * \param dMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rows the rows count
 * \param cols the columns count
 * \param diags the diagonals count
 * \param x the x vector
 * \param beta the beta scalar
 */
void spgpuSdiaspmv (spgpuHandle_t handle,
	__device float *z,
	const __device float *y, 
	float alpha, 
	const __device float* dM,
	const __device int* offsets,
	int dMPitch,
	int rows, 
	int cols,
	int diags,
	const __device float *x, 
	float beta);
	

/** 
* \fn void spgpuDdiaspmv (spgpuHandle_t handle,__device double *z,const __device double *y, double alpha, const __device double* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device double *x, double beta)
 * Computes double precision z = alpha*A*x + beta*y, with A stored in Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param dM The DIA non zero values allocation pointer (stored column-major)
 * \param offsets The DIA diagonals offsets vector
 * \param dMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rows the rows count
 * \param cols the columns count
 * \param diags the diagonals count
 * \param x the x vector
 * \param beta the beta scalar
 */
void spgpuDdiaspmv (spgpuHandle_t handle,
	__device double *z,
	const __device double *y, 
	double alpha, 
	const __device double* dM,
	const __device int* offsets,
	int dMPitch,
	int rows, 
	int cols,
	int diags,
	const __device double *x, 
	double beta);


/** 
* \fn void spgpuCdiaspmv (spgpuHandle_t handle,__device cuFloatComplex *z,const __device cuFloatComplex *y, cuFloatComplex alpha, const __device cuFloatComplex* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device cuFloatComplex *x, cuFloatComplex beta)
 * Computes single precision complex z = alpha*A*x + beta*y, with A stored in Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param dM The DIA non zero values allocation pointer (stored column-major)
 * \param offsets The DIA diagonals offsets vector
 * \param dMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rows the rows count
 * \param cols the columns count
 * \param diags the diagonals count
 * \param x the x vector
 * \param beta the beta scalar
 */
void spgpuCdiaspmv (spgpuHandle_t handle,
	__device cuFloatComplex *z,
	const __device cuFloatComplex *y, 
	cuFloatComplex alpha, 
	const __device cuFloatComplex* dM,
	const __device int* offsets,
	int dMPitch,
	int rows, 
	int cols,
	int diags,
	const __device cuFloatComplex *x, 
	cuFloatComplex beta);
	

/** 
* \fn void spgpuZdiaspmv (spgpuHandle_t handle,__device cuDoubleComplex *z,const __device cuDoubleComplex *y, cuDoubleComplex alpha, const __device cuDoubleComplex* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device cuDoubleComplex *x, cuDoubleComplex beta)
 * Computes double precision complex z = alpha*A*x + beta*y, with A stored in Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param dM The DIA non zero values allocation pointer (stored column-major)
 * \param offsets The DIA diagonals offsets vector
 * \param dMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rows the rows count
 * \param cols the columns count
 * \param diags the diagonals count
 * \param x the x vector
 * \param beta the beta scalar
 */
void spgpuZdiaspmv (spgpuHandle_t handle,
	__device cuDoubleComplex *z,
	const __device cuDoubleComplex *y, 
	cuDoubleComplex alpha, 
	const __device cuDoubleComplex* dM,
	const __device int* offsets,
	int dMPitch,
	int rows, 
	int cols,
	int diags,
	const __device cuDoubleComplex *x, 
	cuDoubleComplex beta);
	
/** @}*/
#ifdef __cplusplus
}
#endif
