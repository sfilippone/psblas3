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

// DIA/HDIA Compressed Matrix Format routines

/// This is the pitch alignment that must be fullfilled by the coefficients and the row pointers allocations.
#define DIA_PITCH_ALIGN_BYTE 128

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
* \fn spgpuShdiaspmm (int count, float *z, int zpitch, const float *y, int ypitch, float alpha, const __device float* dM, const __device int* offsets, int hackSize, const __device int* hackOffsets, int rows, int cols, const __device float *x, int xpitch, float beta)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in Hacked Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param count The cols count
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param zpitch The pitch of the output vector
 * \param y The y input vector
 * \param ypitch The pitch of the y input vector
 * \param alpha The alpha scalar
 * \param dM The stacked HDIA non zero values allocation pointer
 * \param offsets The stacked HDIA diagonals offsets vector
 * \param hackSize The constant size of every hack (must be a multiple of 32) 
 * \param hackOffsets the array of base index offset for every hack of HDIA offsets vector, plus a last value equal to the size of the offsets vector
 * \param rows the rows count
 * \param cols the columns count
 * \param x the x vector
 * \param xpitch The pitch of the x input vector
 * \param beta the beta scalar
 */
void 
spgpuShdiaspmm (spgpuHandle_t handle, 
  int count,
	__device float *z,
	int zpitch,
	const __device float *y,
	int ypitch,
	float alpha, 
	const __device float* dM, 
	const __device int* offsets, 
	int hackSize, 
	const __device int* hackOffsets,
	int rows,
	int cols, 
	const __device float *x, 
	int xpitch,
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
* \fn spgpuDhdiaspmm (int count, double *z, int zpitch, const double *y, int ypitch, double alpha, const __device double* dM, const __device int* offsets, int hackSize, const __device int* hackOffsets, int rows, int cols, const __device float *x, int xpitch, double beta)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in Hacked Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param count The cols count
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param zpitch The pitch of the output vector
 * \param y The y input vector
 * \param ypitch The pitch of the y input vector
 * \param alpha The alpha scalar
 * \param dM The stacked HDIA non zero values allocation pointer
 * \param offsets The stacked HDIA diagonals offsets vector
 * \param hackSize The constant size of every hack (must be a multiple of 32) 
 * \param hackOffsets the array of base index offset for every hack of HDIA offsets vector, plus a last value equal to the size of the offsets vector
 * \param rows the rows count
 * \param cols the columns count
 * \param x the x vector
 * \param xpitch The pitch of the x input vector
 * \param beta the beta scalar
 */
void 
spgpuDhdiaspmm (spgpuHandle_t handle, 
  int count,
	__device double *z,
	int zpitch,
	const __device double *y,
	int ypitch,
	double alpha, 
	const __device double* dM, 
	const __device int* offsets, 
	int hackSize, 
	const __device int* hackOffsets,
	int rows,
	int cols, 
	const __device double *x, 
	int xpitch,
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
* \fn spgpuChdiaspmm (int count, cuFloatComplex *z, int zpitch, const cuFloatComplex *y, int ypitch, cuFloatComplex alpha, const __device cuFloatComplex* dM, const __device int* offsets, int hackSize, const __device int* hackOffsets, int rows, int cols, const __device cuFloatComplex *x, int xpitch, cuFloatComplex beta)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in Hacked Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param count The cols count
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param zpitch The pitch of the output vector
 * \param y The y input vector
 * \param ypitch The pitch of the y input vector
 * \param alpha The alpha scalar
 * \param dM The stacked HDIA non zero values allocation pointer
 * \param offsets The stacked HDIA diagonals offsets vector
 * \param hackSize The constant size of every hack (must be a multiple of 32) 
 * \param hackOffsets the array of base index offset for every hack of HDIA offsets vector, plus a last value equal to the size of the offsets vector
 * \param rows the rows count
 * \param cols the columns count
 * \param x the x vector
 * \param xpitch The pitch of the x input vector
 * \param beta the beta scalar
 */
void 
spgpuChdiaspmm (spgpuHandle_t handle, 
  int count,
	__device cuFloatComplex *z,
	int zpitch,
	const __device cuFloatComplex *y,
	int ypitch,
	cuFloatComplex alpha, 
	const __device cuFloatComplex* dM, 
	const __device int* offsets, 
	int hackSize, 
	const __device int* hackOffsets,
	int rows,
	int cols, 
	const __device cuFloatComplex *x, 
	int xpitch,
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

/** 
* \fn spgpuZhdiaspmm (int count, cuDoubleComplex *z, int zpitch, const cuDoubleComplex *y, int ypitch, cuDoubleComplex alpha, const __device cuDoubleComplex* dM, const __device int* offsets, int hackSize, const __device int* hackOffsets, int rows, int cols, const __device cuDoubleComplex *x, int xpitch, cuDoubleComplex beta)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in Hacked Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param count The cols count
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param zpitch The pitch of the output vector
 * \param y The y input vector
 * \param ypitch The pitch of the y input vector
 * \param alpha The alpha scalar
 * \param dM The stacked HDIA non zero values allocation pointer
 * \param offsets The stacked HDIA diagonals offsets vector
 * \param hackSize The constant size of every hack (must be a multiple of 32) 
 * \param hackOffsets the array of base index offset for every hack of HDIA offsets vector, plus a last value equal to the size of the offsets vector
 * \param rows the rows count
 * \param cols the columns count
 * \param x the x vector
 * \param xpitch The pitch of the x input vector
 * \param beta the beta scalar
 */
void 
spgpuZhdiaspmm (spgpuHandle_t handle, 
  int count,
	__device cuDoubleComplex *z,
	int zpitch,
	const __device cuDoubleComplex *y,
	int ypitch,
	cuDoubleComplex alpha, 
	const __device cuDoubleComplex* dM, 
	const __device int* offsets, 
	int hackSize, 
	const __device int* hackOffsets,
	int rows,
	int cols, 
	const __device cuDoubleComplex *x, 
	int xpitch,
	cuDoubleComplex beta);

		
/** @}*/

#ifdef __cplusplus
}
#endif

