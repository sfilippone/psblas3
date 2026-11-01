#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2013 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
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
 * \fn void spgpuSdiaspmv(spgpuHandle_t handle,__device float *z,const __device float *y, float alpha, const __device float* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device float *x, float beta)
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
void spgpuSdiaspmv(spgpuHandle_t handle,
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
 * \fn void spgpuSdiaspmm(spgpuHandle_t handle,__device float *z, int ldz, const __device float *y, int ldy, float alpha, const __device float* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device float *x, int ldx, int mvectCount, float beta)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output multivector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param ldz The z pitch
 * \param y The y input multivector
 * \param ldy The y pitch
 * \param alpha The alpha scalar
 * \param dM The DIA non zero values allocation pointer (stored column-major)
 * \param offsets The DIA diagonals offsets vector
 * \param dMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rows the rows count
 * \param cols the columns count
 * \param diags the diagonals count
 * \param x the x multivector
 * \param ldx the x pitch
 * \param mvectCount the number of column in x, y and z (must be equals)
 * \param beta the beta scalar
 */
void spgpuSdiaspmm(spgpuHandle_t handle,
	__device float *z,
	int ldz,
	const __device float *y, 
	int ldy,
	float alpha, 
	const __device float* dM,
	const __device int* offsets,
	int dMPitch,
	int rows, 
	int cols,
	int diags,
	const __device float *x, 
	int ldx,
	float beta,
	int mvectCount);

/** 
 * \fn void spgpuDdiaspmv(spgpuHandle_t handle,__device double *z,const __device double *y, double alpha, const __device double* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device double *x, double beta)
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
void spgpuDdiaspmv(spgpuHandle_t handle,
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
 * \fn void spgpuDdiaspmm(spgpuHandle_t handle,__device double *z, int ldz, const __device double *y, int ldy, double alpha, const __device double* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device double *x, int ldx, int mvectCount, double beta)
 * Computes double precision z = alpha*A*x + beta*y, with A stored in Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output multivector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param ldz The z pitch
 * \param y The y input multivector
 * \param ldy The y pitch
 * \param alpha The alpha scalar
 * \param dM The DIA non zero values allocation pointer (stored column-major)
 * \param offsets The DIA diagonals offsets vector
 * \param dMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rows the rows count
 * \param cols the columns count
 * \param diags the diagonals count
 * \param x the x multivector
 * \param ldx the x pitch
 * \param mvectCount the number of column in x, y and z (must be equals)
 * \param beta the beta scalar
 */
void spgpuDdiaspmm(spgpuHandle_t handle,
	__device double *z,
	int ldz,
	const __device double *y, 
	int ldy,
	double alpha, 
	const __device double* dM,
	const __device int* offsets,
	int dMPitch,
	int rows, 
	int cols,
	int diags,
	const __device double *x, 
	int ldx,
	double beta,
	int mvectCount);

/** 
 * \fn void spgpuCdiaspmv(spgpuHandle_t handle,__device cuFloatComplex *z,const __device cuFloatComplex *y, cuFloatComplex alpha, const __device cuFloatComplex* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device cuFloatComplex *x, cuFloatComplex beta)
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
void spgpuCdiaspmv(spgpuHandle_t handle,
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
 * \fn void spgpuCdiaspmm(spgpuHandle_t handle,__device cuFloatComplex *z, int ldz, const __device cuFloatComplex *y, int ldy, cuFloatComplex alpha, const __device cuFloatComplex* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device cuFloatComplex *x, int ldx, int mvectCount, cuFloatComplex beta)
 * Computes single precision complex z = alpha*A*x + beta*y, with A stored in Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output multivector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param ldz The z pitch
 * \param y The y input multivector
 * \param ldy The y pitch
 * \param alpha The alpha scalar
 * \param dM The DIA non zero values allocation pointer (stored column-major)
 * \param offsets The DIA diagonals offsets vector
 * \param dMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rows the rows count
 * \param cols the columns count
 * \param diags the diagonals count
 * \param x the x multivector
 * \param ldx the x pitch
 * \param mvectCount the number of column in x, y and z (must be equals)
 * \param beta the beta scalar
 */
void spgpuCdiaspmm(spgpuHandle_t handle,
	__device cuFloatComplex *z,
	int ldz,
	const __device cuFloatComplex *y, 
	int ldy,
	cuFloatComplex alpha, 
	const __device cuFloatComplex* dM,
	const __device int* offsets,
	int dMPitch,
	int rows, 
	int cols,
	int diags,
	const __device cuFloatComplex *x, 
	int ldx,
	cuFloatComplex beta,
	int mvectCount);

/** 
 * \fn void spgpuZdiaspmv(spgpuHandle_t handle,__device cuDoubleComplex *z,const __device cuDoubleComplex *y, cuDoubleComplex alpha, const __device cuDoubleComplex* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device cuDoubleComplex *x, cuDoubleComplex beta)
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
void spgpuZdiaspmv(spgpuHandle_t handle,
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

/** 
 * \fn void spgpuZdiaspmm(spgpuHandle_t handle,__device cuDoubleComplex *z, int ldz, const __device cuDoubleComplex *y, int ldy, cuDoubleComplex alpha, const __device cuDoubleComplex* dM, const __device int* offsets, int dMPitch, int rows, int cols, int diags, const __device cuDoubleComplex *x, int ldx, int mvectCount, cuDoubleComplex beta)
 * Computes double precision complex z = alpha*A*x + beta*y, with A stored in Diagonal Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output multivector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param ldz The z pitch
 * \param y The y input multivector
 * \param ldy The y pitch
 * \param alpha The alpha scalar
 * \param dM The DIA non zero values allocation pointer (stored column-major)
 * \param offsets The DIA diagonals offsets vector
 * \param dMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rows the rows count
 * \param cols the columns count
 * \param diags the diagonals count
 * \param x the x multivector
 * \param ldx the x pitch
 * \param mvectCount the number of column in x, y and z (must be equals)
 * \param beta the beta scalar
 */
void spgpuZdiaspmm(spgpuHandle_t handle,
	__device cuDoubleComplex *z,
	int ldz,
	const __device cuDoubleComplex *y, 
	int ldy,
	cuDoubleComplex alpha, 
	const __device cuDoubleComplex* dM,
	const __device int* offsets,
	int dMPitch,
	int rows, 
	int cols,
	int diags,
	const __device cuDoubleComplex *x, 
	int ldx,
	cuDoubleComplex beta,
	int mvectCount);
	
/** @}*/
#ifdef __cplusplus
}
#endif
