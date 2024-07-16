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
 
#include "core.h"
#include "cuComplex.h"

/** \addtogroup ellFun ELL/HELL Format
 *  @{
 */
 
#ifdef __cplusplus
extern "C" {
#endif


// HELL Compressed Matrix Format routines

/// This is the pitch alignment that must be fullfilled by the coefficients and the row pointers allocations.
#define HELL_PITCH_ALIGN_BYTE 128

/** 
* \fn void spgpuShellspmv (spgpuHandle_t handle,__device float *z,const __device float *y, float alpha, const __device float* cM, const __device int* rP,int hackSize,const __device int* hackOffsets, const __device int* rS,const __device int* rIdx, int avgNnzPerRow, int rows, const __device float *x, float beta,int baseIndex)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in Hacked ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param cM The HELL non zero values allocation pointer
 * \param rP The HELL column indices allocation pointer
 * \param hackSize The constant size of every hack (must be a multiple of 32).
 * \param hackOffsets the array of base index offset for every hack of HELL non zero values allocation and HELL indices allocation.
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Hell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the HELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param rows the rows count
 * \param x the x vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuShellspmv (spgpuHandle_t handle,
	__device float *z,
	const __device float *y, 
	float alpha, 
	const __device float* cM, 
	const __device int* rP,
	int hackSize,
	const __device int* hackOffsets, 
	const __device int* rS,
	const __device int* rIdx, 
	int avgNnzPerRow,
	int rows, 
	const __device float *x, 
	float beta,
	int baseIndex);



/** 
* \fn void spgpuDhellspmv (spgpuHandle_t handle,__device double *z,const __device double *y, double alpha, const __device double* cM, const __device int* rP,int hackSize,const __device int* hackOffsets, const __device int* rS,const __device int* rIdx, int avgNnzPerRow, int rows, const __device double *x, double beta,int baseIndex)
 * Computes double precision z = alpha*A*x + beta*y, with A stored in Hacked ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param cM The HELL non zero values allocation pointer
 * \param rP The HELL column indices allocation pointer
 * \param hackSize The constant size of every hack (must be a multiple of 32).
 * \param hackOffsets the array of base index offset for every hack of HELL non zero values allocation and HELL indices allocation.
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Hell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the HELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param rows the rows count
 * \param x the x vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuDhellspmv (spgpuHandle_t handle,
	__device double *z,
	const __device double *y, 
	double alpha, 
	const __device double* cM, 
	const __device int* rP,
	int hackSize,
	const __device int* hackOffsets, 
	const __device int* rS,
	const __device int* rIdx, 
	int avgNnzPerRow,
	int rows, 
	const __device double *x, 
	double beta,
	int baseIndex);


/** 
* \fn void spgpuChellspmv (spgpuHandle_t handle,__device cuFloatComplex *z,const __device cuFloatComplex *y, cuFloatComplex alpha, const __device cuFloatComplex* cM, const __device int* rP,int hackSize,const __device int* hackOffsets, const __device int* rS, const __device int* rIdx, int avgNnzPerRow, int rows, const __device cuFloatComplex *x, cuFloatComplex beta, int baseIndex)
 * Computes single precision complex z = alpha*A*x + beta*y, with A stored in Hacked ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param cM The HELL non zero values allocation pointer
 * \param rP The HELL column indices allocation pointer
 * \param hackSize The constant size of every hack (must be a multiple of 32).
 * \param hackOffsets the array of base index offset for every hack of HELL non zero values allocation and HELL indices allocation.
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Hell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the HELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param rows the rows count
 * \param x the x vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuChellspmv (spgpuHandle_t handle,
	__device cuFloatComplex *z,
	const __device cuFloatComplex *y, 
	cuFloatComplex alpha, 
	const __device cuFloatComplex* cM, 
	const __device int* rP,
	int hackSize,
	const __device int* hackOffsets, 
	const __device int* rS,
	const __device int* rIdx, 
	int avgNnzPerRow,
	int rows, 
	const __device cuFloatComplex *x, 
	cuFloatComplex beta,
	int baseIndex);



/** 
* \fn void spgpuZhellspmv (spgpuHandle_t handle,__device cuDoubleComplex *z,const __device cuDoubleComplex *y, cuDoubleComplex alpha, const __device cuDoubleComplex* cM, const __device int* rP,int hackSize,const __device int* hackOffsets, const __device int* rS,const __device int* rIdx, int avgNnzPerRow, int rows, const __device cuDoubleComplex *x, cuDoubleComplex beta, int baseIndex)
 * Computes double precision complex z = alpha*A*x + beta*y, with A stored in Hacked ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param cM The HELL non zero values allocation pointer
 * \param rP The HELL column indices allocation pointer
 * \param hackSize The constant size of every hack (must be a multiple of 32).
 * \param hackOffsets the array of base index offset for every hack of HELL non zero values allocation and HELL indices allocation.
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Hell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the HELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param rows the rows count
 * \param x the x vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuZhellspmv (spgpuHandle_t handle,
	__device cuDoubleComplex *z,
	const __device cuDoubleComplex *y, 
	cuDoubleComplex alpha, 
	const __device cuDoubleComplex* cM, 
	const __device int* rP,
	int hackSize,
	const __device int* hackOffsets, 
	const __device int* rS,
	const __device int* rIdx, 
	int avgNnzPerRow,
	int rows, 
	const __device cuDoubleComplex *x, 
	cuDoubleComplex beta,
	int baseIndex);


/** @}*/

#ifdef __cplusplus
}
#endif
