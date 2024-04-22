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


// ELL/HELL Compressed Matrix Format routines

/// This is the pitch alignment that must be fullfilled by the coefficients and the row pointers allocations.
#define ELL_PITCH_ALIGN_BYTE 128

/** 
* \fn void spgpuSellspmv (spgpuHandle_t handle, __device float *z, const __device float *y, float alpha, const __device float* cM, const __device int* rP, int cMPitch, int rPPitch, const __device int* rS, const __device int* rIdx, int avgNnzPerRow, int maxNnzPerRow, int rows, const __device float *x, float beta,int baseIndex)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Ell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the ELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param maxNnzPerRow Maximum number of non zeroes per row.
 * \param rows the rows count
 * \param x the x vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuSellspmv (spgpuHandle_t handle,
	__device float *z,
	const __device float *y, 
	float alpha, 
	const __device float* cM, 
	const __device int* rP, 
	int cMPitch, 
	int rPPitch, 
	const __device int* rS, 
	const __device int* rIdx, 
	int avgNnzPerRow,
	int maxNnzPerRow,
	int rows, 
	const __device float *x, 
	float beta,
	int baseIndex);

/** 
* \fn void spgpuSellspmm (spgpuHandle_t handle,int count,__device float *z,int zpitch,const __device float *y,int ypitch,float alpha, const __device float* cM, const __device int* rP,int cMPitch,int rPPitch,const __device int* rS,const __device int* rIdx, int avgRowSize,int maxRowSize,int rows, const __device float *x,int xpitch,float beta,int baseIndex)
 * Computes single precision z = alpha*A*x + beta*y, with A stored in ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param count The cols count
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param zpitch The pitch of the output vector
 * \param y The y input vector
 * \param ypitch The pitch of the y input vector
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Hell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the HELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param maxNnzPerRow Maximum number of non zeroes per row.
 * \param rows the rows count
 * \param x the x vector
 * \param xpitch The pitch of the x input vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuSellspmm(spgpuHandle_t handle,
	int count,
	__device float *z,
	int zpitch,
	const __device float *y,
	int ypitch,
	float alpha, 
	const __device float* cM, 
	const __device int* rP,
    int cMPitch,
    int rPPitch,
	const __device int* rS,
	const __device int* rIdx, 
    int avgNnzPerRow,
    int maxNnzPerRow,
	int rows, 
	const __device float *x,
	int xpitch,
	float beta,
	int baseIndex);

/** 
* \fn void spgpuDellspmv (spgpuHandle_t handle,__device double *z,const __device double *y, double alpha, const __device double* cM, const __device int* rP, int cMPitch, int rPPitch, const __device int* rS, const __device int* rIdx, int avgNnzPerRow, int maxNnzPerRow, int rows, const __device double *x, double beta,int baseIndex)
 * Computes double precision z = alpha*A*x + beta*y, with A stored in ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Ell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the ELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param maxNnzPerRow Maximum number of non zeroes per row.
 * \param rows the rows count
 * \param x the x vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuDellspmv (spgpuHandle_t handle,
	__device double *z,
	const __device double *y, 
	double alpha, 
	const __device double* cM, 
	const __device int* rP, 
	int cMPitch, 
	int rPPitch, 
	const __device int* rS, 
	const __device int* rIdx,  
	int avgNnzPerRow,
	int maxNnzPerRow,
	int rows, 
	const __device double *x, 
	double beta,
	int baseIndex);

/** 
* \fn void spgpuDellspmm (	int count,__device double *z,int zpitch,const __device double *y,int ypitch,double alpha, const __device double* cM, const __device int* rP,int cMPitch,int rPPitch,const __device int* rS,const __device int* rIdx, int avgNnzPerRow,int maxNnzPerRow,int rows, const __device double *x,int xpitch,double beta,int baseIndex)
 * Computes double precision z = alpha*A*x + beta*y, with A stored in ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param count The cols count
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param zpitch The pitch of the output vector
 * \param y The y input vector
 * \param ypitch The pitch of the y input vector
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Hell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the HELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param maxNnzPerRow Maximum number of non zeroes per row.
 * \param rows the rows count
 * \param x the x vector
 * \param xpitch The pitch of the x input vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuDellspmm(spgpuHandle_t handle,
	int count,
	__device double *z,
	int zpitch,
	const __device double *y,
	int ypitch,
	double alpha, 
	const __device double* cM, 
	const __device int* rP,
    int cMPitch,
    int rPPitch,
	const __device int* rS,
	const __device int* rIdx, 
    int avgNnzPerRow,
    int maxNnzPerRow,
	int rows, 
	const __device double *x,
	int xpitch,
	double beta,
	int baseIndex);

/** 
* \fn void spgpuCellspmv (spgpuHandle_t handle,__device cuFloatComplex *z,const __device cuFloatComplex *y, cuFloatComplex alpha, const __device cuFloatComplex* cM, const __device int* rP, int cMPitch, int rPPitch, const __device int* rS, const __device int* rIdx, int avgNnzPerRow, int maxNnzPerRow, int rows, const __device cuFloatComplex *x, cuFloatComplex beta, int baseIndex)
 * Computes single precision complex z = alpha*A*x + beta*y, with A stored in ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Ell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the ELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param maxNnzPerRow Maximum number of non zeroes per row.
 * \param rows the rows count
 * \param x the x vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuCellspmv (spgpuHandle_t handle,
	__device cuFloatComplex *z,
	const __device cuFloatComplex *y, 
	cuFloatComplex alpha, 
	const __device cuFloatComplex* cM, 
	const __device int* rP, 
	int cMPitch, 
	int rPPitch, 
	const __device int* rS, 
	const __device int* rIdx,  
	int avgNnzPerRow,
	int maxNnzPerRow,
	int rows, 
	const __device cuFloatComplex *x, 
	cuFloatComplex beta,
	int baseIndex);

/** 
* \fn void spgpuCellspmm (spgpuHandle_t handle,int count,__device cuFloatComplex *z,int zpitch,const __device cuFloatComplex *y,int ypitch,cuFloatComplex alpha, const __device cuFloatComplex* cM, const __device int* rP,int cMPitch,int rPPitch,const __device int* rS,const __device int* rIdx, int avgNnzPerRow,int maxNnzPerRow,int rows, const __device cuFloatComplex *x,int xpitch,cuFloatComplex beta,int baseIndex)
 * Computes single precision complex z = alpha*A*x + beta*y, with A stored in ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param count The cols count
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param zpitch The pitch of the output vector
 * \param y The y input vector
 * \param ypitch The pitch of the y input vector
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Hell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the HELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param maxNnzPerRow Maximum number of non zeroes per row.
 * \param rows the rows count
 * \param x the x vector
 * \param xpitch The pitch of the x input vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuCellspmm(spgpuHandle_t handle,
	int count,
	__device cuFloatComplex *z,
	int zpitch,
	const __device cuFloatComplex *y,
	int ypitch,
	cuFloatComplex alpha, 
	const __device cuFloatComplex* cM, 
	const __device int* rP,
    int cMPitch,
    int rPPitch,
	const __device int* rS,
	const __device int* rIdx, 
    int avgNnzPerRow,
    int maxNnzPerRow,
	int rows, 
	const __device cuFloatComplex *x,
	int xpitch,
	cuFloatComplex beta,
	int baseIndex);

/** 
* \fn void spgpuZellspmv (spgpuHandle_t handle,__device cuDoubleComplex *z,const __device cuDoubleComplex *y, cuDoubleComplex alpha, const __device cuDoubleComplex* cM, const __device int* rP, int cMPitch, int rPPitch, const __device int* rS, const __device int* rIdx, int avgNnzPerRow, int maxNnzPerRow, int rows, const __device cuDoubleComplex *x, cuDoubleComplex beta, int baseIndex)
 * Computes double precision complex z = alpha*A*x + beta*y, with A stored in ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param y The y input vector
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Ell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the ELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param maxNnzPerRow Maximum number of non zeroes per row.
 * \param rows the rows count
 * \param x the x vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuZellspmv (spgpuHandle_t handle,
	__device cuDoubleComplex *z,
	const __device cuDoubleComplex *y, 
	cuDoubleComplex alpha, 
	const __device cuDoubleComplex* cM, 
	const __device int* rP, 
	int cMPitch, 
	int rPPitch, 
	const __device int* rS, 
	const __device int* rIdx,  
	int avgNnzPerRow,
	int maxNnzPerRow,
	int rows, 
	const __device cuDoubleComplex *x, 
	cuDoubleComplex beta,
	int baseIndex);
	
/** 
* \fn void spgpuCellspmm (spgpuHandle_t handle,int count,__device cuDoubleComplex *z,int zpitch,const __device cuDoubleComplex *y,int ypitch,cuDoubleComplex alpha, const __device cuDoubleComplex* cM, const __device int* rP,int cMPitch,int rPPitch,const __device int* rS,const __device int* rIdx, int avgNnzPerRow,int maxNnzPerRow,int rows, const __device cuDoubleComplex *x,int xpitch,cuDoubleComplex beta,int baseIndex)
 * Computes double precision complex z = alpha*A*x + beta*y, with A stored in ELLpack Format on GPU.
 * \param handle The spgpu handle used to call this routine
 * \param count The cols count
 * \param z The output vector of the routine. z could be y, but not y + k (i.e. an overlapping area over y, but starting from a base index different from y).
 * \param zpitch The pitch of the output vector
 * \param y The y input vector
 * \param ypitch The pitch of the y input vector
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param rIdx (optional) An array containing the row index per every row (i.e. the reorder array) of the Hell matrix. Pass NULL if you don't use a reorder array (i.e. the k-th row is stored in the k-th position in the HELL format).
 * \param avgNnzPerRow (optional) Average number of non zeroes per row. Pass 0 if you don't have such information.
 * \param maxNnzPerRow Maximum number of non zeroes per row.
 * \param rows the rows count
 * \param x the x vector
 * \param xpitch The pitch of the x input vector
 * \param beta the beta scalar
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuZellspmm(spgpuHandle_t handle,
	int count,
	__device cuDoubleComplex *z,
	int zpitch,
	const __device cuDoubleComplex *y,
	int ypitch,
	cuDoubleComplex alpha, 
	const __device cuDoubleComplex* cM, 
	const __device int* rP,
    int cMPitch,
    int rPPitch,
	const __device int* rS,
	const __device int* rIdx, 
    int avgNnzPerRow,
    int maxNnzPerRow,
	int rows, 
	const __device cuDoubleComplex *x,
	int xpitch,
	cuDoubleComplex beta,
	int baseIndex);
	
/** 
* \fn void spgpuSellcsput (spgpuHandle_t handle, float alpha, __device float *cM, __device const int* rP, int cMPitch, int rPPitch, __device const int* rS, int nnz, __device int *aI, __device int *aJ, __device float *aVal, int baseIndex)
 * Replaces the values at coordinate (aI[i], aJ[i]) inside A with the value aVal[i] * alpha, with A stored in ELLpack Format on GPU.
 * It assumes that indices from the same row inside rP are sorted by ascending order.
 * Values are single precision floating point numbers.
 * \param handle The spgpu handle used to call this routine
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param nnz the number of triples (aI, aJ, aVal) to process
 * \param aI The row coordinates vector
 * \param aJ The column coordinates vector
 * \param aVal The values vector
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuSellcsput
	(spgpuHandle_t handle, 
	float alpha, 
	__device float *cM, 
	__device const int* rP, 
	int cMPitch, 
	int rPPitch, 
	__device const int* rS,
	int nnz, 
	__device int *aI, 
	__device int *aJ, 
	__device float *aVal, 
	int baseIndex);	

/** 
* \fn void spgpuDellcsput (spgpuHandle_t handle, double alpha, __device double *cM, __device const int* rP, int cMPitch, int rPPitch, __device const int* rS, int nnz, __device int *aI, __device int *aJ, __device double *aVal, int baseIndex)
 * Replaces the values at coordinate (aI[i], aJ[i]) inside A with the value aVal[i] * alpha, with A stored in ELLpack Format on GPU.
 * It assumes that indices from the same row inside rP are sorted by ascending order.
 * Values are double precision floating point numbers.
 * \param handle The spgpu handle used to call this routine
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param nnz the number of triples (aI, aJ, aVal) to process
 * \param aI The row coordinates vector
 * \param aJ The column coordinates vector
 * \param aVal The values vector
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */	
void spgpuDellcsput
	(spgpuHandle_t handle, 
	double alpha, 
	__device double *cM, 
	__device const int* rP, 
	int cMPitch, 
	int rPPitch, 
	__device const int* rS,
	int nnz, 
	__device int *aI, 
	__device int *aJ, 
	__device double *aVal, 
	int baseIndex);	

/** 
* \fn void spgpuCellcsput (spgpuHandle_t handle, cuFloatComplex alpha, __device cuFloatComplex *cM, __device const int* rP, int cMPitch, int rPPitch, __device const int* rS, int nnz, __device int *aI, __device int *aJ, __device cuFloatComplex *aVal, int baseIndex)
 * Replaces the values at coordinate (aI[i], aJ[i]) inside A with the value aVal[i] * alpha, with A stored in ELLpack Format on GPU.
 * It assumes that indices from the same row inside rP are sorted by ascending order.
 * Values are single precision floating point complex numbers.
 * \param handle The spgpu handle used to call this routine
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param nnz the number of triples (aI, aJ, aVal) to process
 * \param aI The row coordinates vector
 * \param aJ The column coordinates vector
 * \param aVal The values vector
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuCellcsput
	(spgpuHandle_t handle, 
	cuFloatComplex alpha, 
	__device cuFloatComplex *cM, 
	__device const int* rP, 
	int cMPitch, 
	int rPPitch, 
	__device const int* rS,
	int nnz, 
	__device int *aI, 
	__device int *aJ, 
	__device cuFloatComplex *aVal, 
	int baseIndex);	

/** 
* \fn void spgpuZellcsput (spgpuHandle_t handle, cuDoubleComplex alpha, __device cuDoubleComplex *cM, __device const int* rP, int cMPitch, int rPPitch, __device const int* rS, int nnz, __device int *aI, __device int *aJ, __device cuDoubleComplex *aVal, int baseIndex)
 * Replaces the values at coordinate (aI[i], aJ[i]) inside A with the value aVal[i] * alpha, with A stored in ELLpack Format on GPU.
 * It assumes that indices from the same row inside rP are sorted by ascending order.
 * Values are double precision floating point complex numbers.
 * \param handle The spgpu handle used to call this routine
 * \param alpha The alpha scalar
 * \param cM The ELL non zero values allocation pointer
 * \param rP The ELL column indices allocation pointer
 * \param cMPitch the pitch (in number of elements) of the allocation containing the matrix non zero values
 * \param rPPitch  the pitch (in number of elements) of the allocation containing the matrix non zero column indices
 * \param rS the array containing the row sized (in non zero elements)
 * \param nnz the number of triples (aI, aJ, aVal) to process
 * \param aI The row coordinates vector
 * \param aJ The column coordinates vector
 * \param aVal The values vector
 * \param baseIndex the ELL format base index used (i.e. 0 for C, 1 for Fortran).
 */
void spgpuZellcsput
	(spgpuHandle_t handle, 
	cuDoubleComplex alpha, 
	__device cuDoubleComplex *cM, 
	__device const int* rP, 
	int cMPitch, 
	int rPPitch, 
	__device const int* rS,
	int nnz, 
	__device int *aI, 
	__device int *aJ, 
	__device cuDoubleComplex *aVal, 
	int baseIndex);		
	
	
	
	
/** @}*/

#ifdef __cplusplus
}
#endif
