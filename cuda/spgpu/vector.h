#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *     Salvatore Filippone - University of Rome Tor Vergata
 *
 */
 
#include "core.h"

/** \addtogroup vecFun Vectors and sparse vectors routines
 *  @{
 */

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * \fn void spgpuIgath(spgpuHandle_t handle, __device int *xValues, int xNnz, const __device int *xIndices, int xBaseIndex, const __device int* y)
 * Integer gather from y to sparse(x). Computes the integer gather from y to xValues (using xIndices).
 * \param handle the spgpu handle used to call this routine
 * \param xValues the destination array for gathered values
 * \param xNnz the number of elements to gather
 * \param xIndices the array of indices for the elements to be gathered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param y the source vector (from which the elements will be gathered)
 */
void spgpuIgath(spgpuHandle_t handle,
	__device int *xValues,
	int xNnz,
	const __device int *xIndices,
	int xBaseIndex,
	const __device int* y);
	
/** 
 * \fn void spgpuIscat(spgpuHandle_t handle, __device int* y, int xNnz, const __device int *xValues, const __device int *xIndices, int xBaseIndex, int beta)
 * Integer scatter from sparse(x) to y. Computes the integer scatter from xValues to y (using xIndices).
 * The scattered element will be, for i in [0, xNnz), y[xIndices[i]] = beta*y[xIndices[i]] + xValues[i] (to be noted that
 * y values will be multiplied with beta just for scattered values).
 * \param handle the spgpu handle used to call this routine
 * \param y the destination vector (to which the elements will be scattered)
 * \param xNnz the number of elements to scatter
 * \param xValues the source array from which the values will be read
 * \param xIndices the array of indices for the elements to be scattered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param beta the beta value
 */
void spgpuIscat(spgpuHandle_t handle,
	__device int* y,
	int xNnz,
	const __device int *xValues,
	const __device int *xIndices,
	int xBaseIndex, int beta);	

/** 
 * \fn float spgpuSdot (spgpuHandle_t handle, int n, __device float* a, __device float* b)
 * Computes single precision dot product of a and b vectors.
 * \param handle The spgpu handle used to call this routine
 * \param n the vectors length
 * \param a the first input vector
 * \param b the second input vector
 * \return the dot product
 */
float spgpuSdot(spgpuHandle_t handle,
	int n,
	__device float* a,
	__device float* b);

/** 
 * \fn void spgpuSmdot (spgpuHandle_t handle, float* y, int n, __device float* a, __device float* b, int count, int pitch)
 * Computes single precision dot product of a and b multivectors.
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector couples from the multivectors a and b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param count the number of vectors in every multivector
 * \param pitch the pitch, in number of elements, of every multivectors (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuSmdot(spgpuHandle_t handle,
	float* y,
	int n,
	__device float* a,
	__device float* b,
	int count,
	int pitch);

/** 
 * \fn void spgpuSmdot (spgpuHandle_t handle, float* y, int n, __device float* a, __device float* b)
 * Computes float precision dot product of a and b vectors (count assumed = 1 for both).
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of a and b
 * \param n the vectors' length
 * \param a the first input vector
 * \param b the second input vector
 */
void spgpuSmsdot(spgpuHandle_t handle,
	float* y,
	int n,
	__device float* a,
	__device float* b);

/** 
 * \fn void spgpuSmvdot (spgpuHandle_t handle, float* y, int n, __device float* a, __device float* b, int countA, int pitchA)
 * Computes float precision dot product of multivector a and vector b (count assumed = 1 for b).
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector in multivector a with the vector b
 * \param n the vectors' length
 * \param a the input multivector
 * \param b the input vector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuSmvdot(spgpuHandle_t handle,
	float* y,
	int n,
	__device float* a,
	__device float* b,
	int countA,
	int pitchA);

/** 
 * \fn void spgpuSmmdot (spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device float* a, __device float* b, int countA, int pitchA, int countB, int pitchB)
 * Computes float precision dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array, made by dot products of every vector in multivector a with the vector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuSmvdot_CFI(spgpuHandle_t handle,
	CFI_cdesc_t* y,
	int n,
	__device float* a,
	__device float* b,
	int countA,
	int pitchA);

/** 
 * \fn void spgpuSmmdot (spgpuHandle_t handle, float* y, int n, __device float* a, __device float* b, int countA, int pitchA, int countB, int pitchB)
 * Computes float precision dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector in multivector a with every vector in multivector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param countB the number of vectors in multivector b
 * \param pitchB the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 */
void spgpuSmmdot(spgpuHandle_t handle,
	float* y,
	int n,
	__device float* a,
	__device float* b,
	int countA,
	int pitchA,
	int countB,
	int pitchB);

/** 
 * \fn void spgpuSmmdot (spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device float* a, __device float* b, int countA, int pitchA, int countB, int pitchB)
 * Computes float precision dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array, made by dot products of every vector in multivector a with every vector in multivector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param countB the number of vectors in multivector b
 * \param pitchB the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 */
void spgpuSmmdot_CFI(spgpuHandle_t handle,
	CFI_cdesc_t* y,
	int n,
	__device float* a,
	__device float* b,
	int countA,
	int pitchA,
	int countB,
	int pitchB);

/** 
 * \fn void spgpuSmlt (spgpuHandle_t handle, __device float* z, int countZ, int pitchZ,
 * 						int n, float alpha, __device float* x, int countX, int pitchY,
 * 											__device float* y, int countY, int pitchY,
 * 						       float beta, bool conjgX, bool conjgY)
 * Computes z = alpha * op(x) * op(y) + beta * z with x, y, z vectors/multivectors possibly using multiple column both in input or output
 * \param handle the spgpu handle used to call this routine
 * \param z the output multivector
 * \param countX the number of vectors in multivector Z
 * \param pitchX the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param n the vectors' length
 * \param alpha the first scalar coefficient
 * \param y the second input multivector
 * \param countY the number of vectors in multivector y
 * \param pitchY the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 * \param x the first input multivector
 * \param countX the number of vectors in multivector x
 * \param pitchX the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param beta the second scalar coefficient
 * \param congjx boolean flag 
 * \param congjy boolean flag 
 */
void spgpuSmlt(spgpuHandle_t handle,
	__device float* z,
	int countZ,
	int pitchZ,
	int n,
	float alpha,
	__device float* x,
	int countX,
	int pitchX,
	__device float* y,
	int countY,
	int pitchY,
	float beta,
	bool conjgX,
	bool conjgY);

/** 
 * \fn void spgpuSabs (spgpuHandle_t handle, __device float* y, int n, float alpha, __device float* x)
 * Computes single precision (y = alpha * abs(x)) for each element in x.
 * \param handle The spgpu handle used to call this routine
 * \param y the resulting vector (could be x)
 * \param n the vectors length
 * \param alpha the alpha value
 * \param x the input vector
 */
void spgpuSabs(spgpuHandle_t handle,
	__device float* y,
	int n,
	float alpha,
	__device float* x);
	
/** 
 * \fn float spgpuSnrm2(spgpuHandle_t handle, int n, __device float* x)
 * Computes the single precision Euclidean vector norm of x. 
 * \param handle the spgpu handle used to call this routine
 * \param n the vector's length
 * \param x the input vector
 * \return the euclidean vector norm
 */
float spgpuSnrm2(spgpuHandle_t handle,
	int n,
	__device float* x);

/** 
 * \fn void spgpuSmnrm2(spgpuHandle_t handle, float *y, int n, __device float *x, int count, int pitch)
 * Computes the float precision Euclidean vector norm for every vector in the multivector x. 
 * \param handle the spgpu handle used to call this routine
 * \param y the array of results
 * \param n the vectors' length in the x multivector
 * \param x the input multivector
 */	
void spgpuSmsnrm2(spgpuHandle_t handle,
	float *y,
	int n,
	__device float *x);

/** 
 * \fn void spgpuSmnrm2(spgpuHandle_t handle, CFI_cdesc_t *y, int n, __device float *x, int count, int pitch)
 * Computes the float precision Euclidean vector norm for every vector in the multivector x. 
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array
 * \param n the vectors' length in the x multivector
 * \param x the input multivector
 * \param count the number of vectors in x
 * \param pitch the multivector's pitch
 */	
void spgpuSmvnrm2(spgpuHandle_t handle,
	CFI_cdesc_t *y,
	int n,
	__device float *x,
	int count,
	int pitch);

/** 
 * \fn void spgpuSscal(spgpuHandle_t handle, __device float *y, int n, float alpha, __device float *x)
 * Computes the single precision y = alpha * x. y could be exactly x (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param y the resulting vector
 * \param n the vectors' length
 * \param alpha the alpha value
 * \param x the input vector
 */	
void spgpuSscal(spgpuHandle_t handle,
	__device float *y,
	int n,
	float alpha,
	__device float *x);

/** 
 * \fn void spgpuSaxpby(spgpuHandle_t handle, __device float *z, int n, float beta, __device float *y, float alpha, __device float* x)
 * Computes the single precision z = beta * y + alpha * x. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param y the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 */
void spgpuSaxpby(spgpuHandle_t handle,
	__device float *z,
	int n,
	float beta,
	__device float *y,
	float alpha,
	__device float* x);

/** 
 * \fn void spgpuDaxpbycz(spgpuHandle_t handle, __device float* w, int n, float gamma, __device float* z, float beta, __device float* y, float alpha, __device float* x);
 * Computes the float precision w = gamma * z + beta * y + alpha * x. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param gamma the gamma value
 * \param x the first input vector
 * \param beta the beta value
 * \param y the second input vector
 * \param alpha the alpha value
 * \param x the third input vector
 */
void spgpuSaxpbycz(spgpuHandle_t handle,
	__device float* w,
	int n,
	float gamma,
	__device float* z,
	float beta,
	__device float* y,
	float alpha,
	__device float* x);

/** 
 * \fn void spgpuDcolspan(spgpuHandle_t handle, int n, __device float* x, int countX, int pitchX, CFI_cdesc_t* coeff, __device float* y, int countY, int pitchY, bool updFlag)
 * Computes the float precision y = x * coeff or y = y + x * coeff depending on the updFlag value. y should be a vector if coeff as rank 1 and must be a multivector if coeff is rank 2.
 * \param handle the spgpu handle used to call this routine
 * \param n the vectors' length
 * \param x the input multivector
 * \param countX the number of vectors in x
 * \param pitchX the x multivector pitch
 * \param coeff the CFI descrictor of the fotran array of coefficients
 * \param y the output (multi)vector
 * \param countY the number of vectors in y
 * \param pitchY the y multivector pitch
 * \param updFlag boolen flag: if true, the result will be added to y, otherwise it will overwrite y
 */
void spgpuScolspan(spgpuHandle_t handle,
	int n,
	__device float* x,
	int countX,
	int pitchX,
	CFI_cdesc_t* coeff,
	__device float* y,
	int countY,
	int pitchY,
	bool updFlag);

void spgpuSupd_xyz(spgpuHandle_t handle,
	int n,
	float alpha,
	float beta,
	float gamma,
	float delta,
	__device float* x,
	__device float *y,
	__device float *z);

void spgpuSxyzw(spgpuHandle_t handle,
	int n,
	float a, float b,
	float c, float d,
	float e, float f,
	__device float* x,
	__device float *y,
	__device float *z,
	__device float *w);
	
/** 
 * \fn void spgpuSmaxpby(spgpuHandle_t handle, __device float *z, int n, float beta, __device float *y, float alpha, __device float* x, int count, int pitch)
 * Computes the single precision z = beta * y + alpha * x of x and y multivectors. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param y the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param count the number of vectors in z, x and y multivectors
 * \param pitch the multivectors pitch
 */
void spgpuSmaxpby(spgpuHandle_t handle,
	__device float *z,
	int n,
	float beta,
	__device float *y,
	float alpha,
	__device float* x,
	int count, int pitch);

/** 
 * \fn void spgpuSaxy(spgpuHandle_t handle, __device float *z, int n, float alpha, __device float *x, __device float* y)
 * Computes the single precision z = alpha * x * y. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param alpha the alpha value
 * \param x the first input vector
 * \param y the second input vector
 */
void spgpuSaxy(spgpuHandle_t handle,
	__device float *z,
	int n,
	float alpha,
	__device float *x,
	__device float *y);

/** 
 * \fn void spgpuSaxypbz(spgpuHandle_t handle, __device float *w, int n, float beta, __device float *z, float alpha, __device float* x, __device float *y)
 * Computes the single precision w = beta * z + alpha * x * y. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param z the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the third input vector
 */
void spgpuSaxypbz(spgpuHandle_t handle,
	__device float *w,
	int n,
	float beta,
	__device float *z,
	float alpha,
	__device float* x,
	__device float *y);

/** 
 * \fn void spgpuSmaxy(spgpuHandle_t handle, __device float *z, int n, float alpha, __device float *x, __device float* y, int count, int pitch)
 * Computes the single precision z = alpha * x * y for z, x and y multivectors. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting multivector
 * \param n the vectors' length in the multivectors
 * \param alpha the alpha value
 * \param x the first input multivector
 * \param y the second input multivector
 * \param count the number of vectors in z, x and y multivectors
 * \param pitch the multivectors pitch
 */
void spgpuSmaxy(spgpuHandle_t handle,
	__device float *z,
	int n,
	float alpha,
	__device float* x,
	__device float *y,
	int count,
	int pitch);
	
/** 
 * \fn void spgpuSmaxypbz(spgpuHandle_t handle, __device float *w, int n, float beta, __device float *z, float alpha, __device float* x, __device float *y, int count, int pitch)
 * Computes the single precision w = beta * z + alpha * x * y. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param z the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the third input vector
 * \param count the number of vectors in w, z, x and y multivectors
 * \param pitch the multivectors' pitch
 */	
void spgpuSmaxypbz(spgpuHandle_t handle,
	__device float *w,
	int n,
	float beta,
	__device float *z,
	float alpha,
	__device float* x,
	__device float *y,
	int count,
	int pitch);	
	
/** 
 * \fn void spgpuSgath(spgpuHandle_t handle, __device float *xValues, int xNnz, const __device int *xIndices, int xBaseIndex, const __device float* y)
 * Single precision gather from y to sparse(x). Computes the single precision gather from y to xValues (using xIndices).
 * \param handle the spgpu handle used to call this routine
 * \param xValues the destination array for gathered values
 * \param xNnz the number of elements to gather
 * \param xIndices the array of indices for the elements to be gathered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param y the source vector (from which the elements will be gathered)
 */
void spgpuSgath(spgpuHandle_t handle,
	__device float *xValues,
	int xNnz,
	const __device int *xIndices,
	int xBaseIndex,
	const __device float* y);
	
/** 
 * \fn void spgpuSscat(spgpuHandle_t handle, __device float* y, int xNnz, const __device float *xValues, const __device int *xIndices, int xBaseIndex, float beta)
 * Single precision scatter from sparse(x) to y. Computes the single precision scatter from xValues to y (using xIndices).
 * The scattered element will be, for i in [0, xNnz), y[xIndices[i]] = beta*y[xIndices[i]] + xValues[i] (to be noted that
 * y values will be multiplied with beta just for scattered values).
 * \param handle the spgpu handle used to call this routine
 * \param y the destination vector (to which the elements will be scattered)
 * \param xNnz the number of elements to scatter
 * \param xValues the source array from which the values will be read
 * \param xIndices the array of indices for the elements to be scattered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param beta the beta value
 */
void spgpuSscat(spgpuHandle_t handle,
	__device float* y,
	int xNnz,
	const __device float *xValues,
	const __device int *xIndices,
	int xBaseIndex, float beta);	

float spgpuSasum(spgpuHandle_t handle,
	int n,
	float* x);

float spgpuSamax(spgpuHandle_t handle,
	int n,
	float* x);

void spgpuSmasum(spgpuHandle_t handle,
	float* y,
	int n,
	float* x,
	int count,
	int pitch);

void spgpuSmamax(spgpuHandle_t handle,
	float* y,
	int n,
	float* x,
	int count,
	int pitch);

/** 
 * \fn double spgpuDdot (spgpuHandle_t handle, int n, __device double* a, __device double* b)
 * Computes double precision dot product of a and b vectors.
 * \param handle The spgpu handle used to call this routine
 * \param n the vectors length
 * \param a the first input vector
 * \param b the second input vector
 * \return the dot product
 */
 double spgpuDdot(spgpuHandle_t handle,
	int n,
	__device double* a,
	__device double* b);

/** 
 * \fn void spgpuDmdot (spgpuHandle_t handle, double* y, int n, __device double* a, __device double* b, int count, int pitch)
 * Computes double precision dot product of a and b multivectors.
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector couples from the multivectors a and b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param count the number of vectors in every multivector
 * \param pitch the pitch, in number of elements, of every multivectors (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuDmdot(spgpuHandle_t handle,
	double* y,
	int n,
	__device double* a,
	__device double* b,
	int count,
	int pitch);

/** 
 * \fn void spgpuDmdot (spgpuHandle_t handle, double* y, int n, __device double* a, __device double* b)
 * Computes double precision dot product of a and b vectors (count assumed = 1 for both).
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of a and b
 * \param n the vectors' length
 * \param a the first input vector
 * \param b the second input vector
 */
void spgpuDmsdot(spgpuHandle_t handle,
	double* y,
	int n,
	__device double* a,
	__device double* b);

/** 
 * \fn void spgpuDmvdot (spgpuHandle_t handle, double* y, int n, __device double* a, __device double* b, int countA, int pitchA)
 * Computes double precision dot product of multivector a and vector b (count assumed = 1 for b).
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector in multivector a with the vector b
 * \param n the vectors' length
 * \param a the input multivector
 * \param b the input vector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuDmvdot(spgpuHandle_t handle,
	double* y,
	int n,
	__device double* a,
	__device double* b,
	int countA,
	int pitchA);

/** 
 * \fn void spgpuDmmdot (spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device double* a, __device double* b, int countA, int pitchA, int countB, int pitchB)
 * Computes double precision dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array, made by dot products of every vector in multivector a with the vector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuDmvdot_CFI(spgpuHandle_t handle,
	CFI_cdesc_t* y,
	int n,
	__device double* a,
	__device double* b,
	int countA,
	int pitchA);

/** 
 * \fn void spgpuDmmdot (spgpuHandle_t handle, double* y, int n, __device double* a, __device double* b, int countA, int pitchA, int countB, int pitchB)
 * Computes double precision dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector in multivector a with every vector in multivector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param countB the number of vectors in multivector b
 * \param pitchB the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 */
void spgpuDmmdot(spgpuHandle_t handle,
	double* y,
	int n,
	__device double* a,
	__device double* b,
	int countA,
	int pitchA,
	int countB,
	int pitchB);

/** 
 * \fn void spgpuDmmdot (spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device double* a, __device double* b, int countA, int pitchA, int countB, int pitchB)
 * Computes double precision dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array, made by dot products of every vector in multivector a with every vector in multivector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param countB the number of vectors in multivector b
 * \param pitchB the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 */
void spgpuDmmdot_CFI(spgpuHandle_t handle,
	CFI_cdesc_t* y,
	int n,
	__device double* a,
	__device double* b,
	int countA,
	int pitchA,
	int countB,
	int pitchB);

/** 
 * \fn void spgpuDmlt (spgpuHandle_t handle, __device double* z, int countZ, int pitchZ,
 * 						int n, double alpha, __device double* x, int countX, int pitchY,
 * 											 __device double* y, int countY, int pitchY,
 * 								double beta, bool conjgX, bool conjgY)
 * Computes z = alpha * op(x) * op(y) + beta * z with x, y, z vectors/multivectors possibly using multiple column both in input or output
 * \param handle the spgpu handle used to call this routine
 * \param z the output multivector
 * \param countX the number of vectors in multivector Z
 * \param pitchX the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param n the vectors' length
 * \param alpha the first scalar coefficient
 * \param y the second input multivector
 * \param countY the number of vectors in multivector y
 * \param pitchY the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 * \param x the first input multivector
 * \param countX the number of vectors in multivector x
 * \param pitchX the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param beta the second scalar coefficient
 * \param congjx boolean flag 
 * \param congjy boolean flag 
 */
void spgpuDmlt(spgpuHandle_t handle,
	__device double* z,
	int countZ,
	int pitchZ,
	int n,
	double alpha,
	__device double* x,
	int countX,
	int pitchX,
	__device double* y,
	int countY,
	int pitchY,
	double beta,
	bool conjgX,
	bool conjgY);

/** 
 * \fn void spgpuDabs (spgpuHandle_t handle, __device double* y, int n, double alpha, __device double* x)
 * Computes double precision (y = alpha * abs(x)) for each element in x.
 * \param handle The spgpu handle used to call this routine
 * \param y the resulting vector (could be x)
 * \param n the vectors length
 * \param alpha the alpha value
 * \param x the input vector
 */
void spgpuDabs(spgpuHandle_t handle,
	__device double* y,
	int n,
	double alpha,
	__device double* x);
	
/** 
 * \fn double spgpuDnrm2(spgpuHandle_t handle, int n, __device double* x)
 * Computes the double precision Euclidean vector norm of x. 
 * \param handle the spgpu handle used to call this routine
 * \param n the vector's length
 * \param x the input vector
 * \return the euclidean vector norm
 */
 double spgpuDnrm2(spgpuHandle_t handle,
	int n,
	__device double* x);

/** 
 * \fn void spgpuDmnrm2(spgpuHandle_t handle, double *y, int n, __device double *x, int count, int pitch)
 * Computes the double precision Euclidean vector norm for every vector in the multivector x. 
 * \param handle the spgpu handle used to call this routine
 * \param y the array of results
 * \param n the vectors' length in the x multivector
 * \param x the input multivector
 */	
void spgpuDmsnrm2(spgpuHandle_t handle,
	double *y,
	int n,
	__device double *x);

/** 
 * \fn void spgpuDmnrm2(spgpuHandle_t handle, CFI_cdesc_t *y, int n, __device double *x, int count, int pitch)
 * Computes the double precision Euclidean vector norm for every vector in the multivector x. 
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array
 * \param n the vectors' length in the x multivector
 * \param x the input multivector
 * \param count the number of vectors in x
 * \param pitch the multivector's pitch
 */	
void spgpuDmvnrm2(spgpuHandle_t handle,
	CFI_cdesc_t *y,
	int n,
	__device double *x,
	int count,
	int pitch);
	
/** 
 * \fn void spgpuDscal(spgpuHandle_t handle, __device double *y, int n, double alpha, __device double *x)
 * Computes the Double precision y = alpha * x. y could be exactly x (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param y the resulting vector
 * \param n the vectors' length
 * \param alpha the alpha value
 * \param x the input vector
 */
void spgpuDscal(spgpuHandle_t handle,
	__device double *y,
	int n,
	double alpha,
	__device double *x);

/** 
 * \fn void spgpuDaxpby(spgpuHandle_t handle, __device double *z, int n, double beta, __device double *y, double alpha, __device double* x)
 * Computes the double precision z = beta * y + alpha * x. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param y the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 */
void spgpuDaxpby(spgpuHandle_t handle,
	__device double *z,
	int n,
	double beta,
	__device double *y,
	double alpha,
	__device double* x);

/** 
 * \fn void spgpuDaxpbycz(spgpuHandle_t handle, __device double* w, int n, double gamma, __device double* z, double beta, __device double* y, double alpha, __device double* x);
 * Computes the double precision w = gamma * z + beta * y + alpha * x. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param gamma the gamma value
 * \param x the first input vector
 * \param beta the beta value
 * \param y the second input vector
 * \param alpha the alpha value
 * \param x the third input vector
 */
void spgpuDaxpbycz(spgpuHandle_t handle,
	__device double* w,
	int n,
	double gamma,
	__device double* z,
	double beta,
	__device double* y,
	double alpha,
	__device double* x);

/** 
 * \fn void spgpuDcolspan(spgpuHandle_t handle, int n, __device double* x, int countX, int pitchX, CFI_cdesc_t* coeff, __device double* y, int countY, int pitchY, bool updFlag)
 * Computes the double precision y = x * coeff or y = y + x * coeff depending on the updFlag value. y should be a vector if coeff as rank 1 and must be a multivector if coeff is rank 2.
 * \param handle the spgpu handle used to call this routine
 * \param n the vectors' length
 * \param x the input multivector
 * \param countX the number of vectors in x
 * \param pitchX the x multivector pitch
 * \param coeff the CFI descrictor of the fotran array of coefficients
 * \param y the output (multi)vector
 * \param countY the number of vectors in y
 * \param pitchY the y multivector pitch
 * \param updFlag boolen flag: if true, the result will be added to y, otherwise it will overwrite y
 */
void spgpuDcolspan(spgpuHandle_t handle,
	int n,
	__device double* x,
	int countX,
	int pitchX,
	CFI_cdesc_t* coeff,
	__device double* y,
	int countY,
	int pitchY,
	bool updFlag);

void spgpuDupd_xyz(spgpuHandle_t handle,
	int n,
	double alpha,
	double beta,
	double gamma,
	double delta,
	__device double* x,
	__device double *y,
	__device double *z);
	
void spgpuDxyzw(spgpuHandle_t handle,
	int n,
	double a, double b,
	double c, double d,
	double e, double f,
	__device double* x,
	__device double *y,
	__device double *z,
	__device double *w);

  /** 
 * \fn void spgpuDmaxpby(spgpuHandle_t handle, __device double *z, int n, double beta, __device double *y, double alpha, __device double* x, int count, int pitch)
 * Computes the double precision z = beta * y + alpha * x of x and y multivectors. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param y the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param count the number of vectors in z, x and y multivectors
 * \param pitch the multivector's pitch
 */
void spgpuDmaxpby(spgpuHandle_t handle,
	__device double *z,
	int n,
	double beta,
	__device double *y,
	double alpha,
	__device double* x,
	int count, int pitch);

/** 
 * \fn void spgpuDaxy(spgpuHandle_t handle, __device double *z, int n, double alpha, __device double *x, __device double* y)
 * Computes the double precision z = alpha * x * y. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length in the x multivector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the first input vector
 */
void spgpuDaxy(spgpuHandle_t handle,
	__device double *z,
	int n,
	double alpha,
	__device double *x,
	__device double *y);

/** 
 * \fn void spgpuDaxypbz(spgpuHandle_t handle, __device double *w, int n, double beta, __device double *z, double alpha, __device double* x, __device double *y)
 * Computes the double precision w = beta * z + alpha * x * y. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param z the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the third input vector
 */
void spgpuDaxypbz(spgpuHandle_t handle,
	__device double *w,
	int n,
	double beta,
	__device double *z,
	double alpha,
	__device double* x,
	__device double *y);

/** 
 * \fn void spgpuDmaxy(spgpuHandle_t handle, __device double *z, int n, double alpha, __device double *x, __device double* y, int count, int pitch)
 * Computes the double precision z = alpha * x * y for z, x and y multivectors. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting multivector
 * \param n the vectors' length in the multivectors
 * \param alpha the alpha value
 * \param x the first input multivector
 * \param y the second input multivector
 * \param count the number of vectors in z, x and y multivectors
 * \param pitch the multivectors pitch
 */	
void spgpuDmaxy(spgpuHandle_t handle,
	__device double *z,
	int n,
	double alpha,
	__device double* x,
	__device double *y,
	int count,
	int pitch);
	
/** 
 * \fn void spgpuDmaxypbz(spgpuHandle_t handle, __device double *w, int n, double beta, __device double *z, double alpha, __device double* x, __device double *y, int count, int pitch)
 * Computes the double precision w = beta * z + alpha * x * y. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param z the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the third input vector
 * \param count the number of vectors in w, z, x and y multivectors
 * \param pitch the multivectors' pitch
 */
void spgpuDmaxypbz(spgpuHandle_t handle,
	__device double *w,
	int n,
	double beta,
	__device double *z,
	double alpha,
	__device double* x,
	__device double *y,
	int count,
	int pitch);

/** 
 * \fn void spgpuDgath(spgpuHandle_t handle, __device double *xValues, int xNnz, const __device int *xIndices, int xBaseIndex, const __device double* y)
 * Double precision gather from y to sparse(x). Computes the double precision gather from y to xValues (using xIndices).
 * \param handle the spgpu handle used to call this routine
 * \param xValues the destination array for gathered values
 * \param xNnz the number of elements to gather
 * \param xIndices the array of indices for the elements to be gathered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param y the source vector (from which the elements will be gathered)
 */
void spgpuDgath(spgpuHandle_t handle,
	__device double *xValues,
	int xNnz,
	const __device int *xIndices,
	int xBaseIndex,
	const __device double* y);
	
/** 
 * \fn void spgpuDscat(spgpuHandle_t handle, __device double* y, int xNnz, const __device double *xValues, const __device int *xIndices, int xBaseIndex, double beta)
 * Double precision scatter from sparse(x) to y. Computes the single precision scatter from xValues to y (using xIndices).
 * The scattered element will be, for i in [0, xNnz), y[xIndices[i]] = beta*y[xIndices[i]] + xValues[i] (to be noted that
 * y values will be multiplied with beta just for scattered values).
 * \param handle the spgpu handle used to call this routine
 * \param y the destination vector (to which the elements will be scattered)
 * \param xNnz the number of elements to scatter
 * \param xValues the source array from which the values will be read
 * \param xIndices the array of indices for the elements to be scattered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param beta the beta value
 */
void spgpuDscat(spgpuHandle_t handle,
	__device double* y,
	int xNnz,
	const __device double *xValues,
	const __device int *xIndices,
	int xBaseIndex, double beta);

double spgpuDasum(spgpuHandle_t handle,
	int n,
	double* x);

double spgpuDamax(spgpuHandle_t handle,
	int n,
	double* x);

void spgpuDmasum(spgpuHandle_t handle,
	double* y,
	int n,
	double* x,
	int count,
	int pitch);

void spgpuDmamax(spgpuHandle_t handle,
	double* y,
	int n,
	double* x,
	int count,
	int pitch);
	
/** 
 * \fn cuFloatComplex spgpuCdot (spgpuHandle_t handle, int n, __device cuFloatComplex* a, __device cuFloatComplex* b)
 * Computes single precision complex dot product of a and b vectors.
 * \param handle The spgpu handle used to call this routine
 * \param n the vectors length
 * \param a the first input vector
 * \param b the second input vector
 * \return the dot product
 */
cuFloatComplex spgpuCdot(spgpuHandle_t handle,
	int n,
	__device cuFloatComplex* a,
	__device cuFloatComplex* b);

/** 
 * \fn void spgpuCmdot (spgpuHandle_t handle, cuFloatComplex* y, int n, __device cuFloatComplex* a, __device cuFloatComplex* b, int count, int pitch)
 * Computes single precision complex dot product of a and b multivectors.
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector couples from the multivectors a and b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param count the number of vectors in every multivector
 * \param pitch the pitch, in number of elements, of every multivectors (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuCmdot(spgpuHandle_t handle,
	cuFloatComplex* y,
	int n,
	__device cuFloatComplex* a,
	__device cuFloatComplex* b,
	int count,
	int pitch);

/** 
 * \fn void spgpuCmdot (spgpuHandle_t handle, cuFloatComplex* y, int n, __device cuFloatComplex* a, __device cuFloatComplex* b)
 * Computes single precision complex dot product of a and b vectors (count assumed = 1 for both).
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of a and b
 * \param n the vectors' length
 * \param a the first input vector
 * \param b the second input vector
 */
void spgpuCmsdot(spgpuHandle_t handle,
	cuFloatComplex* y,
	int n,
	__device cuFloatComplex* a,
	__device cuFloatComplex* b);

/** 
 * \fn void spgpuCmvdot (spgpuHandle_t handle, cuFloatComplex* y, int n, __device cuFloatComplex* a, __device cuFloatComplex* b, int countA, int pitchA)
 * Computes single precision complex dot product of multivector a and vector b (count assumed = 1 for b).
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector in multivector a with the vector b
 * \param n the vectors' length
 * \param a the input multivector
 * \param b the input vector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuCmvdot(spgpuHandle_t handle,
	cuFloatComplex* y,
	int n,
	__device cuFloatComplex* a,
	__device cuFloatComplex* b,
	int countA,
	int pitchA);

/** 
 * \fn void spgpuCmmdot (spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device cuFloatComplex* a, __device cuFloatComplex* b, int countA, int pitchA, int countB, int pitchB)
 * Computes single precision complex dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array, made by dot products of every vector in multivector a with the vector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuCmvdot_CFI(spgpuHandle_t handle,
	CFI_cdesc_t* y,
	int n,
	__device cuFloatComplex* a,
	__device cuFloatComplex* b,
	int countA,
	int pitchA);

/** 
 * \fn void spgpuCmmdot (spgpuHandle_t handle, cuFloatComplex* y, int n, __device cuFloatComplex* a, __device cuFloatComplex* b, int countA, int pitchA, int countB, int pitchB)
 * Computes single precision complex dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector in multivector a with every vector in multivector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param countB the number of vectors in multivector b
 * \param pitchB the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 */
void spgpuCmmdot(spgpuHandle_t handle,
	cuFloatComplex* y,
	int n,
	__device cuFloatComplex* a,
	__device cuFloatComplex* b,
	int countA,
	int pitchA,
	int countB,
	int pitchB);

/** 
 * \fn void spgpuCmmdot (spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device cuFloatComplex* a, __device cuFloatComplex* b, int countA, int pitchA, int countB, int pitchB)
 * Computes single precision complex dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array, made by dot products of every vector in multivector a with every vector in multivector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param countB the number of vectors in multivector b
 * \param pitchB the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 */
void spgpuCmmdot_CFI(spgpuHandle_t handle,
	CFI_cdesc_t* y,
	int n,
	__device cuFloatComplex* a,
	__device cuFloatComplex* b,
	int countA,
	int pitchA,
	int countB,
	int pitchB);

/** 
 * \fn void spgpuCmlt (spgpuHandle_t handle, __device cuFloatComplex* z, int countZ, int pitchZ,
 * 						int n, cuFloatComplex alpha, __device cuFloatComplex* x, int countX, int pitchY,
 * 											 __device cuFloatComplex* y, int countY, int pitchY,
 * 								cuFloatComplex beta, bool conjgX, bool conjgY)
 * Computes z = alpha * op(x) * op(y) + beta * z with x, y, z vectors/multivectors possibly using multiple column both in input or output
 * \param handle the spgpu handle used to call this routine
 * \param z the output multivector
 * \param countX the number of vectors in multivector Z
 * \param pitchX the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param n the vectors' length
 * \param alpha the first scalar coefficient
 * \param y the second input multivector
 * \param countY the number of vectors in multivector y
 * \param pitchY the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 * \param x the first input multivector
 * \param countX the number of vectors in multivector x
 * \param pitchX the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param beta the second scalar coefficient
 * \param congjx boolean flag 
 * \param congjy boolean flag 
 */
void spgpuCmlt(spgpuHandle_t handle,
	__device cuFloatComplex* z,
	int countZ,
	int pitchZ,
	int n,
	cuFloatComplex alpha,
	__device cuFloatComplex* x,
	int countX,
	int pitchX,
	__device cuFloatComplex* y,
	int countY,
	int pitchY,
	cuFloatComplex beta,
	bool conjgX,
	bool conjgY);

/** 
 * \fn void spgpuCabs (spgpuHandle_t handle, __device float* y, int n, float alpha, __device cuFloatComplex* x)
 * Computes single precision complex (y = alpha * abs(x)) for each element in x.
 * \param handle The spgpu handle used to call this routine
 * \param y the resulting vector (could be x)
 * \param n the vectors length
 * \param alpha the alpha value
 * \param x the input vector
 */
void spgpuCabs(spgpuHandle_t handle,
	__device cuFloatComplex* y,
	int n,
        cuFloatComplex alpha,
	__device cuFloatComplex* x);
	
/** 
 * \fn float spgpuCnrm2(spgpuHandle_t handle, int n, __device cuFloatComplex* x)
 * Computes the single precision complex Euclidean vector norm of x. 
 * \param handle the spgpu handle used to call this routine
 * \param n the vector's length
 * \param x the input vector
 * \return the euclidean vector norm
 */
float spgpuCnrm2(spgpuHandle_t handle,
	int n,
	__device cuFloatComplex* x);

/** 
 * \fn void spgpuCmnrm2(spgpuHandle_t handle, float *y, int n, __device cuFloatComplex *x, int count, int pitch)
 * Computes the single precision complex Euclidean vector norm for every vector in the multivector x. 
 * \param handle the spgpu handle used to call this routine
 * \param y the array of results
 * \param n the vectors' length in the x multivector
 * \param x the input multivector
 */	
void spgpuCmsnrm2(spgpuHandle_t handle,
	float *y,
	int n,
	__device cuFloatComplex *x);

/** 
 * \fn void spgpuCmnrm2(spgpuHandle_t handle, CFI_cdesc_t *y, int n, __device cuFloatComplex *x, int count, int pitch)
 * Computes the single precision complex Euclidean vector norm for every vector in the multivector x. 
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array
 * \param n the vectors' length in the x multivector
 * \param x the input multivector
 * \param count the number of vectors in x
 * \param pitch the multivector's pitch
 */	
void spgpuCmvnrm2(spgpuHandle_t handle,
	CFI_cdesc_t *y,
	int n,
	__device cuFloatComplex *x,
	int count,
	int pitch);

/** 
 * \fn void spgpuCscal(spgpuHandle_t handle, __device cuFloatComplex *y, int n, cuFloatComplex alpha, __device cuFloatComplex *x)
 * Computes the single precision complex y = alpha * x. y could be exactly x (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param y the resulting vector
 * \param n the vectors' length
 * \param alpha the alpha value
 * \param x the input vector
 */	
void spgpuCscal(spgpuHandle_t handle,
	__device cuFloatComplex *y,
	int n,
	cuFloatComplex alpha,
	__device cuFloatComplex *x);

/** 
 * \fn void spgpuCaxpby(spgpuHandle_t handle, __device cuFloatComplex *z, int n, cuFloatComplex beta, __device cuFloatComplex *y, cuFloatComplex alpha, __device cuFloatComplex* x)
 * Computes the single precision complex z = beta * y + alpha * x. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param y the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 */
void spgpuCaxpby(spgpuHandle_t handle,
	__device cuFloatComplex *z,
	int n,
	cuFloatComplex beta,
	__device cuFloatComplex *y,
	cuFloatComplex alpha,
	__device cuFloatComplex* x);

/** 
 * \fn void spgpuCcolspan(spgpuHandle_t handle, int n, __device cuFloatComplex* x, int countX, int pitchX, CFI_cdesc_t* coeff, __device cuFloatComplex* y, int countY, int pitchY, bool updFlag)
 * Computes the single precision complex y = x * coeff or y = y + x * coeff depending on the updFlag value. y should be a vector if coeff as rank 1 and must be a multivector if coeff is rank 2.
 * \param handle the spgpu handle used to call this routine
 * \param n the vectors' length
 * \param x the input multivector
 * \param countX the number of vectors in x
 * \param pitchX the x multivector pitch
 * \param coeff the CFI descrictor of the fotran array of coefficients
 * \param y the output (multi)vector
 * \param countY the number of vectors in y
 * \param pitchY the y multivector pitch
 * \param updFlag boolen flag: if true, the result will be added to y, otherwise it will overwrite y
 */
void spgpuCcolspan(spgpuHandle_t handle,
	int n,
	__device cuFloatComplex* x,
	int countX,
	int pitchX,
	CFI_cdesc_t* coeff,
	__device cuFloatComplex* y,
	int countY,
	int pitchY,
	bool updFlag);

void spgpuCupd_xyz(spgpuHandle_t handle,
	int n,
	cuFloatComplex alpha,
	cuFloatComplex beta,
	cuFloatComplex gamma,
	cuFloatComplex delta,
	__device cuFloatComplex* x,
	__device cuFloatComplex *y,
	__device cuFloatComplex *z);

void spgpuCxyzw(spgpuHandle_t handle,
	int n,
	cuFloatComplex a, cuFloatComplex b,
	cuFloatComplex c, cuFloatComplex d,
	cuFloatComplex e, cuFloatComplex f,
	__device cuFloatComplex* x,
	__device cuFloatComplex *y,
	__device cuFloatComplex *z,
	__device cuFloatComplex *w);
	
/** 
 * \fn void spgpuCmaxpby(spgpuHandle_t handle, __device cuFloatComplex *z, int n, cuFloatComplex beta, __device cuFloatComplex *y, cuFloatComplex alpha, __device cuFloatComplex* x, int count, int pitch)
 * Computes the single precision complex z = beta * y + alpha * x of x and y multivectors. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param y the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param count the number of vectors in z, x and y multivectors
 * \param pitch the multivectors pitch
 */

void spgpuCmaxpby(spgpuHandle_t handle,
	__device cuFloatComplex *z,
	int n,
	cuFloatComplex beta,
	__device cuFloatComplex *y,
	cuFloatComplex alpha,
	__device cuFloatComplex* x,
	int count, int pitch);

/** 
 * \fn void spgpuCaxy(spgpuHandle_t handle, __device cuFloatComplex *z, int n, cuFloatComplex alpha, __device cuFloatComplex *x, __device cuFloatComplex* y)
 * Computes the single precision complex z = alpha * x * y. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param alpha the alpha value
 * \param x the first input vector
 * \param y the second input vector
 */
void spgpuCaxy(spgpuHandle_t handle,
	__device cuFloatComplex *z,
	int n,
	cuFloatComplex alpha,
	__device cuFloatComplex *x,
	__device cuFloatComplex *y);

/** 
 * \fn void spgpuCaxypbz(spgpuHandle_t handle, __device cuFloatComplex *w, int n, cuFloatComplex beta, __device cuFloatComplex *z, cuFloatComplex alpha, __device cuFloatComplex* x, __device cuFloatComplex *y)
 * Computes the single precision complex w = beta * z + alpha * x * y. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param z the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the third input vector
 */
void spgpuCaxypbz(spgpuHandle_t handle,
	__device cuFloatComplex *w,
	int n,
	cuFloatComplex beta,
	__device cuFloatComplex *z,
	cuFloatComplex alpha,
	__device cuFloatComplex* x,
	__device cuFloatComplex *y);

/** 
 * \fn void spgpuCmaxy(spgpuHandle_t handle, __device cuFloatComplex *z, int n, cuFloatComplex alpha, __device cuFloatComplex *x, __device cuFloatComplex* y, int count, int pitch)
 * Computes the single precision complex z = alpha * x * y for z, x and y multivectors. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting multivector
 * \param n the vectors' length in the multivectors
 * \param alpha the alpha value
 * \param x the first input multivector
 * \param y the second input multivector
 * \param count the number of vectors in z, x and y multivectors
 * \param pitch the multivectors pitch
 */
void spgpuCmaxy(spgpuHandle_t handle,
	__device cuFloatComplex *z,
	int n,
	cuFloatComplex alpha,
	__device cuFloatComplex* x,
	__device cuFloatComplex *y,
	int count,
	int pitch);
	
/** 
 * \fn void spgpuCmaxypbz(spgpuHandle_t handle, __device cuFloatComplex *w, int n, cuFloatComplex beta, __device cuFloatComplex *z, cuFloatComplex alpha, __device cuFloatComplex* x, __device cuFloatComplex *y, int count, int pitch)
 * Computes the single precision complex w = beta * z + alpha * x * y. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param z the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the third input vector
 * \param count the number of vectors in w, z, x and y multivectors
 * \param pitch the multivectors' pitch
 */
void spgpuCmaxypbz(spgpuHandle_t handle,
	__device cuFloatComplex *w,
	int n,
	cuFloatComplex beta,
	__device cuFloatComplex *z,
	cuFloatComplex alpha,
	__device cuFloatComplex* x,
	__device cuFloatComplex *y,
	int count,
	int pitch);	
	
/** 
 * \fn void spgpuCgath(spgpuHandle_t handle, __device cuFloatComplex *xValues, int xNnz, const __device int *xIndices, int xBaseIndex, const __device cuFloatComplex* y)
 * Single precision complex gather from y to sparse(x). Computes the single precision complex gather from y to xValues (using xIndices).
 * \param handle the spgpu handle used to call this routine
 * \param xValues the destination array for gathered values
 * \param xNnz the number of elements to gather
 * \param xIndices the array of indices for the elements to be gathered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param y the source vector (from which the elements will be gathered)
 */
void spgpuCgath(spgpuHandle_t handle,
	__device cuFloatComplex *xValues,
	int xNnz,
	const __device int *xIndices,
	int xBaseIndex,
	const __device cuFloatComplex* y);
	
/** 
 * \fn void spgpuCscat(spgpuHandle_t handle, __device cuFloatComplex* y, int xNnz, const __device cuFloatComplex *xValues, const __device int *xIndices, int xBaseIndex, cuFloatComplex beta)
 * Single precision complex scatter from sparse(x) to y. Computes the single precision complex scatter from xValues to y (using xIndices).
 * The scattered element will be, for i in [0, xNnz), y[xIndices[i]] = beta*y[xIndices[i]] + xValues[i] (to be noted that
 * y values will be multiplied with beta just for scattered values).
 * \param handle the spgpu handle used to call this routine
 * \param y the destination vector (to which the elements will be scattered)
 * \param xNnz the number of elements to scatter
 * \param xValues the source array from which the values will be read
 * \param xIndices the array of indices for the elements to be scattered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param beta the beta value
 */
void spgpuCscat(spgpuHandle_t handle,
	__device cuFloatComplex* y,
	int xNnz,
	const __device cuFloatComplex *xValues,
	const __device int *xIndices,
	int xBaseIndex, cuFloatComplex beta);	

float spgpuCasum(spgpuHandle_t handle,
	int n,
	cuFloatComplex* x);

float spgpuCamax(spgpuHandle_t handle,
	int n,
	cuFloatComplex* x);

void spgpuCmasum(spgpuHandle_t handle,
	float* y,
	int n,
	cuFloatComplex* x,
	int count,
	int pitch);

void spgpuCmamax(spgpuHandle_t handle,
	float* y,
	int n,
	cuFloatComplex* x,
	int count,
	int pitch);
	
/** 
 * \fn void spgpuZscal(spgpuHandle_t handle, __device cuDoubleComplex *y, int n, cuDoubleComplex alpha, __device cuDoubleComplex *x)
 * Computes the double precision complex y = alpha * x. y could be exactly x (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param y the resulting vector
 * \param n the vectors' length
 * \param alpha the alpha value
 * \param x the input vector
 */
void spgpuZscal(spgpuHandle_t handle,
	__device cuDoubleComplex *y,
	int n,
	cuDoubleComplex alpha,
	__device cuDoubleComplex *x);

/** 
 * \fn cuFloatComplex spgpuZdot (spgpuHandle_t handle, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b)
 * Computes double precision complex dot product of a and b vectors.
 * \param handle The spgpu handle used to call this routine
 * \param n the vectors length
 * \param a the first input vector
 * \param b the second input vector
 * \return the dot product
 */
 cuDoubleComplex spgpuZdot(spgpuHandle_t handle,
	int n,
	__device cuDoubleComplex* a,
	__device cuDoubleComplex* b);

/** 
 * \fn void spgpuZmdot (spgpuHandle_t handle, cuDoubleComplex* y, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b, int count, int pitch)
 * Computes double precision complex dot product of a and b multivectors.
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector couples from the multivectors a and b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param count the number of vectors in every multivector
 * \param pitch the pitch, in number of elements, of every multivectors (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuZmdot(spgpuHandle_t handle,
	cuDoubleComplex* y,
	int n,
	__device cuDoubleComplex* a,
	__device cuDoubleComplex* b,
	int count,
	int pitch);

/** 
 * \fn void spgpuZmdot (spgpuHandle_t handle, cuDoubleComplex* y, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b)
 * Computes single precision complex dot product of a and b vectors (count assumed = 1 for both).
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of a and b
 * \param n the vectors' length
 * \param a the first input vector
 * \param b the second input vector
 */
void spgpuZmsdot(spgpuHandle_t handle,
	cuDoubleComplex* y,
	int n,
	__device cuDoubleComplex* a,
	__device cuDoubleComplex* b);

/** 
 * \fn void spgpuZmvdot (spgpuHandle_t handle, cuDoubleComplex* y, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b, int countA, int pitchA)
 * Computes single precision complex dot product of multivector a and vector b (count assumed = 1 for b).
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector in multivector a with the vector b
 * \param n the vectors' length
 * \param a the input multivector
 * \param b the input vector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuZmvdot(spgpuHandle_t handle,
	cuDoubleComplex* y,
	int n,
	__device cuDoubleComplex* a,
	__device cuDoubleComplex* b,
	int countA,
	int pitchA);

/** 
 * \fn void spgpuZmmdot (spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b, int countA, int pitchA, int countB, int pitchB)
 * Computes single precision complex dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array, made by dot products of every vector in multivector a with the vector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 */
void spgpuZmvdot_CFI(spgpuHandle_t handle,
	CFI_cdesc_t* y,
	int n,
	__device cuDoubleComplex* a,
	__device cuDoubleComplex* b,
	int countA,
	int pitchA);

/** 
 * \fn void spgpuZmmdot (spgpuHandle_t handle, cuDoubleComplex* y, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b, int countA, int pitchA, int countB, int pitchB)
 * Computes single precision complex dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the result, made by dot products of every vector in multivector a with every vector in multivector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param countB the number of vectors in multivector b
 * \param pitchB the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 */
void spgpuZmmdot(spgpuHandle_t handle,
	cuDoubleComplex* y,
	int n,
	__device cuDoubleComplex* a,
	__device cuDoubleComplex* b,
	int countA,
	int pitchA,
	int countB,
	int pitchB);

/** 
 * \fn void spgpuZmmdot (spgpuHandle_t handle, CFI_cdesc_t* y, int n, __device cuDoubleComplex* a, __device cuDoubleComplex* b, int countA, int pitchA, int countB, int pitchB)
 * Computes single precision complex dot product of multivector a and multivector b.
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array, made by dot products of every vector in multivector a with every vector in multivector b
 * \param n the vectors' length
 * \param a the first input multivector
 * \param b the second input multivector
 * \param countA the number of vectors in multivector a
 * \param pitchA the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param countB the number of vectors in multivector b
 * \param pitchB the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 */
void spgpuZmmdot_CFI(spgpuHandle_t handle,
	CFI_cdesc_t* y,
	int n,
	__device cuDoubleComplex* a,
	__device cuDoubleComplex* b,
	int countA,
	int pitchA,
	int countB,
	int pitchB);

/** 
 * \fn void spgpuZmlt (spgpuHandle_t handle, __device cuDoubleComplex* z, int countZ, int pitchZ,
 * 						int n, cuDoubleComplex alpha, __device cuDoubleComplex* x, int countX, int pitchY,
 * 											 __device cuDoubleComplex* y, int countY, int pitchY,
 * 								cuDoubleComplex beta, bool conjgX, bool conjgY)
 * Computes z = alpha * op(x) * op(y) + beta * z with x, y, z vectors/multivectors possibly using multiple column both in input or output
 * \param handle the spgpu handle used to call this routine
 * \param z the output multivector
 * \param countX the number of vectors in multivector Z
 * \param pitchX the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param n the vectors' length
 * \param alpha the first scalar coefficient
 * \param y the second input multivector
 * \param countY the number of vectors in multivector y
 * \param pitchY the pitch, in number of elements, of multivector b (so the second element of the first vector in b will be b[pitch], the third b[2*pitch], etc.).
 * \param x the first input multivector
 * \param countX the number of vectors in multivector x
 * \param pitchX the pitch, in number of elements, of multivector a (so the second element of the first vector in a will be a[pitch], the third a[2*pitch], etc.).
 * \param beta the second scalar coefficient
 * \param congjx boolean flag 
 * \param congjy boolean flag 
 */
void spgpuZmlt(spgpuHandle_t handle,
	__device cuDoubleComplex* z,
	int countZ,
	int pitchZ,
	int n,
	cuDoubleComplex alpha,
	__device cuDoubleComplex* x,
	int countX,
	int pitchX,
	__device cuDoubleComplex* y,
	int countY,
	int pitchY,
	cuDoubleComplex beta,
	bool conjgX,
	bool conjgY);

/**
 * \fn void spgpuZabs (spgpuHandle_t handle, __device double* y, int n, double alpha, __device cuDoubleComplex* x)
 * Computes double precision complex (y = alpha * abs(x)) for each element in x.
 * \param handle The spgpu handle used to call this routine
 * \param y the resulting vector (could be x)
 * \param n the vectors length
 * \param alpha the alpha value
 * \param x the input vector
 */
void spgpuZabs(spgpuHandle_t handle,
	__device cuDoubleComplex* y,
	int n,
	cuDoubleComplex alpha,
	__device cuDoubleComplex* x);

/** 
 * \fn double spgpuZnrm2(spgpuHandle_t handle, int n, __device cuDoubleComplex* x)
 * Computes the double precision complex Euclidean vector norm of x. 
 * \param handle the spgpu handle used to call this routine
 * \param n the vector's length
 * \param x the input vector
 * \return the euclidean vector norm
 */
 double spgpuZnrm2(spgpuHandle_t handle,
	int n,
	__device cuDoubleComplex* x);

/** 
 * \fn void spgpuZmnrm2(spgpuHandle_t handle, double *y, int n, __device cuDoubleComplex *x, int count, int pitch)
 * Computes the single precision complex Euclidean vector norm for every vector in the multivector x. 
 * \param handle the spgpu handle used to call this routine
 * \param y the array of results
 * \param n the vectors' length in the x multivector
 * \param x the input multivector
 */	
void spgpuZmsnrm2(spgpuHandle_t handle,
	double *y,
	int n,
	__device cuDoubleComplex *x);

/** 
 * \fn void spgpuZmnrm2(spgpuHandle_t handle, CFI_cdesc_t *y, int n, __device cuDoubleComplex *x, int count, int pitch)
 * Computes the single precision complex Euclidean vector norm for every vector in the multivector x. 
 * \param handle the spgpu handle used to call this routine
 * \param y the CFI_descrictor of the result fortran assumed-shape array
 * \param n the vectors' length in the x multivector
 * \param x the input multivector
 * \param count the number of vectors in x
 * \param pitch the multivector's pitch
 */	
void spgpuZmvnrm2(spgpuHandle_t handle,
	CFI_cdesc_t *y,
	int n,
	__device cuDoubleComplex *x,
	int count,
	int pitch);

/** 
 * \fn void spgpuZaxpby(spgpuHandle_t handle, __device cuDoubleComplex *z, int n, cuDoubleComplex beta, __device cuDoubleComplex *y, cuDoubleComplex alpha, __device cuDoubleComplex* x)
 * Computes the double precision complex z = beta * y + alpha * x. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param y the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 */
void spgpuZaxpby(spgpuHandle_t handle,
	__device cuDoubleComplex *z,
	int n,
	cuDoubleComplex beta,
	__device cuDoubleComplex *y,
	cuDoubleComplex alpha,
	__device cuDoubleComplex* x);

/** 
 * \fn void spgpuZcolspan(spgpuHandle_t handle, int n, __device cuDoubleComplex* x, int countX, int pitchX, CFI_cdesc_t* coeff, __device cuDoubleComplex* y, int countY, int pitchY, bool updFlag)
 * Computes the single precision complex y = x * coeff or y = y + x * coeff depending on the updFlag value. y should be a vector if coeff as rank 1 and must be a multivector if coeff is rank 2.
 * \param handle the spgpu handle used to call this routine
 * \param n the vectors' length
 * \param x the input multivector
 * \param countX the number of vectors in x
 * \param pitchX the x multivector pitch
 * \param coeff the CFI descrictor of the fotran array of coefficients
 * \param y the output (multi)vector
 * \param countY the number of vectors in y
 * \param pitchY the y multivector pitch
 * \param updFlag boolen flag: if true, the result will be added to y, otherwise it will overwrite y
 */
void spgpuZcolspan(spgpuHandle_t handle,
	int n,
	__device cuDoubleComplex* x,
	int countX,
	int pitchX,
	CFI_cdesc_t* coeff,
	__device cuDoubleComplex* y,
	int countY,
	int pitchY,
	bool updFlag);

void spgpuZupd_xyz(spgpuHandle_t handle,
	int n,
	cuDoubleComplex alpha,
	cuDoubleComplex beta,
	cuDoubleComplex gamma,
	cuDoubleComplex delta,
	__device cuDoubleComplex* x,
	__device cuDoubleComplex *y,
	__device cuDoubleComplex *z);

void spgpuZxyzw(spgpuHandle_t handle,
	int n,
	cuDoubleComplex a, cuDoubleComplex b,
	cuDoubleComplex c, cuDoubleComplex d,
	cuDoubleComplex e, cuDoubleComplex f,
	__device cuDoubleComplex* x,
	__device cuDoubleComplex *y,
	__device cuDoubleComplex *z,
	__device cuDoubleComplex *w);

/** 
 * \fn void spgpuZmaxpby(spgpuHandle_t handle, __device cuDoubleComplex *z, int n, cuDoubleComplex beta, __device cuDoubleComplex *y, cuDoubleComplex alpha, __device cuDoubleComplex* x, int count, int pitch)
 * Computes the double precision complex z = beta * y + alpha * x of x and y multivectors. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param y the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param count the number of vectors in z, x and y multivectors
 * \param pitch the multivector's pitch
 */
void spgpuZmaxpby(spgpuHandle_t handle,
	__device cuDoubleComplex *z,
	int n,
	cuDoubleComplex beta,
	__device cuDoubleComplex *y,
	cuDoubleComplex alpha,
	__device cuDoubleComplex* x,
	int count, int pitch);

/** 
 * \fn void spgpuZaxy(spgpuHandle_t handle, __device cuDoubleComplex *z, int n, cuDoubleComplex alpha, __device cuDoubleComplex *x, __device cuDoubleComplex* y)
 * Computes the double precision complex z = alpha * x * y. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting vector
 * \param n the vectors' length in the x multivector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the first input vector
 */
void spgpuZaxy(spgpuHandle_t handle,
	__device cuDoubleComplex *z,
	int n,
	cuDoubleComplex alpha,
	__device cuDoubleComplex *x,
	__device cuDoubleComplex *y);

/** 
 * \fn void spgpuZaxypbz(spgpuHandle_t handle, __device cuDoubleComplex *w, int n, cuDoubleComplex beta, __device cuDoubleComplex *z, cuDoubleComplex alpha, __device cuDoubleComplex* x, __device cuDoubleComplex *y)
 * Computes the double precision complex w = beta * z + alpha * x * y. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param z the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the third input vector
 */void spgpuZaxypbz(spgpuHandle_t handle,
	__device cuDoubleComplex *w,
	int n,
	cuDoubleComplex beta,
	__device cuDoubleComplex *z,
	cuDoubleComplex alpha,
	__device cuDoubleComplex* x,
	__device cuDoubleComplex *y);

/** 
 * \fn void spgpuZmaxy(spgpuHandle_t handle, __device cuDoubleComplex *z, int n, cuDoubleComplex alpha, __device cuDoubleComplex *x, __device cuDoubleComplex* y, int count, int pitch)
 * Computes the double precision complex z = alpha * x * y for z, x and y multivectors. z could be exactly x or y (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param z the resulting multivector
 * \param n the vectors' length in the multivectors
 * \param alpha the alpha value
 * \param x the first input multivector
 * \param y the second input multivector
 * \param count the number of vectors in z, x and y multivectors
 * \param pitch the multivectors pitch
 */	
void spgpuZmaxy(spgpuHandle_t handle,
	__device cuDoubleComplex *z,
	int n,
	cuDoubleComplex alpha,
	__device cuDoubleComplex* x,
	__device cuDoubleComplex *y,
	int count,
	int pitch);
	
/** 
 * \fn void spgpuZmaxypbz(spgpuHandle_t handle, __device cuDoubleComplex *w, int n, cuDoubleComplex beta, __device cuDoubleComplex *z, cuDoubleComplex alpha, __device cuDoubleComplex* x, __device cuDoubleComplex *y, int count, int pitch)
 * Computes the double precision complex w = beta * z + alpha * x * y. w could be exactly x, y or z (without offset) or another vector.
 * \param handle the spgpu handle used to call this routine
 * \param w the resulting vector
 * \param n the vectors' length
 * \param beta the beta value
 * \param z the first input vector
 * \param alpha the alpha value
 * \param x the second input vector
 * \param y the third input vector
 * \param count the number of vectors in w, z, x and y multivectors
 * \param pitch the multivectors' pitch
 */
void spgpuZmaxypbz(spgpuHandle_t handle,
	__device cuDoubleComplex *w,
	int n,
	cuDoubleComplex beta,
	__device cuDoubleComplex *z,
	cuDoubleComplex alpha,
	__device cuDoubleComplex* x,
	__device cuDoubleComplex *y,
	int count,
	int pitch);

/** 
 * \fn void spgpuZgath(spgpuHandle_t handle, __device cuDoubleComplex *xValues, int xNnz, const __device int *xIndices, int xBaseIndex, const __device cuDoubleComplex* y)
 * Computes the double precision complex gather from y to xValues (using xIndices).
 * \param handle the spgpu handle used to call this routine
 * \param xValues the destination array for gathered values
 * \param xNnz the number of elements to gather
 * \param xIndices the array of indices for the elements to be gathered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param y the source vector (from which the elements will be gathered)
 */
void spgpuZgath(spgpuHandle_t handle,
	__device cuDoubleComplex *xValues,
	int xNnz,
	const __device int *xIndices,
	int xBaseIndex,
	const __device cuDoubleComplex* y);
	
/** 
 * \fn void spgpuZscat(spgpuHandle_t handle, __device cuDoubleComplex* y, int xNnz, const __device cuDoubleComplex *xValues, const __device int *xIndices, int xBaseIndex, cuDoubleComplex beta)
 * Computes the double precision complex scatter from xValues to y (using xIndices).
 * The scattered element will be, for i in [0, xNnz), y[xIndices[i]] = beta*y[xIndices[i]] + xValues[i] (to be noted that
 * y values will be multiplied with beta just for scattered values).
 * \param handle the spgpu handle used to call this routine
 * \param y the destination vector (to which the elements will be scattered)
 * \param xNnz the number of elements to scatter
 * \param xValues the source array from which the values will be read
 * \param xIndices the array of indices for the elements to be scattered
 * \param xBaseIndex the base index used in xIndices (i.e. 0 for C, 1 for Fortran).
 * \param beta the beta value
 */
void spgpuZscat(spgpuHandle_t handle,
		__device cuDoubleComplex* y,
	int xNnz,
	const __device cuDoubleComplex *xValues,
	const __device int *xIndices,
	int xBaseIndex, cuDoubleComplex beta);	
	
double spgpuZasum(spgpuHandle_t handle,
	int n,
	cuDoubleComplex* x);

double spgpuZamax(spgpuHandle_t handle,
	int n,
	cuDoubleComplex* x);

void spgpuZmasum(spgpuHandle_t handle,
	double* y,
	int n,
	cuDoubleComplex* x,
	int count,
	int pitch);

void spgpuZmamax(spgpuHandle_t handle,
	double* y,
	int n,
	cuDoubleComplex* x,
	int count,
	int pitch);	
	
/** 
 * @}
 */

void spgpuIsetscal(spgpuHandle_t handle,
	int first,
	int last,
	int baseIndex,
	int val,
	__device int *y);

void spgpuSsetscal(spgpuHandle_t handle,
	int first,
	int last,
	int baseIndex,
	float val,
	__device float *y);

void spgpuDsetscal(spgpuHandle_t handle,
	int first,
	int last,
	int baseIndex,
	double val,
	__device double *y);

void spgpuCsetscal(spgpuHandle_t handle,
	int first,
	int last,
	int baseIndex,
	cuFloatComplex val,
	__device cuFloatComplex* y);

void spgpuZsetscal(spgpuHandle_t handle,
	int first,
	int last,
	int baseIndex,
	cuDoubleComplex val,
	__device cuDoubleComplex* y);
  
#ifdef __cplusplus
}
#endif