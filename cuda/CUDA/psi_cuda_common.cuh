#pragma once

#define PRE_CONCAT(A, B) A ## B
#define CONCAT(A, B) PRE_CONCAT(A, B)
#define MIN(A,B) ( (A)<(B) ? (A) : (B) )
#define SQUARE(x) ((x)*(x))
#define GET_ADDR(a,ix,iy,nc) a[(nc)*(ix)+(iy)]
#define GET_VAL(a,ix,iy,nc) (GET_ADDR(a,ix,iy,nc))

__device__ __host__ static float zero_float() { return 0.0f; }
__device__ __host__ static cuFloatComplex zero_cuFloatComplex() { return make_cuFloatComplex(0.0, 0.0); }

#if (__CUDA_ARCH__ >= 130) || (!__CUDA_ARCH__)
__device__ __host__ static double zero_double() { return 0.0; }
__device__ __host__ static cuDoubleComplex zero_cuDoubleComplex() { return make_cuDoubleComplex(0.0, 0.0); }
#endif
