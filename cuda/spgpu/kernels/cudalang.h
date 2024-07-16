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
 
// Used to avoid mad.f32 instructions on c.c. 1.*
#if __CUDA_ARCH__ >= 200
#define PREC_FADD(a,b) ((a) + (b))
#define PREC_FMUL(a,b) ((a) * (b))
#else
#define PREC_FADD(a,b) __fadd_rn((a),(b))
#define PREC_FMUL(a,b) __fmul_rn((a),(b))
#endif

#define PREC_DADD(a,b) ((a) + (b))
#define PREC_DMUL(a,b) ((a) * (b))


inline __host__ __device__ double2 make_double2(double s)
{
	return make_double2(s, s);
}

inline __host__ __device__ double2 operator+(double2 a, double2 b)
{
	return make_double2(a.x + b.x, a.y + b.y);
}

inline __host__ __device__ void operator+=(double2 &a, double2 b)
{
	a.x += b.x; a.y += b.y;
}

inline __host__ __device__ double2 operator-(double2 a, double2 b)
{
	return make_double2(a.x - b.x, a.y - b.y);
}

inline __host__ __device__ void operator-=(double2 &a, double2 b)
{
	a.x -= b.x; a.y -= b.y;
}

inline __host__ __device__ double2 operator*(double2 a, double s)
{
	return make_double2(a.x * s, a.y * s);
}

inline __host__ __device__ double2 operator*(double s, double2 a)
{
	return make_double2(a.x * s, a.y * s);
}

inline __host__ __device__ void operator*=(double2 &a, double s)
{
	a.x *= s; a.y *= s;
}
