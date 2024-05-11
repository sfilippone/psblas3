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

#include "cudadebug.h"
#include "cudalang.h"
#include "cuComplex.h"

extern "C"
{
#include "core.h"
#include "ell.h"
  int getGPUSharedMemPerBlock();
  int getGPUMultiProcessors();
  int getGPUMaxThreadsPerMP();
}

#include "debug.h"

#define VALUE_TYPE cuDoubleComplex
#define TYPE_SYMBOL Z
#define TEX_FETCH_TYPE int4
#include "ell_spmm_base.cuh"
