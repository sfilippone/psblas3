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

extern "C"
{
#include "core.h"
#include "hell.h"
  int getGPUSharedMemPerBlock();
  int getGPUMultiProcessors();
  int getGPUMaxThreadsPerMP();
}

#include "debug.h"

#define VALUE_TYPE float
#define TYPE_SYMBOL S
#define TEX_FETCH_TYPE float
#include "hell_spmm_base.cuh"
