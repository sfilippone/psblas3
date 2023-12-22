#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2015
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


/*! \mainpage The spGPU library documentation
 *
 * \section intro_sec Introduction
 *
 * spGPU is a set of custom matrix storages and CUDA kernels for sparse linear algebra computing on GPU. It isn't a replacement for cuBLAS/cuSPARSE that should be used for a full featured linear algebra environment on GPU.\n
 * The main matrix storage used by spGPU is a GPU-friendly ELLpack format, as well as our HELL (Hacked ELLpack) format, an enhanced version of ELLpack with some interesting memory saving properties.\n
 * HELL format provides a better memory storage compared to ELL (it avoids allocation inefficency provided by spikes in row sizes), while providing quite the same performances for sparse matrix-vector multiply routine..
 *
 * \section install_sec How to build spgpu
 * \subsection linuxbuild Linux (and other unix systems)
 * cd spgpu/build/cmake\n
 * sh configure.sh\n
 * make
 * \section cr_sec Copyright
 * Copyright (C) 2010 - 2015\n
 *     Davide Barbieri - University of Rome Tor Vergata\n
 *     Valeria Cardellini - University of Rome Tor Vergata\n
 *     Salvatore Filippone - University of Rome Tor Vergata
 *
 * This program is free software; you can redistribute it and/or\n
 * modify it under the terms of the GNU General Public License\n
 * version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,\n
 * but WITHOUT ANY WARRANTY; without even the implied warranty of\n
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n
 * GNU General Public License for more details.
 */

 
#include "driver_types.h"
#include "cuComplex.h"
 
/** \addtogroup coreFun Core Routines
 *  @{
 */
 
#ifdef __cplusplus
extern "C" {
#endif

/// __host pointers reference host allocations (it's just a placeholder)
#define __host
/// __device pointers reference device allocations (it's just a placeholder)
#define __device

/// The return code for synchronous functions
typedef int spgpuStatus_t;

#define SPGPU_SUCCESS 		0
#define SPGPU_UNSUPPORTED 	1
#define SPGPU_UNSPECIFIED	2
#define SPGPU_OUTOFMEMORY	3

/// Code to identify a primitive type
typedef int spgpuType_t;

#define SPGPU_TYPE_INT			0
#define SPGPU_TYPE_FLOAT		1
#define SPGPU_TYPE_DOUBLE		2
#define SPGPU_TYPE_COMPLEX_FLOAT	3
#define SPGPU_TYPE_COMPLEX_DOUBLE	4

/// this struct should be modified only internally by spgpu
typedef struct spgpuHandleStruct {
	/// the current stream used by every calls on spgpu routines (passing this handle)
	cudaStream_t currentStream;
	/// the default stream, created during the handle creation.
	cudaStream_t defaultStream;
	/// the device associated to this handle
	int device;
	/// the warp size for this device
	int warpSize;
	/// the max threads per block count for this device
	int maxThreadsPerBlock;
	/// the max size for the X coordinate of the grid dimensions
	int maxGridSizeX;
	/// the max size for the Y coordinate of the grid dimensions
	int maxGridSizeY;
	/// the max size for the Z coordinate of the grid dimensions
	int maxGridSizeZ;
        /// Number of SM
        int multiProcessorCount;
	// compute capability
	int capabilityMajor;
	int capabilityMinor;
} SpgpuHandleStruct;

/// A spGPU handle represents a single CUDA device on your platform.
typedef const SpgpuHandleStruct* spgpuHandle_t;

/**
* \fn spgpuStatus_t spgpuCreate(spgpuHandle_t* pHandle, int device)
* Create a spgpu context for a GPU device. Every call to spgpu routines using this
* handle will execute on the same GPU. This is re-entrant, so it will work if used by multiple host threads.
* \param pHandle outputs the handle
* \param device id of the device to be used by this context
*/
spgpuStatus_t spgpuCreate(spgpuHandle_t* pHandle, int device);

/**
* \fn void spgpuDestroy(spgpuHandle_t pHandle)
* Destroy the spgpu context for pHandle.
* \param pHandle the handle previously created with spgpuCreate().
*/
void spgpuDestroy(spgpuHandle_t pHandle);

/**
* \fn void spgpuStreamCreate(spgpuHandle_t pHandle, cudaStream_t* stream)
* Create a cuda stream according to the device of the spgpu handle.
* \param stream outputs the new stream
* \param pHandle the handle used to obtain the device id for the stream
*/
void spgpuStreamCreate(spgpuHandle_t pHandle, cudaStream_t* stream);

/**
* \fn void spgpuStreamDestroy(cudaStream_t stream)
* Destroy a stream, previously created with spgpuStreamCreate().
* \param stream the stream to destroy
*/
void spgpuStreamDestroy(cudaStream_t stream);

/**
* \fn void spgpuSetStream(spgpuHandle_t pHandle, cudaStream_t stream)
* Change the current stream for the handle pHandle.
* \param pHandle the handle to configure.
* \param stream the stream to use for next spgpu routines call. Use 0 to reset to the default stream.
*/
void spgpuSetStream(spgpuHandle_t pHandle, cudaStream_t stream);

/**
* \fn cudaStram_t spgpuGetStream(spgpuHandle_t pHandle)
* Get the current stream from the handle pHandle.
* \param pHandle the handle from which get the stream.
*/
cudaStream_t spgpuGetStream(spgpuHandle_t pHandle);

/**
* \fn size_t spgpuSizeOf(spgpuType_t typeCode)
* Returns the size, in bytes, of the type denoted by typeCode (e.g. 4 for SPGPU_TYPE_FLOAT, 8 for SPGPU_TYPE_DOUBLE).
* \param typeCode outputs the handle
*/
size_t spgpuSizeOf(spgpuType_t typeCode);

/*
typedef struct {
spgpuMatrix

spgpuMatrixType_t MatrixType;
spgpuFillMode_t FillMode;
spgpuDiagType_t DiagType;
int baseIndex;
} spgpuMatrixDesc_t
*/

#define cuFloatComplex_isZero(a) (a.x == 0.0f && a.y == 0.0f)
#define cuDoubleComplex_isZero(a) (a.x == 0.0 && a.y == 0.0)
#define cuFloatComplex_isNotZero(a) (a.x != 0.0f || a.y != 0.0f)
#define cuDoubleComplex_isNotZero(a) (a.x != 0.0 || a.y != 0.0)

#ifdef __cplusplus
}
#endif

/** @}*/

