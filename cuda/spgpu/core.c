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
 
#include "core.h"
#include "stdlib.h"
#include "cuda_runtime.h"

spgpuStatus_t spgpuCreate(spgpuHandle_t* pHandle, int device)
{
	struct cudaDeviceProp deviceProperties;
	cudaError_t err = cudaGetDeviceProperties(&deviceProperties, device);

	SpgpuHandleStruct* handle = (SpgpuHandleStruct*) malloc(sizeof(SpgpuHandleStruct));

	int currentDevice;
	cudaGetDevice(&currentDevice);
	cudaSetDevice(device);
	cudaStreamCreate(&handle->defaultStream);
	handle->currentStream = handle->defaultStream;
	cudaSetDevice(currentDevice);

	handle->device = device;
	handle->warpSize = deviceProperties.warpSize;
	handle->maxThreadsPerBlock = deviceProperties.maxThreadsPerBlock;
	handle->multiProcessorCount = deviceProperties.multiProcessorCount;
	handle->maxGridSizeX = deviceProperties.maxGridSize[0];
	handle->maxGridSizeY = deviceProperties.maxGridSize[1];
	handle->maxGridSizeZ = deviceProperties.maxGridSize[2];
	handle->capabilityMajor = deviceProperties.major;
	handle->capabilityMinor = deviceProperties.minor;
	
	*pHandle = handle;

	if (err == cudaSuccess)
		return SPGPU_SUCCESS;
	else
		return SPGPU_UNSPECIFIED;
}

void spgpuDestroy(spgpuHandle_t pHandle)
{
	cudaStreamDestroy(pHandle->defaultStream);

	free((void*)pHandle);
}

void spgpuStreamCreate(spgpuHandle_t pHandle, cudaStream_t* stream)
{
	int currentDevice;
	cudaGetDevice(&currentDevice);
	cudaSetDevice(pHandle->device);
	cudaStreamCreate(stream);
	cudaSetDevice(currentDevice);
}

void spgpuStreamDestroy(cudaStream_t stream)
{
	cudaStreamDestroy(stream);
}

void spgpuSetStream(spgpuHandle_t pHandle, cudaStream_t stream)
{
	SpgpuHandleStruct* handle = (SpgpuHandleStruct*)pHandle;

	if (stream)
	{
		handle->currentStream = stream;
	}
	else
		handle->currentStream = pHandle->defaultStream;
}

cudaStream_t spgpuGetStream(spgpuHandle_t pHandle)
{
	SpgpuHandleStruct* handle = (SpgpuHandleStruct*)pHandle;
	return handle->currentStream;
}

size_t spgpuSizeOf(spgpuType_t typeCode)
{
	switch (typeCode)
	{
	case SPGPU_TYPE_INT:
		return sizeof(int);
	case SPGPU_TYPE_FLOAT:
		return sizeof(float);
	case SPGPU_TYPE_DOUBLE:
		return sizeof(double);
	case SPGPU_TYPE_COMPLEX_FLOAT:
		return sizeof(cuFloatComplex);
	case SPGPU_TYPE_COMPLEX_DOUBLE:
		return sizeof(cuDoubleComplex);
	default:
		return 0; // error		
	}
}
