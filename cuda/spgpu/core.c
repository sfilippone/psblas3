/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
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

	return (err == cudaSuccess) ? SPGPU_SUCCESS : SPGPU_UNSPECIFIED;
}

void spgpuDestroy(spgpuHandle_t pHandle)
{
	cudaStreamDestroy(pHandle->defaultStream);
	free((void*) pHandle);
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
	SpgpuHandleStruct* handle = (SpgpuHandleStruct*) pHandle;
	handle->currentStream = (stream) ? stream : pHandle->defaultStream;
}

cudaStream_t spgpuGetStream(spgpuHandle_t pHandle)
{
	SpgpuHandleStruct* handle = (SpgpuHandleStruct*) pHandle;
	return handle->currentStream;
}

size_t spgpuSizeOf(spgpuType_t typeCode)
{
	switch (typeCode)
	{
		case SPGPU_TYPE_INT:			return sizeof(int);
		case SPGPU_TYPE_FLOAT:			return sizeof(float);
		case SPGPU_TYPE_DOUBLE:			return sizeof(double);
		case SPGPU_TYPE_COMPLEX_FLOAT:	return sizeof(cuFloatComplex);
		case SPGPU_TYPE_COMPLEX_DOUBLE:	return sizeof(cuDoubleComplex);
		default:						return 0; // error		
	}
}

int resetCFIbuffer(CFIBuffer* buffer)
{
	if(buffer->devBuffer) cudaFree(buffer->devBuffer);
	buffer->devBuffer = NULL;
	buffer->devSize = 0;

	if(buffer->hostBuffer) free(buffer->hostBuffer);
	buffer->hostBuffer = NULL;
	buffer->hostSize = 0;

	return SPGPU_SUCCESS;
}

CFIBuffer dotBuffer, coeffBuffer;

int resetCFIstaticBuffers()
{
	int dotRes = resetCFIbuffer(&(dotBuffer));
	int colspanRes = resetCFIbuffer(&(coeffBuffer));
	return (dotRes == SPGPU_SUCCESS) && (colspanRes == SPGPU_SUCCESS) ? SPGPU_SUCCESS : SPGPU_UNSPECIFIED;
}

void CFIdevBufferAlloc(CFIBuffer* buffer, size_t byteCount)
{
	if(buffer->devSize >= byteCount)
		return;

	cudaFree(buffer->devBuffer);
	cudaMalloc(&(buffer->devBuffer), byteCount);
	buffer->devSize = byteCount;
}

void CFIhostBufferAlloc(CFIBuffer* buffer, size_t byteCount)
{
	if(buffer->hostSize >= byteCount)
		return;

	if(buffer->hostBuffer) free(buffer->hostBuffer);
	buffer->hostBuffer = malloc(byteCount);
	buffer->hostSize = byteCount;
}

void CFIBufferLoad(spgpuHandle_t handle, const CFI_cdesc_t* cdesc, CFIBuffer* buffer)
{
	size_t es = cdesc->elem_len;
    size_t m  = (size_t)cdesc->dim[0].extent;
    size_t n  = (cdesc->rank == 2) ? (size_t)cdesc->dim[1].extent : 1;
    size_t bytes = m * n * es;

	CFIdevBufferAlloc(buffer, bytes);

	if(CFI_is_contiguous(cdesc))
	{
		cudaMemcpyAsync(buffer->devBuffer, cdesc->base_addr, bytes, cudaMemcpyHostToDevice, handle->currentStream);
	}
	else
	{
		CFIhostBufferAlloc(buffer, bytes);
		switch(cdesc->rank)
		{
			case 1:
				for(size_t i = 0; i < m; ++i)
					memcpy(((CFIbyte_t*) buffer->hostBuffer) + i*es, &CFI_AT1(CFIbyte_t, cdesc, i), es);
				break;
			case 2:
				for(size_t i = 0; i < m; ++i)
					for(size_t j = 0; j < n; ++j)
						memcpy(((CFIbyte_t*) buffer->hostBuffer) + (i + j*m)*es, &CFI_AT2(unsigned char, cdesc, i, j), es);
				break;
			default:
				printf("Error in CFIBufferLoad: unsupported CFI rank\n");
				return;
		}
		cudaMemcpyAsync(buffer->devBuffer, buffer->hostBuffer, bytes, cudaMemcpyHostToDevice, handle->currentStream);
	}
}

void CFIBufferStore(spgpuHandle_t handle, CFI_cdesc_t* cdesc, CFIBuffer* buffer, size_t bytes)
{
	size_t es = cdesc->elem_len;
    size_t m  = (size_t)cdesc->dim[0].extent;
    size_t n  = (cdesc->rank == 2) ? (size_t)cdesc->dim[1].extent : 1;
    if(bytes != m * n * es) 
		return;
	
	if(CFI_is_contiguous(cdesc))
	{
		cudaMemcpyAsync(cdesc->base_addr, buffer->devBuffer, bytes, cudaMemcpyDeviceToHost, handle->currentStream);
		cudaStreamSynchronize(handle->currentStream);
	}
	else
	{
		CFIhostBufferAlloc(buffer, bytes);
		cudaMemcpyAsync(buffer->hostBuffer, buffer->devBuffer, bytes, cudaMemcpyDeviceToHost, handle->currentStream);
		cudaStreamSynchronize(handle->currentStream);
		switch(cdesc->rank)
		{
			case 1:
				for(size_t i = 0; i < m; ++i)
					memcpy(&CFI_AT1(CFIbyte_t, cdesc, i), ((CFIbyte_t*) buffer->hostBuffer) + i*es, es);
				break;
			case 2:
				for(size_t i = 0; i < m; ++i)
					for(size_t j = 0; j < n; ++j)
						memcpy(&CFI_AT2(CFIbyte_t, cdesc, i, j), ((CFIbyte_t*) buffer->hostBuffer) + (i + j*m)*es, es);
				break;
			default:
				printf("Error in CFIBufferStore: unsupported CFI rank\n");
				return;
		}
	}
}