  /*             Parallel Sparse BLAS   GPU plugin  */
  /*   (C) Copyright 2013 */

  /*                      Salvatore Filippone */
  /*                      Alessandro Fanfarillo */
 
  /* Redistribution and use in source and binary forms, with or without */
  /* modification, are permitted provided that the following conditions */
  /* are met: */
  /*   1. Redistributions of source code must retain the above copyright */
  /*      notice, this list of conditions and the following disclaimer. */
  /*   2. Redistributions in binary form must reproduce the above copyright */
  /*      notice, this list of conditions, and the following disclaimer in the */
  /*      documentation and/or other materials provided with the distribution. */
  /*   3. The name of the PSBLAS group or the names of its contributors may */
  /*      not be used to endorse or promote products derived from this */
  /*      software without specific written permission. */
 
  /* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
  /* ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED */
  /* TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR */
  /* PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS */
  /* BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR */
  /* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF */
  /* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS */
  /* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN */
  /* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) */
  /* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
  /* POSSIBILITY OF SUCH DAMAGE. */
 
  

#ifndef _CUDA_UTIL_H_
#define _CUDA_UTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

#include "cuda_runtime.h"
#include "core.h" 
#include "cuComplex.h"
#include "fcusparse.h"
#include "cublas_v2.h"

int allocRemoteBuffer(void** buffer, int count);
int allocMappedMemory(void **buffer, void **dp, int size);
int registerMappedMemory(void *buffer, void **dp, int size);
int unregisterMappedMemory(void *buffer);
int writeRemoteBuffer(void* hostSrc, void* buffer, int count);
int writeRemoteBufferR2(void* hostSrc, int hpitch, void* buffer, int count, int pitch, int size);
int readRemoteBuffer(void* hostDest, void* buffer, int count);
int readRemoteBufferR2(void* hostDest, int hpitch, void* buffer, int count, int pitch, int size);
int freeRemoteBuffer(void* buffer);
int gpuInit(int dev);
int getDeviceCount();
int getDevice();
int getDeviceHasUVA();
int setDevice(int dev);
int getGPUMultiProcessors();
int getGPUMemoryBusWidth();
int getGPUMemoryClockRate();
int getGPUWarpSize();
int getGPUMaxThreadsPerBlock();
int getGPUMaxThreadsPerMP();
int getGPUMaxRegistersPerBlock();
void cpyGPUNameString(char *cstring);


void cudaSync();
void cudaReset();
void gpuClose();


spgpuHandle_t psb_cudaGetHandle(); 
void psb_cudaCreateHandle();
void psb_cudaDestroyHandle();
cudaStream_t psb_cudaGetStream();
void  psb_cudaSetStream(cudaStream_t stream);

cublasHandle_t psb_cudaGetCublasHandle(); 
void psb_cudaCreateCublasHandle();
void psb_cudaDestroyCublasHandle();


int allocateInt(void **, int);
int allocateMultiInt(void **, int, int);
int writeInt(void *, int *, int);
int writeMultiInt(void *, int* , int , int );
int readInt(void *, int *, int);
int readMultiInt(void*, int*, int, int );
int writeIntFirst(int,void *, int *, int,int);
int readIntFirst(int,void *, int *, int,int);
void freeInt(void *);

int allocateFloat(void **, int);
int allocateMultiFloat(void **, int, int);
int writeFloat(void *, float *, int);
int writeMultiFloat(void *, float* , int , int );
int readFloat(void *, float*, int);
int readMultiFloat(void*, float*, int, int );
int writeFloatFirst(int, void *, float*, int, int);
int readFloatFirst(int, void *, float*, int, int);
void freeFloat(void *);

int allocateDouble(void **, int);
int allocateMultiDouble(void **, int, int);
int writeDouble(void *, double*, int);
int writeMultiDouble(void *, double* , int , int );
int readDouble(void *, double*, int);
int readMultiDouble(void*, double*, int, int );
int writeDoubleFirst(int, void *, double*, int, int);
int readDoubleFirst(int, void *, double*, int, int);
void freeDouble(void *);

int allocateFloatComplex(void **, int);
int allocateMultiFloatComplex(void **, int, int);
int writeFloatComplex(void *, cuFloatComplex*, int);
int writeMultiFloatComplex(void *, cuFloatComplex* , int , int );
int readFloatComplex(void *, cuFloatComplex*, int);
int readMultiFloatComplex(void*, cuFloatComplex*, int, int );
int writeFloatComplexFirst(int, void *, cuFloatComplex*, int, int);
int readFloatComplexFirst(int, void *, cuFloatComplex*, int, int);
void freeFloatComplex(void *);

int allocateDoubleComplex(void **, int);
int allocateMultiDoubleComplex(void **, int, int);
int writeDoubleComplex(void *, cuDoubleComplex*, int);
int writeMultiDoubleComplex(void *, cuDoubleComplex* , int , int );
int readDoubleComplex(void *, cuDoubleComplex*, int);
int readMultiDoubleComplex(void*, cuDoubleComplex*, int, int );
int writeDoubleComplexFirst(int, void *, cuDoubleComplex*, int, int);
int readDoubleComplexFirst(int, void *, cuDoubleComplex*, int, int);
void freeDoubleComplex(void *);

double etime();


#endif
