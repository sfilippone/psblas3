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
 
  

#include "cuda_util.h"


static int hasUVA=-1;
static struct cudaDeviceProp *prop=NULL;
static spgpuHandle_t psb_cuda_handle = NULL;
static cublasHandle_t psb_cublas_handle = NULL;
#if defined(TRACK_CUDA_MALLOC)
static long long total_cuda_mem = 0;
#endif

int allocRemoteBuffer(void** buffer, int count)
{
  cudaError_t err = cudaMalloc(buffer, count);
#if defined(TRACK_CUDA_MALLOC)
  total_cuda_mem += count;
  fprintf(stderr,"Tracking CUDA allocRemoteBuffer for %ld bytes total  %ld  address %p\n",
	  count, total_cuda_mem, *buffer);
#endif
    
  if (err == cudaSuccess)
    {
      return SPGPU_SUCCESS;
    }
  else
    { 
      fprintf(stderr,"CUDA allocRemoteBuffer for %d bytes Error: %s \n",
	      count, cudaGetErrorString(err));
      if(err == cudaErrorMemoryAllocation)
	return SPGPU_OUTOFMEMORY;
      else
	return SPGPU_UNSPECIFIED;
    }
}

int hostRegisterMapped(void *pointer, long size) 
{
  cudaError_t err = cudaHostRegister(pointer, size, cudaHostRegisterMapped);

 if (err == cudaSuccess)
    {
      return SPGPU_SUCCESS;
    }
  else
    { 
      fprintf(stderr,"CUDA hostRegisterMapped Error: %s\n", cudaGetErrorString(err));
      if(err == cudaErrorMemoryAllocation)
	return SPGPU_OUTOFMEMORY;
      else
	return SPGPU_UNSPECIFIED;
    }
}

int getDevicePointer(void **d_p, void * h_p)
{
  cudaError_t err = cudaHostGetDevicePointer(d_p,h_p,0);

 if (err == cudaSuccess)
    {
      return SPGPU_SUCCESS;
    }
  else
    { 
      fprintf(stderr,"CUDA getDevicePointer Error: %s\n", cudaGetErrorString(err));
      if(err == cudaErrorMemoryAllocation)
	return SPGPU_OUTOFMEMORY;
      else
	return SPGPU_UNSPECIFIED;
    }
}

int registerMappedMemory(void *buffer, void **dp, int size)
{
  //cudaError_t err = cudaHostAlloc(buffer,size,cudaHostAllocMapped);
  cudaError_t err = cudaHostRegister(buffer, size, cudaHostRegisterMapped);
  if (err == cudaSuccess) err = cudaHostGetDevicePointer(dp,buffer,0);

  if (err == cudaSuccess)
    {
      err = cudaHostGetDevicePointer(dp,buffer,0);
      if (err == cudaSuccess)
	{
	  return SPGPU_SUCCESS;
	}
      else
	{
	  fprintf(stderr,"CUDA registerMappedMemory Error: %s\n", cudaGetErrorString(err));
	  return SPGPU_UNSPECIFIED;
	}
    }
  else
    { 
      fprintf(stderr,"CUDA registerMappedMemory Error: %s\n", cudaGetErrorString(err));
      if(err == cudaErrorMemoryAllocation)
	return SPGPU_OUTOFMEMORY;
      else
	return SPGPU_UNSPECIFIED;
    }
}

int allocMappedMemory(void **buffer, void **dp, int size)
{
  cudaError_t err = cudaHostAlloc(buffer,size,cudaHostAllocMapped);
  if (err == 0) err = cudaHostGetDevicePointer(dp,*buffer,0);

  if (err == cudaSuccess)
    {
      return SPGPU_SUCCESS;
    }
  else
    { 
      fprintf(stderr,"CUDA allocMappedMemory Error: %s\n", cudaGetErrorString(err));
      if(err == cudaErrorMemoryAllocation)
	return SPGPU_OUTOFMEMORY;
      else
	return SPGPU_UNSPECIFIED;
    }
}

int unregisterMappedMemory(void *buffer)
{
  //cudaError_t err = cudaHostAlloc(buffer,size,cudaHostAllocMapped);
  cudaError_t err = cudaHostUnregister(buffer);

  if (err == cudaSuccess)
    {
      return SPGPU_SUCCESS;
    }
  else
    { 
      fprintf(stderr,"CUDA unregisterMappedMemory Error: %s\n", cudaGetErrorString(err));
      if(err == cudaErrorMemoryAllocation)
	return SPGPU_OUTOFMEMORY;
      else
	return SPGPU_UNSPECIFIED;
    }
}

int writeRemoteBuffer(void* hostSrc, void* buffer, int count)
{
  cudaError_t err = cudaMemcpy(buffer, hostSrc, count, cudaMemcpyHostToDevice);

  if (err == cudaSuccess)
    return SPGPU_SUCCESS;	
  else {
    fprintf(stderr,"CUDA Error writeRemoteBuffer: %s  %p %p %d\n", 
	    cudaGetErrorString(err),buffer, hostSrc, count);
    return SPGPU_UNSPECIFIED;
  }
}

// TODO
int writeRemoteBufferR2(void* hostSrc, void* buffer, int count, int pitch, int size)
{
  cudaError_t err = cudaMemcpy2D(buffer, pitch, hostSrc, size, size, count, cudaMemcpyHostToDevice);

  if (err == cudaSuccess)
    return SPGPU_SUCCESS;	
  else {
    fprintf(stderr,"CUDA Error writeRemoteBuffer: %s  %p %p %d %d %d\n", 
	    cudaGetErrorString(err),buffer, hostSrc, count, pitch, size);
    return SPGPU_UNSPECIFIED;
  }
}

int readRemoteBuffer(void* hostDest, void* buffer, int count)
{

  
  cudaError_t err1;
  cudaError_t err;
#if 0 
  {
    err1      =cudaGetLastError();
    fprintf(stderr,"CUDA Error prior to readRemoteBuffer: %s %d\n", 
	  cudaGetErrorString(err1),err1);
  }
  
#endif
  err = cudaMemcpy(hostDest, buffer, count, cudaMemcpyDeviceToHost);

  if (err == cudaSuccess)
    return SPGPU_SUCCESS;	
  else {
    fprintf(stderr,"CUDA Error readRemoteBuffer: %s %p  %p %d %d\n", 
	    cudaGetErrorString(err),hostDest,buffer,count,err);
    return SPGPU_UNSPECIFIED;
  }
}

// TODO sistemare pitch e size (si possono gestire senza realloc su fortran)
int readRemoteBufferR2(void* hostDest, void* buffer, int count, int pitch, int size)
{
  cudaError_t err = cudaMemcpy2D(hostDest, size, buffer, pitch, size, count, cudaMemcpyDeviceToHost);

  if (err == cudaSuccess)
    return SPGPU_SUCCESS;	
  else {
    fprintf(stderr,"CUDA Error readRemoteBuffer: %s %p  %p %d %d\n", 
	    cudaGetErrorString(err),hostDest,buffer,count,err);
    return SPGPU_UNSPECIFIED;
  }
}

int freeRemoteBuffer(void* buffer)
{
  cudaError_t err = cudaFree(buffer);
  if (err == cudaSuccess)
    return SPGPU_SUCCESS;	
  else {
    fprintf(stderr,"CUDA Error freeRemoteBuffer: %s  %p\n", cudaGetErrorString(err),buffer);
    return SPGPU_UNSPECIFIED;
  }
}

int gpuInit(int dev)
{

  int count,err;  
  
  if ((err=cudaSetDeviceFlags(cudaDeviceMapHost))!=cudaSuccess) 
    fprintf(stderr,"Error On SetDeviceFlags: %d '%s'\n",err,cudaGetErrorString(err));
  if ((prop=(struct cudaDeviceProp *) malloc(sizeof(struct cudaDeviceProp)))==NULL) {
    fprintf(stderr,"CUDA Error gpuInit3: not malloced prop\n");
    return SPGPU_UNSPECIFIED;
  }
  err = setDevice(dev);
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA Error gpuInit2: %s\n", cudaGetErrorString(err));
    return SPGPU_UNSPECIFIED;
  }
  if (!psb_cublas_handle)
    psb_cudaCreateCublasHandle();
  hasUVA=getDeviceHasUVA();
  FcusparseCreate();
  return err;
  
}

void gpuClose()
{
  cudaStream_t st1, st2;
  if (! psb_cuda_handle)
    st1=spgpuGetStream(psb_cuda_handle);
  if (! psb_cublas_handle)
    cublasGetStream(psb_cublas_handle,&st2);
  FcusparseDestroy();
  psb_cudaDestroyHandle();
  if (st1 != st2) 
    psb_cudaDestroyCublasHandle();
  free(prop);
  prop=NULL;
  hasUVA=-1;
}


int setDevice(int dev)
{
  int count,err,idev;
  
  err = cudaGetDeviceCount(&count);
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA Error setDevice: %s\n", cudaGetErrorString(err));
    return SPGPU_UNSPECIFIED;
  }
  
  if ((0<=dev)&&(dev<count))
    idev = dev;
  else
    idev = 0;
  err = cudaSetDevice(idev);
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA Error gpuInit2: %s\n", cudaGetErrorString(err));
    return SPGPU_UNSPECIFIED;
  }
  err = cudaGetDeviceProperties(prop,idev);
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA Error gpuInit4: %s\n", cudaGetErrorString(err));
    return SPGPU_UNSPECIFIED;
  }
  
  return SPGPU_SUCCESS;
}

int getDevice()
{ int count;
  
  cudaGetDevice(&count);     
  return(count);
}

int getDeviceHasUVA()
{ int count=0;
  if (prop!=NULL) 
    count = prop->unifiedAddressing;
  return(count);
}

int getGPUMultiProcessors()
{ int count=0;
  if (prop!=NULL) 
    count = prop->multiProcessorCount;
  return(count);
}


int getGPUMemoryBusWidth()
{ int count=0;
#if CUDART_VERSION >= 5000
  if (prop!=NULL) 
    count = prop->memoryBusWidth;
#endif
  return(count);
}
int getGPUMemoryClockRate()
{ int count=0;
#if CUDART_VERSION >= 5000
  if (prop!=NULL) 
    count = prop->memoryClockRate;
#endif
  return(count);
}
int getGPUWarpSize()
{ int count=0;
  if (prop!=NULL) 
    count = prop->warpSize;
  return(count);
}
int getGPUMaxThreadsPerBlock()
{ int count=0;
  if (prop!=NULL) 
    count = prop->maxThreadsPerBlock;
  return(count);
}
int getGPUMaxThreadsPerMP()
{ int count=0;
  if (prop!=NULL) 
    count = prop->maxThreadsPerMultiProcessor;
  return(count);
}
int getGPUMaxRegistersPerBlock()
{ int count=0;
  if (prop!=NULL) 
    count = prop->regsPerBlock;
  return(count);
}

void cpyGPUNameString(char *cstring)
{
  *cstring='\0';
  if (prop!=NULL) 
    strcpy(cstring,prop->name);

}

int DeviceHasUVA()
{ 
  return(hasUVA == 1);
}
  

int getDeviceCount()
{ int count;
  cudaError_t err;
  err = cudaGetDeviceCount(&count);
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA Error getDeviceCount: %s\n", cudaGetErrorString(err));
    return SPGPU_UNSPECIFIED;
  }
  return(count);
}

void cudaSync()
{
  cudaError_t err;
  err = cudaDeviceSynchronize();
#if 0
  if (err == cudaSuccess)
    return SPGPU_SUCCESS;	
  else {
    fprintf(stderr,"CUDA Error cudaSync: %s\n", cudaGetErrorString(err));
    return SPGPU_UNSPECIFIED;
  }
#else
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA Error cudaSync: %s\n", cudaGetErrorString(err));
  }
  return ;
#endif
}

void cudaReset()
{
  cudaError_t err;
  err = cudaDeviceReset();
#if 0
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA Error Reset: %s\n", cudaGetErrorString(err));
    return SPGPU_UNSPECIFIED;
  }
#else
  if (err != cudaSuccess) {
    fprintf(stderr,"CUDA Error Reset: %s\n", cudaGetErrorString(err));
  }
  return ;
#endif
}


spgpuHandle_t psb_cudaGetHandle()
{
  return psb_cuda_handle;
}

void psb_cudaCreateHandle()
{
  if (!psb_cuda_handle)
    spgpuCreate(&psb_cuda_handle, getDevice());
  
}

void psb_cudaDestroyHandle()
{
  if (!psb_cuda_handle)
    spgpuDestroy(psb_cuda_handle);
  psb_cuda_handle = NULL; 
}

cudaStream_t psb_cudaGetStream()
{
  return spgpuGetStream(psb_cuda_handle);
}

void  psb_cudaSetStream(cudaStream_t stream)
{
  spgpuSetStream(psb_cuda_handle, stream);
  return ;
}



cublasHandle_t psb_cudaGetCublasHandle()
{
  if (!psb_cublas_handle)
    psb_cudaCreateCublasHandle();
  return psb_cublas_handle;
}
void psb_cudaCreateCublasHandle()
{  if (!psb_cublas_handle)
    cublasCreate(&psb_cublas_handle);
}
void psb_cudaDestroyCublasHandle()
{
  if (!psb_cublas_handle)
    cublasDestroy(psb_cublas_handle);
  psb_cublas_handle=NULL;
}

/* Simple memory tools */ 

int allocateInt(void **d_int, int n)
{
  return allocRemoteBuffer((void **)(d_int), n*sizeof(int));
}

int writeInt(void *d_int, int* h_int, int n)
{
  int i,j;
  int *di;
  i = writeRemoteBuffer((void*)h_int, (void*)d_int, n*sizeof(int));
  return i;
}

int readInt(void* d_int, int* h_int, int n)
{ int i;
  i = readRemoteBuffer((void *) h_int, (void *) d_int, n*sizeof(int));
  //cudaSync();
  return(i);
}

int writeIntFirst(int first, void *d_int, int* h_int, int n, int IndexBase)
{
  int i,j;
  int *di=(int *) d_int;
  di = &(di[first-IndexBase]);
  i = writeRemoteBuffer((void*)h_int, (void*)di, n*sizeof(int));
  return i;
}

int readIntFirst(int first,void* d_int, int* h_int, int n, int IndexBase)
{ int i;
  int *di=(int *) d_int;
  di = &(di[first-IndexBase]);
  i = readRemoteBuffer((void *) h_int, (void *) di, n*sizeof(int));
  //cudaSync();
  return(i);
}

int allocateMultiInt(void **d_int, int m, int n)
{
  return allocRemoteBuffer((void **)(d_int), m*n*sizeof(int));
}

int writeMultiInt(void *d_int, int* h_int, int m, int n)
{
  int i,j;
  int *di;
  i = writeRemoteBuffer((void*)h_int, (void*)d_int, m*n*sizeof(int));
  return i;
}

int readMultiInt(void* d_int, int* h_int, int m, int n)
{ int i;
  i = readRemoteBuffer((void *) h_int, (void *) d_int, m*n*sizeof(int));
  //cudaSync();
  return(i);
}

void freeInt(void *d_int)
{
  //printf("Before freeInt\n");
  freeRemoteBuffer(d_int);
}




int allocateFloat(void **d_float, int n)
{
  return allocRemoteBuffer((void **)(d_float), n*sizeof(float));
}

int writeFloat(void *d_float, float* h_float, int n)
{
  int i;

  i = writeRemoteBuffer((void*)h_float, (void*)d_float, n*sizeof(float));

  return i;
}

int readFloat(void* d_float, float* h_float, int n)
{ int i;
  i = readRemoteBuffer((void *) h_float, (void *) d_float, n*sizeof(float));

  return(i);
}

int writeFloatFirst(int df, void *d_float, float* h_float, int n, int IndexBase)
{
  int i;

  float *dv=(float *) d_float;
  dv = &dv[df-IndexBase];
  i = writeRemoteBuffer((void*)h_float, (void*)dv, n*sizeof(float));

  return i;
}

int readFloatFirst(int df, void* d_float, float* h_float, int n, int IndexBase)
{ int i;
  float *dv=(float *) d_float;
  dv = &dv[df-IndexBase];    
  //fprintf(stderr,"readFloatFirst: %d %p %p %p %d \n",df,d_float,dv,h_float,n);
  i = readRemoteBuffer((void *) h_float, (void *) dv, n*sizeof(float));

  return(i);
}


int allocateMultiFloat(void **d_float, int m, int n)
{
  return allocRemoteBuffer((void **)(d_float), m*n*sizeof(float));
}

int writeMultiFloat(void *d_float, float* h_float, int m, int n)
{
  int i,j;
  i = writeRemoteBuffer((void*)h_float, (void*)d_float, m*n*sizeof(float));
  return i;
}

int readMultiFloat(void* d_float, float* h_float, int m, int n)
{ int i;
  i = readRemoteBuffer((void *) h_float, (void *) d_float, m*n*sizeof(float));
  //cudaSync();
  return(i);
}

void freeFloat(void *d_float)
{
  freeRemoteBuffer(d_float);
}



int allocateDouble(void **d_double, int n)
{
  return allocRemoteBuffer((void **)(d_double), n*sizeof(double));
}

int writeDouble(void *d_double, double* h_double, int n)
{
  int i;

  i = writeRemoteBuffer((void*)h_double, (void*)d_double, n*sizeof(double));

  return i;
}

int readDouble(void* d_double, double* h_double, int n)
{ int i;
  i = readRemoteBuffer((void *) h_double, (void *) d_double, n*sizeof(double));

  return(i);
}

int writeDoubleFirst(int df, void *d_double, double* h_double, int n, int IndexBase)
{
  int i;

  double *dv=(double *) d_double;
  dv = &dv[df-IndexBase];
  i = writeRemoteBuffer((void*)h_double, (void*)dv, n*sizeof(double));

  return i;
}

int readDoubleFirst(int df, void* d_double, double* h_double, int n, int IndexBase)
{ int i;
  double *dv=(double *) d_double;
  dv = &dv[df-IndexBase];    
  //fprintf(stderr,"readDoubleFirst: %d %p %p %p %d \n",df,d_double,dv,h_double,n);
  i = readRemoteBuffer((void *) h_double, (void *) dv, n*sizeof(double));

  return(i);
}

int allocateMultiDouble(void **d_double, int m, int n)
{
  return allocRemoteBuffer((void **)(d_double), m*n*sizeof(double));
}

int writeMultiDouble(void *d_double, double* h_double, int m, int n)
{
  int i,j;
  i = writeRemoteBuffer((void*)h_double, (void*)d_double, m*n*sizeof(double));
  return i;
}

int readMultiDouble(void* d_double, double* h_double, int m, int n)
{ int i;
  i = readRemoteBuffer((void *) h_double, (void *) d_double, m*n*sizeof(double));
  //cudaSync();
  return(i);
}

void freeDouble(void *d_double)
{
  freeRemoteBuffer(d_double);
}



int allocateFloatComplex(void **d_FloatComplex, int n)
{
  return allocRemoteBuffer((void **)(d_FloatComplex), n*sizeof(cuFloatComplex));
}

int writeFloatComplex(void *d_FloatComplex, cuFloatComplex* h_FloatComplex, int n)
{
  int i;

  i = writeRemoteBuffer((void*)h_FloatComplex, (void*)d_FloatComplex, n*sizeof(cuFloatComplex));

  return i;
}

int readFloatComplex(void* d_FloatComplex, cuFloatComplex* h_FloatComplex, int n)
{ int i;
  i = readRemoteBuffer((void *) h_FloatComplex, (void *) d_FloatComplex, n*sizeof(cuFloatComplex));

  return(i);
}

int allocateMultiFloatComplex(void **d_FloatComplex, int m, int n)
{
  return allocRemoteBuffer((void **)(d_FloatComplex), m*n*sizeof(cuFloatComplex));
}

int writeMultiFloatComplex(void *d_FloatComplex, cuFloatComplex* h_FloatComplex, int m, int n)
{
  int i,j;
  i = writeRemoteBuffer((void*)h_FloatComplex, (void*)d_FloatComplex, m*n*sizeof(cuFloatComplex));
  return i;
}

int readMultiFloatComplex(void* d_FloatComplex, cuFloatComplex* h_FloatComplex, int m, int n)
{ int i;
  i = readRemoteBuffer((void *) h_FloatComplex, (void *) d_FloatComplex, m*n*sizeof(cuFloatComplex));
  //cudaSync();
  return(i);
}

int writeFloatComplexFirst(int df, void *d_floatComplex,
			   cuFloatComplex* h_floatComplex, int n, int IndexBase)
{
  int i;
  
  cuFloatComplex *dv=(cuFloatComplex *) d_floatComplex;
  dv = &dv[df-IndexBase];
  i = writeRemoteBuffer((void*)h_floatComplex, (void*)dv, n*sizeof(cuFloatComplex));

  return i;
}

int readFloatComplexFirst(int df, void* d_floatComplex, cuFloatComplex* h_floatComplex,
			  int n, int IndexBase)
{ int i;
  cuFloatComplex *dv=(cuFloatComplex *) d_floatComplex;
  dv = &dv[df-IndexBase];    
  i = readRemoteBuffer((void *) h_floatComplex, (void *) dv, n*sizeof(cuFloatComplex));

  return(i);
}

void freeFloatComplex(void *d_FloatComplex)
{
  freeRemoteBuffer(d_FloatComplex);
}




int allocateDoubleComplex(void **d_DoubleComplex, int n)
{
  return allocRemoteBuffer((void **)(d_DoubleComplex), n*sizeof(cuDoubleComplex));
}

int writeDoubleComplex(void *d_DoubleComplex, cuDoubleComplex* h_DoubleComplex, int n)
{
  int i;

  i = writeRemoteBuffer((void*)h_DoubleComplex, (void*)d_DoubleComplex, n*sizeof(cuDoubleComplex));

  return i;
}

int readDoubleComplex(void* d_DoubleComplex, cuDoubleComplex* h_DoubleComplex, int n)
{ int i;
  i = readRemoteBuffer((void *) h_DoubleComplex, (void *) d_DoubleComplex, n*sizeof(cuDoubleComplex));

  return(i);
}

int writeDoubleComplexFirst(int df, void *d_doubleComplex,
			   cuDoubleComplex* h_doubleComplex, int n, int IndexBase)
{
  int i;
  
  cuDoubleComplex *dv=(cuDoubleComplex *) d_doubleComplex;
  dv = &dv[df-IndexBase];
  i = writeRemoteBuffer((void*)h_doubleComplex, (void*)dv, n*sizeof(cuDoubleComplex));

  return i;
}

int readDoubleComplexFirst(int df, void* d_doubleComplex, cuDoubleComplex* h_doubleComplex,
			  int n, int IndexBase)
{ int i;
  cuDoubleComplex *dv=(cuDoubleComplex *) d_doubleComplex;
  dv = &dv[df-IndexBase];    
  i = readRemoteBuffer((void *) h_doubleComplex, (void *) dv, n*sizeof(cuDoubleComplex));

  return(i);
}

int allocateMultiDoubleComplex(void **d_DoubleComplex, int m, int n)
{
  return allocRemoteBuffer((void **)(d_DoubleComplex), m*n*sizeof(cuDoubleComplex));
}

int writeMultiDoubleComplex(void *d_DoubleComplex, cuDoubleComplex* h_DoubleComplex, int m, int n)
{
  int i,j;
  i = writeRemoteBuffer((void*)h_DoubleComplex, (void*)d_DoubleComplex, m*n*sizeof(cuDoubleComplex));
  return i;
}

int readMultiDoubleComplex(void* d_DoubleComplex, cuDoubleComplex* h_DoubleComplex, int m, int n)
{ int i;
  i = readRemoteBuffer((void *) h_DoubleComplex, (void *) d_DoubleComplex, m*n*sizeof(cuDoubleComplex));
  //cudaSync();
  return(i);
}

void freeDoubleComplex(void *d_DoubleComplex)
{
  freeRemoteBuffer(d_DoubleComplex);
}



double etime() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    exit(-1);
  }
  temp = ((double)tt.tv_sec) + ((double)tt.tv_usec)*1.0e-6;
  return(temp);
}


