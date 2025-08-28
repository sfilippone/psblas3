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
 
#include <sys/time.h>
#include "dnsdev.h"

#define PASS_RS  0

#define IMIN(a,b) ((a)<(b) ? (a) : (b))

DnsDeviceParams getDnsDeviceParams(unsigned int rows, unsigned int columns,
				   unsigned int elementType, unsigned int firstIndex)
{
  DnsDeviceParams params;

  if (elementType == SPGPU_TYPE_DOUBLE)
    {
      params.pitch = ((rows + ELL_PITCH_ALIGN_D - 1)/ELL_PITCH_ALIGN_D)*ELL_PITCH_ALIGN_D;
    }
  else
    {
      params.pitch = ((rows + ELL_PITCH_ALIGN_S - 1)/ELL_PITCH_ALIGN_S)*ELL_PITCH_ALIGN_S;
    }
  //For complex?
  params.elementType = elementType;
  params.rows = rows;
  params.columns = columns;
  params.firstIndex = firstIndex;

  return params;

}
//new
int allocDnsDevice(void ** remoteMatrix, DnsDeviceParams* params)
{
  struct DnsDevice *tmp = (struct DnsDevice *)malloc(sizeof(struct DnsDevice));
  *remoteMatrix = (void *)tmp;
  tmp->rows = params->rows;
  tmp->columns = params->columns;
  tmp->cMPitch = params->pitch;
  tmp->pitch= tmp->cMPitch;
  tmp->allocsize = (int)tmp->columns * tmp->pitch;
  tmp->baseIndex = params->firstIndex;
  //fprintf(stderr,"allocDnsDevice: %d %d %d \n",tmp->pitch, params->maxRowSize, params->avgRowSize);
  if (params->elementType == SPGPU_TYPE_FLOAT)
    allocRemoteBuffer((void **)&(tmp->cM), tmp->allocsize*sizeof(float));
  else if (params->elementType == SPGPU_TYPE_DOUBLE)
    allocRemoteBuffer((void **)&(tmp->cM), tmp->allocsize*sizeof(double));
  else if (params->elementType == SPGPU_TYPE_COMPLEX_FLOAT)
    allocRemoteBuffer((void **)&(tmp->cM), tmp->allocsize*sizeof(cuFloatComplex));
  else if (params->elementType == SPGPU_TYPE_COMPLEX_DOUBLE)
    allocRemoteBuffer((void **)&(tmp->cM), tmp->allocsize*sizeof(cuDoubleComplex));
  else
    return SPGPU_UNSUPPORTED; // Unsupported params
  //fprintf(stderr,"From allocDnsDevice: %d %d %d %p %p %p\n",tmp->maxRowSize,
  //	  tmp->avgRowSize,tmp->allocsize,tmp->rS,tmp->rP,tmp->cM);

  return SPGPU_SUCCESS;
}

void freeDnsDevice(void* remoteMatrix)
{
  struct DnsDevice *devMat = (struct DnsDevice *) remoteMatrix;  
  //fprintf(stderr,"freeDnsDevice\n");
  if (devMat != NULL) {
    freeRemoteBuffer(devMat->cM);
    free(remoteMatrix);
  }
}

//new
int FallocDnsDevice(void** deviceMat, unsigned int rows,
		    unsigned int columns, unsigned int elementType, 
		    unsigned int firstIndex)
{ int i;
  DnsDeviceParams p;

  p = getDnsDeviceParams(rows, columns, elementType, firstIndex);
  i = allocDnsDevice(deviceMat, &p);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","FallocDnsDevice",i);
  }
  return(i);
}


int spmvDnsDeviceFloat(char transa, int m, int n, int k, float *alpha,
			void *deviceMat, void* deviceX, float *beta, void* deviceY)
{
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  int status;

  cublasHandle_t handle=psb_cudaGetCublasHandle();
  cublasOperation_t trans=((transa == 'N')? CUBLAS_OP_N:((transa=='T')? CUBLAS_OP_T:CUBLAS_OP_C));
  /* Note: the M,N,K choices according to TRANS have already been handled in the caller */  
  if (n == 1) {
    status = cublasSgemv(handle, trans, m,k,
			 alpha, devMat->cM,devMat->pitch, x->v_,1,
			 beta,  y->v_,1);
  } else {
    status = cublasSgemm(handle, trans, CUBLAS_OP_N, m,n,k,
			 alpha, devMat->cM,devMat->pitch, x->v_,x->pitch_,
			 beta,  y->v_,y->pitch_);
  }    
  
  if (status == CUBLAS_STATUS_SUCCESS)  
    return SPGPU_SUCCESS;
  else
    return SPGPU_UNSUPPORTED;
}

int spmvDnsDeviceDouble(char transa, int m, int n, int k, double *alpha,
			void *deviceMat, void* deviceX, double *beta, void* deviceY)
{
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  int status;

  cublasHandle_t handle=psb_cudaGetCublasHandle();
  cublasOperation_t trans=((transa == 'N')? CUBLAS_OP_N:((transa=='T')? CUBLAS_OP_T:CUBLAS_OP_C));
  /* Note: the M,N,K choices according to TRANS have already been handled in the caller */  
  if (n == 1) {
    status = cublasDgemv(handle, trans, m,k,
			 alpha, devMat->cM,devMat->pitch, x->v_,1,
			 beta,  y->v_,1);
  } else {
    status = cublasDgemm(handle, trans, CUBLAS_OP_N, m,n,k,
			 alpha, devMat->cM,devMat->pitch, x->v_,x->pitch_,
			 beta,  y->v_,y->pitch_);
  }    
  
  if (status == CUBLAS_STATUS_SUCCESS)  
    return SPGPU_SUCCESS;
  else
    return SPGPU_UNSUPPORTED;
}

int spmvDnsDeviceFloatComplex(char transa, int m, int n, int k, float complex *alpha,
			void *deviceMat, void* deviceX, float complex *beta, void* deviceY)
{
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  int status;

  cublasHandle_t handle=psb_cudaGetCublasHandle();
  cublasOperation_t trans=((transa == 'N')? CUBLAS_OP_N:((transa=='T')? CUBLAS_OP_T:CUBLAS_OP_C));
  /* Note: the M,N,K choices according to TRANS have already been handled in the caller */  
  if (n == 1) {
    status = cublasCgemv(handle, trans, m,k,
			 (const cuComplex *) alpha, devMat->cM,devMat->pitch, x->v_,1,
			 (const cuComplex *) beta,  y->v_,1);
  } else {
    status = cublasCgemm(handle, trans, CUBLAS_OP_N, m,n,k,
			 (const cuComplex *) alpha, devMat->cM,devMat->pitch, x->v_,x->pitch_,
			 (const cuComplex *) beta,  y->v_,y->pitch_);
  }    
  
  if (status == CUBLAS_STATUS_SUCCESS)  
    return SPGPU_SUCCESS;
  else
    return SPGPU_UNSUPPORTED;
}

int spmvDnsDeviceDoubleComplex(char transa, int m, int n, int k, double complex *alpha,
			void *deviceMat, void* deviceX, double complex *beta, void* deviceY)
{
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  int status;

  cublasHandle_t handle=psb_cudaGetCublasHandle();
  cublasOperation_t trans=((transa == 'N')? CUBLAS_OP_N:((transa=='T')? CUBLAS_OP_T:CUBLAS_OP_C));
  /* Note: the M,N,K choices according to TRANS have already been handled in the caller */  
  if (n == 1) {
    status = cublasZgemv(handle, trans, m,k,
			 (const cuDoubleComplex *) alpha, devMat->cM,devMat->pitch, x->v_,1,
			 (const cuDoubleComplex *) beta,  y->v_,1);
  } else {
    status = cublasZgemm(handle, trans, CUBLAS_OP_N, m,n,k,
			 (const cuDoubleComplex *) alpha, devMat->cM,devMat->pitch, x->v_,x->pitch_,
			 (const cuDoubleComplex *) beta,  y->v_,y->pitch_);
  }    
  
  if (status == CUBLAS_STATUS_SUCCESS)  
    return SPGPU_SUCCESS;
  else
    return SPGPU_UNSUPPORTED;
}


int writeDnsDeviceFloat(void* deviceMat, float* val, int lda, int nc)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  int pitch=devMat->pitch; 
  i = cublasSetMatrix(lda,nc,sizeof(float), (void*) val,lda, (void *)devMat->cM, pitch);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","writeDnsDeviceFloat",i);
  }
  return SPGPU_SUCCESS;
}

int writeDnsDeviceDouble(void* deviceMat, double* val, int lda, int nc)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  int pitch=devMat->pitch; 
  i = cublasSetMatrix(lda,nc,sizeof(double), (void*) val,lda, (void *)devMat->cM, pitch);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","writeDnsDeviceDouble",i);
  }
  return SPGPU_SUCCESS;
}


int writeDnsDeviceFloatComplex(void* deviceMat, float complex* val, int lda, int nc)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  int pitch=devMat->pitch; 
  i = cublasSetMatrix(lda,nc,sizeof(cuFloatComplex), (void*) val,lda, (void *)devMat->cM, pitch);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","writeDnsDeviceFloatComplex",i);
  }
  return SPGPU_SUCCESS;
}

int writeDnsDeviceDoubleComplex(void* deviceMat, double complex* val, int lda, int nc)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  int pitch=devMat->pitch; 
  i = cublasSetMatrix(lda,nc,sizeof(cuDoubleComplex), (void*) val,lda, (void *)devMat->cM, pitch);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","writeDnsDeviceDoubleComplex",i);
  }
  return SPGPU_SUCCESS;
}


int readDnsDeviceFloat(void* deviceMat, float* val, int lda, int nc)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  int pitch=devMat->pitch; 
  i = cublasGetMatrix(lda,nc,sizeof(float), (void*) val,lda, (void *)devMat->cM, pitch);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readDnsDeviceFloat",i);
  }
  return SPGPU_SUCCESS;
}

int readDnsDeviceDouble(void* deviceMat, double* val, int lda, int nc)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  int pitch=devMat->pitch; 
  i = cublasGetMatrix(lda,nc,sizeof(double), (void*) val,lda, (void *)devMat->cM, pitch);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readDnsDeviceDouble",i);
  }
  return SPGPU_SUCCESS;
}


int readDnsDeviceFloatComplex(void* deviceMat, float complex* val, int lda, int nc)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  int pitch=devMat->pitch; 
  i = cublasGetMatrix(lda,nc,sizeof(cuFloatComplex), (void*) val,lda, (void *)devMat->cM, pitch);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readDnsDeviceFloatComplex",i);
  }
  return SPGPU_SUCCESS;
}

int readDnsDeviceDoubleComplex(void* deviceMat, double complex* val, int lda, int nc)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  int pitch=devMat->pitch; 
  i = cublasGetMatrix(lda,nc,sizeof(cuDoubleComplex), (void*) val,lda, (void *)devMat->cM, pitch);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readDnsDeviceDoubleComplex",i);
  }
  return SPGPU_SUCCESS;
}


int getDnsDevicePitch(void* deviceMat)
{ int i;
  struct DnsDevice *devMat = (struct DnsDevice *) deviceMat;
  i = devMat->pitch; 
  return(i);
}

