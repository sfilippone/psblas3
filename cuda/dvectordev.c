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
 
  

#include <stdio.h>
#include <complex.h>
//#include "utils.h"
//#include "common.h"
#include "dvectordev.h"


int registerMappedDouble(void  *buff, void **d_p, int n, double dummy)
{
  return registerMappedMemory(buff,d_p,n*sizeof(double));
}

int writeMultiVecDeviceDouble(void* deviceVec, double* hostVec)
{ int i;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  // Ex updateFromHost vector function
  i = writeRemoteBuffer((void*) hostVec, (void *)devVec->v_, devVec->pitch_*devVec->count_*sizeof(double));
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","FallocMultiVecDevice",i);
  }
  return(i);
}

int writeMultiVecDeviceDoubleR2(void* deviceVec, double* hostVec, int ld)
{ int i;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  i = writeRemoteBufferR2((void*) hostVec, (void *)devVec->v_, devVec->count_, devVec->pitch_*sizeof(double), devVec->size_*sizeof(double));
//   i = writeMultiVecDeviceDouble(deviceVec, (void *) hostVec);
  fprintf(stderr,"From routine : %s : %p %p\n","writeMultiVecDeviceDoubleR2",devVec->v_,devVec->v_+devVec->pitch_);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","writeMultiVecDeviceDoubleR2",i);
  }
  return(i);
}

int readMultiVecDeviceDouble(void* deviceVec, double* hostVec)
{ int i,j;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  i = readRemoteBuffer((void *) hostVec, (void *)devVec->v_, 
		       devVec->pitch_*devVec->count_*sizeof(double));
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readMultiVecDeviceDouble",i);
  }
  return(i);
}

int readMultiVecDeviceDoubleR2(void* deviceVec, double* hostVec, int ld)
{ int i;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  i = readRemoteBufferR2((void *) hostVec, (void *)devVec->v_, devVec->count_, devVec->pitch_*sizeof(double), devVec->size_*sizeof(double));
//   i = readMultiVecDeviceDouble(deviceVec, hostVec);
  fprintf(stderr,"From routine : %s : %p \n","readMultiVecDeviceDoubleR2",devVec->v_);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readMultiVecDeviceDoubleR2",i);
  }
  return(i);
}

int setscalMultiVecDeviceDouble(double val, int first, int last, 
				int indexBase, void* devMultiVecX) 
{ int i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  spgpuHandle_t handle=psb_cudaGetHandle();

  spgpuDsetscal(handle, first, last, indexBase, val, (double *) devVecX->v_);
  
  return(i);
}


int geinsMultiVecDeviceDouble(int n, void* devMultiVecIrl, void* devMultiVecVal, 
			      int dupl, int indexBase, void* devMultiVecX)
{ int j=0, i=0,nmin=0,nmax=0;
  int pitch = 0;
  double beta;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecIrl = (struct MultiVectDevice *) devMultiVecIrl;
  struct MultiVectDevice *devVecVal = (struct MultiVectDevice *) devMultiVecVal;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecIrl->pitch_;
  if ((n > devVecIrl->size_) || (n>devVecVal->size_ )) 
    return SPGPU_UNSUPPORTED;

  //fprintf(stderr,"geins: %d %d  %p %p %p\n",dupl,n,devVecIrl->v_,devVecVal->v_,devVecX->v_);

  if (dupl == INS_OVERWRITE) 
    beta = 0.0;
  else if (dupl == INS_ADD) 
    beta = 1.0;
  else
    beta = 0.0;
 
  spgpuDscat(handle, (double *) devVecX->v_, n, (double*)devVecVal->v_,
	     (int*)devVecIrl->v_, indexBase, beta);
  
  return(i);
}


int igathMultiVecDeviceDoubleVecIdx(void* deviceVec, int vectorId, int n,
				    int first, void* deviceIdx, int hfirst,
				    void* host_values, int indexBase)
{
  int i, *idx;
  struct MultiVectDevice *devIdx = (struct MultiVectDevice *) deviceIdx;

  i= igathMultiVecDeviceDouble(deviceVec, vectorId, n,
			       first, (void*) devIdx->v_, hfirst, host_values, indexBase);
  return(i);
}

int igathMultiVecDeviceDouble(void* deviceVec, int vectorId, int n,
			      int first, void* indexes, int hfirst, void* host_values, int indexBase)
{
  int i, *idx =(int *) indexes;;
  double *hv = (double *) host_values;;  
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  spgpuHandle_t handle=psb_cudaGetHandle();
  
  i=0; 
  hv  = &(hv[hfirst-indexBase]);
  idx = &(idx[first-indexBase]);
  spgpuDgath(handle,hv, n, idx,indexBase, (double *) devVec->v_+vectorId*devVec->pitch_);
  return(i);
}

int iscatMultiVecDeviceDoubleVecIdx(void* deviceVec, int vectorId, int n, int first, void *deviceIdx,
				    int hfirst, void* host_values, int indexBase, double beta)
{  
  int i, *idx;
  struct MultiVectDevice *devIdx = (struct MultiVectDevice *) deviceIdx;
  i= iscatMultiVecDeviceDouble(deviceVec, vectorId, n, first, 
			       (void*) devIdx->v_,  hfirst,host_values, indexBase, beta);
  return(i);
}

int iscatMultiVecDeviceDouble(void* deviceVec, int vectorId, int n, int first, void *indexes,
			      int hfirst, void* host_values, int indexBase, double beta)
{ int i=0;
  double *hv  = (double *) host_values;
  int *idx=(int *) indexes;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  spgpuHandle_t handle=psb_cudaGetHandle();

  idx = &(idx[first-indexBase]);
  hv  = &(hv[hfirst-indexBase]);
  spgpuDscat(handle, (double *) devVec->v_, n, hv, idx, indexBase, beta);
  return SPGPU_SUCCESS;
  
}

int scalMultiVecDeviceDouble(double alpha, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;
  // Note: inner kernel can handle aliased input/output
  spgpuDscal(handle, (double *)devVecA->v_, devVecA->pitch_, 
	     alpha,  (double *)devVecA->v_);
  return(i);
}

int nrm2MultiVecDeviceDouble(double* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuDmnrm2(handle, y_res, n,(double *)devVecA->v_, devVecA->count_, devVecA->pitch_);
  return(i);
}

int amaxMultiVecDeviceDouble(double* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuDmamax(handle, y_res, n,(double *)devVecA->v_, devVecA->count_, devVecA->pitch_);
  return(i);
}

int asumMultiVecDeviceDouble(double* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuDmasum(handle, y_res, n,(double *)devVecA->v_, devVecA->count_, devVecA->pitch_);

  return(i);
}

int dotMultiVecDeviceDouble(double* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) devMultiVecB;
  spgpuHandle_t handle=psb_cudaGetHandle();

  spgpuDmdot(handle, y_res, n, (double*)devVecA->v_, (double*)devVecB->v_,devVecA->count_,devVecB->pitch_);
  return(i);
}

int axpbyMultiVecDeviceDouble(int n,double alpha, void* devMultiVecX, 
			      double beta, void* devMultiVecY)
{ int j=0, i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecY = (struct MultiVectDevice *) devMultiVecY;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecY->pitch_;
  if ((n > devVecY->size_) || (n>devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  for(j=0;j<devVecY->count_;j++)
    fprintf(stderr,"CUDA ENTERED %d %d %d %d \n",j, n, pitch, devVecY->size_);
    spgpuDaxpby(handle,(double*)devVecY->v_+pitch*j, n, beta, 
		(double*)devVecY->v_+pitch*j, alpha,(double*) devVecX->v_+pitch*j);
  return(i);
}
 
int abgdxyzMultiVecDeviceDouble(int n,double alpha,double beta, double gamma, double delta, 
				void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ)
{ int j=0, i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecY = (struct MultiVectDevice *) devMultiVecY;
  struct MultiVectDevice *devVecZ = (struct MultiVectDevice *) devMultiVecZ;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecY->pitch_;
  if ((n > devVecY->size_) || (n>devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuDabgdxyz(handle,n, alpha,beta,gamma,delta, 
	      (double*)devVecX->v_,(double*) devVecY->v_,(double*) devVecZ->v_);
  return(i);
}
 
int xyzwMultiVecDeviceDouble(int n,double a, double b, double c, double d, double e, double f, 
			     void* devMultiVecX, void* devMultiVecY,
			     void* devMultiVecZ, void* devMultiVecW)
{ int j=0, i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecY = (struct MultiVectDevice *) devMultiVecY;
  struct MultiVectDevice *devVecZ = (struct MultiVectDevice *) devMultiVecZ;
  struct MultiVectDevice *devVecW = (struct MultiVectDevice *) devMultiVecW;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecY->pitch_;
  if ((n > devVecY->size_) || (n>devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuDxyzw(handle,n, a,b,c,d,e,f,
	      (double*)devVecX->v_,(double*) devVecY->v_,(double*) devVecZ->v_,(double*) devVecW->v_);
  return(i);
}
 
int axyMultiVecDeviceDouble(int n, double alpha, void *deviceVecA, void *deviceVecB)
{ int i = 0; 
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;
  spgpuHandle_t handle=psb_cudaGetHandle();
  if ((n > devVecA->size_) || (n>devVecB->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuDmaxy(handle, (double*)devVecB->v_, n, alpha, (double*)devVecA->v_, 
	     (double*)devVecB->v_, devVecA->count_, devVecA->pitch_);

  return(i);
}

int axybzMultiVecDeviceDouble(int n, double alpha, void *deviceVecA,
			      void *deviceVecB, double beta, void *deviceVecZ)
{ int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;
  struct MultiVectDevice *devVecZ = (struct MultiVectDevice *) deviceVecZ;
  spgpuHandle_t handle=psb_cudaGetHandle();

  if ((n > devVecA->size_) || (n>devVecB->size_ ) || (n>devVecZ->size_ )) 
    return SPGPU_UNSUPPORTED;
  spgpuDmaxypbz(handle, (double*)devVecZ->v_, n, beta, (double*)devVecZ->v_, 
	       alpha, (double*) devVecA->v_, (double*) devVecB->v_,
	       devVecB->count_, devVecB->pitch_);
  return(i);
}

int absMultiVecDeviceDouble2(int n, double alpha, void *deviceVecA,
			      void *deviceVecB)
{ int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;

  spgpuHandle_t handle=psb_cudaGetHandle();

  if ((n > devVecA->size_) || (n>devVecB->size_ ))
    return SPGPU_UNSUPPORTED;

  spgpuDabs(handle, (double*)devVecB->v_, n, alpha, (double*)devVecA->v_);

  return(i);
}
 
int absMultiVecDeviceDouble(int n, double alpha, void *deviceVecA)
{ int i = 0; 
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  spgpuHandle_t handle=psb_cudaGetHandle();
  if (n > devVecA->size_)
    return SPGPU_UNSUPPORTED;

  spgpuDabs(handle, (double*)devVecA->v_, n, alpha, (double*)devVecA->v_);

  return(i);
}

