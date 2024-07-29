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
#include "zvectordev.h"


int registerMappedDoubleComplex(void  *buff, void **d_p, int n, cuDoubleComplex dummy)
{
  return registerMappedMemory(buff,d_p,n*sizeof(cuDoubleComplex));
}

int writeMultiVecDeviceDoubleComplex(void* deviceVec, cuDoubleComplex* hostVec)
{ int i;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  // Ex updateFromHost vector function
  i = writeRemoteBuffer((void*) hostVec, (void *)devVec->v_, 
			devVec->pitch_*devVec->count_*sizeof(cuDoubleComplex));
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","FallocMultiVecDevice",i);
  }
  return(i);
}

int writeMultiVecDeviceDoubleComplexR2(void* deviceVec, cuDoubleComplex* hostVec, int ld)
{ int i;
  i = writeMultiVecDeviceDoubleComplex(deviceVec, (void *) hostVec);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","writeMultiVecDeviceDoubleComplexR2",i);
  }
  return(i);
}

int readMultiVecDeviceDoubleComplex(void* deviceVec, cuDoubleComplex* hostVec)
{ int i,j;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  i = readRemoteBuffer((void *) hostVec, (void *)devVec->v_, 
		       devVec->pitch_*devVec->count_*sizeof(cuDoubleComplex));
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readMultiVecDeviceDoubleComplex",i);
  }
  return(i);
}

int readMultiVecDeviceDoubleComplexR2(void* deviceVec, cuDoubleComplex* hostVec, int ld)
{ int i;
  i = readMultiVecDeviceDoubleComplex(deviceVec, hostVec);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readMultiVecDeviceDoubleComplexR2",i);
  }
  return(i);
}

int setscalMultiVecDeviceDoubleComplex(cuDoubleComplex val, int first, int last, 
				int indexBase, void* devMultiVecX) 
{ int i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  spgpuHandle_t handle=psb_cudaGetHandle();

  spgpuZsetscal(handle, first, last, indexBase, val, (cuDoubleComplex *) devVecX->v_);
  
  return(i);
}

int geinsMultiVecDeviceDoubleComplex(int n, void* devMultiVecIrl, void* devMultiVecVal, 
				     int dupl, int indexBase, void* devMultiVecX)
{ int j=0, i=0,nmin=0,nmax=0;
  int pitch = 0;
  cuDoubleComplex beta;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecIrl = (struct MultiVectDevice *) devMultiVecIrl;
  struct MultiVectDevice *devVecVal = (struct MultiVectDevice *) devMultiVecVal;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecIrl->pitch_;
  if ((n > devVecIrl->size_) || (n>devVecVal->size_ )) 
    return SPGPU_UNSUPPORTED;

  //fprintf(stderr,"geins: %d %d  %p %p %p\n",dupl,n,devVecIrl->v_,devVecVal->v_,devVecX->v_);
  if (dupl == INS_OVERWRITE) 
    beta =  make_cuDoubleComplex(0.0, 0.0);
  else if (dupl == INS_ADD) 
    beta =  make_cuDoubleComplex(1.0, 0.0);
  else
    beta = make_cuDoubleComplex(0.0, 0.0);
 
  spgpuZscat(handle, (cuDoubleComplex *) devVecX->v_, n, (cuDoubleComplex*)devVecVal->v_,
	     (int*)devVecIrl->v_, indexBase, beta);
  
  return(i);
}


int igathMultiVecDeviceDoubleComplexVecIdx(void* deviceVec, int vectorId, int n,
				    int first, void* deviceIdx, int hfirst,
				    void* host_values, int indexBase)
{
  int i, *idx;
  struct MultiVectDevice *devIdx = (struct MultiVectDevice *) deviceIdx;

  i= igathMultiVecDeviceDoubleComplex(deviceVec, vectorId, n,
				      first, (void*) devIdx->v_, 
				      hfirst, host_values, indexBase);
  return(i);
}

int igathMultiVecDeviceDoubleComplex(void* deviceVec, int vectorId, int n,
				     int first, void* indexes, int hfirst, 
				     void* host_values, int indexBase)
{
  int i, *idx =(int *) indexes;;
  cuDoubleComplex *hv = (cuDoubleComplex *) host_values;;  
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  spgpuHandle_t handle=psb_cudaGetHandle();
  
  i=0; 
  hv  = &(hv[hfirst-indexBase]);
  idx = &(idx[first-indexBase]);
  spgpuZgath(handle,hv, n, idx,indexBase, 
	     (cuDoubleComplex *) devVec->v_+vectorId*devVec->pitch_);
  return(i);
}

int iscatMultiVecDeviceDoubleComplexVecIdx(void* deviceVec, int vectorId, int n, 
					   int first, void *deviceIdx,
					   int hfirst, void* host_values, 
					   int indexBase, cuDoubleComplex beta)
{  
  int i, *idx;
  struct MultiVectDevice *devIdx = (struct MultiVectDevice *) deviceIdx;
  i= iscatMultiVecDeviceDoubleComplex(deviceVec, vectorId, n, first, 
			       (void*) devIdx->v_,  hfirst,host_values, indexBase, beta);
  return(i);
}

int iscatMultiVecDeviceDoubleComplex(void* deviceVec, int vectorId, int n, 
				     int first, void *indexes,
				     int hfirst, void* host_values,
				     int indexBase, cuDoubleComplex beta)
{ int i=0;
  cuDoubleComplex *hv  = (cuDoubleComplex *) host_values;
  int *idx=(int *) indexes;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  spgpuHandle_t handle=psb_cudaGetHandle();

  idx = &(idx[first-indexBase]);
  hv  = &(hv[hfirst-indexBase]);
  spgpuZscat(handle, (cuDoubleComplex *) devVec->v_, n, hv, idx, indexBase, beta);
  return SPGPU_SUCCESS;
  
}


int nrm2MultiVecDeviceDoubleComplex(cuDoubleComplex* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuZmnrm2(handle, y_res, n,(cuDoubleComplex *)devVecA->v_, devVecA->count_, devVecA->pitch_);
  return(i);
}

int amaxMultiVecDeviceDoubleComplex(cuDoubleComplex* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuZmamax(handle, y_res, n,(cuDoubleComplex *)devVecA->v_, 
	      devVecA->count_, devVecA->pitch_);
  return(i);
}

int asumMultiVecDeviceDoubleComplex(cuDoubleComplex* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuZmasum(handle, y_res, n,(cuDoubleComplex *)devVecA->v_, 
	      devVecA->count_, devVecA->pitch_);

  return(i);
}

int scalMultiVecDeviceDoubleComplex(cuDoubleComplex alpha, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;
  // Note: inner kernel can handle aliased input/output
  spgpuZscal(handle, (cuDoubleComplex *)devVecA->v_, devVecA->pitch_, 
	     alpha,  (cuDoubleComplex *)devVecA->v_);
  return(i);
}

int dotMultiVecDeviceDoubleComplex(cuDoubleComplex* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) devMultiVecB;
  spgpuHandle_t handle=psb_cudaGetHandle();

  spgpuZmdot(handle, y_res, n, (cuDoubleComplex*)devVecA->v_,
	     (cuDoubleComplex*)devVecB->v_,devVecA->count_,devVecB->pitch_);
  return(i);
}

int upd_xyzMultiVecDeviceDoubleComplex(int n,cuDoubleComplex  alpha,
				       cuDoubleComplex  beta, cuDoubleComplex  gamma, cuDoubleComplex  delta, 
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

  spgpuZupd_xyz(handle,n, alpha,beta,gamma,delta, 
	      (cuDoubleComplex *)devVecX->v_,(cuDoubleComplex *) devVecY->v_,(cuDoubleComplex *) devVecZ->v_);
  return(i);
}

int xyzwMultiVecDeviceDoubleComplex(int n,cuDoubleComplex  a,   cuDoubleComplex  b,
				    cuDoubleComplex  c, cuDoubleComplex  d,
				    cuDoubleComplex  e, cuDoubleComplex  f, 
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

  spgpuZxyzw(handle,n, a,b,c,d,e,f,
	      (cuDoubleComplex *)devVecX->v_,(cuDoubleComplex *) devVecY->v_,
	      (cuDoubleComplex *) devVecZ->v_,(cuDoubleComplex *) devVecW->v_);
  return(i);
}

int axpbyMultiVecDeviceDoubleComplex(int n,cuDoubleComplex alpha, void* devMultiVecX, 
			      cuDoubleComplex beta, void* devMultiVecY)
{ int j=0, i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecY = (struct MultiVectDevice *) devMultiVecY;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecY->pitch_;
  if ((n > devVecY->size_) || (n>devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  for(j=0;j<devVecY->count_;j++)
    spgpuZaxpby(handle,(cuDoubleComplex*)devVecY->v_+pitch*j, n, beta, 
		(cuDoubleComplex*)devVecY->v_+pitch*j, alpha,
		(cuDoubleComplex*) devVecX->v_+pitch*j);
  return(i);
}

int axyMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, 
				   void *deviceVecA, void *deviceVecB)
{ int i = 0; 
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;
  spgpuHandle_t handle=psb_cudaGetHandle();
  if ((n > devVecA->size_) || (n>devVecB->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuZmaxy(handle, (cuDoubleComplex*)devVecB->v_, n, alpha,
	     (cuDoubleComplex*)devVecA->v_, 
	     (cuDoubleComplex*)devVecB->v_, devVecA->count_, devVecA->pitch_);

  return(i);
}

int axybzMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA,
			      void *deviceVecB, cuDoubleComplex beta, void *deviceVecZ)
{ int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;
  struct MultiVectDevice *devVecZ = (struct MultiVectDevice *) deviceVecZ;
  spgpuHandle_t handle=psb_cudaGetHandle();

  if ((n > devVecA->size_) || (n>devVecB->size_ ) || (n>devVecZ->size_ )) 
    return SPGPU_UNSUPPORTED;
  spgpuZmaxypbz(handle, (cuDoubleComplex*)devVecZ->v_, n, beta, 
		(cuDoubleComplex*)devVecZ->v_, 
	       alpha, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_,
	       devVecB->count_, devVecB->pitch_);
  return(i);
}


int absMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void *deviceVecA,
			      void *deviceVecB)
{ int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;

  spgpuHandle_t handle=psb_cudaGetHandle();

  if ((n > devVecA->size_) || (n>devVecB->size_ ))
    return SPGPU_UNSUPPORTED;

  spgpuZabs(handle, (cuDoubleComplex*)devVecB->v_, n,
	    alpha, (cuDoubleComplex*)devVecA->v_); 

  return(i);
}
 
int absMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA)
{ int i = 0; 
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  spgpuHandle_t handle=psb_cudaGetHandle();
  if (n > devVecA->size_)
    return SPGPU_UNSUPPORTED;

  spgpuZabs(handle, (cuDoubleComplex*)devVecA->v_, n, 
	    alpha, (cuDoubleComplex*)devVecA->v_);

  return(i);
}

