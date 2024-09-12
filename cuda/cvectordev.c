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
#include "cvectordev.h"


int registerMappedFloatComplex(void  *buff, void **d_p, int n, cuFloatComplex dummy)
{
  return registerMappedMemory(buff,d_p,n*sizeof(cuFloatComplex));
}

int writeMultiVecDeviceFloatComplex(void* deviceVec, cuFloatComplex* hostVec)
{ int i;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  // Ex updateFromHost vector function
  i = writeRemoteBuffer((void*) hostVec, (void *)devVec->v_, 
			devVec->pitch_*devVec->count_*sizeof(cuFloatComplex));
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","FallocMultiVecDevice",i);
  }
  return(i);
}

int writeMultiVecDeviceFloatComplexR2(void* deviceVec, cuFloatComplex* hostVec, int ld)
{ int i;
  i = writeMultiVecDeviceFloatComplex(deviceVec, (void *) hostVec);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","writeMultiVecDeviceFloatComplexR2",i);
  }
  return(i);
}

int readMultiVecDeviceFloatComplex(void* deviceVec, cuFloatComplex* hostVec)
{ int i,j;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  i = readRemoteBuffer((void *) hostVec, (void *)devVec->v_, 
		       devVec->pitch_*devVec->count_*sizeof(cuFloatComplex));
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readMultiVecDeviceFloat",i);
  }
  return(i);
}

int readMultiVecDeviceFloatComplexR2(void* deviceVec, cuFloatComplex* hostVec, int ld)
{ int i;
  i = readMultiVecDeviceFloatComplex(deviceVec, hostVec);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readMultiVecDeviceFloatComplexR2",i);
  }
  return(i);
}

int setscalMultiVecDeviceFloatComplex(cuFloatComplex val, int first, int last, 
				int indexBase, void* devMultiVecX) 
{ int i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  spgpuHandle_t handle=psb_cudaGetHandle();

  spgpuCsetscal(handle, first, last, indexBase, val, (cuFloatComplex *) devVecX->v_);
  
  return(i);
}

int geinsMultiVecDeviceFloatComplex(int n, void* devMultiVecIrl, void* devMultiVecVal, 
			      int dupl, int indexBase, void* devMultiVecX)
{ int j=0, i=0,nmin=0,nmax=0;
  int pitch = 0;
  cuFloatComplex beta;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecIrl = (struct MultiVectDevice *) devMultiVecIrl;
  struct MultiVectDevice *devVecVal = (struct MultiVectDevice *) devMultiVecVal;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecIrl->pitch_;
  if ((n > devVecIrl->size_) || (n>devVecVal->size_ )) 
    return SPGPU_UNSUPPORTED;

  //fprintf(stderr,"geins: %d %d  %p %p %p\n",dupl,n,devVecIrl->v_,devVecVal->v_,devVecX->v_);

  if (dupl == INS_OVERWRITE) 
    beta =  make_cuFloatComplex(0.0, 0.0);
  else if (dupl == INS_ADD) 
    beta =  make_cuFloatComplex(1.0, 0.0);
  else
    beta = make_cuFloatComplex(0.0, 0.0);
 
  spgpuCscat(handle, (cuFloatComplex *) devVecX->v_, n, (cuFloatComplex*)devVecVal->v_,
	     (int*)devVecIrl->v_, indexBase, beta);
  
  return(i);
}


int igathMultiVecDeviceFloatComplexVecIdx(void* deviceVec, int vectorId, int n,
				    int first, void* deviceIdx, int hfirst,
				    void* host_values, int indexBase)
{
  int i, *idx;
  struct MultiVectDevice *devIdx = (struct MultiVectDevice *) deviceIdx;

  i= igathMultiVecDeviceFloatComplex(deviceVec, vectorId, n,
			       first, (void*) devIdx->v_, hfirst, host_values, indexBase);
  return(i);
}

int igathMultiVecDeviceFloatComplex(void* deviceVec, int vectorId, int n,
				    int first, void* indexes, int hfirst, 
				    void* host_values, int indexBase)
{
  int i, *idx =(int *) indexes;;
  cuFloatComplex *hv = (cuFloatComplex *) host_values;;  
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  spgpuHandle_t handle=psb_cudaGetHandle();
  
  i=0; 
  hv  = &(hv[hfirst-indexBase]);
  idx = &(idx[first-indexBase]);
  spgpuCgath(handle,hv, n, idx,indexBase,
	     (cuFloatComplex *) devVec->v_+vectorId*devVec->pitch_);
  return(i);
}

int iscatMultiVecDeviceFloatComplexVecIdx(void* deviceVec, int vectorId, int n, 
				   int first, void *deviceIdx,
				   int hfirst, void* host_values, 
				   int indexBase, cuFloatComplex beta)
{  
  int i, *idx;
  struct MultiVectDevice *devIdx = (struct MultiVectDevice *) deviceIdx;
  i= iscatMultiVecDeviceFloatComplex(deviceVec, vectorId, n, first, 
				     (void*) devIdx->v_,  hfirst,host_values, 
				     indexBase, beta);
  return(i);
}

int iscatMultiVecDeviceFloatComplex(void* deviceVec, int vectorId, int n, 
				    int first, void *indexes,
				    int hfirst, void* host_values, 
				    int indexBase, cuFloatComplex beta)
{ int i=0;
  cuFloatComplex *hv  = (cuFloatComplex *) host_values;
  int *idx=(int *) indexes;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  spgpuHandle_t handle=psb_cudaGetHandle();

  idx = &(idx[first-indexBase]);
  hv  = &(hv[hfirst-indexBase]);
  spgpuCscat(handle, (cuFloatComplex *) devVec->v_, n, hv, idx, indexBase, beta);
  return SPGPU_SUCCESS;
  
}


int nrm2MultiVecDeviceFloatComplex(cuFloatComplex* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuCmnrm2(handle, y_res, n,(cuFloatComplex *)devVecA->v_, 
	      devVecA->count_, devVecA->pitch_);
  return(i);
}

int amaxMultiVecDeviceFloatComplex(cuFloatComplex* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuCmamax(handle, y_res, n,(cuFloatComplex *)devVecA->v_, 
	      devVecA->count_, devVecA->pitch_);
  return(i);
}

int asumMultiVecDeviceFloatComplex(cuFloatComplex* y_res, int n, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;

  spgpuCmasum(handle, y_res, n,(cuFloatComplex *)devVecA->v_, 
	      devVecA->count_, devVecA->pitch_);

  return(i);
}

int scalMultiVecDeviceFloatComplex(cuFloatComplex alpha, void* devMultiVecA)
{ int i=0;
  spgpuHandle_t handle=psb_cudaGetHandle();
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;
  // Note: inner kernel can handle aliased input/output
  spgpuCscal(handle, (cuFloatComplex *)devVecA->v_, devVecA->pitch_, 
	     alpha,  (cuFloatComplex *)devVecA->v_);
  return(i);
}

int dotMultiVecDeviceFloatComplex(cuFloatComplex* y_res, int n, 
				  void* devMultiVecA, void* devMultiVecB)
{int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) devMultiVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) devMultiVecB;
  spgpuHandle_t handle=psb_cudaGetHandle();

  spgpuCmdot(handle, y_res, n, (cuFloatComplex*)devVecA->v_, 
	     (cuFloatComplex*)devVecB->v_,devVecA->count_,devVecB->pitch_);
  return(i);
}

int axpbyMultiVecDeviceFloatComplex(int n,cuFloatComplex alpha, void* devMultiVecX, 
				    cuFloatComplex beta, void* devMultiVecY)
{ int j=0, i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecY = (struct MultiVectDevice *) devMultiVecY;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecY->pitch_;
  if ((n > devVecY->size_) || (n>devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  for(j=0;j<devVecY->count_;j++)
    spgpuCaxpby(handle,(cuFloatComplex*)devVecY->v_+pitch*j, n, beta, 
		(cuFloatComplex*)devVecY->v_+pitch*j, alpha,
		(cuFloatComplex*) devVecX->v_+pitch*j);
  return(i);
}

int upd_xyzMultiVecDeviceFloatComplex(int n,cuFloatComplex  alpha,cuFloatComplex  beta,
				      cuFloatComplex  gamma, cuFloatComplex  delta, 
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

  spgpuCupd_xyz(handle,n, alpha,beta,gamma,delta, 
		(cuFloatComplex *)devVecX->v_,(cuFloatComplex *) devVecY->v_,(cuFloatComplex *) devVecZ->v_);
  return(i);
}

int xyzwMultiVecDeviceFloatComplex(int n,cuFloatComplex  a,cuFloatComplex  b,
				   cuFloatComplex c, cuFloatComplex  d,
				   cuFloatComplex e, cuFloatComplex  f, 
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

  spgpuCxyzw(handle,n, a,b,c,d,e,f,
	     (cuFloatComplex *)devVecX->v_,(cuFloatComplex *) devVecY->v_,
	     (cuFloatComplex *) devVecZ->v_,(cuFloatComplex *) devVecW->v_);
  return(i);
}

int axyMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, 
				  void *deviceVecA, void *deviceVecB)
{ int i = 0; 
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;
  spgpuHandle_t handle=psb_cudaGetHandle();
  if ((n > devVecA->size_) || (n>devVecB->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuCmaxy(handle, (cuFloatComplex*)devVecB->v_, n, alpha,
	     (cuFloatComplex*)devVecA->v_, 
	     (cuFloatComplex*)devVecB->v_, devVecA->count_, devVecA->pitch_);

  return(i);
}

int axybzMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, void *deviceVecA,
				    void *deviceVecB, cuFloatComplex beta, 
				    void *deviceVecZ)
{ int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;
  struct MultiVectDevice *devVecZ = (struct MultiVectDevice *) deviceVecZ;
  spgpuHandle_t handle=psb_cudaGetHandle();

  if ((n > devVecA->size_) || (n>devVecB->size_ ) || (n>devVecZ->size_ )) 
    return SPGPU_UNSUPPORTED;
  spgpuCmaxypbz(handle, (cuFloatComplex*)devVecZ->v_, n, beta, 
		(cuFloatComplex*)devVecZ->v_, 
	       alpha, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_,
	       devVecB->count_, devVecB->pitch_);
  return(i);
}


int absMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void *deviceVecA,
			      void *deviceVecB)
{ int i=0;
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  struct MultiVectDevice *devVecB = (struct MultiVectDevice *) deviceVecB;

  spgpuHandle_t handle=psb_cudaGetHandle();

  if ((n > devVecA->size_) || (n>devVecB->size_ ))
    return SPGPU_UNSUPPORTED;

  spgpuCabs(handle, (cuFloatComplex*)devVecB->v_, n, 
	    alpha, (cuFloatComplex*)devVecA->v_);

  return(i);
}
 
int absMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, void *deviceVecA)
{ int i = 0; 
  struct MultiVectDevice *devVecA = (struct MultiVectDevice *) deviceVecA;
  spgpuHandle_t handle=psb_cudaGetHandle();
  if (n > devVecA->size_)
    return SPGPU_UNSUPPORTED;

  spgpuCabs(handle, (cuFloatComplex*)devVecA->v_, n, 
	    alpha, (cuFloatComplex*)devVecA->v_);

  return(i);
}


