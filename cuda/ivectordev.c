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
#include "ivectordev.h"


int registerMappedInt(void  *buff, void **d_p, int n, int dummy)
{
  return registerMappedMemory(buff,d_p,n*sizeof(int));
}

int writeMultiVecDeviceInt(void* deviceVec, int* hostVec)
{ int i;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  i = writeRemoteBuffer((void*) hostVec, (void *)devVec->v_, 
			devVec->pitch_*devVec->count_*sizeof(int));
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","FallocMultiVecDevice",i);
  }

  return(i);
}

int writeMultiVecDeviceIntR2(void* deviceVec, int* hostVec, int ld)
{ int i;
  i = writeMultiVecDeviceInt(deviceVec, (void *) hostVec);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","writeMultiVecDeviceIntR2",i);
  }
  return(i);
}

int readMultiVecDeviceInt(void* deviceVec, int* hostVec)
{ int i,j;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  i = readRemoteBuffer((void *) hostVec, (void *)devVec->v_, 
		       devVec->pitch_*devVec->count_*sizeof(int));
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readMultiVecDeviceInt",i);
  }
  return(i);
}

int readMultiVecDeviceIntR2(void* deviceVec, int* hostVec, int ld)
{ int i;
  i = readMultiVecDeviceInt(deviceVec, hostVec);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readMultiVecDeviceIntR2",i);
  }
  return(i);
}


int setscalMultiVecDeviceInt(int val, int first, int last, 
				int indexBase, void* devMultiVecX) 
{ int i=0;
  int pitch = 0;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  spgpuHandle_t handle=psb_cudaGetHandle();

  spgpuIsetscal(handle, first, last, indexBase, val, (int *) devVecX->v_);
  
  return(i);
}

int geinsMultiVecDeviceInt(int n, void* devMultiVecIrl, void* devMultiVecVal, 
			      int dupl, int indexBase, void* devMultiVecX)
{ int j=0, i=0,nmin=0,nmax=0;
  int pitch = 0;
  int beta;
  struct MultiVectDevice *devVecX = (struct MultiVectDevice *) devMultiVecX;
  struct MultiVectDevice *devVecIrl = (struct MultiVectDevice *) devMultiVecIrl;
  struct MultiVectDevice *devVecVal = (struct MultiVectDevice *) devMultiVecVal;
  spgpuHandle_t handle=psb_cudaGetHandle();
  pitch = devVecIrl->pitch_;
  if ((n > devVecIrl->size_) || (n>devVecVal->size_ )) 
    return SPGPU_UNSUPPORTED;

  //fprintf(stderr,"geins: %d %d  %p %p %p\n",dupl,n,devVecIrl->v_,devVecVal->v_,devVecX->v_);

  if (dupl == INS_OVERWRITE) 
    beta = 0;
  else if (dupl == INS_ADD) 
    beta = 1;
  else
    beta = 0;
 
  spgpuIscat(handle, (int *) devVecX->v_, n, (int *)devVecVal->v_,
	     (int*)devVecIrl->v_, indexBase, beta);
  
  return(i);
}


int igathMultiVecDeviceIntVecIdx(void* deviceVec, int vectorId, int n,
				    int first, void* deviceIdx, int hfirst,
				    void* host_values, int indexBase)
{
  int i, *idx;
  struct MultiVectDevice *devIdx = (struct MultiVectDevice *) deviceIdx;

  i= igathMultiVecDeviceInt(deviceVec, vectorId, n,
			       first, (void*) devIdx->v_, hfirst, host_values, indexBase);
  return(i);
}

int igathMultiVecDeviceInt(void* deviceVec, int vectorId, int n,
			      int first, void* indexes, int hfirst, void* host_values, int indexBase)
{
  int i, *idx =(int *) indexes;;
  int *hv = (int *) host_values;;  
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  spgpuHandle_t handle=psb_cudaGetHandle();
  
  i=0; 
  hv  = &(hv[hfirst-indexBase]);
  idx = &(idx[first-indexBase]);
  spgpuIgath(handle,hv, n, idx,indexBase, (int *) devVec->v_+vectorId*devVec->pitch_);
  return(i);
}

int iscatMultiVecDeviceIntVecIdx(void* deviceVec, int vectorId, int n, int first, void *deviceIdx,
				    int hfirst, void* host_values, int indexBase, int beta)
{  
  int i, *idx;
  struct MultiVectDevice *devIdx = (struct MultiVectDevice *) deviceIdx;
  i= iscatMultiVecDeviceInt(deviceVec, vectorId, n, first, 
			       (void*) devIdx->v_,  hfirst,host_values, indexBase, beta);
  return(i);
}

int iscatMultiVecDeviceInt(void* deviceVec, int vectorId, int n, int first, void *indexes,
			      int hfirst, void* host_values, int indexBase, int beta)
{ int i=0;
  int *hv  = (int *) host_values;
  int *idx=(int *) indexes;
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  spgpuHandle_t handle=psb_cudaGetHandle();

  idx = &(idx[first-indexBase]);
  hv  = &(hv[hfirst-indexBase]);
  spgpuIscat(handle, (int *) devVec->v_, n, hv, idx, indexBase, beta);
  return SPGPU_SUCCESS;
  
}
