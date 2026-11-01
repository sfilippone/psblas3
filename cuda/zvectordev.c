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
/*      software without specific prior written permission. */

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

int registerMappedDoubleComplex(void* buff, void** d_p, int n, cuDoubleComplex dummy)
{
  return registerMappedMemory(buff, d_p, ((size_t) n) * sizeof(cuDoubleComplex));
}

int checkMultiVecDeviceDoubleComplex(void* deviceMultiVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceMultiVec;
  fprintf(stderr, "checkMultiVecDeviceInt Size: %d  Pointer %p\n", devVec->size_, devVec->v_);
  return SPGPU_SUCCESS;
}

int writeMultiVecDeviceDoubleComplex(void* deviceVec, cuDoubleComplex* hostVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  // Ex updateFromHost vector function
  int i = writeRemoteBuffer((void*) hostVec, (void*) devVec->v_, ((size_t) devVec->pitch_) * devVec->count_ * sizeof(cuDoubleComplex));
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "FallocMultiVecDevice", i);

  return i;
}

int writeMultiVecDeviceDoubleComplexR2(void* deviceVec, cuDoubleComplex* hostVec, int ld)
{
  int i = writeMultiVecDeviceDoubleComplex(deviceVec, (void *) hostVec);
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeMultiVecDeviceDoubleComplexR2", i);

  return i;
}

int readMultiVecDeviceDoubleComplex(void* deviceVec, cuDoubleComplex* hostVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  int i = readRemoteBuffer((void*) hostVec, (void*) devVec->v_, ((size_t) devVec->pitch_) * devVec->count_ * sizeof(cuDoubleComplex));
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "readMultiVecDeviceDoubleComplex", i);

  return i;
}

int readMultiVecDeviceDoubleComplexR2(void* deviceVec, cuDoubleComplex* hostVec, int ld)
{
  int i = readMultiVecDeviceDoubleComplex(deviceVec, hostVec);
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "readMultiVecDeviceDoubleComplexR2", i);

  return i;
}

int setscalMultiVecDeviceDoubleComplex(cuDoubleComplex val, int first, int last, int indexBase, void* devMultiVecX) 
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  spgpuZsetscal(handle, first, last, indexBase, val, (cuDoubleComplex*) devVecX->v_);
  
  return SPGPU_SUCCESS;
}

int geinsMultiVecDeviceDoubleComplex(int n, void* devMultiVecIrl, void* devMultiVecVal, int dupl, int indexBase, void* devMultiVecX)
{
  cuDoubleComplex beta;
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecIrl = (struct MultiVectDevice*) devMultiVecIrl;
  struct MultiVectDevice* devVecVal = (struct MultiVectDevice*) devMultiVecVal;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecIrl->pitch_;

  if((n > devVecIrl->size_) || (n > devVecVal->size_)) 
    return SPGPU_UNSUPPORTED;

  //fprintf(stderr, "geins: %d %d  %p %p %p\n", dupl, n, devVecIrl->v_, devVecVal->v_, devVecX->v_);

  if(dupl == INS_OVERWRITE) 
    beta = make_cuDoubleComplex(0.0, 0.0);
  else if(dupl == INS_ADD) 
    beta = make_cuDoubleComplex(1.0, 0.0);
  else
    beta = make_cuDoubleComplex(0.0, 0.0);

  spgpuZscat(handle, (cuDoubleComplex*) devVecX->v_, n, (cuDoubleComplex*) devVecVal->v_, (int*) devVecIrl->v_, indexBase, beta);
  
  return SPGPU_SUCCESS;
}

int igathMultiVecDeviceDoubleComplex(void* deviceVec, int vectorId, int n, int first, void* indexes, int hfirst, void* host_values, int indexBase)
{
  int* idx = (int*) indexes;
  cuDoubleComplex* hv = (cuDoubleComplex*) host_values; 
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  hv  = &(hv[hfirst - indexBase]);
  idx = &(idx[first - indexBase]);
  spgpuZgath(handle, hv, n, idx, indexBase, (cuDoubleComplex*) devVec->v_ + (vectorId * devVec->pitch_));
  return SPGPU_SUCCESS;
}

int igathMultiVecDeviceDoubleComplexVecIdx(void* deviceVec, int vectorId, int n, int first, void* deviceIdx, int hfirst,  void* host_values, int indexBase)
{
  struct MultiVectDevice* devIdx = (struct MultiVectDevice*) deviceIdx;
  return igathMultiVecDeviceDoubleComplex(deviceVec, vectorId, n, first, (void*) devIdx->v_, hfirst, host_values, indexBase);
}

int iscatMultiVecDeviceDoubleComplex(void* deviceVec, int vectorId, int n, int first, void *indexes, int hfirst, void* host_values, int indexBase, cuDoubleComplex beta)
{
  int* idx = (int*) indexes;
  cuDoubleComplex* hv  = (cuDoubleComplex*) host_values;
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  spgpuHandle_t handle = psb_cudaGetHandle();

  idx = &(idx[first - indexBase]);
  hv  = &(hv[hfirst - indexBase]);
  spgpuZscat(handle, (cuDoubleComplex*) devVec->v_, n, hv, idx, indexBase, beta);
  return SPGPU_SUCCESS;
}

int iscatMultiVecDeviceDoubleComplexVecIdx(void* deviceVec, int vectorId, int n, int first, void *deviceIdx, int hfirst, void* host_values, int indexBase, cuDoubleComplex beta)
{
  struct MultiVectDevice* devIdx = (struct MultiVectDevice*) deviceIdx;
  return  iscatMultiVecDeviceDoubleComplex(deviceVec, vectorId, n, first, (void*) devIdx->v_, hfirst, host_values, indexBase, beta);
}

//AXPY-like operations
int axpbyMultiVecDeviceDoubleComplex1(int n, cuDoubleComplex alpha, void* devMultiVecX, cuDoubleComplex beta, void* devMultiVecY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  if(devVecX->count_ != devVecY->count_)        //case v_full
    for(int j = 0; j < devVecY->count_; ++j)
      spgpuZaxpby(handle, (cuDoubleComplex*) devVecY->v_ + (devVecY->pitch_ * j), n, beta, (cuDoubleComplex*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (cuDoubleComplex*) devVecX->v_);
  else                                          //case m_full or cudaVect
    for(int j = 0; j < devVecY->count_; ++j)
      spgpuZaxpby(handle, (cuDoubleComplex*) devVecY->v_ + (devVecY->pitch_ * j), n, beta, (cuDoubleComplex*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (cuDoubleComplex*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void* devMultiVecX, int idxX, cuDoubleComplex beta, void* devMultiVecY) 
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = ((cuDoubleComplex*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuDoubleComplex* yColAdd = (cuDoubleComplex*) devVecY->v_;
  
  spgpuZaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);

  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceDoubleComplex3(int n, cuDoubleComplex alpha, void* devMultiVecX, cuDoubleComplex beta, void* devMultiVecY, int idxY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
  
  cuDoubleComplex* xColAdd = (cuDoubleComplex*) devVecX->v_;
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  
  spgpuZaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceDoubleComplex4(int n, cuDoubleComplex alpha, void* devMultiVecX, int idxX, cuDoubleComplex beta, void* devMultiVecY, int idxY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
    
  cuDoubleComplex* xColAdd = ((cuDoubleComplex*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  
  spgpuZaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceDoubleComplex5(int n, cuDoubleComplex alpha, void* devMultiVecX, cuDoubleComplex beta, void* devMultiVecY, void* devMultiVecZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
  
  if((devVecX->count_ != devVecZ->count_) || (devVecX->count_ != devVecZ->count_))
    return SPGPU_UNSUPPORTED;
  
  for(int j = 0; j < devVecZ->count_; ++j)
    spgpuZaxpby(handle, (cuDoubleComplex*) devVecZ->v_ + (devVecZ->pitch_ * j), n, beta, (cuDoubleComplex*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (cuDoubleComplex*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceDoubleComplex1(int n, cuDoubleComplex alpha, void* devMultiVecX, cuDoubleComplex beta, void* devMultiVecY, cuDoubleComplex gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = (cuDoubleComplex*) devVecX->v_;
  cuDoubleComplex* yColAdd = (cuDoubleComplex*) devVecY->v_;
  cuDoubleComplex* zColAdd = ((cuDoubleComplex*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuZaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void* devMultiVecX, cuDoubleComplex beta, void* devMultiVecY, int idxY, cuDoubleComplex gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = (cuDoubleComplex*) devVecX->v_;
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuDoubleComplex* zColAdd = ((cuDoubleComplex*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuZaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceDoubleComplex3(int n, cuDoubleComplex alpha, void* devMultiVecX, int idxX, cuDoubleComplex beta, void* devMultiVecY, int idxY, cuDoubleComplex gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = ((cuDoubleComplex*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuDoubleComplex* zColAdd = ((cuDoubleComplex*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuZaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceDoubleComplex4(int n, cuDoubleComplex alpha, void* devMultiVecX, cuDoubleComplex beta, void* devMultiVecY, cuDoubleComplex gamma, void* devMultiVecZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
  
  if((devVecX->count_ != devVecZ->count_) || (devVecX->count_ != devVecZ->count_))
    return SPGPU_UNSUPPORTED;

  for(int j = 0; j < devVecZ->count_; ++j)
    spgpuZaxpbycz(handle, (cuDoubleComplex*) devVecZ->v_ + (devVecZ->pitch_ * j), n, gamma, (cuDoubleComplex*) devVecZ->v_ + (devVecZ->pitch_ * j), beta, (cuDoubleComplex*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (cuDoubleComplex*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceDoubleComplex5(int n, cuDoubleComplex alpha, void* devMultiVecX, int idxX, cuDoubleComplex beta, void* devMultiVecY, int idxY, cuDoubleComplex gamma, void* devMultiVecZ, int idxZ, void* devMultiVecW, int idxW)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  struct MultiVectDevice* devVecW = (struct MultiVectDevice*) devMultiVecW;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecW->size_) || (n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
  
  cuDoubleComplex* xColAdd = ((cuDoubleComplex*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuDoubleComplex* zColAdd = ((cuDoubleComplex*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  cuDoubleComplex* wColAdd = ((cuDoubleComplex*) devVecW->v_) + (IDX_F2C(idxW) * devVecW->pitch_);

  spgpuZaxpbycz(handle, wColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int colspanDeviceDoubleComplex(int n, void *devMultiVecX, CFI_cdesc_t *coeff, void *devMultiVecY, bool updFlag)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_)) 
    return SPGPU_UNSUPPORTED;
  
  spgpuZcolspan(handle, n, (cuDoubleComplex*) devVecX->v_, devVecX->count_, devVecX->pitch_, coeff, (cuDoubleComplex*) devVecY->v_, devVecY->count_, devVecY->pitch_, updFlag);

  return SPGPU_SUCCESS;
}

int upd_xyzMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, cuDoubleComplex beta, cuDoubleComplex gamma, cuDoubleComplex delta, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecY->pitch_;
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuZupd_xyz(handle, n, alpha, beta, gamma, delta, (cuDoubleComplex*) devVecX->v_, (cuDoubleComplex*) devVecY->v_, (cuDoubleComplex*) devVecZ->v_);
  return SPGPU_SUCCESS;
}

int xyzwMultiVecDeviceDoubleComplex(int n, cuDoubleComplex a, cuDoubleComplex b, cuDoubleComplex c, cuDoubleComplex d, cuDoubleComplex e, cuDoubleComplex f, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ, void* devMultiVecW)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  struct MultiVectDevice* devVecW = (struct MultiVectDevice*) devMultiVecW;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecY->pitch_;
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuZxyzw(handle, n, a, b, c, d, e, f, (cuDoubleComplex*) devVecX->v_, (cuDoubleComplex*) devVecY->v_, (cuDoubleComplex*) devVecZ->v_, (cuDoubleComplex*) devVecW->v_);
  return SPGPU_SUCCESS;
}

//Scaling and norms
int scalMultiVecDeviceDoubleComplex(cuDoubleComplex alpha, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  // Note: inner kernel can handle aliased input/output
  spgpuZscal(handle, (cuDoubleComplex*) devVecA->v_, devVecA->pitch_, alpha, (cuDoubleComplex*) devVecA->v_);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceDoubleComplexS(double* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuZmsnrm2(handle, y_res, n, (cuDoubleComplex*) devVecA->v_);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceDoubleComplexI(double* y_res, int n, void* devMultiVecA, int idx)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  cuDoubleComplex* colAdd = ((cuDoubleComplex*) devVecA->v_) + (IDX_F2C(idx) * devVecA->pitch_);
  spgpuZmsnrm2(handle, y_res, n, colAdd);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceDoubleComplexV(CFI_cdesc_t* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuZmvnrm2(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int amaxMultiVecDeviceDoubleComplex(double* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuZmamax(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int asumMultiVecDeviceDoubleComplex(double* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuZmasum(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

//DOT-like operations
int dotMultiVecDeviceDoubleComplexS(cuDoubleComplex* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuZmsdot(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceDoubleComplexV(cuDoubleComplex* y_res, int n, int mA, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuZmvdot(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceDoubleComplexM(cuDoubleComplex* y_res, int n, int mA, int mB, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuZmmdot(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
  
  // int opMask = ((devVecA->count_ > 1) << 1) | (devVecB->count_ > 1);
  // switch(opMask)
  // {
  //   case 0: /*vv*/ *y_res = spgpuZdot(handle, n, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_);                                                                      break;
  //   case 1: /*vm*/ spgpuZmvdot(handle, y_res, n, (cuDoubleComplex*) devVecB->v_, (cuDoubleComplex*) devVecA->v_, devVecB->count_, devVecB->pitch_);                                    break;
  //   case 2: /*mv*/ spgpuZmvdot(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_);                                    break;
  //   case 3: /*mm*/ spgpuZmmdot(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);  break;
  // }
}

//DOT-like operations with CFI
int dotMultiVecDeviceDoubleComplexV_CFI(CFI_cdesc_t* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuZmmdot_CFI(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceDoubleComplexM_CFI(CFI_cdesc_t* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuZmmdot_CFI(handle, y_res, n, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

//Elementwise operations
int mltMultiVecDeviceDoubleComplex1(int n, cuDoubleComplex alpha, void* devMultiVecX, void* devMultiVecY, cuDoubleComplex beta, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_))
    return SPGPU_UNSUPPORTED;
  
  cuDoubleComplex* xColAdd = (cuDoubleComplex*) devVecX->v_;
  cuDoubleComplex* yColAdd = (cuDoubleComplex*) devVecY->v_;

  spgpuZmlt(handle, yColAdd, devVecY->count_, devVecY->pitch_, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void* devMultiVecX, void* devMultiVecY, cuDoubleComplex beta, int idxY, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_))
    return SPGPU_UNSUPPORTED;
  
  cuDoubleComplex* xColAdd = (cuDoubleComplex*) devVecX->v_;
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);

  spgpuZmlt(handle, yColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 
 
int mltMultiVecDeviceDoubleComplex3(int n, cuDoubleComplex alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, cuDoubleComplex beta, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_)) 
    return SPGPU_UNSUPPORTED;
  
  cuDoubleComplex* xColAdd = ((cuDoubleComplex *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuDoubleComplex* yColAdd = ((cuDoubleComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);

  spgpuZmlt(handle, yColAdd, 1, 0, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceDoubleComplex4(int n, cuDoubleComplex alpha, void* devMultiVecX, void* devMultiVecY, cuDoubleComplex beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))  
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = (cuDoubleComplex *) devVecX->v_;
  cuDoubleComplex* yColAdd = (cuDoubleComplex *) devVecY->v_;
  cuDoubleComplex* zColAdd = (cuDoubleComplex *) devVecZ->v_;

  spgpuZmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, yColAdd, devVecX->count_, devVecX->pitch_, xColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceDoubleComplex5(int n, cuDoubleComplex alpha, void* devMultiVecX, void* devMultiVecY, cuDoubleComplex beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))  
    return SPGPU_UNSUPPORTED;
  
  cuDoubleComplex* xColAdd = (cuDoubleComplex *) devVecX->v_;
  cuDoubleComplex* yColAdd = (cuDoubleComplex *) devVecY->v_;
  cuDoubleComplex* zColAdd = ((cuDoubleComplex *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuZmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceDoubleComplex6(int n, cuDoubleComplex alpha, void* devMultiVecX, void* devMultiVecY, int idxY, cuDoubleComplex beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = (cuDoubleComplex *) devVecX->v_;
  cuDoubleComplex* yColAdd = ((cuDoubleComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuDoubleComplex* zColAdd = ((cuDoubleComplex *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuZmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceDoubleComplex7(int n, cuDoubleComplex alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, cuDoubleComplex beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = ((cuDoubleComplex *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuDoubleComplex* yColAdd = ((cuDoubleComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuDoubleComplex* zColAdd = ((cuDoubleComplex *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuZmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceDoubleComplex8(int n, cuDoubleComplex alpha, void* devMultiVecX, void* devMultiVecY, int idxY, cuDoubleComplex beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_)) 
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = (cuDoubleComplex *) devVecX->v_;
  cuDoubleComplex* yColAdd = ((cuDoubleComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuDoubleComplex* zColAdd = (cuDoubleComplex *) devVecZ->v_;

  spgpuZmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceDoubleComplex9(int n, cuDoubleComplex alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, cuDoubleComplex beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_)) 
    return SPGPU_UNSUPPORTED;

  cuDoubleComplex* xColAdd = ((cuDoubleComplex *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuDoubleComplex* yColAdd = ((cuDoubleComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuDoubleComplex* zColAdd = (cuDoubleComplex *) devVecZ->v_;

  spgpuZmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int axyMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA, void *deviceVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if((n > devVecA->size_) || (n > devVecB->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuZmaxy(handle, (cuDoubleComplex*) devVecB->v_, n, alpha, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_);

  return SPGPU_SUCCESS;
}

int axybzMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA, void *deviceVecB, cuDoubleComplex beta, void *deviceVecZ)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) deviceVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if((n > devVecA->size_) || (n > devVecB->size_ ) || (n > devVecZ->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuZmaxypbz(handle, (cuDoubleComplex*) devVecZ->v_, n, beta, (cuDoubleComplex*) devVecZ->v_, alpha, (cuDoubleComplex*) devVecA->v_, (cuDoubleComplex*) devVecB->v_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

int absMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void *deviceVecA, void *deviceVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecA->size_) || (n > devVecB->size_ ))
    return SPGPU_UNSUPPORTED;

  spgpuZabs(handle, (cuDoubleComplex*) devVecB->v_, n, alpha, (cuDoubleComplex*) devVecA->v_);
  return SPGPU_SUCCESS;
}
 
int absMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if(n > devVecA->size_)
    return SPGPU_UNSUPPORTED;

  spgpuZabs(handle, (cuDoubleComplex*)devVecA->v_, n, alpha, (cuDoubleComplex*)devVecA->v_);
  return SPGPU_SUCCESS;
}