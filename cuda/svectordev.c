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
#include "svectordev.h"

int registerMappedFloat(void* buff, void** d_p, int n, float dummy)
{
  return registerMappedMemory(buff, d_p, ((size_t) n) * sizeof(float));
}

int checkMultiVecDeviceFloat(void* deviceMultiVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceMultiVec;
  fprintf(stderr, "checkMultiVecDeviceInt Size: %d  Pointer %p\n", devVec->size_, devVec->v_);
  return SPGPU_SUCCESS;
}

int writeMultiVecDeviceFloat(void* deviceVec, float* hostVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  // Ex updateFromHost vector function
  int i = writeRemoteBuffer((void*) hostVec, (void*) devVec->v_, ((size_t) devVec->pitch_) * devVec->count_ * sizeof(float));
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "FallocMultiVecDevice", i);

  return i;
}

int writeMultiVecDeviceFloatR2(void* deviceVec, float* hostVec, int ld)
{
  int i = writeMultiVecDeviceFloat(deviceVec, (void *) hostVec);
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeMultiVecDeviceFloatR2", i);

  return i;
}

int readMultiVecDeviceFloat(void* deviceVec, float* hostVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  int i = readRemoteBuffer((void*) hostVec, (void*) devVec->v_, ((size_t) devVec->pitch_) * devVec->count_ * sizeof(float));
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "readMultiVecDeviceFloat", i);

  return i;
}

int readMultiVecDeviceFloatR2(void* deviceVec, float* hostVec, int ld)
{
  int i = readMultiVecDeviceFloat(deviceVec, hostVec);
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "readMultiVecDeviceFloatR2", i);

  return i;
}

int setscalMultiVecDeviceFloat(float val, int first, int last, int indexBase, void* devMultiVecX) 
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  spgpuSsetscal(handle, first, last, indexBase, val, (float*) devVecX->v_);
  
  return SPGPU_SUCCESS;
}

int geinsMultiVecDeviceFloat(int n, void* devMultiVecIrl, void* devMultiVecVal, int dupl, int indexBase, void* devMultiVecX)
{
  float beta;
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecIrl = (struct MultiVectDevice*) devMultiVecIrl;
  struct MultiVectDevice* devVecVal = (struct MultiVectDevice*) devMultiVecVal;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecIrl->pitch_;

  if((n > devVecIrl->size_) || (n > devVecVal->size_)) 
    return SPGPU_UNSUPPORTED;

  //fprintf(stderr, "geins: %d %d  %p %p %p\n", dupl, n, devVecIrl->v_, devVecVal->v_, devVecX->v_);

  if(dupl == INS_OVERWRITE) 
    beta = 0.0;
  else if(dupl == INS_ADD) 
    beta = 1.0;
  else
    beta = 0.0;

  spgpuSscat(handle, (float*) devVecX->v_, n, (float*) devVecVal->v_, (int*) devVecIrl->v_, indexBase, beta);
  
  return SPGPU_SUCCESS;
}

int igathMultiVecDeviceFloat(void* deviceVec, int vectorId, int n, int first, void* indexes, int hfirst, void* host_values, int indexBase)
{
  int* idx = (int*) indexes;
  float* hv = (float*) host_values; 
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  hv  = &(hv[hfirst - indexBase]);
  idx = &(idx[first - indexBase]);
  spgpuSgath(handle, hv, n, idx, indexBase, (float*) devVec->v_ + (vectorId * devVec->pitch_));
  return SPGPU_SUCCESS;
}

int igathMultiVecDeviceFloatVecIdx(void* deviceVec, int vectorId, int n, int first, void* deviceIdx, int hfirst,  void* host_values, int indexBase)
{
  struct MultiVectDevice *devIdx = (struct MultiVectDevice*) deviceIdx;
  return igathMultiVecDeviceFloat(deviceVec, vectorId, n, first, (void*) devIdx->v_, hfirst, host_values, indexBase);
}

int iscatMultiVecDeviceFloat(void* deviceVec, int vectorId, int n, int first, void *indexes, int hfirst, void* host_values, int indexBase, float beta)
{
  int* idx = (int*) indexes;
  float* hv  = (float*) host_values;
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  spgpuHandle_t handle = psb_cudaGetHandle();

  idx = &(idx[first - indexBase]);
  hv  = &(hv[hfirst - indexBase]);
  spgpuSscat(handle, (float*) devVec->v_, n, hv, idx, indexBase, beta);
  return SPGPU_SUCCESS;
}

int iscatMultiVecDeviceFloatVecIdx(void* deviceVec, int vectorId, int n, int first, void *deviceIdx, int hfirst, void* host_values, int indexBase, float beta)
{
  struct MultiVectDevice *devIdx = (struct MultiVectDevice*) deviceIdx;
  return  iscatMultiVecDeviceFloat(deviceVec, vectorId, n, first, (void*) devIdx->v_, hfirst, host_values, indexBase, beta);
}

//AXPY-like operations
int axpbyMultiVecDeviceFloat1(int n, float alpha, void* devMultiVecX, float beta, void* devMultiVecY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  if(devVecX->count_ != devVecY->count_)        //case v_full
    for(int j = 0; j < devVecY->count_; ++j)
      spgpuSaxpby(handle, (float*) devVecY->v_ + (devVecY->pitch_ * j), n, beta, (float*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (float*) devVecX->v_);
  else                                          //case m_full or cudaVect
    for(int j = 0; j < devVecY->count_; ++j)
      spgpuSaxpby(handle, (float*) devVecY->v_ + (devVecY->pitch_ * j), n, beta, (float*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (float*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceFloat2(int n, float alpha, void* devMultiVecX, int idxX, float beta, void* devMultiVecY) 
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  float* xColAdd = ((float*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  float* yColAdd = (float*) devVecY->v_;
  
  spgpuSaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);

  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceFloat3(int n, float alpha, void* devMultiVecX, float beta, void* devMultiVecY, int idxY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
  
  float* xColAdd = (float*) devVecX->v_;
  float* yColAdd = ((float*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  
  spgpuSaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceFloat4(int n, float alpha, void* devMultiVecX, int idxX, float beta, void* devMultiVecY, int idxY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
    
  float* xColAdd = ((float*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  float* yColAdd = ((float*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  
  spgpuSaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceFloat5(int n, float alpha, void* devMultiVecX, float beta, void* devMultiVecY, void* devMultiVecZ)
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
    spgpuSaxpby(handle, (float*) devVecZ->v_ + (devVecZ->pitch_ * j), n, beta, (float*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (float*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloat1(int n, float alpha, void* devMultiVecX, float beta, void* devMultiVecY, float gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  float* xColAdd = (float*) devVecX->v_;
  float* yColAdd = (float*) devVecY->v_;
  float* zColAdd = ((float*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuSaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloat2(int n, float alpha, void* devMultiVecX, float beta, void* devMultiVecY, int idxY, float gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  float* xColAdd = (float*) devVecX->v_;
  float* yColAdd = ((float*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  float* zColAdd = ((float*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuSaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloat3(int n, float alpha, void* devMultiVecX, int idxX, float beta, void* devMultiVecY, int idxY, float gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  float* xColAdd = ((float*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  float* yColAdd = ((float*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  float* zColAdd = ((float*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuSaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloat4(int n, float alpha, void* devMultiVecX, float beta, void* devMultiVecY, float gamma, void* devMultiVecZ)
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
    spgpuSaxpbycz(handle, (float*) devVecZ->v_ + (devVecZ->pitch_ * j), n, gamma, (float*) devVecZ->v_ + (devVecZ->pitch_ * j), beta, (float*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (float*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloat5(int n, float alpha, void* devMultiVecX, int idxX, float beta, void* devMultiVecY, int idxY, float gamma, void* devMultiVecZ, int idxZ, void* devMultiVecW, int idxW)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  struct MultiVectDevice* devVecW = (struct MultiVectDevice*) devMultiVecW;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecW->size_) || (n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
  
  float* xColAdd = ((float*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  float* yColAdd = ((float*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  float* zColAdd = ((float*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  float* wColAdd = ((float*) devVecW->v_) + (IDX_F2C(idxW) * devVecW->pitch_);

  spgpuSaxpbycz(handle, wColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int colspanDeviceFloat(int n, void *devMultiVecX, CFI_cdesc_t *coeff, void *devMultiVecY, bool updFlag)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_)) 
    return SPGPU_UNSUPPORTED;
  
  spgpuScolspan(handle, n, (float*) devVecX->v_, devVecX->count_, devVecX->pitch_, coeff, (float*) devVecY->v_, devVecY->count_, devVecY->pitch_, updFlag);

  return SPGPU_SUCCESS;
}

int upd_xyzMultiVecDeviceFloat(int n, float alpha, float beta, float gamma, float delta, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecY->pitch_;
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuSupd_xyz(handle, n, alpha, beta, gamma, delta, (float*) devVecX->v_, (float*) devVecY->v_, (float*) devVecZ->v_);
  return SPGPU_SUCCESS;
}

int xyzwMultiVecDeviceFloat(int n, float a, float b, float c, float d, float e, float f, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ, void* devMultiVecW)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  struct MultiVectDevice* devVecW = (struct MultiVectDevice*) devMultiVecW;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecY->pitch_;
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuSxyzw(handle, n, a, b, c, d, e, f, (float*) devVecX->v_, (float*) devVecY->v_, (float*) devVecZ->v_, (float*) devVecW->v_);
  return SPGPU_SUCCESS;
}

//Scaling and norms
int scalMultiVecDeviceFloat(float alpha, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  // Note: inner kernel can handle aliased input/output
  spgpuSscal(handle, (float*) devVecA->v_, devVecA->pitch_, alpha, (float*) devVecA->v_);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceFloatS(float* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuSmsnrm2(handle, y_res, n, (float*) devVecA->v_);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceFloatI(float* y_res, int n, void* devMultiVecA, int idx)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  float* colAdd = ((float*) devVecA->v_) + (IDX_F2C(idx) * devVecA->pitch_);
  spgpuSmsnrm2(handle, y_res, n, colAdd);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceFloatV(CFI_cdesc_t* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuSmvnrm2(handle, y_res, n, (float*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int amaxMultiVecDeviceFloat(float* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuSmamax(handle, y_res, n, (float*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int asumMultiVecDeviceFloat(float* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuSmasum(handle, y_res, n, (float*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

//DOT-like operations
int dotMultiVecDeviceFloatS(float* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuSmsdot(handle, y_res, n, (float*) devVecA->v_, (float*) devVecB->v_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceFloatV(float* y_res, int n, int mA, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuSmvdot(handle, y_res, n, (float*) devVecA->v_, (float*) devVecB->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceFloatM(float* y_res, int n, int mA, int mB, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuSmmdot(handle, y_res, n, (float*) devVecA->v_, (float*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
  
  // int opMask = ((devVecA->count_ > 1) << 1) | (devVecB->count_ > 1);
  // switch(opMask)
  // {
  //   case 0: /*vv*/ *y_res = spgpuSdot(handle, n, (float*) devVecA->v_, (float*) devVecB->v_);                                                                      break;
  //   case 1: /*vm*/ spgpuSmvdot(handle, y_res, n, (float*) devVecB->v_, (float*) devVecA->v_, devVecB->count_, devVecB->pitch_);                                    break;
  //   case 2: /*mv*/ spgpuSmvdot(handle, y_res, n, (float*) devVecA->v_, (float*) devVecB->v_, devVecA->count_, devVecA->pitch_);                                    break;
  //   case 3: /*mm*/ spgpuSmmdot(handle, y_res, n, (float*) devVecA->v_, (float*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);  break;
  // }
}

//DOT-like operations with CFI
int dotMultiVecDeviceFloatV_CFI(CFI_cdesc_t* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuSmmdot_CFI(handle, y_res, n, (float*) devVecA->v_, (float*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceFloatM_CFI(CFI_cdesc_t* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuSmmdot_CFI(handle, y_res, n, (float*) devVecA->v_, (float*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

//Elementwise operations
int mltMultiVecDeviceFloat1(int n, float alpha, void* devMultiVecX, void* devMultiVecY, float beta, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_))
    return SPGPU_UNSUPPORTED;
  
  float* xColAdd = (float*) devVecX->v_;
  float* yColAdd = (float*) devVecY->v_;

  spgpuSmlt(handle, yColAdd, devVecY->count_, devVecY->pitch_, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloat2(int n, float alpha, void* devMultiVecX, void* devMultiVecY, float beta, int idxY, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_))
    return SPGPU_UNSUPPORTED;
  
  float* xColAdd = (float*) devVecX->v_;
  float* yColAdd = ((float*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);

  spgpuSmlt(handle, yColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 
 
int mltMultiVecDeviceFloat3(int n, float alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, float beta, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_)) 
    return SPGPU_UNSUPPORTED;
  
  float* xColAdd = ((float *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  float* yColAdd = ((float *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);

  spgpuSmlt(handle, yColAdd, 1, 0, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloat4(int n, float alpha, void* devMultiVecX, void* devMultiVecY, float beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))  
    return SPGPU_UNSUPPORTED;

  float* xColAdd = (float *) devVecX->v_;
  float* yColAdd = (float *) devVecY->v_;
  float* zColAdd = (float *) devVecZ->v_;

  spgpuSmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, yColAdd, devVecX->count_, devVecX->pitch_, xColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloat5(int n, float alpha, void* devMultiVecX, void* devMultiVecY, float beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))  
    return SPGPU_UNSUPPORTED;
  
  float* xColAdd = (float *) devVecX->v_;
  float* yColAdd = (float *) devVecY->v_;
  float* zColAdd = ((float *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuSmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloat6(int n, float alpha, void* devMultiVecX, void* devMultiVecY, int idxY, float beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))
    return SPGPU_UNSUPPORTED;

  float* xColAdd = (float *) devVecX->v_;
  float* yColAdd = ((float *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  float* zColAdd = ((float *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuSmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloat7(int n, float alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, float beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))
    return SPGPU_UNSUPPORTED;

  float* xColAdd = ((float *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  float* yColAdd = ((float *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  float* zColAdd = ((float *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuSmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloat8(int n, float alpha, void* devMultiVecX, void* devMultiVecY, int idxY, float beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_)) 
    return SPGPU_UNSUPPORTED;

  float* xColAdd = (float *) devVecX->v_;
  float* yColAdd = ((float *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  float* zColAdd = (float *) devVecZ->v_;

  spgpuSmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloat9(int n, float alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, float beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_)) 
    return SPGPU_UNSUPPORTED;

  float* xColAdd = ((float *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  float* yColAdd = ((float *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  float* zColAdd = (float *) devVecZ->v_;

  spgpuSmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
}

int axyMultiVecDeviceFloat(int n, float alpha, void *deviceVecA, void *deviceVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if((n > devVecA->size_) || (n > devVecB->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuSmaxy(handle, (float*) devVecB->v_, n, alpha, (float*) devVecA->v_, (float*) devVecB->v_, devVecA->count_, devVecA->pitch_);

  return SPGPU_SUCCESS;
}

int axybzMultiVecDeviceFloat(int n, float alpha, void *deviceVecA, void *deviceVecB, float beta, void *deviceVecZ)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) deviceVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if((n > devVecA->size_) || (n > devVecB->size_ ) || (n > devVecZ->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuSmaxypbz(handle, (float*) devVecZ->v_, n, beta, (float*) devVecZ->v_, alpha, (float*) devVecA->v_, (float*) devVecB->v_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

int absMultiVecDeviceFloat2(int n, float alpha, void *deviceVecA, void *deviceVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecA->size_) || (n > devVecB->size_ ))
    return SPGPU_UNSUPPORTED;

  spgpuSabs(handle, (float*) devVecB->v_, n, alpha, (float*) devVecA->v_);
  return SPGPU_SUCCESS;
}
 
int absMultiVecDeviceFloat(int n, float alpha, void *deviceVecA)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if(n > devVecA->size_)
    return SPGPU_UNSUPPORTED;

  spgpuSabs(handle, (float*)devVecA->v_, n, alpha, (float*)devVecA->v_);
  return SPGPU_SUCCESS;
}