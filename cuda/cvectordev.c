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
#include "cvectordev.h"

int registerMappedFloatComplex(void* buff, void** d_p, int n, cuFloatComplex dummy)
{
  return registerMappedMemory(buff, d_p, ((size_t) n) * sizeof(cuFloatComplex));
}

int checkMultiVecDeviceFloatComplex(void* deviceMultiVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceMultiVec;
  fprintf(stderr, "checkMultiVecDeviceInt Size: %d  Pointer %p\n", devVec->size_, devVec->v_);
  return SPGPU_SUCCESS;
}

int writeMultiVecDeviceFloatComplex(void* deviceVec, cuFloatComplex* hostVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  // Ex updateFromHost vector function
  int i = writeRemoteBuffer((void*) hostVec, (void*) devVec->v_, ((size_t) devVec->pitch_) * devVec->count_ * sizeof(cuFloatComplex));
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "FallocMultiVecDevice", i);

  return i;
}

int writeMultiVecDeviceFloatComplexR2(void* deviceVec, cuFloatComplex* hostVec, int ld)
{
  int i = writeMultiVecDeviceFloatComplex(deviceVec, (void *) hostVec);
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeMultiVecDeviceFloatComplexR2", i);

  return i;
}

int readMultiVecDeviceFloatComplex(void* deviceVec, cuFloatComplex* hostVec)
{
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  int i = readRemoteBuffer((void*) hostVec, (void*) devVec->v_, ((size_t) devVec->pitch_) * devVec->count_ * sizeof(cuFloatComplex));
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "readMultiVecDeviceFloatComplex", i);

  return i;
}

int readMultiVecDeviceFloatComplexR2(void* deviceVec, cuFloatComplex* hostVec, int ld)
{
  int i = readMultiVecDeviceFloatComplex(deviceVec, hostVec);
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "readMultiVecDeviceFloatComplexR2", i);

  return i;
}

int setscalMultiVecDeviceFloatComplex(cuFloatComplex val, int first, int last, int indexBase, void* devMultiVecX) 
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  spgpuCsetscal(handle, first, last, indexBase, val, (cuFloatComplex*) devVecX->v_);
  
  return SPGPU_SUCCESS;
}

int geinsMultiVecDeviceFloatComplex(int n, void* devMultiVecIrl, void* devMultiVecVal, int dupl, int indexBase, void* devMultiVecX)
{
  cuFloatComplex beta;
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecIrl = (struct MultiVectDevice*) devMultiVecIrl;
  struct MultiVectDevice* devVecVal = (struct MultiVectDevice*) devMultiVecVal;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecIrl->pitch_;

  if((n > devVecIrl->size_) || (n > devVecVal->size_)) 
    return SPGPU_UNSUPPORTED;

  //fprintf(stderr, "geins: %d %d  %p %p %p\n", dupl, n, devVecIrl->v_, devVecVal->v_, devVecX->v_);

  if(dupl == INS_OVERWRITE) 
    beta = make_cuFloatComplex(0.0, 0.0);
  else if(dupl == INS_ADD) 
    beta = make_cuFloatComplex(1.0, 0.0);
  else
    beta = make_cuFloatComplex(0.0, 0.0);

  spgpuCscat(handle, (cuFloatComplex*) devVecX->v_, n, (cuFloatComplex*) devVecVal->v_, (int*) devVecIrl->v_, indexBase, beta);
  
  return SPGPU_SUCCESS;
}

int igathMultiVecDeviceFloatComplex(void* deviceVec, int vectorId, int n, int first, void* indexes, int hfirst, void* host_values, int indexBase)
{
  int* idx = (int*) indexes;
  cuFloatComplex* hv = (cuFloatComplex*) host_values; 
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  hv  = &(hv[hfirst - indexBase]);
  idx = &(idx[first - indexBase]);
  spgpuCgath(handle, hv, n, idx, indexBase, (cuFloatComplex*) devVec->v_ + (vectorId * devVec->pitch_));
  return SPGPU_SUCCESS;
}

int igathMultiVecDeviceFloatComplexVecIdx(void* deviceVec, int vectorId, int n, int first, void* deviceIdx, int hfirst,  void* host_values, int indexBase)
{
  struct MultiVectDevice* devIdx = (struct MultiVectDevice*) deviceIdx;
  return igathMultiVecDeviceFloatComplex(deviceVec, vectorId, n, first, (void*) devIdx->v_, hfirst, host_values, indexBase);
}

int iscatMultiVecDeviceFloatComplex(void* deviceVec, int vectorId, int n, int first, void *indexes, int hfirst, void* host_values, int indexBase, cuFloatComplex beta)
{
  int* idx = (int*) indexes;
  cuFloatComplex* hv  = (cuFloatComplex*) host_values;
  struct MultiVectDevice* devVec = (struct MultiVectDevice*) deviceVec;
  spgpuHandle_t handle = psb_cudaGetHandle();

  idx = &(idx[first - indexBase]);
  hv  = &(hv[hfirst - indexBase]);
  spgpuCscat(handle, (cuFloatComplex*) devVec->v_, n, hv, idx, indexBase, beta);
  return SPGPU_SUCCESS;
}

int iscatMultiVecDeviceFloatComplexVecIdx(void* deviceVec, int vectorId, int n, int first, void *deviceIdx, int hfirst, void* host_values, int indexBase, cuFloatComplex beta)
{
  struct MultiVectDevice* devIdx = (struct MultiVectDevice*) deviceIdx;
  return  iscatMultiVecDeviceFloatComplex(deviceVec, vectorId, n, first, (void*) devIdx->v_, hfirst, host_values, indexBase, beta);
}

//AXPY-like operations
int axpbyMultiVecDeviceFloatComplex1(int n, cuFloatComplex alpha, void* devMultiVecX, cuFloatComplex beta, void* devMultiVecY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  if(devVecX->count_ != devVecY->count_)        //case v_full
    for(int j = 0; j < devVecY->count_; ++j)
      spgpuCaxpby(handle, (cuFloatComplex*) devVecY->v_ + (devVecY->pitch_ * j), n, beta, (cuFloatComplex*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (cuFloatComplex*) devVecX->v_);
  else                                          //case m_full or cudaVect
    for(int j = 0; j < devVecY->count_; ++j)
      spgpuCaxpby(handle, (cuFloatComplex*) devVecY->v_ + (devVecY->pitch_ * j), n, beta, (cuFloatComplex*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (cuFloatComplex*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void* devMultiVecX, int idxX, cuFloatComplex beta, void* devMultiVecY) 
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = ((cuFloatComplex*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuFloatComplex* yColAdd = (cuFloatComplex*) devVecY->v_;
  
  spgpuCaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);

  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceFloatComplex3(int n, cuFloatComplex alpha, void* devMultiVecX, cuFloatComplex beta, void* devMultiVecY, int idxY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
  
  cuFloatComplex* xColAdd = (cuFloatComplex*) devVecX->v_;
  cuFloatComplex* yColAdd = ((cuFloatComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  
  spgpuCaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceFloatComplex4(int n, cuFloatComplex alpha, void* devMultiVecX, int idxX, cuFloatComplex beta, void* devMultiVecY, int idxY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
    
  cuFloatComplex* xColAdd = ((cuFloatComplex*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuFloatComplex* yColAdd = ((cuFloatComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  
  spgpuCaxpby(handle, yColAdd, n, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyMultiVecDeviceFloatComplex5(int n, cuFloatComplex alpha, void* devMultiVecX, cuFloatComplex beta, void* devMultiVecY, void* devMultiVecZ)
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
    spgpuCaxpby(handle, (cuFloatComplex*) devVecZ->v_ + (devVecZ->pitch_ * j), n, beta, (cuFloatComplex*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (cuFloatComplex*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloatComplex1(int n, cuFloatComplex alpha, void* devMultiVecX, cuFloatComplex beta, void* devMultiVecY, cuFloatComplex gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = (cuFloatComplex*) devVecX->v_;
  cuFloatComplex* yColAdd = (cuFloatComplex*) devVecY->v_;
  cuFloatComplex* zColAdd = ((cuFloatComplex*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuCaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void* devMultiVecX, cuFloatComplex beta, void* devMultiVecY, int idxY, cuFloatComplex gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = (cuFloatComplex*) devVecX->v_;
  cuFloatComplex* yColAdd = ((cuFloatComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuFloatComplex* zColAdd = ((cuFloatComplex*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuCaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloatComplex3(int n, cuFloatComplex alpha, void* devMultiVecX, int idxX, cuFloatComplex beta, void* devMultiVecY, int idxY, cuFloatComplex gamma, void* devMultiVecZ, int idxZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = ((cuFloatComplex*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuFloatComplex* yColAdd = ((cuFloatComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuFloatComplex* zColAdd = ((cuFloatComplex*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  
  spgpuCaxpbycz(handle, zColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloatComplex4(int n, cuFloatComplex alpha, void* devMultiVecX, cuFloatComplex beta, void* devMultiVecY, cuFloatComplex gamma, void* devMultiVecZ)
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
    spgpuCaxpbycz(handle, (cuFloatComplex*) devVecZ->v_ + (devVecZ->pitch_ * j), n, gamma, (cuFloatComplex*) devVecZ->v_ + (devVecZ->pitch_ * j), beta, (cuFloatComplex*) devVecY->v_ + (devVecY->pitch_ * j), alpha, (cuFloatComplex*) devVecX->v_ + (devVecX->pitch_ * j));
  
  return SPGPU_SUCCESS;
}

int axpbyczMultiVecDeviceFloatComplex5(int n, cuFloatComplex alpha, void* devMultiVecX, int idxX, cuFloatComplex beta, void* devMultiVecY, int idxY, cuFloatComplex gamma, void* devMultiVecZ, int idxZ, void* devMultiVecW, int idxW)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  struct MultiVectDevice* devVecW = (struct MultiVectDevice*) devMultiVecW;
  spgpuHandle_t handle = psb_cudaGetHandle();
  
  if((n > devVecW->size_) || (n > devVecZ->size_) || (n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;
  
  cuFloatComplex* xColAdd = ((cuFloatComplex*) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuFloatComplex* yColAdd = ((cuFloatComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuFloatComplex* zColAdd = ((cuFloatComplex*) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);
  cuFloatComplex* wColAdd = ((cuFloatComplex*) devVecW->v_) + (IDX_F2C(idxW) * devVecW->pitch_);

  spgpuCaxpbycz(handle, wColAdd, n, gamma, zColAdd, beta, yColAdd, alpha, xColAdd);
  
  return SPGPU_SUCCESS;
}

int colspanDeviceFloatComplex(int n, void *devMultiVecX, CFI_cdesc_t *coeff, void *devMultiVecY, bool updFlag)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecY->size_) || (n > devVecX->size_)) 
    return SPGPU_UNSUPPORTED;
  
  spgpuCcolspan(handle, n, (cuFloatComplex*) devVecX->v_, devVecX->count_, devVecX->pitch_, coeff, (cuFloatComplex*) devVecY->v_, devVecY->count_, devVecY->pitch_, updFlag);

  return SPGPU_SUCCESS;
}

int upd_xyzMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, cuFloatComplex beta, cuFloatComplex gamma, cuFloatComplex delta, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecY->pitch_;
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuCupd_xyz(handle, n, alpha, beta, gamma, delta, (cuFloatComplex*) devVecX->v_, (cuFloatComplex*) devVecY->v_, (cuFloatComplex*) devVecZ->v_);
  return SPGPU_SUCCESS;
}

int xyzwMultiVecDeviceFloatComplex(int n, cuFloatComplex a, cuFloatComplex b, cuFloatComplex c, cuFloatComplex d, cuFloatComplex e, cuFloatComplex f, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ, void* devMultiVecW)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  struct MultiVectDevice* devVecW = (struct MultiVectDevice*) devMultiVecW;
  spgpuHandle_t handle = psb_cudaGetHandle();
  int pitch = devVecY->pitch_;
  if((n > devVecY->size_) || (n > devVecX->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuCxyzw(handle, n, a, b, c, d, e, f, (cuFloatComplex*) devVecX->v_, (cuFloatComplex*) devVecY->v_, (cuFloatComplex*) devVecZ->v_, (cuFloatComplex*) devVecW->v_);
  return SPGPU_SUCCESS;
}

//Scaling and norms
int scalMultiVecDeviceFloatComplex(cuFloatComplex alpha, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  // Note: inner kernel can handle aliased input/output
  spgpuCscal(handle, (cuFloatComplex*) devVecA->v_, devVecA->pitch_, alpha, (cuFloatComplex*) devVecA->v_);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceFloatComplexS(float* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuCmsnrm2(handle, y_res, n, (cuFloatComplex*) devVecA->v_);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceFloatComplexI(float* y_res, int n, void* devMultiVecA, int idx)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  cuFloatComplex* colAdd = ((cuFloatComplex*) devVecA->v_) + (IDX_F2C(idx) * devVecA->pitch_);
  spgpuCmsnrm2(handle, y_res, n, colAdd);
  return SPGPU_SUCCESS;
}

int nrm2MultiVecDeviceFloatComplexV(CFI_cdesc_t* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuCmvnrm2(handle, y_res, n, (cuFloatComplex*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int amaxMultiVecDeviceFloatComplex(float* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuCmamax(handle, y_res, n, (cuFloatComplex*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int asumMultiVecDeviceFloatComplex(float* y_res, int n, void* devMultiVecA)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;

  spgpuCmasum(handle, y_res, n, (cuFloatComplex*) devVecA->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

//DOT-like operations
int dotMultiVecDeviceFloatComplexS(cuFloatComplex* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuCmsdot(handle, y_res, n, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceFloatComplexV(cuFloatComplex* y_res, int n, int mA, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuCmvdot(handle, y_res, n, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceFloatComplexM(cuFloatComplex* y_res, int n, int mA, int mB, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuCmmdot(handle, y_res, n, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
  
  // int opMask = ((devVecA->count_ > 1) << 1) | (devVecB->count_ > 1);
  // switch(opMask)
  // {
  //   case 0: /*vv*/ *y_res = spgpuCdot(handle, n, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_);                                                                      break;
  //   case 1: /*vm*/ spgpuCmvdot(handle, y_res, n, (cuFloatComplex*) devVecB->v_, (cuFloatComplex*) devVecA->v_, devVecB->count_, devVecB->pitch_);                                    break;
  //   case 2: /*mv*/ spgpuCmvdot(handle, y_res, n, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_);                                    break;
  //   case 3: /*mm*/ spgpuCmmdot(handle, y_res, n, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);  break;
  // }
}

//DOT-like operations with CFI
int dotMultiVecDeviceFloatComplexV_CFI(CFI_cdesc_t* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuCmmdot_CFI(handle, y_res, n, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

int dotMultiVecDeviceFloatComplexM_CFI(CFI_cdesc_t* y_res, int n, void* devMultiVecA, void* devMultiVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) devMultiVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) devMultiVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  spgpuCmmdot_CFI(handle, y_res, n, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

//Elementwise operations
int mltMultiVecDeviceFloatComplex1(int n, cuFloatComplex alpha, void* devMultiVecX, void* devMultiVecY, cuFloatComplex beta, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_))
    return SPGPU_UNSUPPORTED;
  
  cuFloatComplex* xColAdd = (cuFloatComplex*) devVecX->v_;
  cuFloatComplex* yColAdd = (cuFloatComplex*) devVecY->v_;

  spgpuCmlt(handle, yColAdd, devVecY->count_, devVecY->pitch_, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void* devMultiVecX, void* devMultiVecY, cuFloatComplex beta, int idxY, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_))
    return SPGPU_UNSUPPORTED;
  
  cuFloatComplex* xColAdd = (cuFloatComplex*) devVecX->v_;
  cuFloatComplex* yColAdd = ((cuFloatComplex*) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);

  spgpuCmlt(handle, yColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 
 
int mltMultiVecDeviceFloatComplex3(int n, cuFloatComplex alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, cuFloatComplex beta, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_)) 
    return SPGPU_UNSUPPORTED;
  
  cuFloatComplex* xColAdd = ((cuFloatComplex *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuFloatComplex* yColAdd = ((cuFloatComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);

  spgpuCmlt(handle, yColAdd, 1, 0, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloatComplex4(int n, cuFloatComplex alpha, void* devMultiVecX, void* devMultiVecY, cuFloatComplex beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))  
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = (cuFloatComplex *) devVecX->v_;
  cuFloatComplex* yColAdd = (cuFloatComplex *) devVecY->v_;
  cuFloatComplex* zColAdd = (cuFloatComplex *) devVecZ->v_;

  spgpuCmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, yColAdd, devVecX->count_, devVecX->pitch_, xColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloatComplex5(int n, cuFloatComplex alpha, void* devMultiVecX, void* devMultiVecY, cuFloatComplex beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))  
    return SPGPU_UNSUPPORTED;
  
  cuFloatComplex* xColAdd = (cuFloatComplex *) devVecX->v_;
  cuFloatComplex* yColAdd = (cuFloatComplex *) devVecY->v_;
  cuFloatComplex* zColAdd = ((cuFloatComplex *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuCmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, devVecY->count_, devVecY->pitch_, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloatComplex6(int n, cuFloatComplex alpha, void* devMultiVecX, void* devMultiVecY, int idxY, cuFloatComplex beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = (cuFloatComplex *) devVecX->v_;
  cuFloatComplex* yColAdd = ((cuFloatComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuFloatComplex* zColAdd = ((cuFloatComplex *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuCmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloatComplex7(int n, cuFloatComplex alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, cuFloatComplex beta, void* devMultiVecZ, int idxZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_))
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = ((cuFloatComplex *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuFloatComplex* yColAdd = ((cuFloatComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuFloatComplex* zColAdd = ((cuFloatComplex *) devVecZ->v_) + (IDX_F2C(idxZ) * devVecZ->pitch_);

  spgpuCmlt(handle, zColAdd, 1, 0, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloatComplex8(int n, cuFloatComplex alpha, void* devMultiVecX, void* devMultiVecY, int idxY, cuFloatComplex beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_)) 
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = (cuFloatComplex *) devVecX->v_;
  cuFloatComplex* yColAdd = ((cuFloatComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuFloatComplex* zColAdd = (cuFloatComplex *) devVecZ->v_;

  spgpuCmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, xColAdd, devVecX->count_, devVecX->pitch_, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int mltMultiVecDeviceFloatComplex9(int n, cuFloatComplex alpha, void* devMultiVecX, int idxX, void* devMultiVecY, int idxY, cuFloatComplex beta, void* devMultiVecZ, bool conjgX, bool conjgY)
{
  struct MultiVectDevice* devVecX = (struct MultiVectDevice*) devMultiVecX;
  struct MultiVectDevice* devVecY = (struct MultiVectDevice*) devMultiVecY;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) devMultiVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecX->size_) || (n > devVecY->size_) || (n > devVecZ->size_)) 
    return SPGPU_UNSUPPORTED;

  cuFloatComplex* xColAdd = ((cuFloatComplex *) devVecX->v_) + (IDX_F2C(idxX) * devVecX->pitch_);
  cuFloatComplex* yColAdd = ((cuFloatComplex *) devVecY->v_) + (IDX_F2C(idxY) * devVecY->pitch_);
  cuFloatComplex* zColAdd = (cuFloatComplex *) devVecZ->v_;

  spgpuCmlt(handle, zColAdd, devVecZ->count_, devVecZ->pitch_, n, alpha, xColAdd, 1, 0, yColAdd, 1, 0, beta, conjgX, conjgY);

  return SPGPU_SUCCESS;
} 

int axyMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, void *deviceVecA, void *deviceVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if((n > devVecA->size_) || (n > devVecB->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuCmaxy(handle, (cuFloatComplex*) devVecB->v_, n, alpha, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_, devVecA->count_, devVecA->pitch_);

  return SPGPU_SUCCESS;
}

int axybzMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, void *deviceVecA, void *deviceVecB, cuFloatComplex beta, void *deviceVecZ)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  struct MultiVectDevice* devVecZ = (struct MultiVectDevice*) deviceVecZ;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if((n > devVecA->size_) || (n > devVecB->size_ ) || (n > devVecZ->size_ )) 
    return SPGPU_UNSUPPORTED;

  spgpuCmaxypbz(handle, (cuFloatComplex*) devVecZ->v_, n, beta, (cuFloatComplex*) devVecZ->v_, alpha, (cuFloatComplex*) devVecA->v_, (cuFloatComplex*) devVecB->v_, devVecB->count_, devVecB->pitch_);
  return SPGPU_SUCCESS;
}

int absMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void *deviceVecA, void *deviceVecB)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  struct MultiVectDevice* devVecB = (struct MultiVectDevice*) deviceVecB;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if((n > devVecA->size_) || (n > devVecB->size_ ))
    return SPGPU_UNSUPPORTED;

  spgpuCabs(handle, (cuFloatComplex*) devVecB->v_, n, alpha, (cuFloatComplex*) devVecA->v_);
  return SPGPU_SUCCESS;
}
 
int absMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, void *deviceVecA)
{
  struct MultiVectDevice* devVecA = (struct MultiVectDevice*) deviceVecA;
  spgpuHandle_t handle = psb_cudaGetHandle();
  if(n > devVecA->size_)
    return SPGPU_UNSUPPORTED;

  spgpuCabs(handle, (cuFloatComplex*)devVecA->v_, n, alpha, (cuFloatComplex*)devVecA->v_);
  return SPGPU_SUCCESS;
}