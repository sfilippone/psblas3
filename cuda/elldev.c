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
 
#include <sys/time.h>
#include "elldev.h"

#define PASS_RS  0

EllDeviceParams getEllDeviceParams(unsigned int rows, unsigned int maxRowSize, unsigned int nnzeros, unsigned int columns, unsigned int elementType, unsigned int firstIndex)
{
  EllDeviceParams params;

  if (elementType == SPGPU_TYPE_DOUBLE)
    params.pitch = ((rows + ELL_PITCH_ALIGN_D - 1) / ELL_PITCH_ALIGN_D) * ELL_PITCH_ALIGN_D;
  else
    params.pitch = ((rows + ELL_PITCH_ALIGN_S - 1) / ELL_PITCH_ALIGN_S) * ELL_PITCH_ALIGN_S;
  
  //For complex?
  params.elementType = elementType;
  params.rows = rows;
  params.maxRowSize = maxRowSize;
  params.avgRowSize = (nnzeros+rows-1)/rows;
  params.columns = columns;
  params.firstIndex = firstIndex;

  //params.pitch = computeEllAllocPitch(rows);
  return params;
}

int allocEllDevice(void ** remoteMatrix, EllDeviceParams* params)
{
  struct EllDevice* tmp = (struct EllDevice*) malloc(sizeof(struct EllDevice));
  *remoteMatrix = (void*) tmp;
  tmp->rows = params->rows;
  tmp->cMPitch = computeEllAllocPitch(tmp->rows);
  tmp->rPPitch = tmp->cMPitch;
  tmp->pitch = tmp->cMPitch;
  tmp->maxRowSize = params->maxRowSize;
  tmp->avgRowSize = params->avgRowSize;
  tmp->allocsize = (int)tmp->maxRowSize * tmp->pitch;
  //tmp->allocsize = (int)params->maxRowSize * tmp->cMPitch;
  allocRemoteBuffer((void**) &(tmp->rS), ((size_t) tmp->rows) * sizeof(int));
  allocRemoteBuffer((void**) &(tmp->diag), ((size_t) tmp->rows) * sizeof(int));
  allocRemoteBuffer((void**) &(tmp->rP), ((size_t) tmp->allocsize) * sizeof(int));
  tmp->columns = params->columns;
  tmp->baseIndex = params->firstIndex;
  tmp->dataType = params->elementType;
  //fprintf(stderr, "allocEllDevice: %d %d %d \n", tmp->pitch, params->maxRowSize, params->avgRowSize);

  switch(params->elementType)
  {
    case SPGPU_TYPE_FLOAT:          allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(float));            break;
    case SPGPU_TYPE_DOUBLE:         allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(double));           break;
    case SPGPU_TYPE_COMPLEX_FLOAT:  allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(cuFloatComplex));   break;
    case SPGPU_TYPE_COMPLEX_DOUBLE: allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(cuDoubleComplex));  break;
    default: return SPGPU_UNSUPPORTED; // Unsupported params
  }
  
  //fprintf(stderr, "From allocEllDevice: %d %d %d %p %p %p\n", tmp->maxRowSize, 
  //	  tmp->avgRowSize, tmp->allocsize, tmp->rS, tmp->rP, tmp->cM);

  return SPGPU_SUCCESS;
}

void zeroEllDevice(void *remoteMatrix)
{
  struct EllDevice* tmp = (struct EllDevice*) remoteMatrix;
  
  switch(tmp->dataType)
  {
    case SPGPU_TYPE_FLOAT:          cudaMemset((void*) tmp->cM, 0, ((size_t) tmp->allocsize) * sizeof(float));            break;
    case SPGPU_TYPE_DOUBLE:         cudaMemset((void*) tmp->cM, 0, ((size_t) tmp->allocsize) * sizeof(double));           break;
    case SPGPU_TYPE_COMPLEX_FLOAT:  cudaMemset((void*) tmp->cM, 0, ((size_t) tmp->allocsize) * sizeof(cuFloatComplex));   break;
    case SPGPU_TYPE_COMPLEX_DOUBLE: cudaMemset((void*) tmp->cM, 0, ((size_t) tmp->allocsize) * sizeof(cuDoubleComplex));  break;
    default: return; // Unsupported params
  }
  return;
}

void freeEllDevice(void* remoteMatrix)
{
  struct EllDevice* devMat = (struct EllDevice*) remoteMatrix;  
  //fprintf(stderr, "freeEllDevice\n");
  if (devMat != NULL) 
  {
    freeRemoteBuffer(devMat->rS);
    freeRemoteBuffer(devMat->rP);
    freeRemoteBuffer(devMat->cM);
    free(remoteMatrix);
  }
}

int FallocEllDevice(void** deviceMat, unsigned int rows, unsigned int maxRowSize, unsigned int nnzeros, unsigned int columns, unsigned int elementType, unsigned int firstIndex)
{
  EllDeviceParams p = getEllDeviceParams(rows, maxRowSize, nnzeros, columns, elementType, firstIndex);
  int i = allocEllDevice(deviceMat, &p);
  if (i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "FallocEllDevice", i);
  
  return(i);
}


void sspmdmm_gpu(float *z, int s, int vPitch, float *y, float alpha, float* cM, int* rP, int* rS, int avgRowSize, int maxRowSize, int rows, int pitch, float *x, float beta, int firstIndex)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  int* m_rS = (PASS_RS) ? rS : NULL;

  for(int i = 0; i < s; ++i)
  {
	  spgpuSellspmv (handle, (float*) z, (float*) y, alpha, (float*) cM, rP, pitch, pitch, m_rS, NULL, avgRowSize, maxRowSize, rows, (float*) x, beta, firstIndex);
    
    z += vPitch;
    y += vPitch;
    x += vPitch;		
  }
}

int spmvEllDeviceFloatF(void *deviceMat, float alpha, void* deviceX, float beta, void* deviceY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  /*spgpuSellspmv (handle, (double*) y->v_, (double*)y->v_, alpha, (double*) devMat->cM, devMat->rP, devMat->cMPitch, devMat->rPPitch, devMat->rS, devMat->rows, (double*)x->v_, beta, devMat->baseIndex);*/
  /* fprintf(stderr, "From spmvEllDouble: mat %d %d %d %d y %d %d \n", */
  /* 	  devMat->avgRowSize, devMat->maxRowSize, devMat->rows, */
  /* 	  devMat->pitch, y->count_, y->pitch_); */

  sspmdmm_gpu((float*) y->v_, y->count_, y->pitch_, (float*) y->v_, alpha, (float*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, (float*) x->v_, beta, devMat->baseIndex);
  
  return SPGPU_SUCCESS;
}

int spmvEllDeviceFloatMV(void *deviceMat, float alpha, void* deviceX, float beta, void* deviceY, int idxY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  float* xColAdd = (float*) x->v_;
  float* yColAdd = ((float*) y->v_) + (IDX_F2C(idxY) * y->pitch_);

  sspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, alpha, (float*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, beta, devMat->baseIndex);

  return SPGPU_SUCCESS;
}

int spmvEllDeviceFloatVM(void *deviceMat, float alpha, void* deviceX, int idxX, float beta, void* deviceY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  float* xColAdd = ((float*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  float* yColAdd = (float*) y->v_;

  sspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, alpha, (float*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, beta, devMat->baseIndex);

  return SPGPU_SUCCESS;
}

int spmvEllDeviceFloatMMi(void *deviceMat, float alpha, void* deviceX, int idxX, float beta, void* deviceY, int idxY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  float* xColAdd = ((float*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  float* yColAdd = ((float*) y->v_) + (IDX_F2C(idxY) * y->pitch_);

  sspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, alpha, (float*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, beta, devMat->baseIndex);

  return SPGPU_SUCCESS;
}


void dspmdmm_gpu(double *z, int s, int vPitch, double *y, double alpha, double* cM, int* rP, int* rS, int avgRowSize, int maxRowSize, int rows, int pitch, double *x, double beta, int firstIndex)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  int* m_rS = (PASS_RS) ? rS : NULL;

  for(int i = 0; i < s; ++i)
  {
    spgpuDellspmv(handle, (double*) z, (double*) y, alpha, (double*) cM, rP, pitch, pitch, m_rS, NULL, avgRowSize, maxRowSize, rows, (double*) x, beta, firstIndex);

    z += vPitch;
    y += vPitch;
    x += vPitch;		
  }
}

int spmvEllDeviceDoubleF(void *deviceMat, double alpha, void* deviceX, double beta, void* deviceY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  /*spgpuDellspmv (handle, (double*) y->v_, (double*)y->v_, alpha, (double*) devMat->cM, devMat->rP, devMat->cMPitch, devMat->rPPitch, devMat->rS, devMat->rows, (double*)x->v_, beta, devMat->baseIndex);*/
  /* fprintf(stderr, "From spmvEllDouble: mat %d %d %d %d y %d %d \n", */
  /* 	  devMat->avgRowSize, devMat->maxRowSize, devMat->rows, */
  /* 	  devMat->pitch, y->count_, y->pitch_); */

  dspmdmm_gpu((double*) y->v_, y->count_, y->pitch_, (double*) y->v_, alpha, (double*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, (double*) x->v_, beta, devMat->baseIndex);
  
  return SPGPU_SUCCESS;
}

int spmvEllDeviceDoubleMV(void *deviceMat, double alpha, void* deviceX, double beta, void* deviceY, int idxY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  double* xColAdd = (double*) x->v_;
  double* yColAdd = ((double*) y->v_) + (IDX_F2C(idxY) * y->pitch_);

  dspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, alpha, (double*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, beta, devMat->baseIndex);

  return SPGPU_SUCCESS;
}

int spmvEllDeviceDoubleVM(void *deviceMat, double alpha, void* deviceX, int idxX, double beta, void* deviceY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  double* xColAdd = ((double*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  double* yColAdd = (double*) y->v_;

  dspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, alpha, (double*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, beta, devMat->baseIndex);

  return SPGPU_SUCCESS;
}

int spmvEllDeviceDoubleMMi(void *deviceMat, double alpha, void* deviceX, int idxX, double beta, void* deviceY, int idxY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  double* xColAdd = ((double*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  double* yColAdd = ((double*) y->v_) + (IDX_F2C(idxY) * y->pitch_);

  dspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, alpha, (double*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, beta, devMat->baseIndex);

  return SPGPU_SUCCESS;
}


void cspmdmm_gpu(cuFloatComplex* z, int s, int vPitch, cuFloatComplex* y, cuFloatComplex alpha, cuFloatComplex* cM, int* rP, int* rS, int avgRowSize, int maxRowSize, int rows, int pitch, cuFloatComplex* x, cuFloatComplex beta, int firstIndex)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  int* m_rS = (PASS_RS) ? rS : NULL;

  for(int i = 0; i < s; ++i)
  {
    spgpuCellspmv(handle, (cuFloatComplex*) z, (cuFloatComplex*) y, alpha, (cuFloatComplex*) cM, rP, pitch, pitch, m_rS, NULL, avgRowSize, maxRowSize, rows, (cuFloatComplex*) x, beta, firstIndex);

    z += vPitch;
    y += vPitch;
    x += vPitch;		
  }
}

int spmvEllDeviceFloatComplexF(void *deviceMat, float complex alpha, void* deviceX, float complex beta, void* deviceY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  cuFloatComplex a = make_cuFloatComplex(crealf(alpha), cimagf(alpha));
  cuFloatComplex b = make_cuFloatComplex(crealf(beta), cimagf(beta));

  /*spgpuCellspmv (handle, (cuFloatComplex*) y->v_, (cuFloatComplex*)y->v_, alpha, (cuFloatComplex*) devMat->cM, devMat->rP, devMat->cMPitch, devMat->rPPitch, devMat->rS, devMat->rows, (cuFloatComplex*)x->v_, beta, devMat->baseIndex);*/
  /* fprintf(stderr, "From spmvEllcuFloatComplex: mat %d %d %d %d y %d %d \n", */
  /* 	  devMat->avgRowSize, devMat->maxRowSize, devMat->rows, */
  /* 	  devMat->pitch, y->count_, y->pitch_); */

  cspmdmm_gpu((cuFloatComplex*) y->v_, y->count_, y->pitch_, (cuFloatComplex*) y->v_, a, (cuFloatComplex*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, (cuFloatComplex*) x->v_, b, devMat->baseIndex);
  
  return SPGPU_SUCCESS;
}

int spmvEllDeviceFloatComplexMV(void *deviceMat, float complex alpha, void* deviceX, float complex beta, void* deviceY, int idxY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  cuFloatComplex a = make_cuFloatComplex(crealf(alpha), cimagf(alpha));
  cuFloatComplex b = make_cuFloatComplex(crealf(beta), cimagf(beta));

  cuFloatComplex* xColAdd = (cuFloatComplex*) x->v_;
  cuFloatComplex* yColAdd = ((cuFloatComplex*) y->v_) + (IDX_F2C(idxY) * y->pitch_);

  cspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, a, (cuFloatComplex*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, b, devMat->baseIndex);

  return SPGPU_SUCCESS;
}

int spmvEllDeviceFloatComplexVM(void *deviceMat, float complex alpha, void* deviceX, int idxX, float complex beta, void* deviceY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  cuFloatComplex a = make_cuFloatComplex(crealf(alpha), cimagf(alpha));
  cuFloatComplex b = make_cuFloatComplex(crealf(beta), cimagf(beta));

  cuFloatComplex* xColAdd = ((cuFloatComplex*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  cuFloatComplex* yColAdd = (cuFloatComplex*) y->v_;

  cspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, a, (cuFloatComplex*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, b, devMat->baseIndex);

  return SPGPU_SUCCESS;
}

int spmvEllDeviceFloatComplexMMi(void *deviceMat, float complex alpha, void* deviceX, int idxX, float complex beta, void* deviceY, int idxY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  cuFloatComplex a = make_cuFloatComplex(crealf(alpha), cimagf(alpha));
  cuFloatComplex b = make_cuFloatComplex(crealf(beta), cimagf(beta));

  cuFloatComplex* xColAdd = ((cuFloatComplex*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  cuFloatComplex* yColAdd = ((cuFloatComplex*) y->v_) + (IDX_F2C(idxY) * y->pitch_);

  cspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, a, (cuFloatComplex*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, b, devMat->baseIndex);

  return SPGPU_SUCCESS;
}


void zspmdmm_gpu(cuDoubleComplex* z, int s, int vPitch, cuDoubleComplex* y, cuDoubleComplex alpha, cuDoubleComplex* cM, int* rP, int* rS, int avgRowSize, int maxRowSize, int rows, int pitch, cuDoubleComplex* x, cuDoubleComplex beta, int firstIndex)
{
  spgpuHandle_t handle = psb_cudaGetHandle();
  int* m_rS = (PASS_RS) ? rS : NULL;

  for(int i = 0; i < s; ++i)
  {
    spgpuZellspmv(handle, (cuDoubleComplex*) z, (cuDoubleComplex*) y, alpha, (cuDoubleComplex*) cM, rP, pitch, pitch, m_rS, NULL, avgRowSize, maxRowSize, rows, (cuDoubleComplex*) x, beta, firstIndex);

    z += vPitch;
    y += vPitch;
    x += vPitch;		
  }
}

int spmvEllDeviceDoubleComplexF(void *deviceMat, double complex alpha, void* deviceX, double complex beta, void* deviceY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  cuDoubleComplex a = make_cuDoubleComplex(creal(alpha), cimag(alpha));
  cuDoubleComplex b = make_cuDoubleComplex(creal(beta), cimag(beta));

  /*spgpuZellspmv (handle, (cuDoubleComplex*) y->v_, (cuDoubleComplex*)y->v_, alpha, (cuDoubleComplex*) devMat->cM, devMat->rP, devMat->cMPitch, devMat->rPPitch, devMat->rS, devMat->rows, (cuDoubleComplex*)x->v_, beta, devMat->baseIndex);*/
  /* fprintf(stderr, "From spmvEllcuDoubleComplex: mat %d %d %d %d y %d %d \n", */
  /* 	  devMat->avgRowSize, devMat->maxRowSize, devMat->rows, */
  /* 	  devMat->pitch, y->count_, y->pitch_); */

  zspmdmm_gpu((cuDoubleComplex*) y->v_, y->count_, y->pitch_, (cuDoubleComplex*) y->v_, a, (cuDoubleComplex*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, (cuDoubleComplex*) x->v_, b, devMat->baseIndex);
  
  return SPGPU_SUCCESS;
}

int spmvEllDeviceDoubleComplexMV(void *deviceMat, double complex alpha, void* deviceX, double complex beta, void* deviceY, int idxY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  cuDoubleComplex a = make_cuDoubleComplex(creal(alpha), cimag(alpha));
  cuDoubleComplex b = make_cuDoubleComplex(creal(beta), cimag(beta));

  cuDoubleComplex* xColAdd = (cuDoubleComplex*) x->v_;
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) y->v_) + (IDX_F2C(idxY) * y->pitch_);

  zspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, a, (cuDoubleComplex*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, b, devMat->baseIndex);

  return SPGPU_SUCCESS;
}

int spmvEllDeviceDoubleComplexVM(void *deviceMat, double complex alpha, void* deviceX, int idxX, double complex beta, void* deviceY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  cuDoubleComplex a = make_cuDoubleComplex(creal(alpha), cimag(alpha));
  cuDoubleComplex b = make_cuDoubleComplex(creal(beta), cimag(beta));

  cuDoubleComplex* xColAdd = ((cuDoubleComplex*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  cuDoubleComplex* yColAdd = (cuDoubleComplex*) y->v_;

  zspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, a, (cuDoubleComplex*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, b, devMat->baseIndex);

  return SPGPU_SUCCESS;
}

int spmvEllDeviceDoubleComplexMMi(void *deviceMat, double complex alpha, void* deviceX, int idxX, double complex beta, void* deviceY, int idxY)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;

  cuDoubleComplex a = make_cuDoubleComplex(creal(alpha), cimag(alpha));
  cuDoubleComplex b = make_cuDoubleComplex(creal(beta), cimag(beta));

  cuDoubleComplex* xColAdd = ((cuDoubleComplex*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) y->v_) + (IDX_F2C(idxY) * y->pitch_);

  zspmdmm_gpu(yColAdd, 1, y->pitch_, yColAdd, a, (cuDoubleComplex*) devMat->cM, devMat->rP, devMat->rS, devMat->avgRowSize, devMat->maxRowSize, devMat->rows, devMat->pitch, xColAdd, b, devMat->baseIndex);

  return SPGPU_SUCCESS;
}


int writeEllDeviceFloat(void* deviceMat, float* val, int* ja, int ldj, int* irn, int *idiag)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  // Ex updateFromHost function
  int i = writeRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(float));
  if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));

  /* i = writeEllDevice(deviceMat, (void *) val, ja, irn);
  if (i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceFloat", i);
  }*/
  return SPGPU_SUCCESS;
}

int writeEllDeviceDouble(void* deviceMat, double* val, int* ja, int ldj, int* irn, int *idiag)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  // Ex updateFromHost function
  int i = writeRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(double));
  if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));

  /* i = writeEllDevice(deviceMat, (void *) val, ja, irn);
  if (i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}

int writeEllDeviceFloatComplex(void* deviceMat, float complex* val, int* ja, int ldj, int* irn, int *idiag)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  // Ex updateFromHost function
  int i;
  i = writeRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(cuFloatComplex));
  i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));

  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);
  if (i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}

int writeEllDeviceDoubleComplex(void* deviceMat, double complex* val, int* ja, int ldj, int* irn, int *idiag)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  // Ex updateFromHost function
  int i;
  i = writeRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(cuDoubleComplex));
  i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));

  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);
  if (i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}


int readEllDeviceFloat(void* deviceMat, float* val, int* ja, int ldj, int* irn, int *idiag)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  int i;
  i = readRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(float));
  i = readRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = readRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));

  /*i = readEllDevice(deviceMat, (void*) val, ja, irn);
  if (i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "readEllDeviceFloat", i);
  }*/
  return SPGPU_SUCCESS;
}

int readEllDeviceDouble(void* deviceMat, double* val, int* ja, int ldj, int* irn, int *idiag)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  int i;
  i = readRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(double));
  i = readRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = readRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));

  /*if (i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "readEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}

int readEllDeviceFloatComplex(void* deviceMat, float complex* val, int* ja, int ldj, int* irn, int *idiag)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  int i;
  i = readRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(cuFloatComplex));
  i = readRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = readRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));

  /*if (i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "readEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}

int readEllDeviceDoubleComplex(void* deviceMat, double complex* val, int* ja, int ldj, int* irn, int *idiag)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  int i;
  i = readRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(cuDoubleComplex));
  i = readRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = readRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));

  /*if (i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "readEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}


int getEllDevicePitch(void* deviceMat)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  int i = devMat->pitch; //old
  //i = getPitchEllDevice(deviceMat);
  return(i);
}

int getEllDeviceMaxRowSize(void* deviceMat)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  int i = devMat->maxRowSize;
  return(i);
}


// New copying interface
int psiCopyCooToElgFloat(int nr, int nc, int nza, int hacksz, int ldv, int nzm, int *irn, int *idisp, int *ja, float *val, void *deviceMat)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  float *devVal;
  int *devIdisp, *devJa;
  spgpuHandle_t handle = psb_cudaGetHandle();

  allocRemoteBuffer((void**) &(devIdisp), ((size_t) (nr+1)) * sizeof(int));
  allocRemoteBuffer((void**) &(devJa), ((size_t) (nza)) * sizeof(int));
  allocRemoteBuffer((void**) &(devVal), ((size_t) (nza)) * sizeof(float));

  int i = writeRemoteBuffer((void*) val, (void*) devVal, ((size_t) nza) * sizeof(float));
  if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devJa, ((size_t) nza) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) idisp, (void*) devIdisp, ((size_t) devMat->rows + 1) * sizeof(int));

  if(i == 0) psi_cuda_s_CopyCooToElg(handle, nr, nc, nza, devMat->baseIndex, hacksz, ldv, nzm, 
                                    (int*) devMat->rS, devIdisp, devJa, devVal, 
                                    (int*) devMat->diag, (int*) devMat->rP, (float*) devMat->cM);
  // Ex updateFromHost function
  //i = writeRemoteBuffer((void*) val, (void*) devMat->cM, devMat->allocsize*sizeof(float));
  //if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, devMat->allocsize*sizeof(int));
  //if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, devMat->rows*sizeof(int));

  freeRemoteBuffer(devIdisp);
  freeRemoteBuffer(devJa);
  freeRemoteBuffer(devVal);
  
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);*/
  if (i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceFloat", i);
  
  return SPGPU_SUCCESS;
}

int psiCopyCooToElgDouble(int nr, int nc, int nza, int hacksz, int ldv, int nzm, int *irn, int *idisp, int *ja, double *val, void *deviceMat)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  double *devVal;
  int *devIdisp, *devJa;
  spgpuHandle_t handle = psb_cudaGetHandle();

  allocRemoteBuffer((void**) &(devIdisp), ((size_t) nr+1) * sizeof(int));
  allocRemoteBuffer((void**) &(devJa), ((size_t) nza) * sizeof(int));
  allocRemoteBuffer((void**) &(devVal), ((size_t) nza) * sizeof(double));

  int i = writeRemoteBuffer((void*) val, (void*) devVal, ((size_t) nza) * sizeof(double));
  if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devJa, ((size_t) nza) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) idisp, (void*) devIdisp, ((size_t) devMat->rows + 1) * sizeof(int));

  if(i == 0) psi_cuda_d_CopyCooToElg(handle, nr, nc, nza, devMat->baseIndex, hacksz, ldv, nzm, 
                                    (int*) devMat->rS, devIdisp, devJa, devVal, 
                                    (int*) devMat->diag, (int*) devMat->rP, (double *)devMat->cM);
  // Ex updateFromHost function
  //i = writeRemoteBuffer((void*) val, (void*) devMat->cM, devMat->allocsize*sizeof(double));
  //if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, devMat->allocsize*sizeof(int));
  //if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, devMat->rows*sizeof(int));

  freeRemoteBuffer(devIdisp);
  freeRemoteBuffer(devJa);
  freeRemoteBuffer(devVal);
  
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);*/
  if (i != 0) 
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceDouble", i);
  
  return SPGPU_SUCCESS;
}

int psiCopyCooToElgFloatComplex(int nr, int nc, int nza, int hacksz, int ldv, int nzm, int *irn, int *idisp, int *ja, float complex *val, void *deviceMat)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  float complex *devVal;
  int *devIdisp, *devJa;
  spgpuHandle_t handle = psb_cudaGetHandle();

  allocRemoteBuffer((void**) &(devIdisp), ((size_t) nr+1) * sizeof(int));
  allocRemoteBuffer((void**) &(devJa), ((size_t) nza) * sizeof(int));
  allocRemoteBuffer((void**) &(devVal), ((size_t) nza) * sizeof(cuFloatComplex));

  int i = writeRemoteBuffer((void*) val, (void*) devVal, ((size_t) nza) * sizeof(cuFloatComplex));
  if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devJa, ((size_t) nza) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) idisp, (void*) devIdisp, ((size_t) devMat->rows + 1) * sizeof(int));

  if(i == 0) psi_cuda_c_CopyCooToElg(handle, nr, nc, nza, devMat->baseIndex, hacksz, ldv, nzm, 
                                    (int*) devMat->rS, devIdisp, devJa, devVal, 
                                    (int*) devMat->diag, (int*) devMat->rP, (float complex*) devMat->cM);
  // Ex updateFromHost function
  //i = writeRemoteBuffer((void*) val, (void*) devMat->cM, devMat->allocsize*sizeof(float complex));
  //if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, devMat->allocsize*sizeof(int));
  //if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, devMat->rows*sizeof(int));

  freeRemoteBuffer(devIdisp);
  freeRemoteBuffer(devJa);
  freeRemoteBuffer(devVal);
  
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);*/
  if (i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceFloatComplex", i);
  
  return SPGPU_SUCCESS;
}

int psiCopyCooToElgDoubleComplex(int nr, int nc, int nza, int hacksz, int ldv, int nzm, int *irn, int *idisp, int *ja, double complex *val, void *deviceMat)
{ 
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  double complex *devVal;
  int *devIdisp, *devJa;
  spgpuHandle_t handle = psb_cudaGetHandle();

  allocRemoteBuffer((void**) &(devIdisp), ((size_t) nr+1) * sizeof(int));
  allocRemoteBuffer((void**) &(devJa), ((size_t) nza) * sizeof(int));
  allocRemoteBuffer((void**) &(devVal), ((size_t) nza) * sizeof(cuDoubleComplex));

  int i = writeRemoteBuffer((void*) val, (void*) devVal, ((size_t) nza) * sizeof(cuDoubleComplex));
  if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devJa, ((size_t) nza) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) i = writeRemoteBuffer((void*) idisp, (void*) devIdisp, ((size_t) (devMat->rows + 1)) * sizeof(int));

  if(i == 0) psi_cuda_z_CopyCooToElg(handle, nr, nc, nza, devMat->baseIndex, hacksz, ldv, nzm, 
                                    (int*) devMat->rS, devIdisp, devJa, devVal, 
                                    (int*) devMat->diag, (int*) devMat->rP, (double complex*) devMat->cM);
  // Ex updateFromHost function
  //i = writeRemoteBuffer((void*) val, (void*) devMat->cM, devMat->allocsize*sizeof(double complex));
  //if(i == 0) i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, devMat->allocsize*sizeof(int));
  //if(i == 0) i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, devMat->rows*sizeof(int));

  freeRemoteBuffer(devIdisp);
  freeRemoteBuffer(devJa);
  freeRemoteBuffer(devVal);
  
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);*/
  if (i != 0) 
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceDoubleComplex", i);
  
  return SPGPU_SUCCESS;
}


int dev_csputEllDeviceFloat(void* deviceMat, int nnz, void *ia, void *ja, void *val)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* devVal = (struct MultiVectDevice*) val;
  struct MultiVectDevice* devIa = (struct MultiVectDevice*) ia;
  struct MultiVectDevice* devJa = (struct MultiVectDevice*) ja;
  float alpha = 1.0;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if(nnz <= 0) return SPGPU_SUCCESS; 
  //fprintf(stderr, "Going through csputEllDeviceDouble %d %p %d\n", nnz, devUpdIdx, cnt);
  
  spgpuSellcsput(handle, alpha, (float*) devMat->cM, devMat->rP, devMat->pitch, devMat->pitch,
                devMat->rS, nnz, devIa->v_, devJa->v_, (float*) devVal->v_, 1);

  return SPGPU_SUCCESS; 
}

int dev_csputEllDeviceDouble(void* deviceMat, int nnz, void *ia, void *ja, void *val)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* devVal = (struct MultiVectDevice*) val;
  struct MultiVectDevice* devIa = (struct MultiVectDevice*) ia;
  struct MultiVectDevice* devJa = (struct MultiVectDevice*) ja;
  double alpha = 1.0;
  spgpuHandle_t handle = psb_cudaGetHandle();

  if(nnz <= 0) return SPGPU_SUCCESS; 
  //fprintf(stderr, "Going through csputEllDeviceDouble %d %p %d\n", nnz, devUpdIdx, cnt);
  
  spgpuDellcsput(handle, alpha, (double *) devMat->cM, devMat->rP, devMat->pitch, devMat->pitch,
		            devMat->rS, nnz, devIa->v_, devJa->v_, (double *) devVal->v_, 1);

  return SPGPU_SUCCESS; 
}

int dev_csputEllDeviceFloatComplex(void* deviceMat, int nnz, void *ia, void *ja, void *val)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* devVal = (struct MultiVectDevice*) val;
  struct MultiVectDevice* devIa = (struct MultiVectDevice*) ia;
  struct MultiVectDevice* devJa = (struct MultiVectDevice*) ja;
  cuFloatComplex alpha = make_cuFloatComplex(1.0, 0.0);
  spgpuHandle_t handle = psb_cudaGetHandle();

  if(nnz <= 0) return SPGPU_SUCCESS; 
  //fprintf(stderr, "Going through csputEllDeviceDouble %d %p %d\n", nnz, devUpdIdx, cnt);
  
  spgpuCellcsput(handle, alpha, (cuFloatComplex* ) devMat->cM, devMat->rP, devMat->pitch, devMat->pitch,
		            devMat->rS, nnz, devIa->v_, devJa->v_, (cuFloatComplex* ) devVal->v_, 1);

  return SPGPU_SUCCESS; 
}

int dev_csputEllDeviceDoubleComplex(void* deviceMat, int nnz, void *ia, void *ja, void *val)
{
  struct EllDevice* devMat = (struct EllDevice*) deviceMat;
  struct MultiVectDevice* devVal = (struct MultiVectDevice*) val;
  struct MultiVectDevice* devIa = (struct MultiVectDevice*) ia;
  struct MultiVectDevice* devJa = (struct MultiVectDevice*) ja;
  cuDoubleComplex alpha = make_cuDoubleComplex(1.0, 0.0);
  spgpuHandle_t handle = psb_cudaGetHandle();

  if(nnz <= 0) return SPGPU_SUCCESS; 
  //fprintf(stderr, "Going through csputEllDeviceDouble %d %p %d\n", nnz, devUpdIdx, cnt);
  
  spgpuZellcsput(handle, alpha, (cuDoubleComplex *) devMat->cM, devMat->rP, devMat->pitch, devMat->pitch,
		            devMat->rS, nnz, devIa->v_, devJa->v_, (cuDoubleComplex *) devVal->v_, 1);

  return SPGPU_SUCCESS; 
}
