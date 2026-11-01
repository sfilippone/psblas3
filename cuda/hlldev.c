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

#include "hlldev.h"

HllDeviceParams bldHllDeviceParams(unsigned int hksize, unsigned int rows, unsigned int nzeros, unsigned int allocsize, unsigned int elementType, unsigned int firstIndex)
{
  HllDeviceParams params;

  params.elementType = elementType;
  params.hackSize = hksize;
  //numero di elementi di val
  params.allocsize = allocsize;
  params.rows  = rows;
  params.nzt   = nzeros;
  params.avgNzr = (nzeros+rows-1)/rows;
  params.firstIndex = firstIndex;
  return params;
}

int getHllDeviceParams(HllDevice* mat, int* hksize, int* rows, int* nzeros, int* allocsize, int* hackOffsLength, int* firstIndex, int* avgnzr)
{
  if(mat != NULL) 
  {
    *hackOffsLength = mat->hackOffsLength;
    *hksize         = mat->hackSize;
    *nzeros         = mat->nzt; 
    *allocsize      = mat->allocsize;
    *rows           = mat->rows;
    *avgnzr         = mat->avgNzr;
    *firstIndex     = mat->baseIndex;
    return SPGPU_SUCCESS;
  } 
  else 
  {
    return SPGPU_UNSUPPORTED;
  }
}

int allocHllDevice(void** remoteMatrix, HllDeviceParams* params)
{
  HllDevice* tmp = (HllDevice*) malloc(sizeof(HllDevice));
  int ret = SPGPU_SUCCESS;
  // size_t tt;
  // int ti;
  *remoteMatrix = (void*) tmp;

  tmp->hackSize = params->hackSize;
  tmp->allocsize = params->allocsize;
  tmp->rows   = params->rows;
  tmp->avgNzr = params->avgNzr;
  tmp->nzt    = params->nzt;
  tmp->baseIndex = params->firstIndex;
  //fprintf(stderr, "Allocating HLG with %d avgNzr\n", params->avgNzr);
  tmp->hackOffsLength = (int)(tmp->rows + tmp->hackSize - 1)/tmp->hackSize;

  //printf("hackOffsLength %d\n", tmp->hackOffsLength);
 
  if(ret == SPGPU_SUCCESS)
    ret = allocRemoteBuffer((void**) &(tmp->rP), ((size_t) tmp->allocsize) * sizeof(int));
  
  if(ret == SPGPU_SUCCESS)
    ret = allocRemoteBuffer((void**) &(tmp->rS), ((size_t) tmp->rows) * sizeof(int));

  if(ret == SPGPU_SUCCESS)
    ret = allocRemoteBuffer((void**) &(tmp->diag), ((size_t) tmp->rows) * sizeof(int));

  if(ret == SPGPU_SUCCESS)
    ret = allocRemoteBuffer((void**) &(tmp->hackOffs), (((size_t) tmp->hackOffsLength + 1) * sizeof(int)));


  if(ret != SPGPU_SUCCESS)
    return ret;
  
  switch(params->elementType)
  {
    case SPGPU_TYPE_INT:            ret = allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(int));             break;
    case SPGPU_TYPE_FLOAT:          ret = allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(float));           break;
    case SPGPU_TYPE_DOUBLE:         ret = allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(double));          break;
    case SPGPU_TYPE_COMPLEX_FLOAT:  ret = allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(cuFloatComplex));  break;
    case SPGPU_TYPE_COMPLEX_DOUBLE: ret = allocRemoteBuffer((void**) &(tmp->cM), ((size_t) tmp->allocsize) * sizeof(cuDoubleComplex)); break;
    default: return SPGPU_UNSUPPORTED;
  }

  return ret;
}

void freeHllDevice(void* remoteMatrix)
{
  HllDevice* devMat = (HllDevice*) remoteMatrix;  
  //fprintf(stderr, "freeHllDevice\n");
  if(devMat != NULL)
  {
    freeRemoteBuffer(devMat->rS);
    freeRemoteBuffer(devMat->diag);
    freeRemoteBuffer(devMat->rP);
    freeRemoteBuffer(devMat->cM);
    free(remoteMatrix);
  }
}

int FallocHllDevice(void** deviceMat, unsigned int hksize, unsigned int rows, unsigned int nzeros, unsigned int allocsize, unsigned int elementType, unsigned int firstIndex)
{
  HllDeviceParams p = bldHllDeviceParams(hksize, rows, nzeros, allocsize, elementType, firstIndex);
  int i = allocHllDevice(deviceMat, &p);
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "FallocHllDevice", i);
  
  return(i);
}


int spmvHllDeviceFloatF(void *deviceMat, float alpha, void* deviceX, float beta, void* deviceY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  if((x->count_ == 1) && (y->count_ == 1))
    spgpuShellspmv(handle, (float*) y->v_, (float*) y->v_, alpha, (float*) devMat->cM, 
                  devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                  devMat->avgNzr, devMat->rows, (double*) x->v_, beta, devMat->baseIndex);
  else
    spgpuShellspmm(handle, (float*) y->v_, y->pitch_, (float*) y->v_, y->pitch_, alpha, (float*) devMat->cM, 
                  devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                  devMat->avgNzr, devMat->rows, (float*) x->v_, x->pitch_, beta, x->count_, devMat->baseIndex);

  //cudaSync();
  return SPGPU_SUCCESS;
}

int spmvHllDeviceFloatMV(void *deviceMat, float alpha, void* deviceX, float beta, void* deviceY, int idxY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  float* xColAdd = (float*) x->v_;
  float* yColAdd = ((float*) y->v_) + (IDX_F2C(idxY) * y->pitch_);
  spgpuShellspmv(handle, yColAdd, yColAdd, alpha, (float*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, beta, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmvHllDeviceFloatVM(void *deviceMat, float alpha, void* deviceX, int idxX, float beta, void* deviceY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  float* xColAdd = ((float*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  float* yColAdd = (float*) y->v_;
  spgpuShellspmv(handle, yColAdd, yColAdd, alpha, (float*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, beta, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmvHllDeviceFloatMMi(void *deviceMat, float alpha, void* deviceX, int idxX, float beta, void* deviceY, int idxY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  float* xColAdd = ((float*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  float* yColAdd = ((float*) y->v_) + (IDX_F2C(idxY) * y->pitch_);
  spgpuShellspmv(handle, yColAdd, yColAdd, alpha, (float*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, beta, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}


int spmvHllDeviceDoubleF(void *deviceMat, double alpha, void* deviceX, double beta, void* deviceY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  if((x->count_ == 1) && (y->count_ == 1))
    spgpuDhellspmv(handle, (double*) y->v_, (double*) y->v_, alpha, (double*) devMat->cM, 
                  devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                  devMat->avgNzr, devMat->rows, (double*) x->v_, beta, devMat->baseIndex);
  else
    spgpuDhellspmm(handle, (double*) y->v_, y->pitch_, (double*) y->v_, y->pitch_, alpha, (double*) devMat->cM, 
                  devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                  devMat->avgNzr, devMat->rows, (double*) x->v_, x->pitch_, beta, x->count_, devMat->baseIndex);

  //cudaSync();
  return SPGPU_SUCCESS;
}

int spmvHllDeviceDoubleMV(void *deviceMat, double alpha, void* deviceX, double beta, void* deviceY, int idxY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  double* xColAdd = (double*) x->v_;
  double* yColAdd = ((double*) y->v_) + (IDX_F2C(idxY) * y->pitch_);
  spgpuDhellspmv(handle, yColAdd, yColAdd, alpha, (double*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, beta, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmvHllDeviceDoubleVM(void *deviceMat, double alpha, void* deviceX, int idxX, double beta, void* deviceY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  double* xColAdd = ((double*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  double* yColAdd = (double*) y->v_;
  spgpuDhellspmv(handle, yColAdd, yColAdd, alpha, (double*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, beta, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmvHllDeviceDoubleMMi(void *deviceMat, double alpha, void* deviceX, int idxX, double beta, void* deviceY, int idxY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  double* xColAdd = ((double*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  double* yColAdd = ((double*) y->v_) + (IDX_F2C(idxY) * y->pitch_);
  spgpuDhellspmv(handle, yColAdd, yColAdd, alpha, (double*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, beta, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}


int spmvHllDeviceFloatComplexF(void *deviceMat, float complex alpha, void* deviceX, float complex beta, void* deviceY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  cuFloatComplex a = make_cuFloatComplex(crealf(alpha), cimagf(alpha));
  cuFloatComplex b = make_cuFloatComplex(crealf(beta), cimagf(beta));

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  if((x->count_ == 1) && (y->count_ == 1))
    spgpuChellspmv(handle, (cuFloatComplex*) y->v_, (cuFloatComplex*) y->v_, a, (cuFloatComplex*) devMat->cM, 
                  devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                  devMat->avgNzr, devMat->rows, (cuFloatComplex*) x->v_, b, devMat->baseIndex);
  else
    spgpuChellspmm(handle, (cuFloatComplex*) y->v_, y->pitch_, (cuFloatComplex*) y->v_, y->pitch_, a, (cuFloatComplex*) devMat->cM, 
                  devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                  devMat->avgNzr, devMat->rows, (cuFloatComplex*) x->v_, x->pitch_, b, x->count_, devMat->baseIndex);

  //cudaSync();
  return SPGPU_SUCCESS;
}

int spmvHllDeviceFloatComplexMV(void *deviceMat, float complex alpha, void* deviceX, float complex beta, void* deviceY, int idxY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  cuFloatComplex a = make_cuFloatComplex(crealf(alpha), cimagf(alpha));
  cuFloatComplex b = make_cuFloatComplex(crealf(beta), cimagf(beta));

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  cuFloatComplex* xColAdd = (cuFloatComplex*) x->v_;
  cuFloatComplex* yColAdd = ((cuFloatComplex*) y->v_) + (IDX_F2C(idxY) * y->pitch_);
  spgpuChellspmv(handle, yColAdd, yColAdd, a, (cuFloatComplex*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, b, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmvHllDeviceFloatComplexVM(void *deviceMat, float complex alpha, void* deviceX, int idxX, float complex beta, void* deviceY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  cuFloatComplex a = make_cuFloatComplex(crealf(alpha), cimagf(alpha));
  cuFloatComplex b = make_cuFloatComplex(crealf(beta), cimagf(beta));

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  cuFloatComplex* xColAdd = ((cuFloatComplex*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  cuFloatComplex* yColAdd = (cuFloatComplex*) y->v_;
  spgpuChellspmv(handle, yColAdd, yColAdd, a, (cuFloatComplex*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, b, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmvHllDeviceFloatComplexMMi(void *deviceMat, float complex alpha, void* deviceX, int idxX, float complex beta, void* deviceY, int idxY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  cuFloatComplex a = make_cuFloatComplex(crealf(alpha), cimagf(alpha));
  cuFloatComplex b = make_cuFloatComplex(crealf(beta), cimagf(beta));

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  cuFloatComplex* xColAdd = ((cuFloatComplex*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  cuFloatComplex* yColAdd = ((cuFloatComplex*) y->v_) + (IDX_F2C(idxY) * y->pitch_);
  spgpuChellspmv(handle, yColAdd, yColAdd, a, (cuFloatComplex*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, b, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}


int spmvHllDeviceDoubleComplexF(void *deviceMat, double complex alpha, void* deviceX, double complex beta, void* deviceY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  cuDoubleComplex a = make_cuDoubleComplex(creal(alpha), cimag(alpha));
  cuDoubleComplex b = make_cuDoubleComplex(creal(beta), cimag(beta));

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  if((x->count_ == 1) && (y->count_ == 1))
    spgpuZhellspmv(handle, (cuDoubleComplex*) y->v_, (cuDoubleComplex*) y->v_, a, (cuDoubleComplex*) devMat->cM, 
                  devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                  devMat->avgNzr, devMat->rows, (cuDoubleComplex*) x->v_, b, devMat->baseIndex);
  else
    spgpuZhellspmm(handle, (cuDoubleComplex*) y->v_, y->pitch_, (cuDoubleComplex*) y->v_, y->pitch_, a, (cuDoubleComplex*) devMat->cM, 
                  devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                  devMat->avgNzr, devMat->rows, (cuDoubleComplex*) x->v_, x->pitch_, b, x->count_, devMat->baseIndex);

  //cudaSync();
  return SPGPU_SUCCESS;
}

int spmvHllDeviceDoubleComplexMV(void *deviceMat, double complex alpha, void* deviceX, double complex beta, void* deviceY, int idxY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  cuDoubleComplex a = make_cuDoubleComplex(creal(alpha), cimag(alpha));
  cuDoubleComplex b = make_cuDoubleComplex(creal(beta), cimag(beta));

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  cuDoubleComplex* xColAdd = (cuDoubleComplex*) x->v_;
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) y->v_) + (IDX_F2C(idxY) * y->pitch_);
  spgpuZhellspmv(handle, yColAdd, yColAdd, a, (cuDoubleComplex*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, b, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmvHllDeviceDoubleComplexVM(void *deviceMat, double complex alpha, void* deviceX, int idxX, double complex beta, void* deviceY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  cuDoubleComplex a = make_cuDoubleComplex(creal(alpha), cimag(alpha));
  cuDoubleComplex b = make_cuDoubleComplex(creal(beta), cimag(beta));

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  cuDoubleComplex* xColAdd = ((cuDoubleComplex*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  cuDoubleComplex* yColAdd = (cuDoubleComplex*) y->v_;
  spgpuZhellspmv(handle, yColAdd, yColAdd, a, (cuDoubleComplex*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, b, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmvHllDeviceDoubleComplexMMi(void *deviceMat, double complex alpha, void* deviceX, int idxX, double complex beta, void* deviceY, int idxY)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  struct MultiVectDevice* x = (struct MultiVectDevice*) deviceX;
  struct MultiVectDevice* y = (struct MultiVectDevice*) deviceY;
  spgpuHandle_t handle = psb_cudaGetHandle();

  cuDoubleComplex a = make_cuDoubleComplex(creal(alpha), cimag(alpha));
  cuDoubleComplex b = make_cuDoubleComplex(creal(beta), cimag(beta));

#ifdef VERBOSE
  /*__assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif

  cuDoubleComplex* xColAdd = ((cuDoubleComplex*) x->v_) + (IDX_F2C(idxX) * x->pitch_);
  cuDoubleComplex* yColAdd = ((cuDoubleComplex*) y->v_) + (IDX_F2C(idxY) * y->pitch_);
  spgpuZhellspmv(handle, yColAdd, yColAdd, a, (cuDoubleComplex*) devMat->cM, 
                devMat->rP, devMat->hackSize, devMat->hackOffs, devMat->rS, NULL, 
                devMat->avgNzr, devMat->rows, xColAdd, b, devMat->baseIndex);
  //cudaSync();

  return SPGPU_SUCCESS;
}


int writeHllDeviceFloat(void* deviceMat, float* val, int* ja, int* hkoffs, int* irn, int* idiag)
{
  int i;
  HllDevice* devMat = (HllDevice*) deviceMat;
  // Ex updateFromHost function
  i = writeRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(float));
  i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) hkoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  //i = writeEllDevice(deviceMat, (void*) val, ja, irn);
  /*if(i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceFloat", i);
  }*/
  return SPGPU_SUCCESS;
}

int writeHllDeviceDouble(void* deviceMat, double* val, int* ja, int* hkoffs, int* irn, int* idiag)
{
  int i;
  HllDevice* devMat = (HllDevice*) deviceMat;
  // Ex updateFromHost function
  i = writeRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(double));
  i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) hkoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);
  if(i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}

int writeHllDeviceFloatComplex(void* deviceMat, float complex* val, int* ja, int* hkoffs, int* irn, int* idiag)
{
  int i;
  HllDevice* devMat = (HllDevice*) deviceMat;
  // Ex updateFromHost function
  i = writeRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(cuFloatComplex));
  i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) hkoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);
  if(i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceDouble", i); 
  }*/
  return SPGPU_SUCCESS;
}

int writeHllDeviceDoubleComplex(void* deviceMat, double complex* val, int* ja, int* hkoffs, int* irn, int* idiag)
{ int i;
  HllDevice* devMat = (HllDevice*) deviceMat;
  // Ex updateFromHost function
  i = writeRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(cuDoubleComplex));
  i = writeRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));
  i = writeRemoteBuffer((void*) hkoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);
  if(i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "writeEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}


int readHllDeviceFloat(void* deviceMat, float* val, int* ja, int* hkoffs, int* irn, int* idiag)
{
  int i;
  HllDevice* devMat = (HllDevice*) deviceMat;
  i = readRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(float));
  i = readRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = readRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) hkoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  /*i = readEllDevice(deviceMat, (void*) val, ja, irn);
  if(i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "readEllDeviceFloat", i);
  }*/
  return SPGPU_SUCCESS;
}

int readHllDeviceDouble(void* deviceMat, double* val, int* ja, int* hkoffs, int* irn, int* idiag)
{
  int i;
  HllDevice* devMat = (HllDevice*) deviceMat;
  i = readRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(double));
  i = readRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = readRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) hkoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  /*if(i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "readEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}

int readHllDeviceFloatComplex(void* deviceMat, float complex* val, int* ja, int* hkoffs, int* irn, int* idiag)
{
  int i;
  HllDevice* devMat = (HllDevice*) deviceMat;
  i = readRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(cuFloatComplex));
  i = readRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = readRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) hkoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  /*if(i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "readEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}

int readHllDeviceDoubleComplex(void* deviceMat, double complex* val, int* ja, int* hkoffs, int* irn, int* idiag)
{
  int i;
  HllDevice* devMat = (HllDevice*) deviceMat;
  i = readRemoteBuffer((void*) val, (void*) devMat->cM, ((size_t) devMat->allocsize) * sizeof(cuDoubleComplex));
  i = readRemoteBuffer((void*) ja, (void*) devMat->rP, ((size_t) devMat->allocsize) * sizeof(int));
  i = readRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) idiag, (void*) devMat->diag, ((size_t) devMat->rows) * sizeof(int));
  i = readRemoteBuffer((void*) hkoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  /*if(i != 0) {
    fprintf(stderr, "From routine : %s : %d \n", "readEllDeviceDouble", i);
  }*/
  return SPGPU_SUCCESS;
}

// New copy routines.
int psiCopyCooToHlgFloat(int nr, int nc, int nza, int hacksz, int noffs, int isz, int* irn, int* hoffs, int* idisp, int* ja, float* val, void *deviceMat)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  float* devVal;
  int* devIdisp, *devJa;
  
  //fprintf(stderr, "devMat: %p\n", devMat);
  allocRemoteBuffer((void**) &(devIdisp), ((size_t) nr + 1) * sizeof(int));
  allocRemoteBuffer((void**) &(devJa), ((size_t) nza) * sizeof(int));
  allocRemoteBuffer((void**) &(devVal), ((size_t) nza) * sizeof(float));

  // fprintf(stderr, "Writing: %d %d %d %d %d %d %d\n", nr, devMat->rows, nza, isz, hoffs[noffs], noffs, devMat->hackOffsLength);
  int i = writeRemoteBuffer((void*) val, (void*) devVal, ((size_t) nza) * sizeof(float));
  if(i == 0) 
    i = writeRemoteBuffer((void*) ja, (void*) devJa, ((size_t) nza) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) hoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) idisp, (void*) devIdisp, ((size_t) devMat->rows + 1) * sizeof(int));
  //cudaSync();

  spgpuHandle_t handle = psb_cudaGetHandle();
  psi_cuda_s_CopyCooToHlg(handle, nr, nc, nza, devMat->baseIndex, hacksz, noffs, isz, 
                          (int*) devMat->rS, (int*) devMat->hackOffs, devIdisp, devJa, devVal, 
                          (int*) devMat->diag, (int*) devMat->rP, (float*) devMat->cM);

  freeRemoteBuffer(devIdisp);
  freeRemoteBuffer(devJa);
  freeRemoteBuffer(devVal); 
  
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);*/
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeHllDeviceFloat", i);
  
  return SPGPU_SUCCESS;
}

int psiCopyCooToHlgDouble(int nr, int nc, int nza, int hacksz, int noffs, int isz, int* irn, int* hoffs, int* idisp, int* ja, double *val, void *deviceMat)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  double *devVal;
  int* devIdisp, *devJa;

  //fprintf(stderr, "devMat: %p\n", devMat);
  allocRemoteBuffer((void**) &(devIdisp), ((size_t) nr + 1) * sizeof(int));
  allocRemoteBuffer((void**) &(devJa), ((size_t) nza) * sizeof(int));
  allocRemoteBuffer((void**) &(devVal), ((size_t) nza) * sizeof(double));

  // fprintf(stderr, "Writing: %d %d %d %d %d %d %d\n", nr, devMat->rows, nza, isz, hoffs[noffs], noffs, devMat->hackOffsLength);
  int i = writeRemoteBuffer((void*) val, (void*) devVal, ((size_t) nza) * sizeof(double));
  //fprintf(stderr, "WriteRemoteBuffer   val  %d\n", i);
  if(i == 0) 
    i = writeRemoteBuffer((void*) ja, (void*) devJa, ((size_t) nza) * sizeof(int));
  //fprintf(stderr, "WriteRemoteBuffer   ja  %d\n", i);
  if(i == 0) 
    i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  //fprintf(stderr, "WriteRemoteBuffer   irn  %d\n", i);
  if(i == 0) 
    i = writeRemoteBuffer((void*) hoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  //fprintf(stderr, "WriteRemoteBuffer   hoffs  %d\n", i);
  if(i == 0) 
    i = writeRemoteBuffer((void*) idisp, (void*) devIdisp, ((size_t) devMat->rows + 1) * sizeof(int));
  
  //fprintf(stderr, "WriteRemoteBuffer   idisp  %d\n", i);
  //cudaSync();
  //fprintf(stderr, " hacksz: %d \n", hacksz);
  spgpuHandle_t handle = psb_cudaGetHandle();
  psi_cuda_d_CopyCooToHlg(handle, nr, nc, nza, devMat->baseIndex, hacksz, noffs, isz, 
                          (int*) devMat->rS, (int*) devMat->hackOffs, devIdisp, devJa, devVal, 
                          (int*) devMat->diag, (int*) devMat->rP, (double*) devMat->cM);

  freeRemoteBuffer(devIdisp);
  freeRemoteBuffer(devJa);
  freeRemoteBuffer(devVal); 
  
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);*/
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeHllDeviceDouble", i);
  
  return SPGPU_SUCCESS;
}

int psiCopyCooToHlgFloatComplex(int nr, int nc, int nza, int hacksz, int noffs, int isz, int* irn, int* hoffs, int* idisp, int* ja, float complex *val, void *deviceMat)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  float complex *devVal;
  int* devIdisp, *devJa;
  
  //fprintf(stderr, "devMat: %p\n", devMat);
  allocRemoteBuffer((void**) &(devIdisp), ((size_t) nr + 1) * sizeof(int));
  allocRemoteBuffer((void**) &(devJa), ((size_t) nza) * sizeof(int));
  allocRemoteBuffer((void**) &(devVal), ((size_t) nza) * sizeof(cuFloatComplex));

  // fprintf(stderr, "Writing: %d %d %d %d %d %d %d\n", nr, devMat->rows, nza, isz, hoffs[noffs], noffs, devMat->hackOffsLength);
  int i = writeRemoteBuffer((void*) val, (void*) devVal, ((size_t) nza) * sizeof(cuFloatComplex));
  if(i == 0) 
    i = writeRemoteBuffer((void*) ja, (void*) devJa, ((size_t) nza) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) hoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) idisp, (void*) devIdisp, ((size_t) devMat->rows + 1) * sizeof(int));
  //cudaSync();

  spgpuHandle_t handle = psb_cudaGetHandle();
  psi_cuda_c_CopyCooToHlg(handle, nr, nc, nza, devMat->baseIndex, hacksz, noffs, isz, 
                          (int*) devMat->rS, (int*) devMat->hackOffs, devIdisp, devJa, devVal, 
                          (int*) devMat->diag, (int*) devMat->rP, (float complex *)devMat->cM);

  freeRemoteBuffer(devIdisp);
  freeRemoteBuffer(devJa);
  freeRemoteBuffer(devVal); 
  
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);*/
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeHllDeviceFloatComplex", i);
  
  return SPGPU_SUCCESS;
}

int psiCopyCooToHlgDoubleComplex(int nr, int nc, int nza, int hacksz, int noffs, int isz, int* irn, int* hoffs, int* idisp, int* ja, double complex *val, void *deviceMat)
{
  HllDevice* devMat = (HllDevice*) deviceMat;
  double complex *devVal;
  int* devIdisp, *devJa;
  
  //fprintf(stderr, "devMat: %p\n", devMat);
  allocRemoteBuffer((void**) &(devIdisp), ((size_t) nr + 1) * sizeof(int));
  allocRemoteBuffer((void**) &(devJa), ((size_t) nza) * sizeof(int));
  allocRemoteBuffer((void**) &(devVal), ((size_t) nza) * sizeof(cuDoubleComplex));

  // fprintf(stderr, "Writing: %d %d %d %d %d %d %d\n", nr, devMat->rows, nza, isz, hoffs[noffs], noffs, devMat->hackOffsLength);
  int i = writeRemoteBuffer((void*) val, (void*) devVal, ((size_t) nza) * sizeof(cuDoubleComplex));
  if(i == 0) 
    i = writeRemoteBuffer((void*) ja, (void*) devJa, ((size_t) nza) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) irn, (void*) devMat->rS, ((size_t) devMat->rows) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) hoffs, (void*) devMat->hackOffs, ((size_t) devMat->hackOffsLength + 1) * sizeof(int));
  if(i == 0) 
    i = writeRemoteBuffer((void*) idisp, (void*) devIdisp, ((size_t) devMat->rows + 1) * sizeof(int));
  //cudaSync();

  spgpuHandle_t handle = psb_cudaGetHandle();
  psi_cuda_z_CopyCooToHlg(handle, nr, nc, nza, devMat->baseIndex, hacksz, noffs, isz, 
                          (int*) devMat->rS, (int*) devMat->hackOffs, devIdisp, devJa, devVal, 
                          (int*) devMat->diag, (int*) devMat->rP, (double complex *)devMat->cM);

  freeRemoteBuffer(devIdisp);
  freeRemoteBuffer(devJa);
  freeRemoteBuffer(devVal); 
  
  /*i = writeEllDevice(deviceMat, (void*) val, ja, irn);*/
  if(i != 0)
    fprintf(stderr, "From routine : %s : %d \n", "writeHllDeviceDoubleComplex", i);
  
  return SPGPU_SUCCESS;
}