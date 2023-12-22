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
 

#include "diagdev.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
//new
DiagDeviceParams getDiagDeviceParams(unsigned int rows, unsigned int columns, unsigned int diags, unsigned int elementType)
{
  DiagDeviceParams params;

  params.elementType = elementType;
  //numero di elementi di val
  params.rows = rows;
  params.columns = columns; 
  params.diags = diags;

  return params;

}
//new
int allocDiagDevice(void ** remoteMatrix, DiagDeviceParams* params)
{
  struct DiagDevice *tmp = (struct DiagDevice *)malloc(sizeof(struct DiagDevice));
  int ret=SPGPU_SUCCESS;
  *remoteMatrix = (void *)tmp;

  tmp->rows = params->rows;

  tmp->cols = params->columns;

  tmp->diags = params->diags;
 
  if (ret == SPGPU_SUCCESS)
    ret=allocRemoteBuffer((void **)&(tmp->off), tmp->diags*sizeof(int));
  
  /* tmp->baseIndex = params->firstIndex; */

  if (params->elementType == SPGPU_TYPE_INT)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->rows*tmp->diags*sizeof(int));
    }
  else if (params->elementType == SPGPU_TYPE_FLOAT)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->rows*tmp->diags*sizeof(float));
    }    
  else if (params->elementType == SPGPU_TYPE_DOUBLE)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->rows*tmp->diags*sizeof(double));
    }
  else if (params->elementType == SPGPU_TYPE_COMPLEX_FLOAT)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->rows*tmp->diags*sizeof(cuFloatComplex));
    }
  else if (params->elementType == SPGPU_TYPE_COMPLEX_DOUBLE)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->rows*tmp->diags*sizeof(cuDoubleComplex));
    }
  else
    return SPGPU_UNSUPPORTED; // Unsupported params
  return ret;
}

void freeDiagDevice(void* remoteMatrix)
{
  struct DiagDevice *devMat = (struct DiagDevice *) remoteMatrix;  
  //fprintf(stderr,"freeHllDevice\n");
  if (devMat != NULL) {
    freeRemoteBuffer(devMat->off);
    freeRemoteBuffer(devMat->cM);
    free(remoteMatrix);
  }
}

//new
int FallocDiagDevice(void** deviceMat, unsigned int rows, unsigned int columns,unsigned int diags,unsigned int elementType)
{ int i;
  DiagDeviceParams p;
  
  p = getDiagDeviceParams(rows, columns, diags,elementType);
  i = allocDiagDevice(deviceMat, &p);
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","FallocEllDevice",i);
  }
  return(i);
}

int writeDiagDeviceDouble(void* deviceMat, double* a, int* off, int n)
{ int i,fo,fa;
  char buf_a[255], buf_o[255],tmp[255];
  struct DiagDevice *devMat = (struct DiagDevice *) deviceMat;
  // Ex updateFromHost function
  /* memset(buf_a,'\0',255); */
  /* memset(buf_o,'\0',255); */
  /* memset(tmp,'\0',255); */

  /* strcat(buf_a,"mat_"); */
  /* strcat(buf_o,"off_"); */
  /* sprintf(tmp,"%d_%d.dat",devMat->rows,devMat->cols); */
  /* strcat(buf_a,tmp); */
  /* memset(tmp,'\0',255); */
  /* sprintf(tmp,"%d.dat",devMat->cols); */
  /* strcat(buf_o,tmp); */

  /* fa = open(buf_a, O_CREAT | O_WRONLY | O_TRUNC, 0664); */
  /* fo = open(buf_o, O_CREAT | O_WRONLY | O_TRUNC, 0664); */

  /* i = write(fa, a, sizeof(double)*devMat->cols*devMat->rows); */
  /* i = write(fo, off, sizeof(int)*devMat->cols); */

  /* close(fa); */
  /* close(fo); */

  i = writeRemoteBuffer((void*) a, (void *)devMat->cM, devMat->rows*devMat->diags*sizeof(double));
  i = writeRemoteBuffer((void*) off, (void *)devMat->off, devMat->diags*sizeof(int));

  if(i==0)
    return SPGPU_SUCCESS;
  else
    return SPGPU_UNSUPPORTED;
}

int readDiagDeviceDouble(void* deviceMat, double* a, int* off)
{ int i;
  struct DiagDevice *devMat = (struct DiagDevice *) deviceMat;
  i = readRemoteBuffer((void *) a, (void *)devMat->cM,devMat->rows*devMat->diags*sizeof(double));
  i = readRemoteBuffer((void *) off, (void *)devMat->off, devMat->diags*sizeof(int));
  /*if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readEllDeviceDouble",i);
  }*/
  return SPGPU_SUCCESS;
}

//new
int spmvDiagDeviceDouble(void *deviceMat, double alpha, void* deviceX, 
		       double beta, void* deviceY)
{
  struct DiagDevice *devMat = (struct DiagDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  spgpuHandle_t handle=psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif
  /* spgpuDdiagspmv(handle, (double *)y->v_, (double *)y->v_,alpha,(double *)devMat->cM,devMat->off,devMat->rows,devMat->cols,x->v_,beta,devMat->baseIndex); */

  spgpuDdiaspmv(handle, (double *)y->v_, (double *)y->v_,alpha,(double *)devMat->cM,devMat->off,devMat->rows,devMat->rows,devMat->cols,devMat->diags,x->v_,beta);

  //cudaSync();

  return SPGPU_SUCCESS;
}


int writeDiagDeviceFloat(void* deviceMat, float* a, int* off, int n)
{ int i,fo,fa;
  char buf_a[255], buf_o[255],tmp[255];
  struct DiagDevice *devMat = (struct DiagDevice *) deviceMat;
  // Ex updateFromHost function
  /* memset(buf_a,'\0',255); */
  /* memset(buf_o,'\0',255); */
  /* memset(tmp,'\0',255); */

  /* strcat(buf_a,"mat_"); */
  /* strcat(buf_o,"off_"); */
  /* sprintf(tmp,"%d_%d.dat",devMat->rows,devMat->cols); */
  /* strcat(buf_a,tmp); */
  /* memset(tmp,'\0',255); */
  /* sprintf(tmp,"%d.dat",devMat->cols); */
  /* strcat(buf_o,tmp); */

  /* fa = open(buf_a, O_CREAT | O_WRONLY | O_TRUNC, 0664); */
  /* fo = open(buf_o, O_CREAT | O_WRONLY | O_TRUNC, 0664); */

  /* i = write(fa, a, sizeof(float)*devMat->cols*devMat->rows); */
  /* i = write(fo, off, sizeof(int)*devMat->cols); */

  /* close(fa); */
  /* close(fo); */

  i = writeRemoteBuffer((void*) a, (void *)devMat->cM, devMat->rows*devMat->diags*sizeof(float));
  i = writeRemoteBuffer((void*) off, (void *)devMat->off, devMat->diags*sizeof(int));

  if(i==0)
    return SPGPU_SUCCESS;
  else
    return SPGPU_UNSUPPORTED;
}

int readDiagDeviceFloat(void* deviceMat, float* a, int* off)
{ int i;
  struct DiagDevice *devMat = (struct DiagDevice *) deviceMat;
  i = readRemoteBuffer((void *) a, (void *)devMat->cM,devMat->rows*devMat->diags*sizeof(float));
  i = readRemoteBuffer((void *) off, (void *)devMat->off, devMat->diags*sizeof(int));
  /*if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readEllDeviceFloat",i);
  }*/
  return SPGPU_SUCCESS;
}

//new
int spmvDiagDeviceFloat(void *deviceMat, float alpha, void* deviceX, 
		       float beta, void* deviceY)
{
  struct DiagDevice *devMat = (struct DiagDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  spgpuHandle_t handle=psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif
  /* spgpuDdiagspmv(handle, (float *)y->v_, (float *)y->v_,alpha,(float *)devMat->cM,devMat->off,devMat->rows,devMat->cols,x->v_,beta,devMat->baseIndex); */

  spgpuSdiaspmv(handle, (float *)y->v_, (float *)y->v_,alpha,(float *)devMat->cM,devMat->off,devMat->rows,devMat->rows,devMat->cols,devMat->diags,x->v_,beta);

  //cudaSync();

  return SPGPU_SUCCESS;
}

