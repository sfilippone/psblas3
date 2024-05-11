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
 

#include "hdiagdev.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#define DEBUG 0
void freeHdiagDevice(void* remoteMatrix)
{
  struct HdiagDevice *devMat = (struct HdiagDevice *) remoteMatrix;  
  //fprintf(stderr,"freeHllDevice\n");
  if (devMat != NULL) {
    freeRemoteBuffer(devMat->hackOffsets);
    freeRemoteBuffer(devMat->cM);
    free(remoteMatrix);
  }
}


HdiagDeviceParams getHdiagDeviceParams(unsigned int rows, unsigned int columns, 
				       unsigned int allocationHeight, unsigned int hackSize, 
				       unsigned int hackCount, unsigned int elementType) 
{
  HdiagDeviceParams params;

  params.elementType = elementType;
  //numero di elementi di val
  params.rows = rows;
  params.columns = columns; 
  params.allocationHeight = allocationHeight;
  params.hackSize = hackSize;
  params.hackCount = hackCount;

  return params;

}

int allocHdiagDevice(void **remoteMatrix, HdiagDeviceParams* params)
{
  struct HdiagDevice *tmp = (struct HdiagDevice *)malloc(sizeof(struct HdiagDevice));
  int ret=SPGPU_SUCCESS;
  int *tmpOff = NULL;

  *remoteMatrix = (void *) tmp;
#if DEBUG
  fprintf(stderr,"From alloc: %p\n",*remoteMatrix);
#endif
  
  tmp->rows = params->rows;

  tmp->hackSize = params->hackSize;

  tmp->cols = params->columns;

  tmp->allocationHeight = params->allocationHeight;

  tmp->hackCount = params->hackCount;
  


#if DEBUG
  fprintf(stderr,"hackcount %d  allocationHeight %d\n",tmp->hackCount,tmp->allocationHeight);
#endif

  if (ret == SPGPU_SUCCESS)
    ret=allocRemoteBuffer((void **)&(tmp->hackOffsets), (tmp->hackCount+1)*sizeof(int));
  

  if (ret == SPGPU_SUCCESS)
    ret=allocRemoteBuffer((void **)&(tmp->hdiaOffsets), tmp->allocationHeight*sizeof(int));
  
  /* tmp->baseIndex = params->firstIndex; */

  if (params->elementType == SPGPU_TYPE_INT)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->hackSize*tmp->allocationHeight*sizeof(int));
    }
  else if (params->elementType == SPGPU_TYPE_FLOAT)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->hackSize*tmp->allocationHeight*sizeof(float));
    }    
  else if (params->elementType == SPGPU_TYPE_DOUBLE)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->hackSize*tmp->allocationHeight*sizeof(double));
    }
  else if (params->elementType == SPGPU_TYPE_COMPLEX_FLOAT)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->hackSize*tmp->allocationHeight*sizeof(cuFloatComplex));
    }
  else if (params->elementType == SPGPU_TYPE_COMPLEX_DOUBLE)
    {
      if (ret == SPGPU_SUCCESS)
	ret=allocRemoteBuffer((void **)&(tmp->cM), tmp->hackSize*tmp->allocationHeight*sizeof(cuDoubleComplex));
    }
  else
    return SPGPU_UNSUPPORTED; // Unsupported params
  return ret;
}

int FallocHdiagDevice(void** deviceMat, unsigned int rows, unsigned int cols, 
		      unsigned int allocationHeight, unsigned int hackSize,
		      unsigned int hackCount, unsigned int elementType)
{ int i=0;
  HdiagDeviceParams p;
 
  p = getHdiagDeviceParams(rows, cols, allocationHeight, hackSize, hackCount,elementType);
  
  i = allocHdiagDevice(deviceMat, &p);
#if DEBUG
  fprintf(stderr," Falloc  %p \n",*deviceMat);
#endif

  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","FallocEllDevice",i);
  }
  return(i);

}

int writeHdiagDeviceDouble(void* deviceMat, double* val, int* hdiaOffsets, int *hackOffsets)
{ int i=0,fo,fa,j,k,p;
  char buf_a[255], buf_o[255],tmp[255];
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  
  i=SPGPU_SUCCESS; 
  
  
#if DEBUG
  fprintf(stderr," Write  %p \n",devMat);
  
  fprintf(stderr,"HDIAG writing to device memory: allocationHeight %d hackCount %d\n",
	  devMat->allocationHeight,devMat->hackCount);
  fprintf(stderr,"HackOffsets: ");
  for (j=0; j<devMat->hackCount+1; j++)
    fprintf(stderr," %d",hackOffsets[j]);
  fprintf(stderr,"\n");
  fprintf(stderr,"diaOffsets: ");
  for (j=0; j<devMat->allocationHeight; j++)
    fprintf(stderr," %d",hdiaOffsets[j]);
  fprintf(stderr,"\n");
#if 1
  fprintf(stderr,"values: \n");
  p=0;
  for (j=0; j<devMat->hackCount; j++){
    fprintf(stderr,"Hack no: %d\n",j+1);
    for (k=0; k<devMat->hackSize*(devMat->allocationHeight/devMat->hackCount); k++){
      fprintf(stderr," %d %lf\n",p+1,val[p]); p++;
    }
  }
  fprintf(stderr,"\n");
#endif
#endif
  
  
  if(i== SPGPU_SUCCESS)
    i = writeRemoteBuffer((void *) hackOffsets,(void *) devMat->hackOffsets,
			  (devMat->hackCount+1)*sizeof(int));
  
  if(i== SPGPU_SUCCESS)
    i = writeRemoteBuffer((void*) hdiaOffsets, (void *)devMat->hdiaOffsets, 
			  devMat->allocationHeight*sizeof(int));
  if(i== SPGPU_SUCCESS)
    i = writeRemoteBuffer((void*) val, (void *)devMat->cM, 
			  devMat->allocationHeight*devMat->hackSize*sizeof(double));
  if (i!=0) 
    fprintf(stderr,"Error in writeHdiagDeviceDouble %d\n",i);

#if DEBUG
  fprintf(stderr," EndWrite  %p \n",devMat);
#endif
  
  if(i==0)
    return SPGPU_SUCCESS;
  else
    return SPGPU_UNSUPPORTED;
}



long long int sizeofHdiagDeviceDouble(void* deviceMat)
{ int i=0,fo,fa;
  int *hoff=NULL,*hackoff=NULL;
  long long int memsize=0;
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  
  
  memsize += (devMat->hackCount+1)*sizeof(int);
  memsize += devMat->allocationHeight*sizeof(int);
  memsize += devMat->allocationHeight*devMat->hackSize*sizeof(double);
  return(memsize);
}



int readHdiagDeviceDouble(void* deviceMat, double* a, int* off)
{ int i;
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  /* i = readRemoteBuffer((void *) a, (void *)devMat->cM,devMat->rows*devMat->diags*sizeof(double)); */
  /* i = readRemoteBuffer((void *) off, (void *)devMat->off, devMat->diags*sizeof(int)); */


  /*if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readEllDeviceDouble",i);
  }*/
  return SPGPU_SUCCESS;
}

int spmvHdiagDeviceDouble(void *deviceMat, double alpha, void* deviceX, 
			  double beta, void* deviceY)
{
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  spgpuHandle_t handle=psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif
#if DEBUG
  fprintf(stderr," First  %p \n",devMat);
  fprintf(stderr,"%d %d %d  %p %p %p\n",devMat->rows,devMat->cols, devMat->hackSize,
	  devMat->hackOffsets, devMat->hdiaOffsets, devMat->cM);
#endif 
  spgpuDhdiaspmv (handle, (double*)y->v_, (double *)y->v_, alpha,
		  (double *)devMat->cM,devMat->hdiaOffsets, 
		  devMat->hackSize, devMat->hackOffsets, devMat->rows,devMat->cols,
		  x->v_, beta);
  
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmmHdiagDeviceDouble(void *deviceMat, double alpha, void* deviceX, 
			  double beta, void* deviceY)
{
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  spgpuHandle_t handle=psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif
  spgpuDhdiaspmm (handle, y->count_, (double*)y->v_, y->pitch_,
      (double *)y->v_, y->pitch_, alpha, (double *)devMat->cM,
      devMat->hdiaOffsets, devMat->hackSize, devMat->hackOffsets, 
      devMat->rows, devMat->cols, (double *)x->v_, x->pitch_, beta);

  return SPGPU_SUCCESS;
}

int writeHdiagDeviceFloat(void* deviceMat, float* val, int* hdiaOffsets, int *hackOffsets)
{ int i=0,fo,fa,j,k,p;
  char buf_a[255], buf_o[255],tmp[255];
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  
  i=SPGPU_SUCCESS; 
  
  
#if DEBUG
  fprintf(stderr," Write  %p \n",devMat);
  
  fprintf(stderr,"HDIAG writing to device memory: allocationHeight %d hackCount %d\n",
	  devMat->allocationHeight,devMat->hackCount);
  fprintf(stderr,"HackOffsets: ");
  for (j=0; j<devMat->hackCount+1; j++)
    fprintf(stderr," %d",hackOffsets[j]);
  fprintf(stderr,"\n");
  fprintf(stderr,"diaOffsets: ");
  for (j=0; j<devMat->allocationHeight; j++)
    fprintf(stderr," %d",hdiaOffsets[j]);
  fprintf(stderr,"\n");
#if 1
  fprintf(stderr,"values: \n");
  p=0;
  for (j=0; j<devMat->hackCount; j++){
    fprintf(stderr,"Hack no: %d\n",j+1);
    for (k=0; k<devMat->hackSize*(devMat->allocationHeight/devMat->hackCount); k++){
      fprintf(stderr," %d %lf\n",p+1,val[p]); p++;
    }
  }
  fprintf(stderr,"\n");
#endif
#endif
  
  
  if(i== SPGPU_SUCCESS)
    i = writeRemoteBuffer((void *) hackOffsets,(void *) devMat->hackOffsets,
			  (devMat->hackCount+1)*sizeof(int));
  
  if(i== SPGPU_SUCCESS)
    i = writeRemoteBuffer((void*) hdiaOffsets, (void *)devMat->hdiaOffsets, 
			  devMat->allocationHeight*sizeof(int));
  if(i== SPGPU_SUCCESS)
    i = writeRemoteBuffer((void*) val, (void *)devMat->cM, 
			  devMat->allocationHeight*devMat->hackSize*sizeof(float));
  if (i!=0) 
    fprintf(stderr,"Error in writeHdiagDeviceFloat %d\n",i);

#if DEBUG
  fprintf(stderr," EndWrite  %p \n",devMat);
#endif
  
  if(i==0)
    return SPGPU_SUCCESS;
  else
    return SPGPU_UNSUPPORTED;
}



long long int sizeofHdiagDeviceFloat(void* deviceMat)
{ int i=0,fo,fa;
  int *hoff=NULL,*hackoff=NULL;
  long long int memsize=0;
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  
  
  memsize += (devMat->hackCount+1)*sizeof(int);
  memsize += devMat->allocationHeight*sizeof(int);
  memsize += devMat->allocationHeight*devMat->hackSize*sizeof(float);
  
  return(memsize);
}



int readHdiagDeviceFloat(void* deviceMat, float* a, int* off)
{ int i;
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  /* i = readRemoteBuffer((void *) a, (void *)devMat->cM,devMat->rows*devMat->diags*sizeof(float)); */
  /* i = readRemoteBuffer((void *) off, (void *)devMat->off, devMat->diags*sizeof(int)); */


  /*if (i != 0) {
    fprintf(stderr,"From routine : %s : %d \n","readEllDeviceFloat",i);
  }*/
  return SPGPU_SUCCESS;
}

int spmvHdiagDeviceFloat(void *deviceMat, float alpha, void* deviceX, 
			  float beta, void* deviceY)
{
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  spgpuHandle_t handle=psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif
#if DEBUG
  fprintf(stderr," First  %p \n",devMat);
  fprintf(stderr,"%d %d %d  %p %p %p\n",devMat->rows,devMat->cols, devMat->hackSize,
	  devMat->hackOffsets, devMat->hdiaOffsets, devMat->cM);
#endif 
  spgpuShdiaspmv (handle, (float*)y->v_, (float *)y->v_, alpha,
		  (float *)devMat->cM,devMat->hdiaOffsets, 
		  devMat->hackSize, devMat->hackOffsets, devMat->rows,devMat->cols,
		  x->v_, beta);
  
  //cudaSync();

  return SPGPU_SUCCESS;
}

int spmmHdiagDeviceFloat(void *deviceMat, float alpha, void* deviceX, 
			  float beta, void* deviceY)
{
  struct HdiagDevice *devMat = (struct HdiagDevice *) deviceMat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  spgpuHandle_t handle=psb_cudaGetHandle();

#ifdef VERBOSE
  /*__assert(x->count_ == x->count_, "ERROR: x and y don't share the same number of vectors");*/
  /*__assert(x->size_ >= devMat->columns, "ERROR: x vector's size is not >= to matrix size (columns)");*/
  /*__assert(y->size_ >= devMat->rows, "ERROR: y vector's size is not >= to matrix size (rows)");*/
#endif
  spgpuShdiaspmm (handle, y->count_, (float*)y->v_, y->pitch_,
      (float *)y->v_, y->pitch_, alpha, (float *)devMat->cM,
      devMat->hdiaOffsets, devMat->hackSize, devMat->hackOffsets, 
      devMat->rows, devMat->cols, (float *)x->v_, x->pitch_, beta);

  return SPGPU_SUCCESS;
}
