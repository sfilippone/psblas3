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
#if defined(HAVE_SPGPU)
#include "cuComplex.h"
#include "vectordev.h"
#include "cuda_runtime.h"
#include "core.h"

//new
MultiVectorDeviceParams getMultiVectorDeviceParams(unsigned int count, unsigned int size, 
						   unsigned int elementType)
{
  struct MultiVectorDeviceParams params;

  if (count == 1)
    params.pitch = size;
  else
    if (elementType == SPGPU_TYPE_INT)
      {
	//fprintf(stderr,"Getting parms for  a DOUBLE vector\n");
	params.pitch = (((size*sizeof(int) + 255)/256)*256)/sizeof(int);
      }
    else if (elementType == SPGPU_TYPE_DOUBLE)
      {
	//fprintf(stderr,"Getting parms for  a DOUBLE vector\n");
	params.pitch = (((size*sizeof(double) + 255)/256)*256)/sizeof(double);
      }
    else if (elementType == SPGPU_TYPE_FLOAT)
      {
	params.pitch = (((size*sizeof(float) + 255)/256)*256)/sizeof(float);
      }
    else if (elementType == SPGPU_TYPE_COMPLEX_FLOAT)
      {
	params.pitch = (((size*sizeof(cuFloatComplex) + 255)/256)*256)/sizeof(cuFloatComplex);
      }
    else if (elementType == SPGPU_TYPE_COMPLEX_DOUBLE)
      {
	params.pitch = (((size*sizeof(cuDoubleComplex) + 255)/256)*256)/sizeof(cuDoubleComplex);
      }
    else
      params.pitch = 0;

  params.elementType = elementType;
	
  params.count = count;
  params.size = size;

  return params;

}
//new
int allocMultiVecDevice(void ** remoteMultiVec, struct MultiVectorDeviceParams *params)
{
  if (params->pitch == 0)
    return SPGPU_UNSUPPORTED; // Unsupported params
  
  struct MultiVectDevice *tmp = (struct MultiVectDevice *)malloc(sizeof(struct MultiVectDevice));
  *remoteMultiVec = (void *)tmp;
  tmp->size_ = params->size;
  tmp->count_ = params->count;

  if (params->elementType == SPGPU_TYPE_INT)
    {
      if (params->count == 1)
	tmp->pitch_ = params->size;
      else
	tmp->pitch_ = (((params->size*sizeof(int) + 255)/256)*256)/sizeof(int);
      //fprintf(stderr,"Allocating  an INT vector %ld\n",tmp->pitch_*tmp->count_*sizeof(double));
      
      return allocRemoteBuffer((void **)&(tmp->v_), tmp->pitch_*params->count*sizeof(int));
    }
  else if (params->elementType == SPGPU_TYPE_FLOAT)
    {
      if (params->count == 1)
	tmp->pitch_ = params->size;
      else
	tmp->pitch_ = (((params->size*sizeof(float) + 255)/256)*256)/sizeof(float);

      return allocRemoteBuffer((void **)&(tmp->v_), tmp->pitch_*params->count*sizeof(float));
    }
  else if (params->elementType == SPGPU_TYPE_DOUBLE)
    {
      
      if (params->count == 1)
	tmp->pitch_ = params->size;
      else
	tmp->pitch_ = (int)(((params->size*sizeof(double) + 255)/256)*256)/sizeof(double);
      //fprintf(stderr,"Allocating  a DOUBLE vector %ld\n",tmp->pitch_*tmp->count_*sizeof(double));
 
      return allocRemoteBuffer((void **)&(tmp->v_), tmp->pitch_*tmp->count_*sizeof(double));
    }
  else if (params->elementType == SPGPU_TYPE_COMPLEX_FLOAT)
    {
      if (params->count == 1)
	tmp->pitch_ = params->size;
      else
	tmp->pitch_ = (int)(((params->size*sizeof(cuFloatComplex) + 255)/256)*256)/sizeof(cuFloatComplex);
      return allocRemoteBuffer((void **)&(tmp->v_), tmp->pitch_*tmp->count_*sizeof(cuFloatComplex));
    }
  else if (params->elementType == SPGPU_TYPE_COMPLEX_DOUBLE)
    {
      if (params->count == 1)
	tmp->pitch_ = params->size;
      else
	tmp->pitch_ = (int)(((params->size*sizeof(cuDoubleComplex) + 255)/256)*256)/sizeof(cuDoubleComplex);
      return allocRemoteBuffer((void **)&(tmp->v_), tmp->pitch_*tmp->count_*sizeof(cuDoubleComplex));
    }
  else
    return SPGPU_UNSUPPORTED; // Unsupported params
  return SPGPU_SUCCESS; // Success
}


int unregisterMapped(void *buff)
{
  return unregisterMappedMemory(buff);
}

void freeMultiVecDevice(void* deviceVec)
{
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  // fprintf(stderr,"freeMultiVecDevice\n");
  if (devVec != NULL) {
    //fprintf(stderr,"Before freeMultiVecDevice% ld\n",devVec->pitch_*devVec->count_*sizeof(double));
    freeRemoteBuffer(devVec->v_);
    free(deviceVec);
  }
}

int FallocMultiVecDevice(void** deviceMultiVec, unsigned int count,
			 unsigned int size, unsigned int elementType)
{ int i;
  struct MultiVectorDeviceParams p;

  p = getMultiVectorDeviceParams(count, size, elementType);
  i = allocMultiVecDevice(deviceMultiVec, &p);
  //cudaSync();
  if (i != 0) {
    fprintf(stderr,"From routine : %s : %d, %d %d \n","FallocMultiVecDevice",i, count, size);
  }
  return(i);
}

int getMultiVecDeviceSize(void* deviceVec)
{ int i;
  struct MultiVectDevice *dev = (struct MultiVectDevice *) deviceVec;
  i = dev->size_;
  return(i);
}

int getMultiVecDeviceCount(void* deviceVec)
{ int i;
  struct MultiVectDevice *dev = (struct MultiVectDevice *) deviceVec;
  i = dev->count_;
  return(i);
}

int getMultiVecDevicePitch(void* deviceVec)
{ int i;
  struct MultiVectDevice *dev = (struct MultiVectDevice *) deviceVec;
  i = dev->pitch_;
  return(i);
}

#endif

