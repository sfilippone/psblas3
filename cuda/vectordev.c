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
#include "cuComplex.h"
#include "vectordev.h"
#include "cuda_runtime.h"
#include "core.h"

int unregisterMapped(void *buff)
{
  return unregisterMappedMemory(buff);
}

//Getters
MultiVectorDeviceParams getMultiVectorDeviceParams(unsigned int count, unsigned int size, unsigned int elementType)
{
  const size_t alignment = 256;
  size_t elementSize = 0;
  switch(elementType)
  {
    case SPGPU_TYPE_INT:            elementSize = sizeof(int);              break;
    case SPGPU_TYPE_DOUBLE:         elementSize = sizeof(double);           break;
    case SPGPU_TYPE_FLOAT:          elementSize = sizeof(float);            break;
    case SPGPU_TYPE_COMPLEX_FLOAT:  elementSize = sizeof(cuFloatComplex);   break;
    case SPGPU_TYPE_COMPLEX_DOUBLE: elementSize = sizeof(cuDoubleComplex);  break;
    default:                        elementSize = 0;                        break;
  }

  struct MultiVectorDeviceParams params;
  params.pitch = (count == 1) ? size : (((size * elementSize + alignment - 1) / alignment) * alignment) / elementSize;
  params.elementType = elementType;
  params.count = count;
  params.size = size;
  return params;
}

int getMultiVecDeviceSize(void* deviceVec)
{
  struct MultiVectDevice *dev = (struct MultiVectDevice *) deviceVec;
  return dev->size_;
}

int getMultiVecDeviceCount(void* deviceVec)
{
  struct MultiVectDevice *dev = (struct MultiVectDevice *) deviceVec;
  return dev->count_;
}

int getMultiVecDevicePitch(void* deviceVec)
{
  struct MultiVectDevice *dev = (struct MultiVectDevice *) deviceVec;
  return dev->pitch_;
}

//Alloc/free
int FallocMultiVecDevice(void** deviceMultiVec, unsigned int count, unsigned int size, unsigned int elementType)
{
  struct MultiVectorDeviceParams p = getMultiVectorDeviceParams(count, size, elementType);
  int i = allocMultiVecDevice(deviceMultiVec, &p);
  //cudaSync();
  if(i != SPGPU_SUCCESS)
    fprintf(stderr, "From routine : %s : %d, %d %d \n", "FallocMultiVecDevice", i, count, size);

  return i;
}

int allocMultiVecDevice(void** remoteMultiVec, struct MultiVectorDeviceParams *params)
{
  if(params->pitch == 0)
    return SPGPU_UNSUPPORTED; // Unsupported params
  
  const size_t alignment = 256;
  struct MultiVectDevice *tmp = (struct MultiVectDevice *) malloc(sizeof(struct MultiVectDevice));
  *remoteMultiVec = (void*) tmp;
  tmp->size_ = params->size;
  tmp->count_ = params->count;

  size_t elementSize = 0;
  switch(params->elementType)
  {
    case SPGPU_TYPE_INT:            elementSize = sizeof(int);              break;
    case SPGPU_TYPE_DOUBLE:         elementSize = sizeof(double);           break;
    case SPGPU_TYPE_FLOAT:          elementSize = sizeof(float);            break;
    case SPGPU_TYPE_COMPLEX_FLOAT:  elementSize = sizeof(cuFloatComplex);   break;
    case SPGPU_TYPE_COMPLEX_DOUBLE: elementSize = sizeof(cuDoubleComplex);  break;
    default:                        return SPGPU_UNSUPPORTED;   // Unsupported params
  }

  tmp->pitch_ = (params->count == 1) ? params->size : (((params->size * elementSize + alignment - 1) / alignment) * alignment) / elementSize;
  return allocRemoteBuffer((void**) &(tmp->v_), ((size_t) tmp->pitch_) * tmp->count_ * elementSize);
}

void freeMultiVecDevice(void* deviceVec)
{
  struct MultiVectDevice *devVec = (struct MultiVectDevice *) deviceVec;
  // fprintf(stderr,"freeMultiVecDevice\n");
  if(devVec != NULL)
  {
    //fprintf(stderr,"Before freeMultiVecDevice% ld\n",devVec->pitch_*devVec->count_*sizeof(double));
    freeRemoteBuffer(devVec->v_);
    free(deviceVec);
  }
}
