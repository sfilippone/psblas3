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
 

#pragma once
//#include "utils.h"
#include "cuda_runtime.h"
//#include "common.h"
//#include "cintrf.h"
#include "cuda_util.h"
#include <complex.h>

struct MultiVectDevice
{
  // number of vectors
  int count_;

  //number of elements for a single vector
  int size_;

  //pithc in number of elements
  int pitch_;

  // Vectors in device memory (single allocation)
  void *v_;
};

typedef struct MultiVectorDeviceParams
{			
	// number on vectors
	unsigned int count; //1 for a simple vector

	// The resulting allocation will be pitch*s*(size of the elementType)
	unsigned int elementType;
	
	// Pitch (in number of elements)
	unsigned int pitch;

	// Size of a single vector (in number of elements).
	unsigned int size; 
} MultiVectorDeviceParams;


#define INS_OVERWRITE 0
#define INS_ADD       1


int unregisterMapped(void *);

MultiVectorDeviceParams getMultiVectorDeviceParams(unsigned int count,
						   unsigned int size,  
						   unsigned int elementType);

int FallocMultiVecDevice(void** deviceMultiVec, unsigned count, 
			 unsigned int size, unsigned int elementType);
void freeMultiVecDevice(void* deviceVec);
int allocMultiVecDevice(void ** remoteMultiVec, struct MultiVectorDeviceParams *params);
int getMultiVecDeviceSize(void* deviceVec);
int getMultiVecDeviceCount(void* deviceVec);
int getMultiVecDevicePitch(void* deviceVec);

