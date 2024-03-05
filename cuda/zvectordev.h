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
#include <complex.h>
#include "cuComplex.h"
#include "vectordev.h"
#include "cuda_runtime.h"
#include "core.h"

int registerMappedDoubleComplex(void *, void **, int, cuDoubleComplex);
int writeMultiVecDeviceDoubleComplex(void* deviceMultiVec, cuDoubleComplex* hostMultiVec);
int writeMultiVecDeviceDoubleComplexR2(void* deviceMultiVec, 
				       cuDoubleComplex* hostMultiVec, int ld);
int readMultiVecDeviceDoubleComplex(void* deviceMultiVec, cuDoubleComplex* hostMultiVec);
int readMultiVecDeviceDoubleComplexR2(void* deviceMultiVec, 
				      cuDoubleComplex* hostMultiVec, int ld);
int setscalMultiVecDeviceDoubleComplex(cuDoubleComplex val, int first, int last, 
				       int indexBase, void* devVecX); 

int geinsMultiVecDeviceDoubleComplex(int n, void* devVecIrl, void* devVecVal, 
				     int dupl, int indexBase, void* devVecX); 

int igathMultiVecDeviceDoubleComplexVecIdx(void* deviceVec, int vectorId, int n,
					   int first, void* deviceIdx, int hfirst,
					   void* host_values, int indexBase);
int igathMultiVecDeviceDoubleComplex(void* deviceVec, int vectorId, int n,
				     int first, void* indexes, 
				     int hfirst, void* host_values, 
				     int indexBase);
int iscatMultiVecDeviceDoubleComplexVecIdx(void* deviceVec, int vectorId, 
					   int n, int first, 
					   void *deviceIdx, int hfirst, 
					   void* host_values, 
					   int indexBase, cuDoubleComplex beta);
int iscatMultiVecDeviceDoubleComplex(void* deviceVec, int vectorId, int n, 
				     int first, void *indexes,
				     int hfirst, void* host_values, 
				     int indexBase, cuDoubleComplex beta);

int scalMultiVecDeviceDoubleComplex(cuDoubleComplex alpha, void* devMultiVecA);
int nrm2MultiVecDeviceDoubleComplex(cuDoubleComplex* y_res, int n, void* devVecA);
int amaxMultiVecDeviceDoubleComplex(cuDoubleComplex* y_res, int n, void* devVecA);
int asumMultiVecDeviceDoubleComplex(cuDoubleComplex* y_res, int n, void* devVecA);
int dotMultiVecDeviceDoubleComplex(cuDoubleComplex* y_res, int n, 
				   void* devVecA, void* devVecB);

int axpbyMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void* devVecX, 
				     cuDoubleComplex beta, void* devVecY);
int abgdxyzMultiVecDeviceDoubleComplex(int n,cuDoubleComplex  alpha,
				       cuDoubleComplex  beta, cuDoubleComplex  gamma, cuDoubleComplex  delta, 
				       void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ);
int xyzwMultiVecDeviceDoubleComplex(int n,cuDoubleComplex  a, cuDoubleComplex  b,
				    cuDoubleComplex  c, cuDoubleComplex  d,
				    cuDoubleComplex  e, cuDoubleComplex  f, 
				    void* devMultiVecX, void* devMultiVecY,
				    void* devMultiVecZ, void* devMultiVecW);
int axyMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, 
				   void *deviceVecA, void *deviceVecB);
int axybzMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA,
				     void *deviceVecB, cuDoubleComplex beta, 
				     void *deviceVecZ);
int absMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA);
int absMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, 
				    void *deviceVecA, void *deviceVecB);

