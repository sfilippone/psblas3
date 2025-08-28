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
#include "vectordev.h"
#include "cuda_runtime.h"
#include "core.h"
#include "vector.h"

int registerMappedFloat(void *, void **, int, float);
int writeMultiVecDeviceFloat(void* deviceMultiVec, float* hostMultiVec);
int writeMultiVecDeviceFloatR2(void* deviceMultiVec, float* hostMultiVec, int ld);
int readMultiVecDeviceFloat(void* deviceMultiVec, float* hostMultiVec);
int readMultiVecDeviceFloatR2(void* deviceMultiVec, float* hostMultiVec, int ld);

int setscalMultiVecDeviceFloat(float val, int first, int last, 
			       int indexBase, void* devVecX); 

int geinsMultiVecDeviceFloat(int n, void* devVecIrl, void* devVecVal, 
			      int dupl, int indexBase, void* devVecX); 

int igathMultiVecDeviceFloatVecIdx(void* deviceVec, int vectorId, int n,
				    int first, void* deviceIdx, int hfirst,
				    void* host_values, int indexBase);
int igathMultiVecDeviceFloat(void* deviceVec, int vectorId, int n,
			      int first, void* indexes, int hfirst, void* host_values, 
			      int indexBase);
int iscatMultiVecDeviceFloatVecIdx(void* deviceVec, int vectorId, int n, int first, 
				    void *deviceIdx, int hfirst, void* host_values, 
				    int indexBase, float beta);
int iscatMultiVecDeviceFloat(void* deviceVec, int vectorId, int n, int first, void *indexes,
			      int hfirst, void* host_values, int indexBase, float beta);

int scalMultiVecDeviceFloat(float alpha, void* devMultiVecA);
int nrm2MultiVecDeviceFloat(float* y_res, int n, void* devVecA);
int amaxMultiVecDeviceFloat(float* y_res, int n, void* devVecA);
int asumMultiVecDeviceFloat(float* y_res, int n, void* devVecA);
int dotMultiVecDeviceFloat(float* y_res, int n, void* devVecA, void* devVecB);

int axpbyMultiVecDeviceFloat(int n, float alpha, void* devVecX, float beta, void* devVecY);
int upd_xyzMultiVecDeviceFloat(int n,float alpha,float beta, float gamma, float delta, 
			       void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ);
int xyzwMultiVecDeviceFloat(int n,float a,float b, float c, float d, float e, float f, 
			    void* devMultiVecX, void* devMultiVecY,
			    void* devMultiVecZ, void* devMultiVecW);
int axyMultiVecDeviceFloat(int n, float alpha, void *deviceVecA, void *deviceVecB);
int axybzMultiVecDeviceFloat(int n, float alpha, void *deviceVecA,
			      void *deviceVecB, float beta, void *deviceVecZ);
int absMultiVecDeviceFloat(int n, float alpha, void *deviceVecA);
int absMultiVecDeviceFloat2(int n, float alpha, void *deviceVecA, void *deviceVecB);


