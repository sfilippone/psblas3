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
#if defined(HAVE_SPGPU)
//#include "utils.h"
#include "vectordev.h"
#include "cuda_runtime.h"
#include "core.h"

int registerMappedDouble(void *, void **, int, double);
int writeMultiVecDeviceDouble(void* deviceMultiVec, double* hostMultiVec);
int writeMultiVecDeviceDoubleR2(void* deviceMultiVec, double* hostMultiVec, int ld);
int readMultiVecDeviceDouble(void* deviceMultiVec, double* hostMultiVec);
int readMultiVecDeviceDoubleR2(void* deviceMultiVec, double* hostMultiVec, int ld);

int setscalMultiVecDeviceDouble(double val, int first, int last, 
				int indexBase, void* devVecX); 

int geinsMultiVecDeviceDouble(int n, void* devVecIrl, void* devVecVal, 
			      int dupl, int indexBase, void* devVecX); 

int igathMultiVecDeviceDoubleVecIdx(void* deviceVec, int vectorId, int n,
				    int first, void* deviceIdx, int hfirst,
				    void* host_values, int indexBase);
int igathMultiVecDeviceDouble(void* deviceVec, int vectorId, int n,
			      int first, void* indexes, int hfirst, void* host_values, 
			      int indexBase);
int iscatMultiVecDeviceDoubleVecIdx(void* deviceVec, int vectorId, int n, int first, 
				    void *deviceIdx, int hfirst, void* host_values, 
				    int indexBase, double beta);
int iscatMultiVecDeviceDouble(void* deviceVec, int vectorId, int n, int first, void *indexes,
			      int hfirst, void* host_values, int indexBase, double beta);

int scalMultiVecDeviceDouble(double alpha, void* devMultiVecA);
int nrm2MultiVecDeviceDouble(double* y_res, int n, void* devVecA);
int amaxMultiVecDeviceDouble(double* y_res, int n, void* devVecA);
int asumMultiVecDeviceDouble(double* y_res, int n, void* devVecA);
int dotMultiVecDeviceDouble(double* y_res, int n, void* devVecA, void* devVecB);

int axpbyMultiVecDeviceDouble(int n, double alpha, void* devVecX, double beta, void* devVecY);
int axyMultiVecDeviceDouble(int n, double alpha, void *deviceVecA, void *deviceVecB);
int axybzMultiVecDeviceDouble(int n, double alpha, void *deviceVecA,
			      void *deviceVecB, double beta, void *deviceVecZ);
int absMultiVecDeviceDouble(int n, double alpha, void *deviceVecA);
int absMultiVecDeviceDouble2(int n, double alpha, void *deviceVecA, void *deviceVecB);


#endif
