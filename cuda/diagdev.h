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
 

#ifndef _DIAGDEV_H_
#define _DIAGDEV_H_

#ifdef HAVE_SPGPU
#include "cintrf.h"
#include "dia.h"

struct DiagDevice
{
  // Compressed matrix
  void *cM; //it can be float or double

  // offset (same size of cM)
  int *off;

  int rows;

  int cols;

  int diags;

};

typedef struct DiagDeviceParams
{			
  
  unsigned int elementType;
  
  // Number of rows.
  // Used to allocate rS array
  unsigned int rows;
  //unsigned int hackOffsLength;
  
  // Number of columns.
  // Used for error-checking
  unsigned int columns; 

  unsigned int diags;

} DiagDeviceParams;
DiagDeviceParams getDiagDeviceParams(unsigned int rows, unsigned int columns,
				     unsigned int elementType, unsigned int firstIndex);
int FallocDiagDevice(void** deviceMat, unsigned int rows, unsigned int cols, 
		     unsigned int elementType, unsigned int firstIndex);
int allocDiagDevice(void ** remoteMatrix, DiagDeviceParams* params);
void freeDiagDevice(void* remoteMatrix);

int readDiagDeviceDouble(void* deviceMat, double* a, int* off);
int writeDiagDeviceDouble(void* deviceMat, double* a, int* off, int n);
int spmvDiagDeviceDouble(void *deviceMat, double alpha, void* deviceX, 
			 double beta, void* deviceY);

int readDiagDeviceFloat(void* deviceMat, float* a, int* off);
int writeDiagDeviceFloat(void* deviceMat, float* a, int* off, int n);
int spmvDiagDeviceFloat(void *deviceMat, float alpha, void* deviceX, 
			float beta, void* deviceY);



#else
#define CINTRF_UNSUPPORTED   -1
#endif

#endif
