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
 
#ifndef DCSGA_
#define DCSGA_
  

#include <stdio.h>
#include <stdlib.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include "cintrf.h"

#include "dcusparse.h"
#include "fcusparse.h"
#include "fcusparse_dat.h"

#define  MAX_NNZ_PER_WG  4096

int d_CSGADeviceFree(d_Cmat *Matrix);
int d_CSGADeviceAlloc(d_Cmat *Matrix,int nr, int nc, int nz);
void d_CSGAComputeRowBlocks(int totalRows, int* irp, int* numBlocks, int *rowBlocks);
int d_CSGAHost2Device(d_Cmat *Matrix,int nr, int nc, int nz,
		      int *irp, int *ja, double *val, int numBlocks, int *rowBlocks);
int d_CSGADevice2Host(d_Cmat *Matrix,int nr, int nc, int nz,
		      int *irp, int *ja, double *val, int *numBlocks, int *rowBlocks);
int d_spmvCSGADevice(d_Cmat *Matrix, double alpha, void* deviceX, 
		     double beta, void* deviceY, int *rb);

int dCSGAMV(spgpuHandle_t handle, double beta, double* y,  double alpha,
	    const double* as, const int* ja, const int* irp,
	    int m, int n,  int ncol, int  numBlocks,
	    const int* rowBlocks,   const double  *x,
	    int baseIndex, int *rb);

#endif

