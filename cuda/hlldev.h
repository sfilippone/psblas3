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
 

#ifndef _HLLDEV_H_
#define _HLLDEV_H_

#include "cintrf.h"
#include "hell.h"


typedef struct hlldevice
{
  // Compressed matrix
  void *cM; //it can be float or double

  // row pointers (same size of cM)
  int *rP;
  
  // row size and diagonal position
  int *rS;
  int *diag;

  int *hackOffs;

  int rows;
  int avgNzr;
  int hackOffsLength;
  int nzt;
  int hackSize; //must be multiple of 32

  //matrix size (uncompressed)
  //int rows;
  //int columns;

  //allocation size
  int allocsize;

  /*(i.e. 0 for C, 1 for Fortran)*/
  int baseIndex;
} HllDevice;

typedef struct hlldeviceparams
{			
  
  unsigned int elementType;

  unsigned int hackSize;
  
  // Number of rows.
  // Used to allocate rS array
  unsigned int rows;
  unsigned int avgNzr;
  unsigned int nzt;
  //unsigned int hackOffsLength;
  
  // Number of columns.
  // Used for error-checking
  //  unsigned int columns; 

  unsigned int allocsize;
  
  // First index (e.g 0 or 1)
  unsigned int firstIndex;

} HllDeviceParams;


HllDeviceParams bldHllDeviceParams(unsigned int hksize, unsigned int rows, unsigned int nzeros,
				   unsigned int allocsize, 
				   unsigned int elementType, unsigned int firstIndex);
int getHllDeviceParams(HllDevice* mat, int *hksize, int *rows, int *nzeros,
		       int *allocsize, int *hackOffsLength,  int *firstIndex, int *avgnzr);
int FallocHllDevice(void** deviceMat,unsigned int hksize, unsigned int rows,  unsigned int nzeros,
		    unsigned int allocsize, unsigned int elementType, unsigned int firstIndex);
int allocHllDevice(void ** remoteMatrix, HllDeviceParams* params);
void freeHllDevice(void* remoteMatrix);
int writeHllDeviceFloat(void* deviceMat, float* val, int* ja, int *hkoffs, int* irn, int *idiag);
int writeHllDeviceDouble(void* deviceMat, double* val, int* ja, int *hkoffs, int* irn, int *idiag);
int writeHllDeviceFloatComplex(void* deviceMat, float complex* val, 
			       int* ja, int *hkoffs, int* irn, int *idiag);
int writeHllDeviceDoubleComplex(void* deviceMat, double complex* val, 
				int* ja, int *hkoffs, int* irn, int *idiag);
int readHllDeviceFloat(void* deviceMat, float* val, int* ja, int *hkoffs, int* irn, int *idiag);
int readHllDeviceDouble(void* deviceMat, double* val, int* ja, int *hkoffs, int* irn, int *idiag);
int readHllDeviceFloatComplex(void* deviceMat, float complex* val,
			      int* ja, int *hkoffs, int* irn, int *idiag);
int readHllDeviceDoubleComplex(void* deviceMat, double complex* val, 
			       int* ja, int *hkoffs, int* irn, int *idiag);


int psiCopyCooToHlgFloat(int nr, int nc, int nza, int hacksz, int noffs, int isz,
			 int *irn, int *hoffs,  int *idisp, int *ja, 
			 float *val, void *deviceMat);
int psiCopyCooToHlgDouble(int nr, int nc, int nza, int hacksz, int noffs, int isz,
			  int *irn, int *hoffs,  int *idisp, int *ja, 
			  double *val, void *deviceMat);
int psiCopyCooToHlgFloatComplex(int nr, int nc, int nza, int hacksz, 
				int noffs, int isz, int *irn, 
				int *hoffs,  int *idisp, int *ja, 
				float complex *val, void *deviceMat);
int psiCopyCooToHlgDoubleComplex(int nr, int nc, int nza, int hacksz, 
				 int noffs, int isz,  int *irn, 
				 int *hoffs,  int *idisp, int *ja, 
				 double complex *val, void *deviceMat);

int psi_cuda_s_CopyCooToHlg(spgpuHandle_t handle,int nr, int nc, int nza, 
			    int baseIdx, int hacksz, int noffs, int isz,
			    int *irn, int *hoffs,  int *idisp,
			    int *ja, float *val,
			    int *idiag, int *rP, float *cM);
int psi_cuda_d_CopyCooToHlg(spgpuHandle_t handle,int nr, int nc, int nza, 
			    int baseIdx, int hacksz, int noffs, int isz,
			    int *irn, int *hoffs,  int *idisp,
			    int *ja, double *val,
			    int *idiag, int *rP, double *cM);
int psi_cuda_c_CopyCooToHlg(spgpuHandle_t handle,int nr, int nc, int nza, 
			    int baseIdx, int hacksz, int noffs, int isz,
			    int *irn, int *hoffs,  int *idisp,
			    int *ja, float complex *val,
			    int *idiag, int *rP, float complex *cM);
int psi_cuda_z_CopyCooToHlg(spgpuHandle_t handle,int nr, int nc, int nza, 
			    int baseIdx, int hacksz, int noffs, int isz,
			    int *irn, int *hoffs,  int *idisp,
			    int *ja, double complex *val,
			    int *idiag, int *rP, double complex *cM);


#endif
