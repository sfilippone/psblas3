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
 
  

#ifndef _ELLDEV_H_
#define _ELLDEV_H_

#include "cintrf.h"
#include "vectordev.h"
#include "cuComplex.h"
#include "ell.h"

struct EllDevice
{
  // Compressed matrix
  void *cM; //it can be float or double

  // row pointers (same size of cM)
  int *rP;
  int *diag;
  // row size
  int *rS;

  //matrix size (uncompressed)
  int rows;
  int columns;

  int pitch; //old

  int cMPitch;
  
  int rPPitch;

  int maxRowSize;
  int avgRowSize;
  
  //allocation size (in elements)
  int allocsize;

  /*(i.e. 0 for C, 1 for Fortran)*/
  int baseIndex;
  /* real/complex, single/double */
  int  dataType;
  
};

typedef struct EllDeviceParams
{			
	// The resulting allocation for cM and rP will be pitch*maxRowSize*(size of the elementType)
	unsigned int elementType;
	
	// Pitch (in number of elements)
	unsigned int pitch;

	// Number of rows.
	// Used to allocate rS array
	unsigned int rows; 
		
	// Number of columns.
	// Used for error-checking
	unsigned int columns; 
	
	// Largest row size
        unsigned int maxRowSize;
        unsigned int avgRowSize;
	
	// First index (e.g 0 or 1)
	unsigned int firstIndex;
} EllDeviceParams;

int computeEllAllocPitch(int rowsCount);
int FallocEllDevice(void** deviceMat, unsigned int rows, unsigned int maxRowSize, 
		    unsigned int nnzeros,
		    unsigned int columns, unsigned int elementType, 
		    unsigned int firstIndex);
int allocEllDevice(void ** remoteMatrix, EllDeviceParams* params);
void freeEllDevice(void* remoteMatrix);

int writeEllDeviceFloat(void* deviceMat, float* val, int* ja, int ldj, int* irn, int *idiag);
int writeEllDeviceDouble(void* deviceMat, double* val, int* ja, int ldj, int* irn, int *idiag);
int writeEllDeviceFloatComplex(void* deviceMat, float complex* val, int* ja, int ldj, int* irn, int *idiag);
int writeEllDeviceDoubleComplex(void* deviceMat, double complex* val, int* ja, int ldj, int* irn, int *idiag);

int readEllDeviceFloat(void* deviceMat, float* val, int* ja, int ldj, int* irn, int *idiag);
int readEllDeviceDouble(void* deviceMat, double* val, int* ja, int ldj, int* irn, int *idiag);
int readEllDeviceFloatComplex(void* deviceMat, float complex* val, int* ja, int ldj, int* irn, int *idiag);
int readEllDeviceDoubleComplex(void* deviceMat, double complex* val, int* ja, int ldj, int* irn, int *idiag);

int spmvEllDeviceFloat(void *deviceMat, float alpha, void* deviceX, 
		       float beta, void* deviceY);
int spmvEllDeviceDouble(void *deviceMat, double alpha, void* deviceX, 
			double beta, void* deviceY);
int spmvEllDeviceFloatComplex(void *deviceMat, float complex alpha, void* deviceX,
			      float complex beta, void* deviceY);
int spmvEllDeviceDoubleComplex(void *deviceMat, double complex alpha, void* deviceX,
			       double complex beta, void* deviceY);



int psiCopyCooToElgFloat(int nr, int nc, int nza, int hacksz, int ldv, int nzm, int *irn,
			  int *idisp, int *ja, float *val, void *deviceMat);

int psiCopyCooToElgDouble(int nr, int nc, int nza, int hacksz, int ldv, int nzm, int *irn,
			  int *idisp, int *ja, double *val, void *deviceMat);

int psiCopyCooToElgFloatComplex(int nr, int nc, int nza, int hacksz, int ldv, int nzm, int *irn,
			  int *idisp, int *ja, float complex *val, void *deviceMat);

int psiCopyCooToElgDoubleComplex(int nr, int nc, int nza, int hacksz, int ldv, int nzm, int *irn,
			  int *idisp, int *ja, double complex *val, void *deviceMat);


void psi_cuda_s_CopyCooToElg(spgpuHandle_t handle, int nr, int nc, int nza, int baseIdx,
			     int hacksz, int ldv, int nzm,
			     int *rS,int *devIdisp, int *devJa, float *devVal,
			     int *idiag, int *rP, float *cM);

void psi_cuda_d_CopyCooToElg(spgpuHandle_t handle, int nr, int nc, int nza, int baseIdx,
			     int hacksz, int ldv, int nzm,
			     int *rS,int *devIdisp, int *devJa, double *devVal,
			     int *idiag, int *rP, double *cM);

void psi_cuda_c_CopyCooToElg(spgpuHandle_t handle, int nr, int nc, int nza, int baseIdx,
			     int hacksz, int ldv, int nzm,
			     int *rS,int *devIdisp, int *devJa, float complex *devVal,
			     int *idiag, int *rP, float complex *cM);

void psi_cuda_z_CopyCooToElg(spgpuHandle_t handle, int nr, int nc, int nza, int baseIdx,
			     int hacksz, int ldv, int nzm,
			     int *rS,int *devIdisp, int *devJa, double complex *devVal,
			     int *idiag, int *rP, double complex *cM);


int dev_csputEllDeviceFloat(void* deviceMat, int nnz,
			    void *ia, void *ja, void *val);
int dev_csputEllDeviceDouble(void* deviceMat, int nnz,
			     void *ia, void *ja, void *val);
int dev_csputEllDeviceFloatComplex(void* deviceMat,  int nnz,
				   void *ia, void *ja, void *val);
int dev_csputEllDeviceDoubleComplex(void* deviceMat, int nnz,
				    void *ia, void *ja, void *val);

void zeroEllDevice(void* deviceMat);

int getEllDevicePitch(void* deviceMat);

// sparse Ell matrix-vector product
//int spmvEllDeviceFloat(void *deviceMat, float* alpha, void* deviceX, float* beta, void* deviceY);
//int spmvEllDeviceDouble(void *deviceMat, double* alpha, void* deviceX, double* beta, void* deviceY);

#endif
