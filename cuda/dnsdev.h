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
 
  

#ifndef _DNSDEV_H_
#define _DNSDEV_H_

#include "cintrf.h"
#include "cuComplex.h"
#include "cublas_v2.h"


struct DnsDevice
{
   // Compressed matrix
  void *cM; //it can be float or double


  //matrix size (uncompressed)
  int rows;
  int columns;

  int pitch; //old

  int cMPitch;
  
  //allocation size (in elements)
  int allocsize;

  /*(i.e. 0 for C, 1 for Fortran)*/
  int baseIndex;
};

typedef struct DnsDeviceParams
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
	
	// First index (e.g 0 or 1)
	unsigned int firstIndex;
} DnsDeviceParams;

int FallocDnsDevice(void** deviceMat, unsigned int rows, 
		    unsigned int columns, unsigned int elementType, 
		    unsigned int firstIndex);
int allocDnsDevice(void ** remoteMatrix, DnsDeviceParams* params);
void freeDnsDevice(void* remoteMatrix);

int writeDnsDeviceFloat(void* deviceMat, float* val, int lda, int nc);
int writeDnsDeviceDouble(void* deviceMat, double* val, int lda, int nc);
int writeDnsDeviceFloatComplex(void* deviceMat, float complex* val, int lda, int nc);
int writeDnsDeviceDoubleComplex(void* deviceMat, double complex* val, int lda, int nc);

int readDnsDeviceFloat(void* deviceMat, float* val, int lda, int nc);
int readDnsDeviceDouble(void* deviceMat, double* val, int lda, int nc);
int readDnsDeviceFloatComplex(void* deviceMat, float complex* val, int lda, int nc);
int readDnsDeviceDoubleComplex(void* deviceMat, double complex* val, int lda, int nc);

int spmvDnsDeviceFloat(char transa, int m, int n, int k,
		       float *alpha, void *deviceMat, void* deviceX, 
		       float *beta, void* deviceY);
int spmvDnsDeviceDouble(char transa, int m, int n, int k,
			double *alpha, void *deviceMat, void* deviceX, 
			double *beta, void* deviceY);
int spmvDnsDeviceFloatComplex(char transa, int m, int n, int k,
		       float complex *alpha, void *deviceMat, void* deviceX, 
		       float complex *beta, void* deviceY);
int spmvDnsDeviceDoubleComplex(char transa, int m, int n, int k,
			double complex *alpha, void *deviceMat, void* deviceX, 
			double complex *beta, void* deviceY);

int getDnsDevicePitch(void* deviceMat);

// sparse Dns matrix-vector product
//int spmvDnsDeviceFloat(void *deviceMat, float* alpha, void* deviceX, float* beta, void* deviceY);
//int spmvDnsDeviceDouble(void *deviceMat, double* alpha, void* deviceX, double* beta, void* deviceY);

#endif
