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
#ifndef FCUSPARSE_DAT_
#define FCUSPARSE_DAT_
  

typedef struct T_CSRGDeviceMat
{			
#if CUDA_SHORT_VERSION <= 10
  cusparseMatDescr_t descr;
  cusparseSolveAnalysisInfo_t triang;
#elif CUDA_VERSION <  11030
  cusparseMatDescr_t descr;
  csrsv2Info_t   triang;
  size_t        mvbsize, svbsize;
  void      *mvbuffer, *svbuffer;
#else
  cusparseSpMatDescr_t *spmvDescr;
  cusparseSpSVDescr_t  *spsvDescr;
  size_t        mvbsize, svbsize;
  void      *mvbuffer, *svbuffer;
#endif
  int                 m, n, nz;
  TYPE               *val;
  int                *irp;
  int                *ja;
  int                *rowBlocks;
  int                 numBlocks;
} T_CSRGDeviceMat;

/* Interoperability: type coming from Fortran side to distinguish D/S/C/Z. */
typedef struct T_Cmat
{
  T_CSRGDeviceMat *mat;
} T_Cmat;

#if CUDA_SHORT_VERSION <= 10
typedef struct T_HYBGDeviceMat
{			
  cusparseMatDescr_t descr;
  cusparseSolveAnalysisInfo_t triang;
  cusparseHybMat_t hybA;
  int                m, n, nz;
  TYPE              *val;
  int               *irp;
  int               *ja;
} T_HYBGDeviceMat;


/* Interoperability: type coming from Fortran side to distinguish D/S/C/Z. */
typedef struct T_Hmat
{
  T_HYBGDeviceMat *mat;
} T_Hmat;
#endif

int T_spmvCSRGDevice(T_Cmat *Mat, TYPE alpha, void *deviceX, 
		     TYPE beta, void *deviceY);
int T_spsvCSRGDevice(T_Cmat *Mat, TYPE alpha, void *deviceX, 
		     TYPE beta, void *deviceY);
int T_CSRGDeviceAlloc(T_Cmat *Mat,int nr, int nc, int nz);
int T_CSRGDeviceFree(T_Cmat *Mat);


int T_CSRGHost2Device(T_Cmat *Mat, int m, int n, int nz,
		      int *irp, int *ja, TYPE *val);
int T_CSRGDevice2Host(T_Cmat *Mat, int m, int n, int nz,
		      int *irp, int *ja, TYPE *val);

int T_CSRGDeviceGetParms(T_Cmat *Mat,int *nr, int *nc, int *nz);

#if CUDA_SHORT_VERSION <= 10
int T_CSRGDeviceSetMatType(T_Cmat *Mat, int type);
int T_CSRGDeviceSetMatFillMode(T_Cmat *Mat, int type);
int T_CSRGDeviceSetMatDiagType(T_Cmat *Mat, int type);
int T_CSRGDeviceSetMatIndexBase(T_Cmat *Mat, int type);
int T_CSRGDeviceCsrsmAnalysis(T_Cmat *Mat);
#elif CUDA_VERSION <  11030
int T_CSRGDeviceSetMatType(T_Cmat *Mat, int type);
int T_CSRGDeviceSetMatFillMode(T_Cmat *Mat, int type);
int T_CSRGDeviceSetMatDiagType(T_Cmat *Mat, int type);
int T_CSRGDeviceSetMatIndexBase(T_Cmat *Mat, int type);
#else

int T_CSRGCreateSpMVDescr(T_CSRGDeviceMat *cMat);
int T_CSRGIsNullSvBuffer(T_CSRGDeviceMat *cMat);
int T_CSRGIsNullSvDescr(T_CSRGDeviceMat *cMat);
int T_CSRGIsNullMvDescr(T_CSRGDeviceMat *cMat);
#endif



#if CUDA_SHORT_VERSION <= 10


int T_HYBGDeviceFree(T_Hmat *Matrix);
int T_spmvHYBGDevice(T_Hmat *Matrix, TYPE alpha, void *deviceX,
		     TYPE beta, void *deviceY);
int T_HYBGDeviceAlloc(T_Hmat *Matrix,int nr, int nc, int nz);
int T_HYBGDeviceSetMatDiagType(T_Hmat *Matrix, int type);
int T_HYBGDeviceSetMatIndexBase(T_Hmat *Matrix, int type);
int T_HYBGDeviceSetMatType(T_Hmat *Matrix, int type);
int T_HYBGDeviceSetMatFillMode(T_Hmat *Matrix, int type);
int T_HYBGDeviceHybsmAnalysis(T_Hmat *Matrix);
int T_spsvHYBGDevice(T_Hmat *Matrix, TYPE alpha, void *deviceX,
		     TYPE beta, void *deviceY);
int T_HYBGHost2Device(T_Hmat *Matrix, int m, int n, int nz,
			  int *irp, int *ja, TYPE *val);
#endif

#endif
