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
 
#ifndef DCUSPARSE_
#define DCUSPARSE_
  

#include <stdio.h>
#include <stdlib.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include "cintrf.h"

/*    Double precision real   */ 
#define TYPE 			       double
#define CUSPARSE_BASE_TYPE             CUDA_R_64F
#define T_CSRGDeviceMat		       d_CSRGDeviceMat
#define T_Cmat			       d_Cmat
#define T_spmvCSRGDevice	       d_spmvCSRGDevice
#define T_spsvCSRGDevice	       d_spsvCSRGDevice
#define T_CSRGDeviceAlloc	       d_CSRGDeviceAlloc
#define T_CSRGDeviceFree	       d_CSRGDeviceFree
#define T_CSRGHost2Device	       d_CSRGHost2Device
#define T_CSRGDevice2Host	       d_CSRGDevice2Host
#define T_CSRGDeviceSetMatFillMode     d_CSRGDeviceSetMatFillMode
#define T_CSRGDeviceSetMatDiagType     d_CSRGDeviceSetMatDiagType
#define T_CSRGDeviceGetParms	       d_CSRGDeviceGetParms

#if CUDA_SHORT_VERSION <= 10  
#define T_CSRGDeviceSetMatType	       d_CSRGDeviceSetMatType
#define T_CSRGDeviceSetMatIndexBase    d_CSRGDeviceSetMatIndexBase
#define T_CSRGDeviceCsrsmAnalysis      d_CSRGDeviceCsrsmAnalysis
#define cusparseTcsrmv		       cusparseDcsrmv
#define cusparseTcsrsv_solve	       cusparseDcsrsv_solve
#define cusparseTcsrsv_analysis	       cusparseDcsrsv_analysis
#define T_HYBGDeviceMat		       d_HYBGDeviceMat
#define T_Hmat			       d_Hmat
#define T_HYBGDeviceFree	       d_HYBGDeviceFree
#define T_spmvHYBGDevice	       d_spmvHYBGDevice
#define T_HYBGDeviceAlloc	       d_HYBGDeviceAlloc
#define T_HYBGDeviceSetMatDiagType     d_HYBGDeviceSetMatDiagType
#define T_HYBGDeviceSetMatIndexBase    d_HYBGDeviceSetMatIndexBase
#define T_HYBGDeviceSetMatType	       d_HYBGDeviceSetMatType
#define T_HYBGDeviceSetMatFillMode     d_HYBGDeviceSetMatFillMode
#define T_HYBGDeviceHybsmAnalysis      d_HYBGDeviceHybsmAnalysis
#define T_spsvHYBGDevice	       d_spsvHYBGDevice
#define T_HYBGHost2Device	       d_HYBGHost2Device
#define cusparseThybmv		       cusparseDhybmv
#define cusparseThybsv_solve	       cusparseDhybsv_solve
#define cusparseThybsv_analysis	       cusparseDhybsv_analysis
#define cusparseTcsr2hyb               cusparseDcsr2hyb               

#elif CUDA_VERSION <  11030

#define T_CSRGDeviceSetMatType	       d_CSRGDeviceSetMatType
#define T_CSRGDeviceSetMatIndexBase    d_CSRGDeviceSetMatIndexBase
#define T_CSRGDeviceCsrsv2Analysis     d_CSRGDeviceCsrsv2Analysis
#define cusparseTcsrsv2_bufferSize     cusparseDcsrsv2_bufferSize
#define cusparseTcsrsv2_analysis       cusparseDcsrsv2_analysis
#define cusparseTcsrsv2_solve	       cusparseDcsrsv2_solve
#else

#define T_CSRGIsNullSvBuffer	       d_CSRGIsNullSvBuffer
#define T_CSRGIsNullSvDescr	       d_CSRGIsNullSvDescr
#define T_CSRGIsNullMvDescr	       d_CSRGIsNullMvDescr
#define T_CSRGCreateSpMVDescr	       d_CSRGCreateSpMVDescr

#endif

#include "fcusparse.h"

#endif

