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
 
#ifndef SCUSPARSE_
#define SCUSPARSE_
  

#include <stdio.h>
#include <stdlib.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include "cintrf.h"

/*    Double precision real   */ 
#define TYPE 			       float 
#define CUSPARSE_BASE_TYPE             CUDA_R_32F
#define T_CSRGDeviceMat		       s_CSRGDeviceMat
#define T_Cmat			       s_Cmat
#define T_spmvCSRGDevice	       s_spmvCSRGDevice
#define T_spsvCSRGDevice	       s_spsvCSRGDevice
#define T_CSRGDeviceAlloc	       s_CSRGDeviceAlloc
#define T_CSRGDeviceFree	       s_CSRGDeviceFree
#define T_CSRGHost2Device	       s_CSRGHost2Device
#define T_CSRGDevice2Host	       s_CSRGDevice2Host
#define T_CSRGDeviceSetMatFillMode     s_CSRGDeviceSetMatFillMode
#define T_CSRGDeviceSetMatDiagType     s_CSRGDeviceSetMatDiagType
#define T_CSRGDeviceGetParms	       s_CSRGDeviceGetParms

#if CUDA_SHORT_VERSION <= 10  
#define T_CSRGDeviceSetMatType	       s_CSRGDeviceSetMatType
#define T_CSRGDeviceSetMatIndexBase    s_CSRGDeviceSetMatIndexBase
#define T_CSRGDeviceCsrsmAnalysis      s_CSRGDeviceCsrsmAnalysis
#define cusparseTcsrmv		       cusparseScsrmv
#define cusparseTcsrsv_solve	       cusparseScsrsv_solve
#define cusparseTcsrsv_analysis	       cusparseScsrsv_analysis
#define T_HYBGDeviceMat		       s_HYBGDeviceMat
#define T_Hmat			       s_Hmat
#define T_HYBGDeviceFree	       s_HYBGDeviceFree
#define T_spmvHYBGDevice	       s_spmvHYBGDevice
#define T_HYBGDeviceAlloc	       s_HYBGDeviceAlloc
#define T_HYBGDeviceSetMatDiagType     s_HYBGDeviceSetMatDiagType
#define T_HYBGDeviceSetMatIndexBase    s_HYBGDeviceSetMatIndexBase
#define T_HYBGDeviceSetMatType	       s_HYBGDeviceSetMatType
#define T_HYBGDeviceSetMatFillMode     s_HYBGDeviceSetMatFillMode
#define T_HYBGDeviceHybsmAnalysis      s_HYBGDeviceHybsmAnalysis
#define T_spsvHYBGDevice	       s_spsvHYBGDevice
#define T_HYBGHost2Device	       s_HYBGHost2Device
#define cusparseThybmv		       cusparseShybmv
#define cusparseThybsv_solve	       cusparseShybsv_solve
#define cusparseThybsv_analysis	       cusparseShybsv_analysis
#define cusparseTcsr2hyb               cusparseScsr2hyb               

#elif CUDA_VERSION <  11030

#define T_CSRGDeviceSetMatType	       s_CSRGDeviceSetMatType
#define T_CSRGDeviceSetMatIndexBase    s_CSRGDeviceSetMatIndexBase
#define T_CSRGDeviceCsrsv2Analysis     s_CSRGDeviceCsrsv2Analysis
#define cusparseTcsrsv2_bufferSize     cusparseScsrsv2_bufferSize
#define cusparseTcsrsv2_analysis       cusparseScsrsv2_analysis
#define cusparseTcsrsv2_solve	       cusparseScsrsv2_solve
#else

#define T_CSRGIsNullSvBuffer	       s_CSRGIsNullSvBuffer
#define T_CSRGIsNullSvDescr	       s_CSRGIsNullSvDescr
#define T_CSRGIsNullMvDescr	       s_CSRGIsNullMvDescr
#define T_CSRGCreateSpMVDescr	       s_CSRGCreateSpMVDescr

#endif

#include "fcusparse.h"

#endif
