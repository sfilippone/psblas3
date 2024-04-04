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
 
  

#include <stdio.h>
#include <stdlib.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include "cintrf.h"
#include "fcusparse.h"


/*    Double precision real   */ 
#define TYPE 			       double complex
#define CUSPARSE_BASE_TYPE             CUDA_C_64F
#define T_CSRGDeviceMat		       z_CSRGDeviceMat
#define T_Cmat			       z_Cmat
#define T_spmvCSRGDevice	       z_spmvCSRGDevice
#define T_spsvCSRGDevice	       z_spsvCSRGDevice
#define T_CSRGDeviceAlloc	       z_CSRGDeviceAlloc
#define T_CSRGDeviceFree	       z_CSRGDeviceFree
#define T_CSRGHost2Device	       z_CSRGHost2Device
#define T_CSRGDevice2Host	       z_CSRGDevice2Host
#define T_CSRGDeviceSetMatFillMode     z_CSRGDeviceSetMatFillMode
#define T_CSRGDeviceSetMatDiagType     z_CSRGDeviceSetMatDiagType
#define T_CSRGDeviceGetParms	       z_CSRGDeviceGetParms

#if CUDA_SHORT_VERSION <= 10  
#define T_CSRGDeviceSetMatType	       z_CSRGDeviceSetMatType
#define T_CSRGDeviceSetMatIndexBase    z_CSRGDeviceSetMatIndexBase
#define T_CSRGDeviceCsrsmAnalysis      z_CSRGDeviceCsrsmAnalysis
#define cusparseTcsrmv		       cusparseZcsrmv
#define cusparseTcsrsv_solve	       cusparseZcsrsv_solve
#define cusparseTcsrsv_analysis	       cusparseZcsrsv_analysis
#define T_HYBGDeviceMat		       z_HYBGDeviceMat
#define T_Hmat			       z_Hmat
#define T_HYBGDeviceFree	       z_HYBGDeviceFree
#define T_spmvHYBGDevice	       z_spmvHYBGDevice
#define T_HYBGDeviceAlloc	       z_HYBGDeviceAlloc
#define T_HYBGDeviceSetMatDiagType     z_HYBGDeviceSetMatDiagType
#define T_HYBGDeviceSetMatIndexBase    z_HYBGDeviceSetMatIndexBase
#define T_HYBGDeviceSetMatType	       z_HYBGDeviceSetMatType
#define T_HYBGDeviceSetMatFillMode     z_HYBGDeviceSetMatFillMode
#define T_HYBGDeviceHybsmAnalysis      z_HYBGDeviceHybsmAnalysis
#define T_spsvHYBGDevice	       z_spsvHYBGDevice
#define T_HYBGHost2Device	       z_HYBGHost2Device
#define cusparseThybmv		       cusparseZhybmv
#define cusparseThybsv_solve	       cusparseZhybsv_solve
#define cusparseThybsv_analysis	       cusparseZhybsv_analysis
#define cusparseTcsr2hyb               cusparseZcsr2hyb               

#elif CUDA_VERSION <  11030

#define T_CSRGDeviceSetMatType	       z_CSRGDeviceSetMatType
#define T_CSRGDeviceSetMatIndexBase    z_CSRGDeviceSetMatIndexBase
#define T_CSRGDeviceCsrsv2Analysis     z_CSRGDeviceCsrsv2Analysis
#define cusparseTcsrsv2_bufferSize     cusparseZcsrsv2_bufferSize
#define cusparseTcsrsv2_analysis       cusparseZcsrsv2_analysis
#define cusparseTcsrsv2_solve	       cusparseZcsrsv2_solve
#else

#define T_CSRGIsNullSvBuffer	       z_CSRGIsNullSvBuffer
#define T_CSRGIsNullSvDescr	       z_CSRGIsNullSvDescr
#define T_CSRGIsNullMvDescr	       z_CSRGIsNullMvDescr
#define T_CSRGCreateSpMVDescr	       z_CSRGCreateSpMVDescr

#endif

#include "fcusparse_fct.h"

