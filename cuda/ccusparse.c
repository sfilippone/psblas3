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

#ifdef HAVE_SPGPU
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include "cintrf.h"
#include "fcusparse.h"

/*    Single precision complex   */ 
#define TYPE 			       float complex                     
#define CUSPARSE_BASE_TYPE             CUDA_C_32F
#define T_CSRGDeviceMat		       c_CSRGDeviceMat
#define T_Cmat			       c_Cmat
#define T_spmvCSRGDevice	       c_spmvCSRGDevice
#define T_spsvCSRGDevice	       c_spsvCSRGDevice
#define T_CSRGDeviceAlloc	       c_CSRGDeviceAlloc
#define T_CSRGDeviceFree	       c_CSRGDeviceFree
#define T_CSRGHost2Device	       c_CSRGHost2Device
#define T_CSRGDevice2Host	       c_CSRGDevice2Host
#define T_CSRGDeviceSetMatFillMode     c_CSRGDeviceSetMatFillMode
#define T_CSRGDeviceSetMatDiagType     c_CSRGDeviceSetMatDiagType
#define T_CSRGDeviceGetParms	       c_CSRGDeviceGetParms

#if CUDA_SHORT_VERSION <= 10  

#define T_CSRGDeviceSetMatType	       c_CSRGDeviceSetMatType
#define T_CSRGDeviceSetMatIndexBase    c_CSRGDeviceSetMatIndexBase
#define T_CSRGDeviceCsrsmAnalysis      c_CSRGDeviceCsrsmAnalysis
#define cusparseTcsrmv		       cusparseCcsrmv
#define cusparseTcsrsv_solve	       cusparseCcsrsv_solve
#define cusparseTcsrsv_analysis	       cusparseCcsrsv_analysis

#elif CUDA_VERSION <  11030

#define T_CSRGDeviceSetMatType	       c_CSRGDeviceSetMatType
#define T_CSRGDeviceSetMatIndexBase    c_CSRGDeviceSetMatIndexBase
#define T_CSRGDeviceCsrsv2Analysis     c_CSRGDeviceCsrsv2Analysis
#define cusparseTcsrsv2_bufferSize     cusparseCcsrsv2_bufferSize
#define cusparseTcsrsv2_analysis       cusparseCcsrsv2_analysis
#define cusparseTcsrsv2_solve	       cusparseCcsrsv2_solve

#else

#define T_HYBGDeviceMat		       c_HYBGDeviceMat
#define T_Hmat			       c_Hmat
#define T_HYBGDeviceFree	       c_HYBGDeviceFree
#define T_spmvHYBGDevice	       c_spmvHYBGDevice
#define T_HYBGDeviceAlloc	       c_HYBGDeviceAlloc
#define T_HYBGDeviceSetMatDiagType     c_HYBGDeviceSetMatDiagType
#define T_HYBGDeviceSetMatIndexBase    c_HYBGDeviceSetMatIndexBase
#define T_HYBGDeviceSetMatType	       c_HYBGDeviceSetMatType
#define T_HYBGDeviceSetMatFillMode     c_HYBGDeviceSetMatFillMode
#define T_HYBGDeviceHybsmAnalysis      c_HYBGDeviceHybsmAnalysis
#define T_spsvHYBGDevice	       c_spsvHYBGDevice
#define T_HYBGHost2Device	       c_HYBGHost2Device
#define cusparseThybmv		       cusparseChybmv
#define cusparseThybsv_solve	       cusparseChybsv_solve
#define cusparseThybsv_analysis	       cusparseChybsv_analysis
#define cusparseTcsr2hyb               cusparseCcsr2hyb               
#endif

#include "fcusparse_fct.h"

#endif 
