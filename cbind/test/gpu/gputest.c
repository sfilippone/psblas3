/*----------------------------------------------------------------------------------*/
/*              Parallel Sparse BLAS  v 3.9.0					                    */
/*    (C) Copyright 2017 Salvatore Filippone                                        */
/* 										                                            */
/*  Redistribution and use in source and binary forms, with or without		        */
/*  modification, are permitted provided that the following conditions		        */
/*  are met:									                                    */
/*    1. Redistributions of source code must retain the above copyright		        */
/*       notice, this list of conditions and the following disclaimer.		        */
/*    2. Redistributions in binary form must reproduce the above copyright	        */
/*       notice, this list of conditions, and the following disclaimer in the	    */
/*       documentation and/or other materials provided with the distribution.	    */
/*    3. The name of the PSBLAS group or the names of its contributors may	        */
/*       not be used to endorse or promote products derived from this		        */
/*       software without specific written permission.				                */
/* 										                                            */
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS		        */
/*  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED	    */
/*  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR	    */
/*  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS  */
/*  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR	        */
/*  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF	        */
/*  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS	    */
/*  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN	        */
/*  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)	        */
/*  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE	    */
/*  POSSIBILITY OF SUCH DAMAGE.							                            */
/* 										                                            */
/*  File: gputest.c							                                        */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "psb_base_cbind.h"

int main(int argc, char **argv)
{
    psb_c_ctxt *cctxt;
    psb_c_i_t iam, np;

    // Initialize PSBLAS context
    cctxt = psb_c_new_ctxt();
    psb_c_init(cctxt);
    // Initialze GPU support
#ifdef PSB_HAVE_CUDA
    psb_c_cuda_init(cctxt);
#endif

    psb_c_info(*cctxt, &iam, &np);
    printf("Hello from process %d of %d\n", iam, np);
#ifdef PSB_HAVE_CUDA
    printf("Number of available GPU devices: %d seen by process %d\n", psb_c_cuda_getDeviceCount(), iam);
#endif

    // Exit from GPU support
#ifdef PSB_HAVE_CUDA
    psb_c_cuda_exit();
#endif
    // Finalize PSBLAS context
    psb_c_exit(*cctxt);
    free(cctxt);
}
