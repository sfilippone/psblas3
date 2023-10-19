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
#include "cintrf.h"
#include "fcusparse.h"

static   cusparseHandle_t *cusparse_handle=NULL;


void setHandle(cusparseHandle_t);

int FcusparseCreate() 
{
  int ret=CUSPARSE_STATUS_SUCCESS;
  cusparseHandle_t *handle;
  if (cusparse_handle == NULL) {
    if ((handle = (cusparseHandle_t *)malloc(sizeof(cusparseHandle_t)))==NULL) 
      return((int) CUSPARSE_STATUS_ALLOC_FAILED);
    ret = (int)cusparseCreate(handle);
    if (ret == CUSPARSE_STATUS_SUCCESS)
      cusparse_handle = handle;
  }
  return (ret);
}

int FcusparseDestroy() 
{
  int val;
  val = (int) cusparseDestroy(*cusparse_handle);
  free(cusparse_handle);
  cusparse_handle=NULL;
  return(val);
}
cusparseHandle_t *getHandle()
{
  if (cusparse_handle == NULL)
    FcusparseCreate();
  return(cusparse_handle);
}



#endif 
