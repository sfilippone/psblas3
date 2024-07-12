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
 
#include <sys/time.h>
#if defined(HAVE_RSB)
#include "rsb.h"
#include "rsb_int.h"

int rsbInit()
{
  rsb_err_t errval = RSB_ERR_NO_ERROR;

  if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS))!=RSB_ERR_NO_ERROR)
    {
      printf("Error initializing the library!\n");
      return 1;
    }
  
  return 0;
}

int rsbExit()
{
  rsb_err_t errval = RSB_ERR_NO_ERROR;

  if((errval = rsb_lib_exit(RSB_NULL_INIT_OPTIONS))!=RSB_ERR_NO_ERROR)
    {
      printf("Error finalizing the library!\n");
      return 1;
    }
  
  return 0;
}

int Rsb_double_from_coo(void **rsbMat, double *va, int *ia,int *ja,int nnz,int nr,
			int nc, int br, int bc)
{
  int i=0;
  rsb_err_t errval = RSB_ERR_NO_ERROR;

  *rsbMat = rsb_mtx_alloc_from_coo_const(va,ia,ja,nnz,RSB_NUMERICAL_TYPE_DOUBLE,nr,nc,br,bc,RSB_FLAG_FORTRAN_INDICES_INTERFACE,&errval);

  if((!*rsbMat) || (errval != RSB_ERR_NO_ERROR))
    {
      printf("Error while allocating the matrix!\n");
      return 1;
    }
  return 0;
}

//X is the input and y is the output
int Rsb_double_spmv(void *rsbMat, double *x, double alfa, double *y, double beta,char trans)
{
  rsb_err_t errval = RSB_ERR_NO_ERROR;

  if(trans=='N')
    errval = rsb_spmv(RSB_TRANSPOSITION_N,&alfa,(struct rsb_mtx_t *)rsbMat,x,1,&beta,y,1);
  else
    errval = rsb_spmv(RSB_TRANSPOSITION_T,&alfa,(struct rsb_mtx_t *)rsbMat,x,1,&beta,y,1);
  
  if(errval != RSB_ERR_NO_ERROR)
    {
      printf("Error performing a multiplication!\n");
      return 1;
    }
  
  return 0;
}

//Should it return a long instead of integer?
int Rsb_getNZeros(void *rsbMat)
{
  int res = 0;
  rsb_mtx_get_info((struct rsb_mtx_t *)rsbMat,RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T,(void *)&res);
  return res;
}

void freeRsbMat(void *rsbMat)
{
  rsb_mtx_free(rsbMat);
}

#endif
