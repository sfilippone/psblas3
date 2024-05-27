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
#include "dcsga.h"


int d_CSGADeviceFree(d_Cmat *Matrix)
{
  d_CSRGDeviceMat *cMat= Matrix->mat;
  
  if (cMat!=NULL)  d_CSRGDeviceFree(cMat);
  return(CUSPARSE_STATUS_SUCCESS);
}


int d_CSGADeviceAlloc(d_Cmat *Matrix,int nr, int nc, int nz)
{
  int rc=0;
  d_CSRGDeviceMat *cMat;
  
  if ((rc=d_CSRGDeviceAlloc(Matrix,nr,nc,nz))!=0)
    return(rc);
  
  cMat = Matrix->mat;
  if (nr <= 0) nr = 1;
      
  if ((rc= allocRemoteBuffer(((void **) &(cMat->rowBlocks)), ((nr+1)*sizeof(int)))) != 0)
    return(rc);
  
  return(CUSPARSE_STATUS_SUCCESS);
}

void d_CSGAComputeRowBlocks(int totalRows, int* irp, int* numBlocks, int *rowBlocks){
    rowBlocks[0] = 1; 
    int sum = 0, last_i= 1, ctr=1; 
    for(int i = 1; i < totalRows; i++){
        sum += irp[i]-irp[i-1];
        if(sum == MAX_NNZ_PER_WG){
            last_i = i+1;
            rowBlocks[ctr++] = i+1;
            sum = 0;

        }
        else if( sum > MAX_NNZ_PER_WG){
            if(i - last_i > 1){
                rowBlocks[ctr++] = i-1 +1;
                i--;
            }
            else if(i - last_i == 1){
                rowBlocks[ctr++] = i +1;
            }
            last_i = i+1; 
            sum = 0;
        } 
    }
    
    //printf("%d %d\n",ctr,totalRows);
    *numBlocks = ctr;
    rowBlocks[ctr++] = totalRows;
    return ;
} 

int d_CSGAHost2Device(d_Cmat *Matrix,int nr, int nc, int nz,
		      int *irp, int *ja, double *val, int numBlocks, int *rowBlocks)
{
  int rc=0;
  d_CSRGDeviceMat *cMat= Matrix->mat;

  if (cMat!=NULL) {
    if ((rc=d_CSRGHost2Device(Matrix,nr,nc,nz,irp,ja,val))
	!= SPGPU_SUCCESS)
      return(rc);
    cMat->numBlocks = numBlocks;
    //    fprintf(stderr," CSGAH2D: %d   (%d:%d)  %p\n",numBlocks,
    //	    rowBlocks[0],rowBlocks[1],cMat->rowBlocks);    
    if ((rc=writeRemoteBuffer((void *) rowBlocks,(void *) cMat->rowBlocks, 
			      (numBlocks+1)*sizeof(int)))
	!= SPGPU_SUCCESS) 
      return(rc);
    //fprintf(stderr," CSGAH2D ok\n");
  } else {
    return(-1);
  }
  return(CUSPARSE_STATUS_SUCCESS);
}

int d_CSGADevice2Host(d_Cmat *Matrix,int nr, int nc, int nz,
		      int *irp, int *ja, double *val, int *numBlocks, int *rowBlocks)
{
  int rc=0;
  d_CSRGDeviceMat *cMat= Matrix->mat;
  
  if (cMat!=NULL) {
    if ((rc=d_CSRGDevice2Host(Matrix,nr,nc,nz,irp,ja,val))
	!= SPGPU_SUCCESS)
      return(rc);
    *numBlocks = cMat->numBlocks  ;
    if ((rc=readRemoteBuffer((void *) rowBlocks,(void *) cMat->rowBlocks, 
			     ((*numBlocks)+1)*sizeof(int)))
	!= SPGPU_SUCCESS) 
      return(rc);
    
  }
  return(CUSPARSE_STATUS_SUCCESS);
}


int d_spmvCSGADevice(d_Cmat *Matrix, double alpha, void* deviceX, 
		     double beta, void* deviceY, int *rb)
{
  d_CSRGDeviceMat *devMat = (d_CSRGDeviceMat *) Matrix->mat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY;
  spgpuHandle_t handle=psb_cudaGetHandle();
  int indexBase=1;
  
#if 0 
  fprintf(stderr,"devMat %p  m %d n %d nB %d  rB %p    x  %p   y %p \n",devMat,
	  devMat->m,devMat->n,devMat->numBlocks,devMat->rowBlocks,x->v_,y->v_);
  fprintf(stderr,"x_count %d  y_count %d  xsize  %d ysize %d  \n",x->count_,y->count_,
	  x->size_,y->size_);
#endif
#if  0&&defined(VERBOSE)
  __assert(x->count_ == y->count_, "ERROR: x and y don't share the same number of vectors");
  __assert(x->size_ >= devMat->n, "ERROR: x vector's size is not >= to matrix size (columns)");
  __assert(y->size_ >= devMat->m, "ERROR: y vector's size is not >= to matrix size (rows)");
#endif
  //fprintf(stderr,"Calling dCSGAMV \n");
  dCSGAMV(handle, beta,(double *)y->v_, alpha,
	  (double *)devMat->val, devMat->ja, devMat->irp, 
	  devMat->m,devMat->n,y->count_, devMat->numBlocks, devMat->rowBlocks, 
	  (double *)x->v_, 1, rb);
 
  
}

