#include "stdio.h"
#include "cudalang.h"
#include "cudadebug.h"
extern "C"
{
#include "core.h"
#include "csga.h"
}

#include "debug.h"

//#define MAX_NNZ_PER_WG 6144
#define MAX_NNZ_PER_WG 4096
#define THREAD_BLOCK   1024
#define MAX_GRID_SIZE 65536
#define WARP_SIZE 32


__device__ double warp_reduce(double val){
  for(int offset=warpSize/2; offset>0; offset/=2){
    val += __shfl_down_sync(0xffffffff,val, offset);
  }
  return val;
}

__global__ void dCSGAmvINNER(double alpha, const double* as, const int* ja, const int* irp,
			     const double* multivector, int m, int n, int col_multivector,
			     const int* rowBlocks, double beta, double* resultData, int baseIndex){
  __shared__ double vals[MAX_NNZ_PER_WG];
  __shared__ int cols[MAX_NNZ_PER_WG];
  
  int startRow = rowBlocks[blockIdx.x];
  int stopRow = rowBlocks[blockIdx.x+1];
  long int numRows = stopRow - startRow;
  int nnz = irp[stopRow]-irp[startRow];
  int tid = threadIdx.x; // indice del thread nel blocco

  if (numRows > 1){
    //CSR-Stream
    //printf("csr stream\n");
    
    int localCol;
    
    for (int i = tid; i < nnz; i+= blockDim.x){ 
      localCol = irp[startRow]+i;
      vals[i] = as[localCol];
      //vals[i] *= multivector[ja[localCol]*col_multivector+j];
      cols[i] = ja[localCol];
    }
    //return;
    int firstCol = irp[startRow];
    
    __syncthreads();
    for (int t = tid; t < numRows*col_multivector; t += blockDim.x){
      int localRow = startRow + t/col_multivector;
      int j = t%col_multivector;
      double temp = 0; 
      for (int i = irp[localRow]-firstCol; i < irp[localRow+1]-firstCol; i++){
	temp += vals[i]*multivector[cols[i]*col_multivector + j];
      }
      if (beta == 0.0) {
	resultData[localRow*col_multivector +j] = alpha*temp;
      } else {
	resultData[localRow*col_multivector +j] = alpha*temp + 	beta*resultData[localRow*col_multivector +j];
      }
    }
    
    __syncthreads();    
    
  } else {
    //CSR-Vector
    //printf("csr vector\n");
    int warpId = tid / 32; // Global warp index
    int lane = tid &(0xFFFF); // thread index within the warp
    //one warp per row
    double val; 
    int col;
    double sum[64] = {0};   
    if (nnz < 4096){
      int localCol;
      for (int i = tid; i < nnz; i+= blockDim.x){ 
	localCol = irp[startRow]+i;
	vals[i] = as[localCol];
	cols[i] = ja[localCol];
      }
    }
    //return;
    __syncthreads();
    if (warpId < col_multivector){
      for (int col_m = warpId; col_m < col_multivector; col_m +=32){
	for (int i = irp[startRow] + lane; i < irp[startRow+1]; i +=32){
	  if (nnz < 4096){
	    val = vals[i-irp[startRow]];
	    col = cols[i-irp[startRow]];
	  } else {
	    val = as[i];
	    col = ja[i];
	  }
	  sum[col_m] += val*multivector[col*col_multivector + col_m];     
	}
	sum[col_m] = warp_reduce(sum[col_m]);
	if (lane == 0){
	  if (beta == 0.0) {
	    resultData[startRow*col_multivector + col_m] = alpha*sum[col_m];   
	  } else {
	    resultData[startRow*col_multivector + col_m] = alpha*sum[col_m] +
	      beta*resultData[startRow*col_multivector + col_m];
	  }	  
	}
      }
    }
  }
}


__host__ int dCSGAMV(spgpuHandle_t handle, 
		     double beta,
		     double* y, 
		     double alpha, 
		     const double* as, 
		     const int* ja,
		     const int* irp,
		     int m,
		     int n,
		     int ncol,
		     int  numBlocks,
		     const int* rowBlocks,
		     const double *x,
		     int baseIndex,
		     int *rb)
{ 
  //  fprintf(stderr," dcsgamv  %d   \n",numBlocks);
  int maxBForACall = min(handle->maxGridSizeX, numBlocks);
  //int maxBForACall = 1024;
  int blockX = THREAD_BLOCK;
  int gridX = maxBForACall;
  int rp,rows, blcks, bp, numb;
  dim3 blockSize(blockX);
  dim3 gridSize(gridX);
  //fprintf(stderr," dcsgamv  start %d %d %d  \n",m,n,ncol);

  bp = 0;
  rp = 0;
  numb = numBlocks;  
  while (numb > maxBForACall) {//managing large vectors
    blcks = maxBForACall;
    rp = rb[bp]-baseIndex;
    rows = rb[bp+blcks]-rp;
    fprintf(stderr,"  rp %d  rows %d  bp %d  blcks %d\n",rp,rows,bp,blcks);
    dCSGAmvINNER<<<gridSize,blockSize>>>(alpha,as,ja,irp,x,rows,n,ncol,
					 &(rowBlocks[bp]),beta,&(y[rp]),baseIndex);
    bp   += blcks;
    numb -= blcks;
  }
  blcks = numb;
  rp = rb[bp]-baseIndex;
  rows = rb[bp+blcks]-rp;
  //fprintf(stderr,"  rp %d  rows %d  bp %d  blcks %d\n",rp,rows,bp,blcks);
  dCSGAmvINNER<<<gridSize,blockSize>>>(alpha,as,ja,irp,x,rows,n,ncol,
				       &(rowBlocks[bp]),beta,&(y[rp]),baseIndex);
  rp += rows;
  //fprintf(stderr,"  Final  rows %d  %d  \n",rows,rp);
  return(SPGPU_SUCCESS);
}
