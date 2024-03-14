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
  
int T_spmvCSRGDevice(T_Cmat *Matrix, TYPE alpha, void *deviceX,
		     TYPE beta, void *deviceY)
{
  T_CSRGDeviceMat *cMat=Matrix->mat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY; 
  void *vX, *vY;
  int r,n;
  cusparseHandle_t *my_handle=getHandle();
  TYPE   ealpha=alpha, ebeta=beta;
#if CUDA_SHORT_VERSION <= 10
  /* getAddrMultiVecDevice(deviceX, &vX); */
  /*   getAddrMultiVecDevice(deviceY, &vY);  */
  vX=x->v_;
  vY=y->v_;

  CHECK_CUSPARSE(cusparseTcsrmv(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
				cMat->m,cMat->n,cMat->nz,(const TYPE *) &alpha,cMat->descr,
				cMat->val, cMat->irp, cMat->ja,
				(const TYPE *) vX, (const TYPE *) &beta, (TYPE *) vY));
  
#elif CUDA_VERSION <  11030
  size_t bfsz;
  vX=x->v_;
  vY=y->v_;
#if 1
  CHECK_CUSPARSE(cusparseCsrmvEx_bufferSize(*my_handle,CUSPARSE_ALG_MERGE_PATH,
					    CUSPARSE_OPERATION_NON_TRANSPOSE,
					    cMat->m,cMat->n,cMat->nz,	    
					    (const void *) &ealpha,CUSPARSE_BASE_TYPE,
					    cMat->descr,
					    (const void *) cMat->val,
					    CUSPARSE_BASE_TYPE,
					    (const int *) cMat->irp,
					    (const int *) cMat->ja,
					    (const void *) vX, CUSPARSE_BASE_TYPE,
					    (const void *) &ebeta, CUSPARSE_BASE_TYPE,
					    (void *) vY, CUSPARSE_BASE_TYPE,
					    CUSPARSE_BASE_TYPE, &bfsz));
#else 
  bfsz=cMat->nz;
#endif
  
  if (bfsz > cMat->mvbsize) {
    if (cMat->mvbuffer != NULL) {
      CHECK_CUDA(cudaFree(cMat->mvbuffer));
      cMat->mvbuffer = NULL;
    }
    //CHECK_CUDA(cudaMalloc((void **) &(cMat->mvbuffer), bfsz));
    allocRemoteBuffer((void **) &(cMat->mvbuffer), bfsz);
    cMat->mvbsize = bfsz;
  }
  CHECK_CUSPARSE(cusparseCsrmvEx(*my_handle,
				 CUSPARSE_ALG_MERGE_PATH,
				 CUSPARSE_OPERATION_NON_TRANSPOSE,
				 cMat->m,cMat->n,cMat->nz,	    
				 (const void *) &ealpha,CUSPARSE_BASE_TYPE,
				 cMat->descr,
				 (const void *) cMat->val, CUSPARSE_BASE_TYPE,
				 (const int *) cMat->irp, (const int *) cMat->ja,
				 (const void *) vX, CUSPARSE_BASE_TYPE,
				 (const void *) &ebeta, CUSPARSE_BASE_TYPE,
				 (void *) vY, CUSPARSE_BASE_TYPE,
				 CUSPARSE_BASE_TYPE, (void *) cMat->mvbuffer));

#else
  cusparseDnVecDescr_t vecX, vecY;
  size_t bfsz;

  if (T_CSRGIsNullMvDescr(cMat)) {
    cMat->spmvDescr = (cusparseSpMatDescr_t  *) malloc(sizeof(cusparseSpMatDescr_t  *));
  }
  T_CSRGCreateSpMVDescr(cMat);
  vX=x->v_;
  vY=y->v_;
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecY, cMat->m, vY, CUSPARSE_BASE_TYPE) );
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecX, cMat->n, vX, CUSPARSE_BASE_TYPE) );
  CHECK_CUSPARSE(cusparseSpMV_bufferSize(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
					 &alpha,(*(cMat->spmvDescr)),vecX,&beta,vecY,
					 CUSPARSE_BASE_TYPE,CUSPARSE_SPMV_ALG_DEFAULT,
					 &bfsz));
  if (bfsz > cMat->mvbsize) {
    if (cMat->mvbuffer != NULL) {
      CHECK_CUDA(cudaFree(cMat->mvbuffer));
      cMat->mvbuffer = NULL;
    }
    //CHECK_CUDA(cudaMalloc((void **) &(cMat->mvbuffer), bfsz));
    allocRemoteBuffer((void **) &(cMat->mvbuffer), bfsz);
    
    cMat->mvbsize = bfsz;
  }
  CHECK_CUSPARSE(cusparseSpMV(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
			      &alpha,(*(cMat->spmvDescr)),vecX,&beta,vecY,
			      CUSPARSE_BASE_TYPE,CUSPARSE_SPMV_ALG_DEFAULT,
			      cMat->mvbuffer));
  CHECK_CUSPARSE(cusparseDestroyDnVec(vecX) );
  CHECK_CUSPARSE(cusparseDestroyDnVec(vecY) );
  CHECK_CUSPARSE(cusparseDestroySpMat(*(cMat->spmvDescr)));
#endif
}

int T_spsvCSRGDevice(T_Cmat *Matrix, TYPE alpha, void *deviceX,
		     TYPE beta, void *deviceY)
{
  T_CSRGDeviceMat *cMat=Matrix->mat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY; 
  void *vX, *vY;  
  int r,n;
  cusparseHandle_t *my_handle=getHandle();
#if CUDA_SHORT_VERSION <= 10
  vX=x->v_;
  vY=y->v_;

  return cusparseTcsrsv_solve(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
			      cMat->m,(const TYPE *) &alpha,cMat->descr,
			      cMat->val, cMat->irp, cMat->ja, cMat->triang,
			      (const TYPE *) vX,  (TYPE *) vY);
#elif CUDA_VERSION <  11030
  vX=x->v_;
  vY=y->v_;
  CHECK_CUSPARSE(cusparseTcsrsv2_solve(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
				       cMat->m,cMat->nz,
				       (const TYPE *) &alpha,
				       cMat->descr,
				       cMat->val, cMat->irp, cMat->ja,
				       cMat->triang,
				       (const TYPE *) vX,  (TYPE *) vY,
				       CUSPARSE_SOLVE_POLICY_USE_LEVEL,
				       (void *) cMat->svbuffer));  
#else
  cusparseDnVecDescr_t vecX, vecY;
  size_t bfsz;
  vX=x->v_;
  vY=y->v_;
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecY, cMat->m, vY, CUSPARSE_BASE_TYPE) );
  CHECK_CUSPARSE( cusparseCreateDnVec(&vecX, cMat->n, vX, CUSPARSE_BASE_TYPE) );
  if (T_CSRGIsNullMvDescr(cMat)) {
    cMat->spmvDescr = (cusparseSpMatDescr_t  *) malloc(sizeof(cusparseSpMatDescr_t  *));
  }
  T_CSRGCreateSpMVDescr(cMat);
  //  fprintf(stderr,"Entry to SpSVDevice:   %d   %p\n",
  //	  T_CSRGIsNullSvDescr(cMat),cMat->spsvDescr);
  if (T_CSRGIsNullSvDescr(cMat)) {
    cMat->spsvDescr=(cusparseSpSVDescr_t  *) malloc(sizeof(cusparseSpSVDescr_t  *));
    cMat->svbsize=0;
    CHECK_CUSPARSE( cusparseSpSV_createDescr(cMat->spsvDescr) );
    //fprintf(stderr,"Entry to SpSVDevice:   %d   %p  %d\n",
    //	  T_CSRGIsNullSvDescr(cMat),cMat->spsvDescr,cMat->svbsize);
    CHECK_CUSPARSE(cusparseSpSV_bufferSize(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
					   &alpha,*(cMat->spmvDescr),vecX,vecY,
					   CUSPARSE_BASE_TYPE,
					   CUSPARSE_SPSV_ALG_DEFAULT,
					   *(cMat->spsvDescr),
					   &bfsz));
    if (bfsz > cMat->svbsize) {
      if (cMat->svbuffer != NULL) {
	CHECK_CUDA(cudaFree(cMat->svbuffer));
	cMat->svbuffer = NULL;
      }
      //CHECK_CUDA(cudaMalloc((void **) &(cMat->svbuffer), bfsz));
      allocRemoteBuffer((void **) &(cMat->svbuffer), bfsz);
    
      cMat->svbsize=bfsz;
      CHECK_CUSPARSE(cusparseSpSV_analysis(*my_handle,
					   CUSPARSE_OPERATION_NON_TRANSPOSE,
					   &alpha,
					   *(cMat->spmvDescr),
					   vecX, vecY,
					   CUSPARSE_BASE_TYPE,
					   CUSPARSE_SPSV_ALG_DEFAULT,
					   *(cMat->spsvDescr),
					   cMat->svbuffer));
    }
    if (T_CSRGIsNullSvBuffer(cMat)) {
      fprintf(stderr,"SpSV_SOLVE NULL spsv-buffer\n");
    }
  }
  CHECK_CUSPARSE(cusparseSpSV_solve(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
				    &alpha,*(cMat->spmvDescr),vecX,vecY,
				    CUSPARSE_BASE_TYPE,
				    CUSPARSE_SPSV_ALG_DEFAULT,
				    *(cMat->spsvDescr)));
  CHECK_CUSPARSE(cusparseDestroyDnVec(vecX) );
  CHECK_CUSPARSE(cusparseDestroyDnVec(vecY) );
  CHECK_CUSPARSE(cusparseDestroySpMat(*(cMat->spmvDescr)));
#endif
}

#if CUDA_VERSION >=  11030
T_CSRGCreateSpMVDescr(T_CSRGDeviceMat *cMat)
{
  int64_t tr,tc,tz;
  tr = cMat->m;
  tc = cMat->n;
  tz = cMat->nz;
  CHECK_CUSPARSE(cusparseCreateCsr(cMat->spmvDescr,
				   tr,tc,tz,
				   (void *) cMat->irp,
				   (void *) cMat->ja,
				   (void *) cMat->val,
				   CUSPARSE_INDEX_32I,
				   CUSPARSE_INDEX_32I,
				   CUSPARSE_INDEX_BASE_ONE,
				   CUSPARSE_BASE_TYPE) );
}
#endif
int T_CSRGDeviceAlloc(T_Cmat *Matrix,int nr, int nc, int nz)
{
  T_CSRGDeviceMat *cMat;
  int nr1=nr, nz1=nz, rc;
  cusparseHandle_t *my_handle=getHandle();
  int bfsz;

  if ((nr<0)||(nc<0)||(nz<0)) 
    return((int) CUSPARSE_STATUS_INVALID_VALUE);
  if ((cMat=(T_CSRGDeviceMat *) malloc(sizeof(T_CSRGDeviceMat)))==NULL)
    return((int) CUSPARSE_STATUS_ALLOC_FAILED);
  cMat->m  = nr;
  cMat->n  = nc;
  cMat->nz = nz;
  if (nr1 == 0) nr1 = 1;
  if (nz1 == 0) nz1 = 1;
  if ((rc= allocRemoteBuffer(((void **) &(cMat->irp)), ((nr1+1)*sizeof(int)))) != 0)
    return(rc);
  if ((rc= allocRemoteBuffer(((void **) &(cMat->ja)), ((nz1)*sizeof(int)))) != 0)
    return(rc);
  if ((rc= allocRemoteBuffer(((void **) &(cMat->val)), ((nz1)*sizeof(TYPE)))) != 0)
    return(rc);
#if CUDA_SHORT_VERSION <= 10  
  if ((rc= cusparseCreateMatDescr(&(cMat->descr))) !=0) 
    return(rc);
  if ((rc= cusparseCreateSolveAnalysisInfo(&(cMat->triang))) !=0)
    return(rc);
#elif CUDA_VERSION <  11030
  if ((rc= cusparseCreateMatDescr(&(cMat->descr))) !=0) 
    return(rc);
  CHECK_CUSPARSE(cusparseSetMatType(cMat->descr,CUSPARSE_MATRIX_TYPE_GENERAL));
  CHECK_CUSPARSE(cusparseSetMatDiagType(cMat->descr,CUSPARSE_DIAG_TYPE_NON_UNIT));
  CHECK_CUSPARSE(cusparseSetMatIndexBase(cMat->descr,CUSPARSE_INDEX_BASE_ONE));
  CHECK_CUSPARSE(cusparseCreateCsrsv2Info(&(cMat->triang)));
  if (cMat->nz > 0) {
    CHECK_CUSPARSE(cusparseTcsrsv2_bufferSize(*my_handle,
					      CUSPARSE_OPERATION_NON_TRANSPOSE,
					      cMat->m,cMat->nz, cMat->descr,
					      cMat->val, cMat->irp, cMat->ja,
					      cMat->triang, &bfsz));
  } else {
    bfsz = 0;
  }
    
  /* if (cMat->svbuffer != NULL) { */
  /*   fprintf(stderr,"Calling cudaFree\n"); */
  /*   CHECK_CUDA(cudaFree(cMat->svbuffer)); */
  /*   cMat->svbuffer = NULL; */
  /* } */
  if (bfsz > 0) {
    //CHECK_CUDA(cudaMalloc((void **) &(cMat->svbuffer), bfsz));
    allocRemoteBuffer((void **) &(cMat->svbuffer), bfsz);

  } else {
    cMat->svbuffer=NULL;
  }
  cMat->svbsize=bfsz;
  
  cMat->mvbuffer=NULL;
  cMat->mvbsize = 0;
  

#else

  cMat->spmvDescr=NULL;
  cMat->spsvDescr=NULL;
  cMat->mvbuffer=NULL;
  cMat->svbuffer=NULL;
  cMat->mvbsize=0;
  cMat->svbsize=0;
#endif  
  Matrix->mat = cMat;
  return(CUSPARSE_STATUS_SUCCESS);
}

int T_CSRGDeviceFree(T_Cmat *Matrix)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  
  if (cMat!=NULL) {
    freeRemoteBuffer(cMat->irp);
    freeRemoteBuffer(cMat->ja);
    freeRemoteBuffer(cMat->val);
#if CUDA_SHORT_VERSION <= 10  
    cusparseDestroyMatDescr(cMat->descr);
    cusparseDestroySolveAnalysisInfo(cMat->triang);
#elif CUDA_VERSION <  11030
    cusparseDestroyMatDescr(cMat->descr);
    cusparseDestroyCsrsv2Info(cMat->triang);
#else
    if (!T_CSRGIsNullMvDescr(cMat)) {
      // already destroyed spmvDescr, just free the pointer
      free(cMat->spmvDescr);
      cMat->spmvDescr=NULL;
    }
    if (cMat->mvbuffer!=NULL)
      CHECK_CUDA( cudaFree(cMat->mvbuffer));
    cMat->mvbuffer=NULL;
    cMat->mvbsize=0;   
    if (!T_CSRGIsNullSvDescr(cMat)) {
      CHECK_CUSPARSE(cusparseSpSV_destroyDescr(*(cMat->spsvDescr)));
      free(cMat->spsvDescr);
      cMat->spsvDescr=NULL;
    }
    if (cMat->svbuffer!=NULL)
      CHECK_CUDA( cudaFree(cMat->svbuffer));
    cMat->svbuffer=NULL;
    cMat->svbsize=0;
#endif
    free(cMat);
    Matrix->mat = NULL;
  }
  return(CUSPARSE_STATUS_SUCCESS);
}

int T_CSRGDeviceGetParms(T_Cmat *Matrix,int *nr, int *nc, int *nz)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  
  if (cMat!=NULL) {
    *nr = cMat->m   ;
    *nc = cMat->n   ;
    *nz = cMat->nz  ;
    return(CUSPARSE_STATUS_SUCCESS);
  } else {
    return((int) CUSPARSE_STATUS_ALLOC_FAILED);
  }
}

#if CUDA_SHORT_VERSION <= 10  

int T_CSRGDeviceSetMatType(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  return ((int) cusparseSetMatType(cMat->descr,type));
}

int T_CSRGDeviceSetMatFillMode(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  return ((int) cusparseSetMatFillMode(cMat->descr,type));
}

int T_CSRGDeviceSetMatDiagType(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  return ((int) cusparseSetMatDiagType(cMat->descr,type));
}

int T_CSRGDeviceSetMatIndexBase(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  return ((int) cusparseSetMatIndexBase(cMat->descr,type));
}

int T_CSRGDeviceCsrsmAnalysis(T_Cmat *Matrix)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;  
  int rc, buffersize;
  cusparseHandle_t *my_handle=getHandle();
  cusparseSolveAnalysisInfo_t info;

  rc= (int)  cusparseTcsrsv_analysis(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
				     cMat->m,cMat->nz,cMat->descr,
				     cMat->val, cMat->irp, cMat->ja,
				     cMat->triang);
  if (rc !=0) {
    fprintf(stderr,"From csrsv_analysis: %d\n",rc);
  }
  return(rc);
}

#elif CUDA_VERSION <  11030
int T_CSRGDeviceSetMatType(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  return ((int) cusparseSetMatType(cMat->descr,type));
}

int T_CSRGDeviceSetMatFillMode(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  return ((int) cusparseSetMatFillMode(cMat->descr,type));
}

int T_CSRGDeviceSetMatDiagType(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  return ((int) cusparseSetMatDiagType(cMat->descr,type));
}

int T_CSRGDeviceSetMatIndexBase(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  return ((int) cusparseSetMatIndexBase(cMat->descr,type));
}

#else

int T_CSRGDeviceSetMatFillMode(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  cusparseFillMode_t  mode=type;

  CHECK_CUSPARSE(cusparseSpMatSetAttribute(cMat->spmvDescr,
					   CUSPARSE_SPMAT_FILL_MODE,
					   (const void*) &mode,
					   sizeof(cusparseFillMode_t)));
 return(0);
}

int T_CSRGDeviceSetMatDiagType(T_Cmat *Matrix, int type)
{
  T_CSRGDeviceMat *cMat= Matrix->mat;
  cusparseDiagType_t  cutype=type;
  CHECK_CUSPARSE(cusparseSpMatSetAttribute(cMat->spmvDescr,
					   CUSPARSE_SPMAT_DIAG_TYPE,
					   (const void*) &cutype,
					   sizeof(cusparseDiagType_t)));
  return(0);
}

int T_CSRGIsNullMvDescr(T_CSRGDeviceMat *cMat)
{
  return(cMat->spmvDescr == NULL);
}

int T_CSRGIsNullSvBuffer(T_CSRGDeviceMat *cMat)
{
  return(cMat->svbuffer == NULL);
}
int T_CSRGIsNullSvDescr(T_CSRGDeviceMat *cMat)
{
  return(cMat->spsvDescr == NULL);
}

#endif

int T_CSRGHost2Device(T_Cmat *Matrix, int m, int n, int nz,
		      int *irp, int *ja, TYPE *val) 
{
  int rc;
  T_CSRGDeviceMat *cMat= Matrix->mat;
  cusparseHandle_t *my_handle=getHandle();
  
  if ((rc=writeRemoteBuffer((void *) irp, (void *) cMat->irp, 
			    (m+1)*sizeof(int)))
      != SPGPU_SUCCESS) 
    return(rc);
  
  if ((rc=writeRemoteBuffer((void *) ja,(void *) cMat->ja, 
			    (nz)*sizeof(int)))
      != SPGPU_SUCCESS) 
    return(rc);
  if ((rc=writeRemoteBuffer((void *) val, (void *) cMat->val, 
			    (nz)*sizeof(TYPE)))
      != SPGPU_SUCCESS) 
    return(rc);
#if (CUDA_SHORT_VERSION > 10  ) && (CUDA_VERSION <  11030)
  if (cusparseGetMatType(cMat->descr)== CUSPARSE_MATRIX_TYPE_TRIANGULAR) {
    // Why do we need to set TYPE_GENERAL??? cuSPARSE can be misterious sometimes. 
    cusparseSetMatType(cMat->descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    CHECK_CUSPARSE(cusparseTcsrsv2_analysis(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
					  cMat->m,cMat->nz,  cMat->descr,
					  cMat->val, cMat->irp, cMat->ja,
					  cMat->triang, CUSPARSE_SOLVE_POLICY_USE_LEVEL,
					  cMat->svbuffer));
  }
#else
  //cusparseSetMatType(*(cMat->spmvDescr),CUSPARSE_MATRIX_TYPE_GENERAL);
#endif
  return(CUSPARSE_STATUS_SUCCESS);
}

int T_CSRGDevice2Host(T_Cmat *Matrix, int m, int n, int nz,
		      int *irp, int *ja, TYPE *val) 
{
  int rc;
  T_CSRGDeviceMat *cMat = Matrix->mat;
  
  if ((rc=readRemoteBuffer((void *) irp, (void *) cMat->irp, (m+1)*sizeof(int))) 
      != SPGPU_SUCCESS) 
    return(rc);

  if ((rc=readRemoteBuffer((void *) ja, (void *) cMat->ja, (nz)*sizeof(int))) 
      != SPGPU_SUCCESS) 
    return(rc);
  if ((rc=readRemoteBuffer((void *) val, (void *) cMat->val, (nz)*sizeof(TYPE))) 
      != SPGPU_SUCCESS) 
    return(rc);

  return(CUSPARSE_STATUS_SUCCESS);
}

#if CUDA_SHORT_VERSION <= 10
int T_HYBGDeviceFree(T_Hmat *Matrix)
{
  T_HYBGDeviceMat *hMat= Matrix->mat;
  if (hMat != NULL) {
    cusparseDestroyMatDescr(hMat->descr);
    cusparseDestroySolveAnalysisInfo(hMat->triang);
    cusparseDestroyHybMat(hMat->hybA);
    free(hMat);
  }
  Matrix->mat = NULL;
  return(CUSPARSE_STATUS_SUCCESS);
}

int T_spmvHYBGDevice(T_Hmat *Matrix, TYPE alpha, void *deviceX,
		     TYPE beta, void *deviceY)
{
  T_HYBGDeviceMat *hMat=Matrix->mat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY; 
  void *vX, *vY;
  int r,n,rc;
  cusparseMatrixType_t type;
  cusparseHandle_t *my_handle=getHandle();

  /*getAddrMultiVecDevice(deviceX, &vX);
    getAddrMultiVecDevice(deviceY, &vY); */
  vX=x->v_;
  vY=y->v_;

  /* rc = (int) cusparseGetMatType(hMat->descr); */
  /* fprintf(stderr,"Spmv MatType: %d\n",rc); */
  /* rc = (int) cusparseGetMatDiagType(hMat->descr); */
  /* fprintf(stderr,"Spmv DiagType: %d\n",rc); */
  /* rc = (int) cusparseGetMatFillMode(hMat->descr); */
  /* fprintf(stderr,"Spmv FillMode: %d\n",rc); */
  /* Dirty trick: apparently hybmv does not accept a triangular
     matrix even though it should not make a difference. So 
     we claim it's general anyway */ 
  type =  cusparseGetMatType(hMat->descr);
  rc = cusparseSetMatType(hMat->descr,CUSPARSE_MATRIX_TYPE_GENERAL);
  if (rc == 0) 
    rc = (int) cusparseThybmv(*my_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
			      (const TYPE *) &alpha, hMat->descr, hMat->hybA, 
			      (const TYPE *) vX, (const TYPE *) &beta,
			      (TYPE *) vY);
  if (rc == 0) 
    rc = cusparseSetMatType(hMat->descr,type);
  return(rc);
}

int T_HYBGDeviceAlloc(T_Hmat *Matrix,int nr, int nc, int nz)
{
  T_HYBGDeviceMat *hMat;
  int nr1=nr, nz1=nz, rc;
  if ((nr<0)||(nc<0)||(nz<0)) 
    return((int) CUSPARSE_STATUS_INVALID_VALUE);
  if ((hMat=(T_HYBGDeviceMat *) malloc(sizeof(T_HYBGDeviceMat)))==NULL)
    return((int) CUSPARSE_STATUS_ALLOC_FAILED);
  hMat->m  = nr;
  hMat->n  = nc;
  hMat->nz = nz;

  if ((rc= cusparseCreateMatDescr(&(hMat->descr))) !=0) 
    return(rc);
  if ((rc= cusparseCreateSolveAnalysisInfo(&(hMat->triang))) !=0)
    return(rc);
  if((rc = cusparseCreateHybMat(&(hMat->hybA))) != 0)
    return(rc);
  Matrix->mat = hMat;
  return(CUSPARSE_STATUS_SUCCESS);
}

int T_HYBGDeviceSetMatDiagType(T_Hmat *Matrix, int type)
{
  T_HYBGDeviceMat *hMat= Matrix->mat;
  return ((int) cusparseSetMatDiagType(hMat->descr,type));
}

int T_HYBGDeviceSetMatIndexBase(T_Hmat *Matrix, int type)
{
  T_HYBGDeviceMat *hMat= Matrix->mat;
  return ((int) cusparseSetMatIndexBase(hMat->descr,type));
}

int T_HYBGDeviceSetMatType(T_Hmat *Matrix, int type)
{
  T_HYBGDeviceMat *hMat= Matrix->mat;
  return ((int) cusparseSetMatType(hMat->descr,type));
}

int T_HYBGDeviceSetMatFillMode(T_Hmat *Matrix, int type)
{
  T_HYBGDeviceMat *hMat= Matrix->mat;
  return ((int) cusparseSetMatFillMode(hMat->descr,type));
}

int T_spsvHYBGDevice(T_Hmat *Matrix, TYPE alpha, void *deviceX,
		     TYPE beta, void *deviceY)
{
  //beta??
  T_HYBGDeviceMat *hMat=Matrix->mat;
  struct MultiVectDevice *x = (struct MultiVectDevice *) deviceX;
  struct MultiVectDevice *y = (struct MultiVectDevice *) deviceY; 
  void *vX, *vY;
  int r,n;
  cusparseHandle_t *my_handle=getHandle();
  /*getAddrMultiVecDevice(deviceX, &vX);
    getAddrMultiVecDevice(deviceY, &vY); */
  vX=x->v_;
  vY=y->v_;

  return cusparseThybsv_solve(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
			      (const TYPE *) &alpha, hMat->descr,
			      hMat->hybA, hMat->triang,
			      (const TYPE *) vX,  (TYPE *) vY);
}

int T_HYBGDeviceHybsmAnalysis(T_Hmat *Matrix)
{
  T_HYBGDeviceMat *hMat= Matrix->mat;  
  cusparseSolveAnalysisInfo_t info;
  int rc;
  cusparseHandle_t *my_handle=getHandle();

  /* rc = (int) cusparseGetMatType(hMat->descr); */
  /* fprintf(stderr,"Analysis MatType: %d\n",rc); */
  /* rc = (int) cusparseGetMatDiagType(hMat->descr); */
  /* fprintf(stderr,"Analysis DiagType: %d\n",rc); */
  /* rc = (int) cusparseGetMatFillMode(hMat->descr); */
  /* fprintf(stderr,"Analysis FillMode: %d\n",rc); */
  rc = (int) cusparseThybsv_analysis(*my_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
				     hMat->descr, hMat->hybA, hMat->triang);

  if (rc !=0) {
    fprintf(stderr,"From csrsv_analysis: %d\n",rc);
  }
  return(rc);
}

int T_HYBGHost2Device(T_Hmat *Matrix, int m, int n, int nz,
		      int *irp, int *ja, TYPE *val) 
{
  int rc; double t1,t2;
  int nr1=m, nz1=nz;
  T_HYBGDeviceMat *hMat= Matrix->mat;
  cusparseHandle_t *my_handle=getHandle();

  if (nr1 == 0) nr1 = 1;
  if (nz1 == 0) nz1 = 1;
  if ((rc= allocRemoteBuffer(((void **) &(hMat->irp)), ((nr1+1)*sizeof(int)))) != 0)
    return(rc);
  if ((rc= allocRemoteBuffer(((void **) &(hMat->ja)), ((nz1)*sizeof(int)))) != 0)
    return(rc);
  if ((rc= allocRemoteBuffer(((void **) &(hMat->val)), ((nz1)*sizeof(TYPE)))) != 0)
    return(rc);

  if ((rc=writeRemoteBuffer((void *) irp, (void *) hMat->irp, 
			    (m+1)*sizeof(int)))
      != SPGPU_SUCCESS) 
    return(rc);
  
  if ((rc=writeRemoteBuffer((void *) ja,(void *) hMat->ja, 
			    (nz)*sizeof(int)))
      != SPGPU_SUCCESS) 
    return(rc);
  if ((rc=writeRemoteBuffer((void *) val, (void *) hMat->val, 
			    (nz)*sizeof(TYPE)))
      != SPGPU_SUCCESS) 
    return(rc);
  /* rc = (int) cusparseGetMatType(hMat->descr); */
  /* fprintf(stderr,"Conversion MatType: %d\n",rc); */
  /* rc = (int) cusparseGetMatDiagType(hMat->descr); */
  /* fprintf(stderr,"Conversion DiagType: %d\n",rc); */
  /* rc = (int) cusparseGetMatFillMode(hMat->descr); */
  /* fprintf(stderr,"Conversion FillMode: %d\n",rc); */
  //t1=etime();
  rc = (int) cusparseTcsr2hyb(*my_handle, m, n,
		   hMat->descr, 
		   (const TYPE *)hMat->val,
		   (const int *)hMat->irp, (const int *)hMat->ja, 
		   hMat->hybA,0,
		   CUSPARSE_HYB_PARTITION_AUTO);

  freeRemoteBuffer(hMat->irp);  hMat->irp = NULL;
  freeRemoteBuffer(hMat->ja);   hMat->ja  = NULL;
  freeRemoteBuffer(hMat->val);  hMat->val = NULL;

  //cudaSync();
  //t2 = etime();
  //fprintf(stderr,"Inner call to cusparseTcsr2hyb: %lf\n",(t2-t1));
  if (rc != 0) {
    fprintf(stderr,"From csr2hyb: %d\n",rc);
  }
  return(rc);
}
#endif

