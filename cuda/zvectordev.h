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
/*      software without specific prior written permission. */

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

#pragma once
//#include "utils.h"
#include <complex.h>
#include "cuComplex.h"
#include "vectordev.h"
#include "cuda_runtime.h"
#include "core.h"
#include "vector.h"

int registerMappedDoubleComplex(void*, void**, int, cuDoubleComplex);

int checkMultiVecDeviceDoubleComplex(void* deviceMultiVec);

int writeMultiVecDeviceDoubleComplex(void* deviceMultiVec, cuDoubleComplex* hostMultiVec);
int writeMultiVecDeviceDoubleComplexR2(void* deviceMultiVec, cuDoubleComplex* hostMultiVec, int ld);

int readMultiVecDeviceDoubleComplex(void* deviceMultiVec, cuDoubleComplex* hostMultiVec);
int readMultiVecDeviceDoubleComplexR2(void* deviceMultiVec, cuDoubleComplex* hostMultiVec, int ld);

int setscalMultiVecDeviceDoubleComplex(cuDoubleComplex val, int first, int last, int indexBase, void* devVecX); 

int geinsMultiVecDeviceDoubleComplex(int n, void* devVecIrl, void* devVecVal, int dupl, int indexBase, void* devVecX); 

int igathMultiVecDeviceDoubleComplex(void* deviceVec, int vectorId, int n, int first, void* indexes, int hfirst, void* host_values, int indexBase);
int igathMultiVecDeviceDoubleComplexVecIdx(void* deviceVec, int vectorId, int n, int first, void* deviceIdx, int hfirst, void* host_values, int indexBase);

int iscatMultiVecDeviceDoubleComplex(void* deviceVec, int vectorId, int n, int first, void *indexes, int hfirst, void* host_values, int indexBase, cuDoubleComplex beta);
int iscatMultiVecDeviceDoubleComplexVecIdx(void* deviceVec, int vectorId, int n, int first,  void *deviceIdx, int hfirst, void* host_values, int indexBase, cuDoubleComplex beta);

//AXPY-like operations
int axpbyMultiVecDeviceDoubleComplex1(int n, cuDoubleComplex alpha, void* devVecX, cuDoubleComplex beta, void* devVecY);                     //cudaVect, v_full, m_full
int axpbyMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void* devVecX, int idxX, cuDoubleComplex beta, void* devVecY);           //ext_col
int axpbyMultiVecDeviceDoubleComplex3(int n, cuDoubleComplex alpha, void* devVecX, cuDoubleComplex beta, void* devVecY, int idxY);           //v_idxs
int axpbyMultiVecDeviceDoubleComplex4(int n, cuDoubleComplex alpha, void* devVecX, int idxX, cuDoubleComplex beta, void* devVecY, int idxY); //m_idxs
int axpbyMultiVecDeviceDoubleComplex5(int n, cuDoubleComplex alpha, void* devVecX, cuDoubleComplex beta, void* devVecY, void* devVecZ);      //m_full_out

int axpbyczMultiVecDeviceDoubleComplex1(int n, cuDoubleComplex alpha, void* devVecX, cuDoubleComplex beta, void* devVecY, cuDoubleComplex gamma, void* devVecZ, int idxZ);                                              //axpbycz_vv
int axpbyczMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void* devVecX, cuDoubleComplex beta, void* devVecY, int idxY, cuDoubleComplex gamma, void* devVecZ, int idxZ);                                    //axpbycz_mv
int axpbyczMultiVecDeviceDoubleComplex3(int n, cuDoubleComplex alpha, void* devVecX, int idxX, cuDoubleComplex beta, void* devVecY, int idxY, cuDoubleComplex gamma, void* devVecZ, int idxZ);                          //axpbycz_mm_idxs
int axpbyczMultiVecDeviceDoubleComplex4(int n, cuDoubleComplex alpha, void* devVecX, cuDoubleComplex beta, void* devVecY, cuDoubleComplex gamma, void* devVecZ);                                                        //axpbycz_mm_full
int axpbyczMultiVecDeviceDoubleComplex5(int n, cuDoubleComplex alpha, void* devVecX, int idxX, cuDoubleComplex beta, void* devVecY, int idxY, cuDoubleComplex gamma, void* devVecZ, int idxZ, void* devVecW, int idxW); //axpbycz_mm_fullout

int colspanDeviceDoubleComplex(int n, void* devVecX, CFI_cdesc_t* y_res, void* devVecY, bool updFlag);

int upd_xyzMultiVecDeviceDoubleComplex(int n,cuDoubleComplex alpha,cuDoubleComplex beta, cuDoubleComplex gamma, cuDoubleComplex delta, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ);
int xyzwMultiVecDeviceDoubleComplex(int n, cuDoubleComplex a, cuDoubleComplex b, cuDoubleComplex c, cuDoubleComplex d, cuDoubleComplex e, cuDoubleComplex f, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ, void* devMultiVecW);

//Scaling and norms
int scalMultiVecDeviceDoubleComplex(cuDoubleComplex alpha, void* devMultiVecA);
int nrm2MultiVecDeviceDoubleComplexS(double* y_res, int n, void* devVecA);
int nrm2MultiVecDeviceDoubleComplexI(double* y_res, int n, void* devVecA, int idx);
int nrm2MultiVecDeviceDoubleComplexV(CFI_cdesc_t* y_res, int n, void* devVecA);
int amaxMultiVecDeviceDoubleComplex(double* y_res, int n, void* devVecA);
int asumMultiVecDeviceDoubleComplex(double* y_res, int n, void* devVecA);

//DOT-like operations
int dotMultiVecDeviceDoubleComplexS(cuDoubleComplex* y_res, int n, void* devVecA, void* devVecB);
int dotMultiVecDeviceDoubleComplexV(cuDoubleComplex* y_res, int n, int mA, void* devVecA, void* devVecB);
int dotMultiVecDeviceDoubleComplexM(cuDoubleComplex* y_res, int n, int mA, int mB, void* devVecA, void* devVecB);
//DOT-like operations with CFI
int dotMultiVecDeviceDoubleComplexV_CFI(CFI_cdesc_t* y_res, int n, void* devVecA, void* devVecB);
int dotMultiVecDeviceDoubleComplexM_CFI(CFI_cdesc_t* y_res, int n, void* devVecA, void* devVecB);

//Elementwise operations
int mltMultiVecDeviceDoubleComplex1(int n, cuDoubleComplex alpha, void* devVecX, void* devVecY, cuDoubleComplex beta, bool conjgX, bool conjgY);                                               //cudaVect, v_full, m_full
int mltMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void* devVecX, void* devVecY, cuDoubleComplex beta, int idxY, bool conjgX, bool conjgY);                                     //v_idxs
int mltMultiVecDeviceDoubleComplex3(int n, cuDoubleComplex alpha, void* devVecX, int idxX, void* devVecY, int idxY, cuDoubleComplex beta, bool conjgX, bool conjgY);                           //m_idxs
int mltMultiVecDeviceDoubleComplex4(int n, cuDoubleComplex alpha, void* devVecX, void* devVecY, cuDoubleComplex beta, void* devVecZ, bool conjgX, bool conjgY);                                //vv_full_out, vm_full_out, mm_full_out
int mltMultiVecDeviceDoubleComplex5(int n, cuDoubleComplex alpha, void* devVecX, void* devVecY, cuDoubleComplex beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);                      //vv_idxs_out
int mltMultiVecDeviceDoubleComplex6(int n, cuDoubleComplex alpha, void* devVecX, void* devVecY, int idxY, cuDoubleComplex beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);            //vm_idxs_out
int mltMultiVecDeviceDoubleComplex7(int n, cuDoubleComplex alpha, void* devVecX, int idxX, void* devVecY, int idxY, cuDoubleComplex beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);  //mm_idxs_out
int mltMultiVecDeviceDoubleComplex8(int n, cuDoubleComplex alpha, void* devVecX, void* devVecY, int idxY, cuDoubleComplex beta, void* devVecZ, bool conjgX, bool conjgY);                      //vm_ext
int mltMultiVecDeviceDoubleComplex9(int n, cuDoubleComplex alpha, void* devVecX, int idxX, void* devVecY, int idxY, cuDoubleComplex beta, void* devVecZ, bool conjgX, bool conjgY);            //mm_ext

int axyMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA, void *deviceVecB);
int axybzMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA, void *deviceVecB, cuDoubleComplex beta, void *deviceVecZ);

int absMultiVecDeviceDoubleComplex(int n, cuDoubleComplex alpha, void *deviceVecA);
int absMultiVecDeviceDoubleComplex2(int n, cuDoubleComplex alpha, void *deviceVecA, void *deviceVecB);