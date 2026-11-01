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

int registerMappedFloatComplex(void*, void**, int, cuFloatComplex);

int checkMultiVecDeviceFloatComplex(void* deviceMultiVec);

int writeMultiVecDeviceFloatComplex(void* deviceMultiVec, cuFloatComplex* hostMultiVec);
int writeMultiVecDeviceFloatComplexR2(void* deviceMultiVec, cuFloatComplex* hostMultiVec, int ld);

int readMultiVecDeviceFloatComplex(void* deviceMultiVec, cuFloatComplex* hostMultiVec);
int readMultiVecDeviceFloatComplexR2(void* deviceMultiVec, cuFloatComplex* hostMultiVec, int ld);

int setscalMultiVecDeviceFloatComplex(cuFloatComplex val, int first, int last, int indexBase, void* devVecX); 

int geinsMultiVecDeviceFloatComplex(int n, void* devVecIrl, void* devVecVal, int dupl, int indexBase, void* devVecX); 

int igathMultiVecDeviceFloatComplex(void* deviceVec, int vectorId, int n, int first, void* indexes, int hfirst, void* host_values, int indexBase);
int igathMultiVecDeviceFloatComplexVecIdx(void* deviceVec, int vectorId, int n, int first, void* deviceIdx, int hfirst, void* host_values, int indexBase);

int iscatMultiVecDeviceFloatComplex(void* deviceVec, int vectorId, int n, int first, void *indexes, int hfirst, void* host_values, int indexBase, cuFloatComplex beta);
int iscatMultiVecDeviceFloatComplexVecIdx(void* deviceVec, int vectorId, int n, int first,  void *deviceIdx, int hfirst, void* host_values, int indexBase, cuFloatComplex beta);

//AXPY-like operations
int axpbyMultiVecDeviceFloatComplex1(int n, cuFloatComplex alpha, void* devVecX, cuFloatComplex beta, void* devVecY);                     //cudaVect, v_full, m_full
int axpbyMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void* devVecX, int idxX, cuFloatComplex beta, void* devVecY);           //ext_col
int axpbyMultiVecDeviceFloatComplex3(int n, cuFloatComplex alpha, void* devVecX, cuFloatComplex beta, void* devVecY, int idxY);           //v_idxs
int axpbyMultiVecDeviceFloatComplex4(int n, cuFloatComplex alpha, void* devVecX, int idxX, cuFloatComplex beta, void* devVecY, int idxY); //m_idxs
int axpbyMultiVecDeviceFloatComplex5(int n, cuFloatComplex alpha, void* devVecX, cuFloatComplex beta, void* devVecY, void* devVecZ);      //m_full_out

int axpbyczMultiVecDeviceFloatComplex1(int n, cuFloatComplex alpha, void* devVecX, cuFloatComplex beta, void* devVecY, cuFloatComplex gamma, void* devVecZ, int idxZ);                                              //axpbycz_vv
int axpbyczMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void* devVecX, cuFloatComplex beta, void* devVecY, int idxY, cuFloatComplex gamma, void* devVecZ, int idxZ);                                    //axpbycz_mv
int axpbyczMultiVecDeviceFloatComplex3(int n, cuFloatComplex alpha, void* devVecX, int idxX, cuFloatComplex beta, void* devVecY, int idxY, cuFloatComplex gamma, void* devVecZ, int idxZ);                          //axpbycz_mm_idxs
int axpbyczMultiVecDeviceFloatComplex4(int n, cuFloatComplex alpha, void* devVecX, cuFloatComplex beta, void* devVecY, cuFloatComplex gamma, void* devVecZ);                                                        //axpbycz_mm_full
int axpbyczMultiVecDeviceFloatComplex5(int n, cuFloatComplex alpha, void* devVecX, int idxX, cuFloatComplex beta, void* devVecY, int idxY, cuFloatComplex gamma, void* devVecZ, int idxZ, void* devVecW, int idxW); //axpbycz_mm_fullout

int colspanDeviceFloatComplex(int n, void* devVecX, CFI_cdesc_t* y_res, void* devVecY, bool updFlag);

int upd_xyzMultiVecDeviceFloatComplex(int n,cuFloatComplex alpha,cuFloatComplex beta, cuFloatComplex gamma, cuFloatComplex delta, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ);
int xyzwMultiVecDeviceFloatComplex(int n, cuFloatComplex a, cuFloatComplex b, cuFloatComplex c, cuFloatComplex d, cuFloatComplex e, cuFloatComplex f, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ, void* devMultiVecW);

//Scaling and norms
int scalMultiVecDeviceFloatComplex(cuFloatComplex alpha, void* devMultiVecA);
int nrm2MultiVecDeviceFloatComplexS(float* y_res, int n, void* devVecA);
int nrm2MultiVecDeviceFloatComplexI(float* y_res, int n, void* devVecA, int idx);
int nrm2MultiVecDeviceFloatComplexV(CFI_cdesc_t* y_res, int n, void* devVecA);
int amaxMultiVecDeviceFloatComplex(float* y_res, int n, void* devVecA);
int asumMultiVecDeviceFloatComplex(float* y_res, int n, void* devVecA);

//DOT-like operations
int dotMultiVecDeviceFloatComplexS(cuFloatComplex* y_res, int n, void* devVecA, void* devVecB);
int dotMultiVecDeviceFloatComplexV(cuFloatComplex* y_res, int n, int mA, void* devVecA, void* devVecB);
int dotMultiVecDeviceFloatComplexM(cuFloatComplex* y_res, int n, int mA, int mB, void* devVecA, void* devVecB);
//DOT-like operations with CFI
int dotMultiVecDeviceFloatComplexV_CFI(CFI_cdesc_t* y_res, int n, void* devVecA, void* devVecB);
int dotMultiVecDeviceFloatComplexM_CFI(CFI_cdesc_t* y_res, int n, void* devVecA, void* devVecB);

//Elementwise operations
int mltMultiVecDeviceFloatComplex1(int n, cuFloatComplex alpha, void* devVecX, void* devVecY, cuFloatComplex beta, bool conjgX, bool conjgY);                                               //cudaVect, v_full, m_full
int mltMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void* devVecX, void* devVecY, cuFloatComplex beta, int idxY, bool conjgX, bool conjgY);                                     //v_idxs
int mltMultiVecDeviceFloatComplex3(int n, cuFloatComplex alpha, void* devVecX, int idxX, void* devVecY, int idxY, cuFloatComplex beta, bool conjgX, bool conjgY);                           //m_idxs
int mltMultiVecDeviceFloatComplex4(int n, cuFloatComplex alpha, void* devVecX, void* devVecY, cuFloatComplex beta, void* devVecZ, bool conjgX, bool conjgY);                                //vv_full_out, vm_full_out, mm_full_out
int mltMultiVecDeviceFloatComplex5(int n, cuFloatComplex alpha, void* devVecX, void* devVecY, cuFloatComplex beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);                      //vv_idxs_out
int mltMultiVecDeviceFloatComplex6(int n, cuFloatComplex alpha, void* devVecX, void* devVecY, int idxY, cuFloatComplex beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);            //vm_idxs_out
int mltMultiVecDeviceFloatComplex7(int n, cuFloatComplex alpha, void* devVecX, int idxX, void* devVecY, int idxY, cuFloatComplex beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);  //mm_idxs_out
int mltMultiVecDeviceFloatComplex8(int n, cuFloatComplex alpha, void* devVecX, void* devVecY, int idxY, cuFloatComplex beta, void* devVecZ, bool conjgX, bool conjgY);                      //vm_ext
int mltMultiVecDeviceFloatComplex9(int n, cuFloatComplex alpha, void* devVecX, int idxX, void* devVecY, int idxY, cuFloatComplex beta, void* devVecZ, bool conjgX, bool conjgY);            //mm_ext

int axyMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, void *deviceVecA, void *deviceVecB);
int axybzMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, void *deviceVecA, void *deviceVecB, cuFloatComplex beta, void *deviceVecZ);

int absMultiVecDeviceFloatComplex(int n, cuFloatComplex alpha, void *deviceVecA);
int absMultiVecDeviceFloatComplex2(int n, cuFloatComplex alpha, void *deviceVecA, void *deviceVecB);