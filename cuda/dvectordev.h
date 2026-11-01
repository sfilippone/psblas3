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
/* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;OR BUSINESS */
/* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN */
/* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) */
/* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#pragma once
//#include "utils.h"
#include "vectordev.h"
#include "cuda_runtime.h"
#include "core.h"
#include "vector.h"

int registerMappedDouble(void*, void**, int, double);

int checkMultiVecDeviceDouble(void* deviceMultiVec);

int writeMultiVecDeviceDouble(void* deviceMultiVec, double* hostMultiVec);
int writeMultiVecDeviceDoubleR2(void* deviceMultiVec, double* hostMultiVec, int ld);

int readMultiVecDeviceDouble(void* deviceMultiVec, double* hostMultiVec);
int readMultiVecDeviceDoubleR2(void* deviceMultiVec, double* hostMultiVec, int ld);

int setscalMultiVecDeviceDouble(double val, int first, int last, int indexBase, void* devVecX);

int geinsMultiVecDeviceDouble(int n, void* devVecIrl, void* devVecVal, int dupl, int indexBase, void* devVecX);

int igathMultiVecDeviceDouble(void* deviceVec, int vectorId, int n, int first, void* indexes, int hfirst, void* host_values, int indexBase);
int igathMultiVecDeviceDoubleVecIdx(void* deviceVec, int vectorId, int n, int first, void* deviceIdx, int hfirst, void* host_values, int indexBase);

int iscatMultiVecDeviceDouble(void* deviceVec, int vectorId, int n, int first, void *indexes, int hfirst, void* host_values, int indexBase, double beta);
int iscatMultiVecDeviceDoubleVecIdx(void* deviceVec, int vectorId, int n, int first, void *deviceIdx, int hfirst, void* host_values, int indexBase, double beta);

//AXPY-like operations
int axpbyMultiVecDeviceDouble1(int n, double alpha, void* devVecX, double beta, void* devVecY);                     //cudaVect, v_full, m_full
int axpbyMultiVecDeviceDouble2(int n, double alpha, void* devVecX, int idxX, double beta, void* devVecY);           //ext_col
int axpbyMultiVecDeviceDouble3(int n, double alpha, void* devVecX, double beta, void* devVecY, int idxY);           //v_idxs
int axpbyMultiVecDeviceDouble4(int n, double alpha, void* devVecX, int idxX, double beta, void* devVecY, int idxY); //m_idxs
int axpbyMultiVecDeviceDouble5(int n, double alpha, void* devVecX, double beta, void* devVecY, void* devVecZ);      //m_full_out

int axpbyczMultiVecDeviceDouble1(int n, double alpha, void* devVecX, double beta, void* devVecY, double gamma, void* devVecZ, int idxZ);                                                //axpbycz_vv
int axpbyczMultiVecDeviceDouble2(int n, double alpha, void* devVecX, double beta, void* devVecY, int idxY, double gamma, void* devVecZ, int idxZ);                                      //axpbycz_mv
int axpbyczMultiVecDeviceDouble3(int n, double alpha, void* devVecX, int idxX, double beta, void* devVecY, int idxY, double gamma, void* devVecZ, int idxZ);                            //axpbycz_mm_idxs
int axpbyczMultiVecDeviceDouble4(int n, double alpha, void* devVecX, double beta, void* devVecY, double gamma, void* devVecZ);                                                          //axpbycz_mm_full
int axpbyczMultiVecDeviceDouble5(int n, double alpha, void* devVecX, int idxX, double beta, void* devVecY, int idxY, double gamma, void* devVecZ, int idxZ, void* devVecW, int idxW);   //axpbycz_mm_fullout

int colspanDeviceDouble(int n, void* devVecX, CFI_cdesc_t* y_res, void* devVecY, bool updFlag);

int upd_xyzMultiVecDeviceDouble(int n, double alpha,double beta, double gamma, double delta, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ);
int xyzwMultiVecDeviceDouble(int n, double a, double b, double c, double d, double e, double f, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ, void* devMultiVecW);

//Scaling and norms
int scalMultiVecDeviceDouble(double alpha, void* devMultiVecA);
int nrm2MultiVecDeviceDoubleS(double* y_res, int n, void* devVecA);
int nrm2MultiVecDeviceDoubleI(double* y_res, int n, void* devVecA, int idx);
int nrm2MultiVecDeviceDoubleV(CFI_cdesc_t* y_res, int n, void* devVecA);
int amaxMultiVecDeviceDouble(double* y_res, int n, void* devVecA);
int asumMultiVecDeviceDouble(double* y_res, int n, void* devVecA);

//DOT-like operations
int dotMultiVecDeviceDoubleS(double* y_res, int n, void* devVecA, void* devVecB);
int dotMultiVecDeviceDoubleV(double* y_res, int n, int mA, void* devVecA, void* devVecB);
int dotMultiVecDeviceDoubleM(double* y_res, int n, int mA, int mB, void* devVecA, void* devVecB);
//DOT-like operations with CFI
int dotMultiVecDeviceDoubleV_CFI(CFI_cdesc_t* y_res, int n, void* devVecA, void* devVecB);
int dotMultiVecDeviceDoubleM_CFI(CFI_cdesc_t* y_res, int n, void* devVecA, void* devVecB);

//Elementwise operations
int mltMultiVecDeviceDouble1(int n, double alpha, void* devVecX, void* devVecY, double beta, bool conjgX, bool conjgY);                                                 //cudaVect, v_full, m_full
int mltMultiVecDeviceDouble2(int n, double alpha, void* devVecX, void* devVecY, double beta, int idxY, bool conjgX, bool conjgY);                                       //v_idxs
int mltMultiVecDeviceDouble3(int n, double alpha, void* devVecX, int idxX, void* devVecY, int idxY, double beta, bool conjgX, bool conjgY);                             //m_idxs
int mltMultiVecDeviceDouble4(int n, double alpha, void* devVecX, void* devVecY, double beta, void* devVecZ, bool conjgX, bool conjgY);                                  //vv_full_out, vm_full_out, mm_full_out
int mltMultiVecDeviceDouble5(int n, double alpha, void* devVecX, void* devVecY, double beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);                        //vv_idxs_out
int mltMultiVecDeviceDouble6(int n, double alpha, void* devVecX, void* devVecY, int idxY, double beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);              //vm_idxs_out
int mltMultiVecDeviceDouble7(int n, double alpha, void* devVecX, int idxX, void* devVecY, int idxY, double beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);    //mm_idxs_out
int mltMultiVecDeviceDouble8(int n, double alpha, void* devVecX, void* devVecY, int idxY, double beta, void* devVecZ, bool conjgX, bool conjgY);                        //vm_ext
int mltMultiVecDeviceDouble9(int n, double alpha, void* devVecX, int idxX, void* devVecY, int idxY, double beta, void* devVecZ, bool conjgX, bool conjgY);              //mm_ext

int axyMultiVecDeviceDouble(int n, double alpha, void *deviceVecA, void *deviceVecB);
int axybzMultiVecDeviceDouble(int n, double alpha, void *deviceVecA, void *deviceVecB, double beta, void *deviceVecZ);

int absMultiVecDeviceDouble(int n, double alpha, void *deviceVecA);
int absMultiVecDeviceDouble2(int n, double alpha, void *deviceVecA, void *deviceVecB);