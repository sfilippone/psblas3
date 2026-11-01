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

int registerMappedFloat(void*, void**, int, float);

int checkMultiVecDeviceFloat(void* deviceMultiVec);

int writeMultiVecDeviceFloat(void* deviceMultiVec, float* hostMultiVec);
int writeMultiVecDeviceFloatR2(void* deviceMultiVec, float* hostMultiVec, int ld);

int readMultiVecDeviceFloat(void* deviceMultiVec, float* hostMultiVec);
int readMultiVecDeviceFloatR2(void* deviceMultiVec, float* hostMultiVec, int ld);

int setscalMultiVecDeviceFloat(float val, int first, int last, int indexBase, void* devVecX);

int geinsMultiVecDeviceFloat(int n, void* devVecIrl, void* devVecVal, int dupl, int indexBase, void* devVecX);

int igathMultiVecDeviceFloat(void* deviceVec, int vectorId, int n, int first, void* indexes, int hfirst, void* host_values, int indexBase);
int igathMultiVecDeviceFloatVecIdx(void* deviceVec, int vectorId, int n, int first, void* deviceIdx, int hfirst, void* host_values, int indexBase);

int iscatMultiVecDeviceFloat(void* deviceVec, int vectorId, int n, int first, void *indexes, int hfirst, void* host_values, int indexBase, float beta);
int iscatMultiVecDeviceFloatVecIdx(void* deviceVec, int vectorId, int n, int first, void *deviceIdx, int hfirst, void* host_values, int indexBase, float beta);

//AXPY-like operations
int axpbyMultiVecDeviceFloat1(int n, float alpha, void* devVecX, float beta, void* devVecY);                     //cudaVect, v_full, m_full
int axpbyMultiVecDeviceFloat2(int n, float alpha, void* devVecX, int idxX, float beta, void* devVecY);           //ext_col
int axpbyMultiVecDeviceFloat3(int n, float alpha, void* devVecX, float beta, void* devVecY, int idxY);           //v_idxs
int axpbyMultiVecDeviceFloat4(int n, float alpha, void* devVecX, int idxX, float beta, void* devVecY, int idxY); //m_idxs
int axpbyMultiVecDeviceFloat5(int n, float alpha, void* devVecX, float beta, void* devVecY, void* devVecZ);      //m_full_out

int axpbyczMultiVecDeviceFloat1(int n, float alpha, void* devVecX, float beta, void* devVecY, float gamma, void* devVecZ, int idxZ);                                                //axpbycz_vv
int axpbyczMultiVecDeviceFloat2(int n, float alpha, void* devVecX, float beta, void* devVecY, int idxY, float gamma, void* devVecZ, int idxZ);                                      //axpbycz_mv
int axpbyczMultiVecDeviceFloat3(int n, float alpha, void* devVecX, int idxX, float beta, void* devVecY, int idxY, float gamma, void* devVecZ, int idxZ);                            //axpbycz_mm_idxs
int axpbyczMultiVecDeviceFloat4(int n, float alpha, void* devVecX, float beta, void* devVecY, float gamma, void* devVecZ);                                                          //axpbycz_mm_full
int axpbyczMultiVecDeviceFloat5(int n, float alpha, void* devVecX, int idxX, float beta, void* devVecY, int idxY, float gamma, void* devVecZ, int idxZ, void* devVecW, int idxW);   //axpbycz_mm_fullout

int colspanDeviceFloat(int n, void* devVecX, CFI_cdesc_t* y_res, void* devVecY, bool updFlag);

int upd_xyzMultiVecDeviceFloat(int n, float alpha,float beta, float gamma, float delta, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ);
int xyzwMultiVecDeviceFloat(int n, float a, float b, float c, float d, float e, float f, void* devMultiVecX, void* devMultiVecY, void* devMultiVecZ, void* devMultiVecW);

//Scaling and norms
int scalMultiVecDeviceFloat(float alpha, void* devMultiVecA);
int nrm2MultiVecDeviceFloatS(float* y_res, int n, void* devVecA);
int nrm2MultiVecDeviceFloatI(float* y_res, int n, void* devVecA, int idx);
int nrm2MultiVecDeviceFloatV(CFI_cdesc_t* y_res, int n, void* devVecA);
int amaxMultiVecDeviceFloat(float* y_res, int n, void* devVecA);
int asumMultiVecDeviceFloat(float* y_res, int n, void* devVecA);

//DOT-like operations
int dotMultiVecDeviceFloatS(float* y_res, int n, void* devVecA, void* devVecB);
int dotMultiVecDeviceFloatV(float* y_res, int n, int mA, void* devVecA, void* devVecB);
int dotMultiVecDeviceFloatM(float* y_res, int n, int mA, int mB, void* devVecA, void* devVecB);
//DOT-like operations with CFI
int dotMultiVecDeviceFloatV_CFI(CFI_cdesc_t* y_res, int n, void* devVecA, void* devVecB);
int dotMultiVecDeviceFloatM_CFI(CFI_cdesc_t* y_res, int n, void* devVecA, void* devVecB);

//Elementwise operations
int mltMultiVecDeviceFloat1(int n, float alpha, void* devVecX, void* devVecY, float beta, bool conjgX, bool conjgY);                                                 //cudaVect, v_full, m_full
int mltMultiVecDeviceFloat2(int n, float alpha, void* devVecX, void* devVecY, float beta, int idxY, bool conjgX, bool conjgY);                                       //v_idxs
int mltMultiVecDeviceFloat3(int n, float alpha, void* devVecX, int idxX, void* devVecY, int idxY, float beta, bool conjgX, bool conjgY);                             //m_idxs
int mltMultiVecDeviceFloat4(int n, float alpha, void* devVecX, void* devVecY, float beta, void* devVecZ, bool conjgX, bool conjgY);                                  //vv_full_out, vm_full_out, mm_full_out
int mltMultiVecDeviceFloat5(int n, float alpha, void* devVecX, void* devVecY, float beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);                        //vv_idxs_out
int mltMultiVecDeviceFloat6(int n, float alpha, void* devVecX, void* devVecY, int idxY, float beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);              //vm_idxs_out
int mltMultiVecDeviceFloat7(int n, float alpha, void* devVecX, int idxX, void* devVecY, int idxY, float beta, void* devVecZ, int idxZ, bool conjgX, bool conjgY);    //mm_idxs_out
int mltMultiVecDeviceFloat8(int n, float alpha, void* devVecX, void* devVecY, int idxY, float beta, void* devVecZ, bool conjgX, bool conjgY);                        //vm_ext
int mltMultiVecDeviceFloat9(int n, float alpha, void* devVecX, int idxX, void* devVecY, int idxY, float beta, void* devVecZ, bool conjgX, bool conjgY);              //mm_ext

int axyMultiVecDeviceFloat(int n, float alpha, void *deviceVecA, void *deviceVecB);
int axybzMultiVecDeviceFloat(int n, float alpha, void *deviceVecA, void *deviceVecB, float beta, void *deviceVecZ);

int absMultiVecDeviceFloat(int n, float alpha, void *deviceVecA);
int absMultiVecDeviceFloat2(int n, float alpha, void *deviceVecA, void *deviceVecB);