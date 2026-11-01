!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific prior written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
module psb_d_vectordev_mod
  use psb_base_vectordev_mod
 
  interface registerMapped
    function registerMappedDouble(buf, d_p, n, dummy) result(res) &
        & bind(c, name = 'registerMappedDouble')
      use iso_c_binding
      type(c_ptr), value    :: buf
      type(c_ptr)           :: d_p
      integer(c_int), value :: n
      real(c_double), value :: dummy
      integer(c_int) :: res
    end function registerMappedDouble
  end interface

  interface 
    function checkMultiVecDeviceDouble(deviceVec) result(res) &
        & bind(c, name = 'checkMultiVecDeviceDouble')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      integer(c_int)  :: res
    end function checkMultiVecDeviceDouble
  end interface
  
  interface writeMultiVecDevice 
    function writeMultiVecDeviceDouble(deviceVec, hostVec) result(res) &
        & bind(c, name = 'writeMultiVecDeviceDouble')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      real(c_double)      :: hostVec(*)
      integer(c_int)  :: res
    end function writeMultiVecDeviceDouble

    function writeMultiVecDeviceDoubleR2(deviceVec, hostVec, ld) result(res) &
        & bind(c, name = 'writeMultiVecDeviceDoubleR2')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec
      real(c_double)        :: hostVec(ld, *)
      integer(c_int), value :: ld
      integer(c_int)  :: res
    end function writeMultiVecDeviceDoubleR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceDouble(deviceVec, hostVec) result(res) &
        & bind(c, name = 'readMultiVecDeviceDouble')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      real(c_double)      :: hostVec(*)
      integer(c_int)  :: res
    end function readMultiVecDeviceDouble

    function readMultiVecDeviceDoubleR2(deviceVec, hostVec, ld) result(res) &
        & bind(c, name = 'readMultiVecDeviceDoubleR2')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec
      real(c_double)        :: hostVec(ld, *)
      integer(c_int), value :: ld
      integer(c_int)  :: res
    end function readMultiVecDeviceDoubleR2
  end interface 

  interface allocateDouble
    function allocateDouble(didx, n) result(res) &
        & bind(c, name = 'allocateDouble') 
      use iso_c_binding
      type(c_ptr)           :: didx
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function allocateDouble
    
    function allocateMultiDouble(didx, m, n) result(res) &
        & bind(c, name = 'allocateMultiDouble') 
      use iso_c_binding
      type(c_ptr)           :: didx
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function allocateMultiDouble
  end interface

  interface writeDouble
    function writeDouble(didx, hidx, n) result(res) &
        & bind(c, name = 'writeDouble')
      use iso_c_binding
      type(c_ptr), value    :: didx
      real(c_double)        :: hidx(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function writeDouble

    function writeDoubleFirst(first, didx, hidx, n, indexBase) result(res) &
        & bind(c, name = 'writeDoubleFirst')
      use iso_c_binding
      integer(c_int), value :: first, n, indexBase
      type(c_ptr), value    :: didx
      real(c_double)        :: hidx(*)
      integer(c_int)  :: res
    end function writeDoubleFirst

    function writeMultiDouble(didx, hidx, m, n) result(res) &
        & bind(c, name = 'writeMultiDouble')
      use iso_c_binding
      type(c_ptr), value    :: didx
      real(c_double)        :: hidx(m, *)
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function writeMultiDouble
  end interface
  
  interface readDouble
    function readDouble(didx, hidx, n) result(res) &
        & bind(c, name = 'readDouble')
      use iso_c_binding
      type(c_ptr), value    :: didx
      real(c_double)        :: hidx(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function readDouble

    function readDoubleFirst(first, didx, hidx, n, indexBase) result(res) &
        & bind(c, name = 'readDoubleFirst')
      use iso_c_binding
      integer(c_int), value :: first, n, indexBase
      type(c_ptr), value    :: didx
      real(c_double)        :: hidx(*)
      integer(c_int)  :: res
    end function readDoubleFirst

    function readMultiDouble(didx, hidx, m, n) result(res) &
        & bind(c, name = 'readMultiDouble')
      use iso_c_binding
      type(c_ptr), value    :: didx
      real(c_double)        :: hidx(m, *)
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function readMultiDouble
  end interface
  
  interface
    subroutine freeDouble(didx) &
        & bind(c, name = 'freeDouble')
      use iso_c_binding
      type(c_ptr), value  :: didx
    end subroutine freeDouble
  end interface
  
  interface setScalDevice
    function setScalMultiVecDeviceDouble(val, first, last, indexBase, deviceVecX) result(res) &
        & bind(c, name = 'setscalMultiVecDeviceDouble')
      use iso_c_binding
      real(c_double), value :: val
      integer(c_int), value :: first, last, indexBase
      type(c_ptr),  value   :: deviceVecX
      integer(c_int)  :: res
    end function setScalMultiVecDeviceDouble
  end interface

  interface 
    function geinsMultiVecDeviceDouble(n, deviceVecIrl, deviceVecVal, dupl, indexBase, deviceVecX) result(res) &
        & bind(c, name = 'geinsMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int), value :: n, dupl, indexBase
      type(c_ptr), value    :: deviceVecIrl, deviceVecVal, deviceVecX
      integer(c_int)  :: res
    end function geinsMultiVecDeviceDouble
  end interface

  ! New gather functions
  interface 
    function igathMultiVecDeviceDouble(deviceVec, vectorId, n, first, idx, hfirst, hostVec, indexBase) result(res) &
        & bind(c, name = 'igathMultiVecDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, n, first, hfirst, indexBase
      integer(c_int)  :: res
    end function igathMultiVecDeviceDouble
  end interface

  interface 
    function igathMultiVecDeviceDoubleVecIdx(deviceVec, vectorId, n, first, idx, hfirst, hostVec, indexBase) result(res) &
        & bind(c, name = 'igathMultiVecDeviceDoubleVecIdx')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      integer(c_int)  :: res
    end function igathMultiVecDeviceDoubleVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceDouble(deviceVec, vectorId, first, n, idx, &
                                      & hfirst, hostVec, indexBase, beta) result(res) & 
        & bind(c, name = 'iscatMultiVecDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      real(c_double), value :: beta
      integer(c_int)  :: res
    end function iscatMultiVecDeviceDouble
  end interface

  interface 
    function iscatMultiVecDeviceDoubleVecIdx(deviceVec, vectorId, first, n, idx, &
                                      & hfirst, hostVec, indexBase, beta) result(res) &
        & bind(c, name = 'iscatMultiVecDeviceDoubleVecIdx')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      real(c_double), value :: beta
      integer(c_int)  :: res
    end function iscatMultiVecDeviceDoubleVecIdx
  end interface

  interface scalMultiVecDevice
    function scalMultiVecDeviceDouble(alpha, deviceVecA) result(val) &
        & bind(c, name = 'scalMultiVecDeviceDouble')
      use iso_c_binding
      real(c_double), value :: alpha
      type(c_ptr), value    :: deviceVecA
      integer(c_int)        :: val
    end function scalMultiVecDeviceDouble
  end interface

  interface dotMultiVecDevice
    function dotMultiVecDeviceDoubleS(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleS')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleS
    
    function dotMultiVecDeviceDoubleV(res, n, mA, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleV')
      use iso_c_binding
      real(c_double)        :: res(mA)
      integer(c_int), value :: n, mA
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleV

    function dotMultiVecDeviceDoubleM(res, n, mA, mB, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleM')
      use iso_c_binding
      real(c_double)        :: res(mA, mB)
      integer(c_int), value :: n, mA, mB
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleM

    ! with CFI
    function dotMultiVecDeviceDoubleV_CFI(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleV_CFI')
      use iso_c_binding
      real(c_double)        :: res(:)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleV_CFI

    function dotMultiVecDeviceDoubleM_CFI(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleM_CFI')
      use iso_c_binding
      real(c_double)        :: res(:, :)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleM_CFI
  end interface
  
  interface nrm2MultiVecDevice
    function nrm2MultiVecDeviceDoubleS(res, n, deviceVecA) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceDoubleS')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceDoubleS

    function nrm2MultiVecDeviceDoubleI(res, n, deviceVecA, idx) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceDoubleI')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int), value :: idx
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceDoubleI

    function nrm2MultiVecDeviceDoubleV(res, n, deviceVecA) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceDoubleV')
      use iso_c_binding
      real(c_double)        :: res(:)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceDoubleV
  end interface

  interface amaxMultiVecDevice
    function amaxMultiVecDeviceDouble(res, n, deviceVecA) result(val) &
        & bind(c, name = 'amaxMultiVecDeviceDouble')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function amaxMultiVecDeviceDouble
  end interface

  interface asumMultiVecDevice
    function asumMultiVecDeviceDouble(res, n, deviceVecA) result(val) &
        & bind(c, name = 'asumMultiVecDeviceDouble')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function asumMultiVecDeviceDouble
  end interface

  interface axpbyMultiVecDevice
    function axpbyMultiVecDeviceDouble1(n, alpha, deviceVecA, &
                                          & beta, deviceVecB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDouble1')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDouble1

    function axpbyMultiVecDeviceDouble2(n, alpha, deviceVecA, idxA, &
                                          & beta, deviceVecB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDouble2')
      use iso_c_binding
      integer(c_int), value :: n, idxA
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDouble2

    function axpbyMultiVecDeviceDouble3(n, alpha, deviceVecA, &
                                          & beta, deviceVecB, idxB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDouble3')
      use iso_c_binding
      integer(c_int), value :: n, idxB
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDouble3

    function axpbyMultiVecDeviceDouble4(n, alpha, deviceVecA, idxA, &
                                          & beta, deviceVecB, idxB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDouble4')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDouble4

    function axpbyMultiVecDeviceDouble5(n, alpha, deviceVecA, &
                                          & beta, deviceVecB, deviceVecC) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDouble5')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDouble5

    function axpbyczMultiVecDeviceDouble1(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDouble1')
      use iso_c_binding
      integer(c_int), value :: n, idxC
      real(c_double), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDouble1

    function axpbyczMultiVecDeviceDouble2(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDouble2')
      use iso_c_binding
      integer(c_int), value :: n, idxB, idxC
      real(c_double), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDouble2

    function axpbyczMultiVecDeviceDouble3(n, alpha, deviceVecA, idxA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDouble3')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB, idxC
      real(c_double), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDouble3

    function axpbyczMultiVecDeviceDouble4(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, &
                                            & gamma, deviceVecC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDouble4')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_double), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDouble4

    function axpbyczMultiVecDeviceDouble5(n, alpha, deviceVecA, idxA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC, deviceVecD, idxD) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDouble5')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB, idxC, idxD
      real(c_double), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC, deviceVecD
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDouble5

    function colspan1DDeviceDouble(n, deviceVecA, coeff, deviceVecB, upd_flag) result(res) &
        & bind(c, name = 'colspanDeviceDouble')
      use iso_c_binding
      integer(c_int), value   :: n
      type(c_ptr), value      :: deviceVecA, deviceVecB
      real(c_double)          :: coeff(:)
      logical(c_bool), value  :: upd_flag 
      integer(c_int)  :: res
    end function colspan1DDeviceDouble
    
    function colspan2DDeviceDouble(n, deviceVecA, coeff, deviceVecB, upd_flag) result(res) &
        & bind(c, name = 'colspanDeviceDouble')
      use iso_c_binding
      integer(c_int), value   :: n
      type(c_ptr), value      :: deviceVecA, deviceVecB
      real(c_double)          :: coeff(:, :)
      logical(c_bool), value  :: upd_flag 
      integer(c_int)  :: res
    end function colspan2DDeviceDouble
  end interface

  interface mltMultiVecDevice
    function mltMultiVecDeviceDouble1(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble1')
      use iso_c_binding
      integer(c_int), value   :: n
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble1
    
    function mltMultiVecDeviceDouble2(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble2')
      use iso_c_binding
      integer(c_int), value   :: n, idxB
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble2
    
    function mltMultiVecDeviceDouble3(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble3')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble3
    
    function mltMultiVecDeviceDouble4(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble4')
      use iso_c_binding
      integer(c_int), value   :: n
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble4
    
    function mltMultiVecDeviceDouble5(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble5')
      use iso_c_binding
      integer(c_int), value   :: n, idxC
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble5
    
    function mltMultiVecDeviceDouble6(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble6')
      use iso_c_binding
      integer(c_int), value   :: n, idxB, idxC
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble6
    
    function mltMultiVecDeviceDouble7(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble7')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB, idxC
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble7

    function mltMultiVecDeviceDouble8(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble8')
      use iso_c_binding
      integer(c_int), value   :: n, idxB
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble8

    function mltMultiVecDeviceDouble9(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDouble9')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB
      real(c_double), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDouble9
  end interface

  interface upd_xyzMultiVecDevice
    function upd_xyzMultiVecDeviceDouble(n, alpha, beta, gamma, delta, deviceVecX, &
         & deviceVecY, deviceVecZ) &
         & result(res) bind(c, name = 'upd_xyzMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      real(c_double), value :: alpha, beta, gamma, delta
      type(c_ptr), value  :: deviceVecX, deviceVecY, deviceVecZ
    end function upd_xyzMultiVecDeviceDouble
  end interface

  interface xyzwMultiVecDevice
    function xyzwMultiVecDeviceDouble(n, a, b, c, d, e, f, deviceVecX, deviceVecY, deviceVecZ, deviceVecW) result(res) &
        & bind(c, name = 'xyzwMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_double), value :: a, b, c, d, e, f
      type(c_ptr), value    :: deviceVecX, deviceVecY, deviceVecZ, deviceVecW
      integer(c_int)  :: res
    end function xyzwMultiVecDeviceDouble
  end interface

  interface axyMultiVecDevice
    function axyMultiVecDeviceDouble(n, alpha, deviceVecA, deviceVecB) result(res) &
        & bind(c, name = 'axyMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axyMultiVecDeviceDouble
  end interface

  interface axybzMultiVecDevice
    function axybzMultiVecDeviceDouble(n, alpha, deviceVecA, deviceVecB, beta, deviceVecZ) result(res) &
        & bind(c, name = 'axybzMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecZ
      integer(c_int)  :: res
    end function axybzMultiVecDeviceDouble
  end interface

  interface absMultiVecDevice
    function absMultiVecDeviceDouble(n, alpha, deviceVecA) result(res) &
        & bind(c, name = 'absMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: res
    end function absMultiVecDeviceDouble

    function absMultiVecDeviceDouble2(n, alpha, deviceVecA, deviceVecB) result(res) &
        & bind(c, name = 'absMultiVecDeviceDouble2')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function absMultiVecDeviceDouble2
  end interface

  interface inner_register
    module procedure inner_registerDouble
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterDouble
  end interface

contains
  function inner_registerDouble(buffer, dval) result(res)
    real(c_double), allocatable, target :: buffer(:)
    type(c_ptr)                         :: dval
    integer(c_int)  :: res

    real(c_double)  :: dummy
    res = registerMapped(c_loc(buffer), dval, size(buffer), dummy)        
  end function inner_registerDouble

  subroutine inner_unregisterDouble(buffer)
    real(c_double), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterDouble
end module psb_d_vectordev_mod