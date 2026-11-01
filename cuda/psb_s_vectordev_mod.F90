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
module psb_s_vectordev_mod
  use psb_base_vectordev_mod
 
  interface registerMapped
    function registerMappedFloat(buf, d_p, n, dummy) result(res) &
        & bind(c, name = 'registerMappedFloat')
      use iso_c_binding
      type(c_ptr), value    :: buf
      type(c_ptr)           :: d_p
      integer(c_int), value :: n
      real(c_float), value :: dummy
      integer(c_int) :: res
    end function registerMappedFloat
  end interface

  interface 
    function checkMultiVecDeviceFloat(deviceVec) result(res) &
        & bind(c, name = 'checkMultiVecDeviceFloat')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      integer(c_int)  :: res
    end function checkMultiVecDeviceFloat
  end interface
  
  interface writeMultiVecDevice 
    function writeMultiVecDeviceFloat(deviceVec, hostVec) result(res) &
        & bind(c, name = 'writeMultiVecDeviceFloat')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      real(c_float)      :: hostVec(*)
      integer(c_int)  :: res
    end function writeMultiVecDeviceFloat

    function writeMultiVecDeviceFloatR2(deviceVec, hostVec, ld) result(res) &
        & bind(c, name = 'writeMultiVecDeviceFloatR2')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec
      real(c_float)        :: hostVec(ld, *)
      integer(c_int), value :: ld
      integer(c_int)  :: res
    end function writeMultiVecDeviceFloatR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceFloat(deviceVec, hostVec) result(res) &
        & bind(c, name = 'readMultiVecDeviceFloat')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      real(c_float)      :: hostVec(*)
      integer(c_int)  :: res
    end function readMultiVecDeviceFloat

    function readMultiVecDeviceFloatR2(deviceVec, hostVec, ld) result(res) &
        & bind(c, name = 'readMultiVecDeviceFloatR2')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec
      real(c_float)        :: hostVec(ld, *)
      integer(c_int), value :: ld
      integer(c_int)  :: res
    end function readMultiVecDeviceFloatR2
  end interface 

  interface allocateFloat
    function allocateFloat(didx, n) result(res) &
        & bind(c, name = 'allocateFloat') 
      use iso_c_binding
      type(c_ptr)           :: didx
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function allocateFloat
    
    function allocateMultiFloat(didx, m, n) result(res) &
        & bind(c, name = 'allocateMultiFloat') 
      use iso_c_binding
      type(c_ptr)           :: didx
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function allocateMultiFloat
  end interface

  interface writeFloat
    function writeFloat(didx, hidx, n) result(res) &
        & bind(c, name = 'writeFloat')
      use iso_c_binding
      type(c_ptr), value    :: didx
      real(c_float)        :: hidx(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function writeFloat

    function writeFloatFirst(first, didx, hidx, n, indexBase) result(res) &
        & bind(c, name = 'writeFloatFirst')
      use iso_c_binding
      integer(c_int), value :: first, n, indexBase
      type(c_ptr), value    :: didx
      real(c_float)        :: hidx(*)
      integer(c_int)  :: res
    end function writeFloatFirst

    function writeMultiFloat(didx, hidx, m, n) result(res) &
        & bind(c, name = 'writeMultiFloat')
      use iso_c_binding
      type(c_ptr), value    :: didx
      real(c_float)        :: hidx(m, *)
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function writeMultiFloat
  end interface
  
  interface readFloat
    function readFloat(didx, hidx, n) result(res) &
        & bind(c, name = 'readFloat')
      use iso_c_binding
      type(c_ptr), value    :: didx
      real(c_float)        :: hidx(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function readFloat

    function readFloatFirst(first, didx, hidx, n, indexBase) result(res) &
        & bind(c, name = 'readFloatFirst')
      use iso_c_binding
      integer(c_int), value :: first, n, indexBase
      type(c_ptr), value    :: didx
      real(c_float)        :: hidx(*)
      integer(c_int)  :: res
    end function readFloatFirst

    function readMultiFloat(didx, hidx, m, n) result(res) &
        & bind(c, name = 'readMultiFloat')
      use iso_c_binding
      type(c_ptr), value    :: didx
      real(c_float)        :: hidx(m, *)
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function readMultiFloat
  end interface
  
  interface
    subroutine freeFloat(didx) &
        & bind(c, name = 'freeFloat')
      use iso_c_binding
      type(c_ptr), value  :: didx
    end subroutine freeFloat
  end interface
  
  interface setScalDevice
    function setScalMultiVecDeviceFloat(val, first, last, indexBase, deviceVecX) result(res) &
        & bind(c, name = 'setscalMultiVecDeviceFloat')
      use iso_c_binding
      real(c_float), value :: val
      integer(c_int), value :: first, last, indexBase
      type(c_ptr),  value   :: deviceVecX
      integer(c_int)  :: res
    end function setScalMultiVecDeviceFloat
  end interface

  interface 
    function geinsMultiVecDeviceFloat(n, deviceVecIrl, deviceVecVal, dupl, indexBase, deviceVecX) result(res) &
        & bind(c, name = 'geinsMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int), value :: n, dupl, indexBase
      type(c_ptr), value    :: deviceVecIrl, deviceVecVal, deviceVecX
      integer(c_int)  :: res
    end function geinsMultiVecDeviceFloat
  end interface

  ! New gather functions
  interface 
    function igathMultiVecDeviceFloat(deviceVec, vectorId, n, first, idx, hfirst, hostVec, indexBase) result(res) &
        & bind(c, name = 'igathMultiVecDeviceFloat')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, n, first, hfirst, indexBase
      integer(c_int)  :: res
    end function igathMultiVecDeviceFloat
  end interface

  interface 
    function igathMultiVecDeviceFloatVecIdx(deviceVec, vectorId, n, first, idx, hfirst, hostVec, indexBase) result(res) &
        & bind(c, name = 'igathMultiVecDeviceFloatVecIdx')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      integer(c_int)  :: res
    end function igathMultiVecDeviceFloatVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceFloat(deviceVec, vectorId, first, n, idx, &
                                      & hfirst, hostVec, indexBase, beta) result(res) & 
        & bind(c, name = 'iscatMultiVecDeviceFloat')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      real(c_float), value :: beta
      integer(c_int)  :: res
    end function iscatMultiVecDeviceFloat
  end interface

  interface 
    function iscatMultiVecDeviceFloatVecIdx(deviceVec, vectorId, first, n, idx, &
                                      & hfirst, hostVec, indexBase, beta) result(res) &
        & bind(c, name = 'iscatMultiVecDeviceFloatVecIdx')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      real(c_float), value :: beta
      integer(c_int)  :: res
    end function iscatMultiVecDeviceFloatVecIdx
  end interface

  interface scalMultiVecDevice
    function scalMultiVecDeviceFloat(alpha, deviceVecA) result(val) &
        & bind(c, name = 'scalMultiVecDeviceFloat')
      use iso_c_binding
      real(c_float), value :: alpha
      type(c_ptr), value    :: deviceVecA
      integer(c_int)        :: val
    end function scalMultiVecDeviceFloat
  end interface

  interface dotMultiVecDevice
    function dotMultiVecDeviceFloatS(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatS')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatS
    
    function dotMultiVecDeviceFloatV(res, n, mA, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatV')
      use iso_c_binding
      real(c_float)        :: res(mA)
      integer(c_int), value :: n, mA
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatV

    function dotMultiVecDeviceFloatM(res, n, mA, mB, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatM')
      use iso_c_binding
      real(c_float)        :: res(mA, mB)
      integer(c_int), value :: n, mA, mB
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatM

    ! with CFI
    function dotMultiVecDeviceFloatV_CFI(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatV_CFI')
      use iso_c_binding
      real(c_float)        :: res(:)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatV_CFI

    function dotMultiVecDeviceFloatM_CFI(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatM_CFI')
      use iso_c_binding
      real(c_float)        :: res(:, :)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatM_CFI
  end interface
  
  interface nrm2MultiVecDevice
    function nrm2MultiVecDeviceFloatS(res, n, deviceVecA) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceFloatS')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceFloatS

    function nrm2MultiVecDeviceFloatI(res, n, deviceVecA, idx) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceFloatI')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int), value :: idx
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceFloatI

    function nrm2MultiVecDeviceFloatV(res, n, deviceVecA) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceFloatV')
      use iso_c_binding
      real(c_float)        :: res(:)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceFloatV
  end interface

  interface amaxMultiVecDevice
    function amaxMultiVecDeviceFloat(res, n, deviceVecA) result(val) &
        & bind(c, name = 'amaxMultiVecDeviceFloat')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function amaxMultiVecDeviceFloat
  end interface

  interface asumMultiVecDevice
    function asumMultiVecDeviceFloat(res, n, deviceVecA) result(val) &
        & bind(c, name = 'asumMultiVecDeviceFloat')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function asumMultiVecDeviceFloat
  end interface

  interface axpbyMultiVecDevice
    function axpbyMultiVecDeviceFloat1(n, alpha, deviceVecA, &
                                          & beta, deviceVecB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloat1')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_float), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloat1

    function axpbyMultiVecDeviceFloat2(n, alpha, deviceVecA, idxA, &
                                          & beta, deviceVecB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloat2')
      use iso_c_binding
      integer(c_int), value :: n, idxA
      real(c_float), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloat2

    function axpbyMultiVecDeviceFloat3(n, alpha, deviceVecA, &
                                          & beta, deviceVecB, idxB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloat3')
      use iso_c_binding
      integer(c_int), value :: n, idxB
      real(c_float), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloat3

    function axpbyMultiVecDeviceFloat4(n, alpha, deviceVecA, idxA, &
                                          & beta, deviceVecB, idxB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloat4')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB
      real(c_float), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloat4

    function axpbyMultiVecDeviceFloat5(n, alpha, deviceVecA, &
                                          & beta, deviceVecB, deviceVecC) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloat5')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_float), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloat5

    function axpbyczMultiVecDeviceFloat1(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloat1')
      use iso_c_binding
      integer(c_int), value :: n, idxC
      real(c_float), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloat1

    function axpbyczMultiVecDeviceFloat2(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloat2')
      use iso_c_binding
      integer(c_int), value :: n, idxB, idxC
      real(c_float), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloat2

    function axpbyczMultiVecDeviceFloat3(n, alpha, deviceVecA, idxA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloat3')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB, idxC
      real(c_float), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloat3

    function axpbyczMultiVecDeviceFloat4(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, &
                                            & gamma, deviceVecC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloat4')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_float), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloat4

    function axpbyczMultiVecDeviceFloat5(n, alpha, deviceVecA, idxA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC, deviceVecD, idxD) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloat5')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB, idxC, idxD
      real(c_float), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC, deviceVecD
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloat5

    function colspan1DDeviceFloat(n, deviceVecA, coeff, deviceVecB, upd_flag) result(res) &
        & bind(c, name = 'colspanDeviceFloat')
      use iso_c_binding
      integer(c_int), value   :: n
      type(c_ptr), value      :: deviceVecA, deviceVecB
      real(c_float)          :: coeff(:)
      logical(c_bool), value  :: upd_flag 
      integer(c_int)  :: res
    end function colspan1DDeviceFloat
    
    function colspan2DDeviceFloat(n, deviceVecA, coeff, deviceVecB, upd_flag) result(res) &
        & bind(c, name = 'colspanDeviceFloat')
      use iso_c_binding
      integer(c_int), value   :: n
      type(c_ptr), value      :: deviceVecA, deviceVecB
      real(c_float)          :: coeff(:, :)
      logical(c_bool), value  :: upd_flag 
      integer(c_int)  :: res
    end function colspan2DDeviceFloat
  end interface

  interface mltMultiVecDevice
    function mltMultiVecDeviceFloat1(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat1')
      use iso_c_binding
      integer(c_int), value   :: n
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat1
    
    function mltMultiVecDeviceFloat2(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat2')
      use iso_c_binding
      integer(c_int), value   :: n, idxB
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat2
    
    function mltMultiVecDeviceFloat3(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat3')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat3
    
    function mltMultiVecDeviceFloat4(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat4')
      use iso_c_binding
      integer(c_int), value   :: n
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat4
    
    function mltMultiVecDeviceFloat5(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat5')
      use iso_c_binding
      integer(c_int), value   :: n, idxC
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat5
    
    function mltMultiVecDeviceFloat6(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat6')
      use iso_c_binding
      integer(c_int), value   :: n, idxB, idxC
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat6
    
    function mltMultiVecDeviceFloat7(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat7')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB, idxC
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat7

    function mltMultiVecDeviceFloat8(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat8')
      use iso_c_binding
      integer(c_int), value   :: n, idxB
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat8

    function mltMultiVecDeviceFloat9(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloat9')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB
      real(c_float), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloat9
  end interface

  interface upd_xyzMultiVecDevice
    function upd_xyzMultiVecDeviceFloat(n, alpha, beta, gamma, delta, deviceVecX, &
         & deviceVecY, deviceVecZ) &
         & result(res) bind(c, name = 'upd_xyzMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      real(c_float), value :: alpha, beta, gamma, delta
      type(c_ptr), value  :: deviceVecX, deviceVecY, deviceVecZ
    end function upd_xyzMultiVecDeviceFloat
  end interface

  interface xyzwMultiVecDevice
    function xyzwMultiVecDeviceFloat(n, a, b, c, d, e, f, deviceVecX, deviceVecY, deviceVecZ, deviceVecW) result(res) &
        & bind(c, name = 'xyzwMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_float), value :: a, b, c, d, e, f
      type(c_ptr), value    :: deviceVecX, deviceVecY, deviceVecZ, deviceVecW
      integer(c_int)  :: res
    end function xyzwMultiVecDeviceFloat
  end interface

  interface axyMultiVecDevice
    function axyMultiVecDeviceFloat(n, alpha, deviceVecA, deviceVecB) result(res) &
        & bind(c, name = 'axyMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_float), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axyMultiVecDeviceFloat
  end interface

  interface axybzMultiVecDevice
    function axybzMultiVecDeviceFloat(n, alpha, deviceVecA, deviceVecB, beta, deviceVecZ) result(res) &
        & bind(c, name = 'axybzMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_float), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecZ
      integer(c_int)  :: res
    end function axybzMultiVecDeviceFloat
  end interface

  interface absMultiVecDevice
    function absMultiVecDeviceFloat(n, alpha, deviceVecA) result(res) &
        & bind(c, name = 'absMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_float), value :: alpha
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: res
    end function absMultiVecDeviceFloat

    function absMultiVecDeviceFloat2(n, alpha, deviceVecA, deviceVecB) result(res) &
        & bind(c, name = 'absMultiVecDeviceFloat2')
      use iso_c_binding
      integer(c_int), value :: n
      real(c_float), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function absMultiVecDeviceFloat2
  end interface

  interface inner_register
    module procedure inner_registerFloat
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterFloat
  end interface

contains
  function inner_registerFloat(buffer, dval) result(res)
    real(c_float), allocatable, target :: buffer(:)
    type(c_ptr)                         :: dval
    integer(c_int)  :: res

    real(c_float)  :: dummy
    res = registerMapped(c_loc(buffer), dval, size(buffer), dummy)        
  end function inner_registerFloat

  subroutine inner_unregisterFloat(buffer)
    real(c_float), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterFloat
end module psb_s_vectordev_mod