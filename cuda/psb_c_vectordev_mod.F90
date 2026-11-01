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
module psb_c_vectordev_mod
  use psb_base_vectordev_mod
 
  interface registerMapped
    function registerMappedFloatComplex(buf, d_p, n, dummy) result(res) &
        & bind(c, name = 'registerMappedFloatComplex')
      use iso_c_binding
      type(c_ptr), value    :: buf
      type(c_ptr)           :: d_p
      integer(c_int), value :: n
      complex(c_float_complex), value :: dummy
      integer(c_int) :: res
    end function registerMappedFloatComplex
  end interface

  interface 
    function checkMultiVecDeviceFloatComplex(deviceVec) result(res) &
        & bind(c, name = 'checkMultiVecDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      integer(c_int)  :: res
    end function checkMultiVecDeviceFloatComplex
  end interface
  
  interface writeMultiVecDevice 
    function writeMultiVecDeviceFloatComplex(deviceVec, hostVec) result(res) &
        & bind(c, name = 'writeMultiVecDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      complex(c_float_complex)      :: hostVec(*)
      integer(c_int)  :: res
    end function writeMultiVecDeviceFloatComplex

    function writeMultiVecDeviceFloatComplexR2(deviceVec, hostVec, ld) result(res) &
        & bind(c, name = 'writeMultiVecDeviceFloatComplexR2')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec
      complex(c_float_complex)        :: hostVec(ld, *)
      integer(c_int), value :: ld
      integer(c_int)  :: res
    end function writeMultiVecDeviceFloatComplexR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceFloatComplex(deviceVec, hostVec) result(res) &
        & bind(c, name = 'readMultiVecDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      complex(c_float_complex)      :: hostVec(*)
      integer(c_int)  :: res
    end function readMultiVecDeviceFloatComplex

    function readMultiVecDeviceFloatComplexR2(deviceVec, hostVec, ld) result(res) &
        & bind(c, name = 'readMultiVecDeviceFloatComplexR2')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec
      complex(c_float_complex)        :: hostVec(ld, *)
      integer(c_int), value :: ld
      integer(c_int)  :: res
    end function readMultiVecDeviceFloatComplexR2
  end interface 

  interface allocateFloatComplex
    function allocateFloatComplex(didx, n) result(res) &
        & bind(c, name = 'allocateFloatComplex') 
      use iso_c_binding
      type(c_ptr)           :: didx
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function allocateFloatComplex
    
    function allocateMultiFloatComplex(didx, m, n) result(res) &
        & bind(c, name = 'allocateMultiFloatComplex') 
      use iso_c_binding
      type(c_ptr)           :: didx
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function allocateMultiFloatComplex
  end interface

  interface writeFloatComplex
    function writeFloatComplex(didx, hidx, n) result(res) &
        & bind(c, name = 'writeFloatComplex')
      use iso_c_binding
      type(c_ptr), value    :: didx
      complex(c_float_complex)        :: hidx(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function writeFloatComplex

    function writeFloatComplexFirst(first, didx, hidx, n, indexBase) result(res) &
        & bind(c, name = 'writeFloatComplexFirst')
      use iso_c_binding
      integer(c_int), value :: first, n, indexBase
      type(c_ptr), value    :: didx
      complex(c_float_complex)        :: hidx(*)
      integer(c_int)  :: res
    end function writeFloatComplexFirst

    function writeMultiFloatComplex(didx, hidx, m, n) result(res) &
        & bind(c, name = 'writeMultiFloatComplex')
      use iso_c_binding
      type(c_ptr), value    :: didx
      complex(c_float_complex)        :: hidx(m, *)
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function writeMultiFloatComplex
  end interface
  
  interface readFloatComplex
    function readFloatComplex(didx, hidx, n) result(res) &
        & bind(c, name = 'readFloatComplex')
      use iso_c_binding
      type(c_ptr), value    :: didx
      complex(c_float_complex)        :: hidx(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function readFloatComplex

    function readFloatComplexFirst(first, didx, hidx, n, indexBase) result(res) &
        & bind(c, name = 'readFloatComplexFirst')
      use iso_c_binding
      integer(c_int), value :: first, n, indexBase
      type(c_ptr), value    :: didx
      complex(c_float_complex)        :: hidx(*)
      integer(c_int)  :: res
    end function readFloatComplexFirst

    function readMultiFloatComplex(didx, hidx, m, n) result(res) &
        & bind(c, name = 'readMultiFloatComplex')
      use iso_c_binding
      type(c_ptr), value    :: didx
      complex(c_float_complex)        :: hidx(m, *)
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function readMultiFloatComplex
  end interface
  
  interface
    subroutine freeFloatComplex(didx) &
        & bind(c, name = 'freeFloatComplex')
      use iso_c_binding
      type(c_ptr), value  :: didx
    end subroutine freeFloatComplex
  end interface
  
  interface setScalDevice
    function setScalMultiVecDeviceFloatComplex(val, first, last, indexBase, deviceVecX) result(res) &
        & bind(c, name = 'setscalMultiVecDeviceFloatComplex')
      use iso_c_binding
      complex(c_float_complex), value :: val
      integer(c_int), value :: first, last, indexBase
      type(c_ptr),  value   :: deviceVecX
      integer(c_int)  :: res
    end function setScalMultiVecDeviceFloatComplex
  end interface

  interface 
    function geinsMultiVecDeviceFloatComplex(n, deviceVecIrl, deviceVecVal, dupl, indexBase, deviceVecX) result(res) &
        & bind(c, name = 'geinsMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int), value :: n, dupl, indexBase
      type(c_ptr), value    :: deviceVecIrl, deviceVecVal, deviceVecX
      integer(c_int)  :: res
    end function geinsMultiVecDeviceFloatComplex
  end interface

  ! New gather functions
  interface 
    function igathMultiVecDeviceFloatComplex(deviceVec, vectorId, n, first, idx, hfirst, hostVec, indexBase) result(res) &
        & bind(c, name = 'igathMultiVecDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, n, first, hfirst, indexBase
      integer(c_int)  :: res
    end function igathMultiVecDeviceFloatComplex
  end interface

  interface 
    function igathMultiVecDeviceFloatComplexVecIdx(deviceVec, vectorId, n, first, idx, hfirst, hostVec, indexBase) result(res) &
        & bind(c, name = 'igathMultiVecDeviceFloatComplexVecIdx')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      integer(c_int)  :: res
    end function igathMultiVecDeviceFloatComplexVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceFloatComplex(deviceVec, vectorId, first, n, idx, &
                                      & hfirst, hostVec, indexBase, beta) result(res) & 
        & bind(c, name = 'iscatMultiVecDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      complex(c_float_complex), value :: beta
      integer(c_int)  :: res
    end function iscatMultiVecDeviceFloatComplex
  end interface

  interface 
    function iscatMultiVecDeviceFloatComplexVecIdx(deviceVec, vectorId, first, n, idx, &
                                      & hfirst, hostVec, indexBase, beta) result(res) &
        & bind(c, name = 'iscatMultiVecDeviceFloatComplexVecIdx')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      complex(c_float_complex), value :: beta
      integer(c_int)  :: res
    end function iscatMultiVecDeviceFloatComplexVecIdx
  end interface

  interface scalMultiVecDevice
    function scalMultiVecDeviceFloatComplex(alpha, deviceVecA) result(val) &
        & bind(c, name = 'scalMultiVecDeviceFloatComplex')
      use iso_c_binding
      complex(c_float_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA
      integer(c_int)        :: val
    end function scalMultiVecDeviceFloatComplex
  end interface

  interface dotMultiVecDevice
    function dotMultiVecDeviceFloatComplexS(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatComplexS')
      use iso_c_binding
      complex(c_float_complex)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatComplexS
    
    function dotMultiVecDeviceFloatComplexV(res, n, mA, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatComplexV')
      use iso_c_binding
      complex(c_float_complex)        :: res(mA)
      integer(c_int), value :: n, mA
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatComplexV

    function dotMultiVecDeviceFloatComplexM(res, n, mA, mB, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatComplexM')
      use iso_c_binding
      complex(c_float_complex)        :: res(mA, mB)
      integer(c_int), value :: n, mA, mB
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatComplexM

    ! with CFI
    function dotMultiVecDeviceFloatComplexV_CFI(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatComplexV_CFI')
      use iso_c_binding
      complex(c_float_complex)        :: res(:)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatComplexV_CFI

    function dotMultiVecDeviceFloatComplexM_CFI(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceFloatComplexM_CFI')
      use iso_c_binding
      complex(c_float_complex)        :: res(:, :)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceFloatComplexM_CFI
  end interface
  
  interface nrm2MultiVecDeviceComplex
    function nrm2MultiVecDeviceFloatComplexS(res, n, deviceVecA) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceFloatComplexS')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceFloatComplexS

    function nrm2MultiVecDeviceFloatComplexI(res, n, deviceVecA, idx) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceFloatComplexI')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int), value :: idx
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceFloatComplexI

    function nrm2MultiVecDeviceFloatComplexV(res, n, deviceVecA) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceFloatComplexV')
      use iso_c_binding
      real(c_float)        :: res(:)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceFloatComplexV
  end interface

  interface amaxMultiVecDeviceComplex
    function amaxMultiVecDeviceFloatComplex(res, n, deviceVecA) result(val) &
        & bind(c, name = 'amaxMultiVecDeviceFloatComplex')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function amaxMultiVecDeviceFloatComplex
  end interface

  interface asumMultiVecDeviceComplex
    function asumMultiVecDeviceFloatComplex(res, n, deviceVecA) result(val) &
        & bind(c, name = 'asumMultiVecDeviceFloatComplex')
      use iso_c_binding
      real(c_float)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function asumMultiVecDeviceFloatComplex
  end interface

  interface axpbyMultiVecDevice
    function axpbyMultiVecDeviceFloatComplex1(n, alpha, deviceVecA, &
                                          & beta, deviceVecB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloatComplex1')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloatComplex1

    function axpbyMultiVecDeviceFloatComplex2(n, alpha, deviceVecA, idxA, &
                                          & beta, deviceVecB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloatComplex2')
      use iso_c_binding
      integer(c_int), value :: n, idxA
      complex(c_float_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloatComplex2

    function axpbyMultiVecDeviceFloatComplex3(n, alpha, deviceVecA, &
                                          & beta, deviceVecB, idxB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloatComplex3')
      use iso_c_binding
      integer(c_int), value :: n, idxB
      complex(c_float_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloatComplex3

    function axpbyMultiVecDeviceFloatComplex4(n, alpha, deviceVecA, idxA, &
                                          & beta, deviceVecB, idxB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloatComplex4')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB
      complex(c_float_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloatComplex4

    function axpbyMultiVecDeviceFloatComplex5(n, alpha, deviceVecA, &
                                          & beta, deviceVecB, deviceVecC) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceFloatComplex5')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceFloatComplex5

    function axpbyczMultiVecDeviceFloatComplex1(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloatComplex1')
      use iso_c_binding
      integer(c_int), value :: n, idxC
      complex(c_float_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloatComplex1

    function axpbyczMultiVecDeviceFloatComplex2(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloatComplex2')
      use iso_c_binding
      integer(c_int), value :: n, idxB, idxC
      complex(c_float_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloatComplex2

    function axpbyczMultiVecDeviceFloatComplex3(n, alpha, deviceVecA, idxA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloatComplex3')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB, idxC
      complex(c_float_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloatComplex3

    function axpbyczMultiVecDeviceFloatComplex4(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, &
                                            & gamma, deviceVecC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloatComplex4')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloatComplex4

    function axpbyczMultiVecDeviceFloatComplex5(n, alpha, deviceVecA, idxA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC, deviceVecD, idxD) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceFloatComplex5')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB, idxC, idxD
      complex(c_float_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC, deviceVecD
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceFloatComplex5

    function colspan1DDeviceFloatComplex(n, deviceVecA, coeff, deviceVecB, upd_flag) result(res) &
        & bind(c, name = 'colspanDeviceFloatComplex')
      use iso_c_binding
      integer(c_int), value   :: n
      type(c_ptr), value      :: deviceVecA, deviceVecB
      complex(c_float_complex)          :: coeff(:)
      logical(c_bool), value  :: upd_flag 
      integer(c_int)  :: res
    end function colspan1DDeviceFloatComplex
    
    function colspan2DDeviceFloatComplex(n, deviceVecA, coeff, deviceVecB, upd_flag) result(res) &
        & bind(c, name = 'colspanDeviceFloatComplex')
      use iso_c_binding
      integer(c_int), value   :: n
      type(c_ptr), value      :: deviceVecA, deviceVecB
      complex(c_float_complex)          :: coeff(:, :)
      logical(c_bool), value  :: upd_flag 
      integer(c_int)  :: res
    end function colspan2DDeviceFloatComplex
  end interface

  interface mltMultiVecDevice
    function mltMultiVecDeviceFloatComplex1(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex1')
      use iso_c_binding
      integer(c_int), value   :: n
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex1
    
    function mltMultiVecDeviceFloatComplex2(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex2')
      use iso_c_binding
      integer(c_int), value   :: n, idxB
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex2
    
    function mltMultiVecDeviceFloatComplex3(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex3')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex3
    
    function mltMultiVecDeviceFloatComplex4(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex4')
      use iso_c_binding
      integer(c_int), value   :: n
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex4
    
    function mltMultiVecDeviceFloatComplex5(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex5')
      use iso_c_binding
      integer(c_int), value   :: n, idxC
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex5
    
    function mltMultiVecDeviceFloatComplex6(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex6')
      use iso_c_binding
      integer(c_int), value   :: n, idxB, idxC
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex6
    
    function mltMultiVecDeviceFloatComplex7(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex7')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB, idxC
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex7

    function mltMultiVecDeviceFloatComplex8(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex8')
      use iso_c_binding
      integer(c_int), value   :: n, idxB
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex8

    function mltMultiVecDeviceFloatComplex9(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceFloatComplex9')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB
      complex(c_float_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceFloatComplex9
  end interface

  interface upd_xyzMultiVecDevice
    function upd_xyzMultiVecDeviceFloatComplex(n, alpha, beta, gamma, delta, deviceVecX, &
         & deviceVecY, deviceVecZ) &
         & result(res) bind(c, name = 'upd_xyzMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha, beta, gamma, delta
      type(c_ptr), value  :: deviceVecX, deviceVecY, deviceVecZ
    end function upd_xyzMultiVecDeviceFloatComplex
  end interface

  interface xyzwMultiVecDevice
    function xyzwMultiVecDeviceFloatComplex(n, a, b, c, d, e, f, deviceVecX, deviceVecY, deviceVecZ, deviceVecW) result(res) &
        & bind(c, name = 'xyzwMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_float_complex), value :: a, b, c, d, e, f
      type(c_ptr), value    :: deviceVecX, deviceVecY, deviceVecZ, deviceVecW
      integer(c_int)  :: res
    end function xyzwMultiVecDeviceFloatComplex
  end interface

  interface axyMultiVecDevice
    function axyMultiVecDeviceFloatComplex(n, alpha, deviceVecA, deviceVecB) result(res) &
        & bind(c, name = 'axyMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axyMultiVecDeviceFloatComplex
  end interface

  interface axybzMultiVecDevice
    function axybzMultiVecDeviceFloatComplex(n, alpha, deviceVecA, deviceVecB, beta, deviceVecZ) result(res) &
        & bind(c, name = 'axybzMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecZ
      integer(c_int)  :: res
    end function axybzMultiVecDeviceFloatComplex
  end interface

  interface absMultiVecDevice
    function absMultiVecDeviceFloatComplex(n, alpha, deviceVecA) result(res) &
        & bind(c, name = 'absMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: res
    end function absMultiVecDeviceFloatComplex

    function absMultiVecDeviceFloatComplex2(n, alpha, deviceVecA, deviceVecB) result(res) &
        & bind(c, name = 'absMultiVecDeviceFloatComplex2')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function absMultiVecDeviceFloatComplex2
  end interface

  interface inner_register
    module procedure inner_registerFloatComplex
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterFloatComplex
  end interface

contains
  function inner_registerFloatComplex(buffer, dval) result(res)
    complex(c_float_complex), allocatable, target :: buffer(:)
    type(c_ptr)                         :: dval
    integer(c_int)  :: res

    complex(c_float_complex)  :: dummy
    res = registerMapped(c_loc(buffer), dval, size(buffer), dummy)        
  end function inner_registerFloatComplex

  subroutine inner_unregisterFloatComplex(buffer)
    complex(c_float_complex), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterFloatComplex
end module psb_c_vectordev_mod