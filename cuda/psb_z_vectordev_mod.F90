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
module psb_z_vectordev_mod
  use psb_base_vectordev_mod
 
  interface registerMapped
    function registerMappedDoubleComplex(buf, d_p, n, dummy) result(res) &
        & bind(c, name = 'registerMappedDoubleComplex')
      use iso_c_binding
      type(c_ptr), value    :: buf
      type(c_ptr)           :: d_p
      integer(c_int), value :: n
      complex(c_double_complex), value :: dummy
      integer(c_int) :: res
    end function registerMappedDoubleComplex
  end interface

  interface 
    function checkMultiVecDeviceDoubleComplex(deviceVec) result(res) &
        & bind(c, name = 'checkMultiVecDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      integer(c_int)  :: res
    end function checkMultiVecDeviceDoubleComplex
  end interface
  
  interface writeMultiVecDevice 
    function writeMultiVecDeviceDoubleComplex(deviceVec, hostVec) result(res) &
        & bind(c, name = 'writeMultiVecDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      complex(c_double_complex)      :: hostVec(*)
      integer(c_int)  :: res
    end function writeMultiVecDeviceDoubleComplex

    function writeMultiVecDeviceDoubleComplexR2(deviceVec, hostVec, ld) result(res) &
        & bind(c, name = 'writeMultiVecDeviceDoubleComplexR2')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec
      complex(c_double_complex)        :: hostVec(ld, *)
      integer(c_int), value :: ld
      integer(c_int)  :: res
    end function writeMultiVecDeviceDoubleComplexR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceDoubleComplex(deviceVec, hostVec) result(res) &
        & bind(c, name = 'readMultiVecDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      complex(c_double_complex)      :: hostVec(*)
      integer(c_int)  :: res
    end function readMultiVecDeviceDoubleComplex

    function readMultiVecDeviceDoubleComplexR2(deviceVec, hostVec, ld) result(res) &
        & bind(c, name = 'readMultiVecDeviceDoubleComplexR2')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec
      complex(c_double_complex)        :: hostVec(ld, *)
      integer(c_int), value :: ld
      integer(c_int)  :: res
    end function readMultiVecDeviceDoubleComplexR2
  end interface 

  interface allocateDoubleComplex
    function allocateDoubleComplex(didx, n) result(res) &
        & bind(c, name = 'allocateDoubleComplex') 
      use iso_c_binding
      type(c_ptr)           :: didx
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function allocateDoubleComplex
    
    function allocateMultiDoubleComplex(didx, m, n) result(res) &
        & bind(c, name = 'allocateMultiDoubleComplex') 
      use iso_c_binding
      type(c_ptr)           :: didx
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function allocateMultiDoubleComplex
  end interface

  interface writeDoubleComplex
    function writeDoubleComplex(didx, hidx, n) result(res) &
        & bind(c, name = 'writeDoubleComplex')
      use iso_c_binding
      type(c_ptr), value    :: didx
      complex(c_double_complex)        :: hidx(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function writeDoubleComplex

    function writeDoubleComplexFirst(first, didx, hidx, n, indexBase) result(res) &
        & bind(c, name = 'writeDoubleComplexFirst')
      use iso_c_binding
      integer(c_int), value :: first, n, indexBase
      type(c_ptr), value    :: didx
      complex(c_double_complex)        :: hidx(*)
      integer(c_int)  :: res
    end function writeDoubleComplexFirst

    function writeMultiDoubleComplex(didx, hidx, m, n) result(res) &
        & bind(c, name = 'writeMultiDoubleComplex')
      use iso_c_binding
      type(c_ptr), value    :: didx
      complex(c_double_complex)        :: hidx(m, *)
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function writeMultiDoubleComplex
  end interface
  
  interface readDoubleComplex
    function readDoubleComplex(didx, hidx, n) result(res) &
        & bind(c, name = 'readDoubleComplex')
      use iso_c_binding
      type(c_ptr), value    :: didx
      complex(c_double_complex)        :: hidx(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function readDoubleComplex

    function readDoubleComplexFirst(first, didx, hidx, n, indexBase) result(res) &
        & bind(c, name = 'readDoubleComplexFirst')
      use iso_c_binding
      integer(c_int), value :: first, n, indexBase
      type(c_ptr), value    :: didx
      complex(c_double_complex)        :: hidx(*)
      integer(c_int)  :: res
    end function readDoubleComplexFirst

    function readMultiDoubleComplex(didx, hidx, m, n) result(res) &
        & bind(c, name = 'readMultiDoubleComplex')
      use iso_c_binding
      type(c_ptr), value    :: didx
      complex(c_double_complex)        :: hidx(m, *)
      integer(c_int), value :: m, n
      integer(c_int)  :: res
    end function readMultiDoubleComplex
  end interface
  
  interface
    subroutine freeDoubleComplex(didx) &
        & bind(c, name = 'freeDoubleComplex')
      use iso_c_binding
      type(c_ptr), value  :: didx
    end subroutine freeDoubleComplex
  end interface
  
  interface setScalDevice
    function setScalMultiVecDeviceDoubleComplex(val, first, last, indexBase, deviceVecX) result(res) &
        & bind(c, name = 'setscalMultiVecDeviceDoubleComplex')
      use iso_c_binding
      complex(c_double_complex), value :: val
      integer(c_int), value :: first, last, indexBase
      type(c_ptr),  value   :: deviceVecX
      integer(c_int)  :: res
    end function setScalMultiVecDeviceDoubleComplex
  end interface

  interface 
    function geinsMultiVecDeviceDoubleComplex(n, deviceVecIrl, deviceVecVal, dupl, indexBase, deviceVecX) result(res) &
        & bind(c, name = 'geinsMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int), value :: n, dupl, indexBase
      type(c_ptr), value    :: deviceVecIrl, deviceVecVal, deviceVecX
      integer(c_int)  :: res
    end function geinsMultiVecDeviceDoubleComplex
  end interface

  ! New gather functions
  interface 
    function igathMultiVecDeviceDoubleComplex(deviceVec, vectorId, n, first, idx, hfirst, hostVec, indexBase) result(res) &
        & bind(c, name = 'igathMultiVecDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, n, first, hfirst, indexBase
      integer(c_int)  :: res
    end function igathMultiVecDeviceDoubleComplex
  end interface

  interface 
    function igathMultiVecDeviceDoubleComplexVecIdx(deviceVec, vectorId, n, first, idx, hfirst, hostVec, indexBase) result(res) &
        & bind(c, name = 'igathMultiVecDeviceDoubleComplexVecIdx')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      integer(c_int)  :: res
    end function igathMultiVecDeviceDoubleComplexVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceDoubleComplex(deviceVec, vectorId, first, n, idx, &
                                      & hfirst, hostVec, indexBase, beta) result(res) & 
        & bind(c, name = 'iscatMultiVecDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      complex(c_double_complex), value :: beta
      integer(c_int)  :: res
    end function iscatMultiVecDeviceDoubleComplex
  end interface

  interface 
    function iscatMultiVecDeviceDoubleComplexVecIdx(deviceVec, vectorId, first, n, idx, &
                                      & hfirst, hostVec, indexBase, beta) result(res) &
        & bind(c, name = 'iscatMultiVecDeviceDoubleComplexVecIdx')
      use iso_c_binding
      type(c_ptr), value    :: deviceVec, idx, hostVec
      integer(c_int), value :: vectorId, first, n, hfirst, indexBase
      complex(c_double_complex), value :: beta
      integer(c_int)  :: res
    end function iscatMultiVecDeviceDoubleComplexVecIdx
  end interface

  interface scalMultiVecDevice
    function scalMultiVecDeviceDoubleComplex(alpha, deviceVecA) result(val) &
        & bind(c, name = 'scalMultiVecDeviceDoubleComplex')
      use iso_c_binding
      complex(c_double_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA
      integer(c_int)        :: val
    end function scalMultiVecDeviceDoubleComplex
  end interface

  interface dotMultiVecDevice
    function dotMultiVecDeviceDoubleComplexS(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleComplexS')
      use iso_c_binding
      complex(c_double_complex)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleComplexS
    
    function dotMultiVecDeviceDoubleComplexV(res, n, mA, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleComplexV')
      use iso_c_binding
      complex(c_double_complex)        :: res(mA)
      integer(c_int), value :: n, mA
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleComplexV

    function dotMultiVecDeviceDoubleComplexM(res, n, mA, mB, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleComplexM')
      use iso_c_binding
      complex(c_double_complex)        :: res(mA, mB)
      integer(c_int), value :: n, mA, mB
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleComplexM

    ! with CFI
    function dotMultiVecDeviceDoubleComplexV_CFI(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleComplexV_CFI')
      use iso_c_binding
      complex(c_double_complex)        :: res(:)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleComplexV_CFI

    function dotMultiVecDeviceDoubleComplexM_CFI(res, n, deviceVecA, deviceVecB) result(val) &
        & bind(c, name = 'dotMultiVecDeviceDoubleComplexM_CFI')
      use iso_c_binding
      complex(c_double_complex)        :: res(:, :)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: val
    end function dotMultiVecDeviceDoubleComplexM_CFI
  end interface
  
  interface nrm2MultiVecDeviceComplex
    function nrm2MultiVecDeviceDoubleComplexS(res, n, deviceVecA) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceDoubleComplexS')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceDoubleComplexS

    function nrm2MultiVecDeviceDoubleComplexI(res, n, deviceVecA, idx) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceDoubleComplexI')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int), value :: idx
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceDoubleComplexI

    function nrm2MultiVecDeviceDoubleComplexV(res, n, deviceVecA) result(val) &
        & bind(c, name = 'nrm2MultiVecDeviceDoubleComplexV')
      use iso_c_binding
      real(c_double)        :: res(:)
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function nrm2MultiVecDeviceDoubleComplexV
  end interface

  interface amaxMultiVecDeviceComplex
    function amaxMultiVecDeviceDoubleComplex(res, n, deviceVecA) result(val) &
        & bind(c, name = 'amaxMultiVecDeviceDoubleComplex')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function amaxMultiVecDeviceDoubleComplex
  end interface

  interface asumMultiVecDeviceComplex
    function asumMultiVecDeviceDoubleComplex(res, n, deviceVecA) result(val) &
        & bind(c, name = 'asumMultiVecDeviceDoubleComplex')
      use iso_c_binding
      real(c_double)        :: res
      integer(c_int), value :: n
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: val
    end function asumMultiVecDeviceDoubleComplex
  end interface

  interface axpbyMultiVecDevice
    function axpbyMultiVecDeviceDoubleComplex1(n, alpha, deviceVecA, &
                                          & beta, deviceVecB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDoubleComplex1')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDoubleComplex1

    function axpbyMultiVecDeviceDoubleComplex2(n, alpha, deviceVecA, idxA, &
                                          & beta, deviceVecB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDoubleComplex2')
      use iso_c_binding
      integer(c_int), value :: n, idxA
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDoubleComplex2

    function axpbyMultiVecDeviceDoubleComplex3(n, alpha, deviceVecA, &
                                          & beta, deviceVecB, idxB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDoubleComplex3')
      use iso_c_binding
      integer(c_int), value :: n, idxB
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDoubleComplex3

    function axpbyMultiVecDeviceDoubleComplex4(n, alpha, deviceVecA, idxA, &
                                          & beta, deviceVecB, idxB) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDoubleComplex4')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDoubleComplex4

    function axpbyMultiVecDeviceDoubleComplex5(n, alpha, deviceVecA, &
                                          & beta, deviceVecB, deviceVecC) result(res) &
        & bind(c, name = 'axpbyMultiVecDeviceDoubleComplex5')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyMultiVecDeviceDoubleComplex5

    function axpbyczMultiVecDeviceDoubleComplex1(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDoubleComplex1')
      use iso_c_binding
      integer(c_int), value :: n, idxC
      complex(c_double_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDoubleComplex1

    function axpbyczMultiVecDeviceDoubleComplex2(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDoubleComplex2')
      use iso_c_binding
      integer(c_int), value :: n, idxB, idxC
      complex(c_double_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDoubleComplex2

    function axpbyczMultiVecDeviceDoubleComplex3(n, alpha, deviceVecA, idxA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDoubleComplex3')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB, idxC
      complex(c_double_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDoubleComplex3

    function axpbyczMultiVecDeviceDoubleComplex4(n, alpha, deviceVecA, &
                                            & beta, deviceVecB, &
                                            & gamma, deviceVecC) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDoubleComplex4')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDoubleComplex4

    function axpbyczMultiVecDeviceDoubleComplex5(n, alpha, deviceVecA, idxA, &
                                            & beta, deviceVecB, idxB, &
                                            & gamma, deviceVecC, idxC, deviceVecD, idxD) result(res) &
        & bind(c, name = 'axpbyczMultiVecDeviceDoubleComplex5')
      use iso_c_binding
      integer(c_int), value :: n, idxA, idxB, idxC, idxD
      complex(c_double_complex), value :: alpha, beta, gamma
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecC, deviceVecD
      integer(c_int)  :: res
    end function axpbyczMultiVecDeviceDoubleComplex5

    function colspan1DDeviceDoubleComplex(n, deviceVecA, coeff, deviceVecB, upd_flag) result(res) &
        & bind(c, name = 'colspanDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int), value   :: n
      type(c_ptr), value      :: deviceVecA, deviceVecB
      complex(c_double_complex)          :: coeff(:)
      logical(c_bool), value  :: upd_flag 
      integer(c_int)  :: res
    end function colspan1DDeviceDoubleComplex
    
    function colspan2DDeviceDoubleComplex(n, deviceVecA, coeff, deviceVecB, upd_flag) result(res) &
        & bind(c, name = 'colspanDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int), value   :: n
      type(c_ptr), value      :: deviceVecA, deviceVecB
      complex(c_double_complex)          :: coeff(:, :)
      logical(c_bool), value  :: upd_flag 
      integer(c_int)  :: res
    end function colspan2DDeviceDoubleComplex
  end interface

  interface mltMultiVecDevice
    function mltMultiVecDeviceDoubleComplex1(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex1')
      use iso_c_binding
      integer(c_int), value   :: n
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex1
    
    function mltMultiVecDeviceDoubleComplex2(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex2')
      use iso_c_binding
      integer(c_int), value   :: n, idxB
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex2
    
    function mltMultiVecDeviceDoubleComplex3(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex3')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex3
    
    function mltMultiVecDeviceDoubleComplex4(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex4')
      use iso_c_binding
      integer(c_int), value   :: n
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex4
    
    function mltMultiVecDeviceDoubleComplex5(n, alpha, deviceVecA, deviceVecB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex5')
      use iso_c_binding
      integer(c_int), value   :: n, idxC
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex5
    
    function mltMultiVecDeviceDoubleComplex6(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex6')
      use iso_c_binding
      integer(c_int), value   :: n, idxB, idxC
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex6
    
    function mltMultiVecDeviceDoubleComplex7(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, deviceVecC, idxC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex7')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB, idxC
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex7

    function mltMultiVecDeviceDoubleComplex8(n, alpha, deviceVecA, deviceVecB, idxB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex8')
      use iso_c_binding
      integer(c_int), value   :: n, idxB
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex8

    function mltMultiVecDeviceDoubleComplex9(n, alpha, deviceVecA, idxA, deviceVecB, idxB, &
                                    & beta, deviceVecC, conjgA, conjgB) result(res) &
        & bind(c, name = 'mltMultiVecDeviceDoubleComplex9')
      use iso_c_binding
      integer(c_int), value   :: n, idxA, idxB
      complex(c_double_complex), value   :: alpha, beta
      type(c_ptr), value      :: deviceVecA, deviceVecB, deviceVecC
      logical(c_bool), value  :: conjgA, conjgB
      integer(c_int)  :: res
    end function mltMultiVecDeviceDoubleComplex9
  end interface

  interface upd_xyzMultiVecDevice
    function upd_xyzMultiVecDeviceDoubleComplex(n, alpha, beta, gamma, delta, deviceVecX, &
         & deviceVecY, deviceVecZ) &
         & result(res) bind(c, name = 'upd_xyzMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha, beta, gamma, delta
      type(c_ptr), value  :: deviceVecX, deviceVecY, deviceVecZ
    end function upd_xyzMultiVecDeviceDoubleComplex
  end interface

  interface xyzwMultiVecDevice
    function xyzwMultiVecDeviceDoubleComplex(n, a, b, c, d, e, f, deviceVecX, deviceVecY, deviceVecZ, deviceVecW) result(res) &
        & bind(c, name = 'xyzwMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_double_complex), value :: a, b, c, d, e, f
      type(c_ptr), value    :: deviceVecX, deviceVecY, deviceVecZ, deviceVecW
      integer(c_int)  :: res
    end function xyzwMultiVecDeviceDoubleComplex
  end interface

  interface axyMultiVecDevice
    function axyMultiVecDeviceDoubleComplex(n, alpha, deviceVecA, deviceVecB) result(res) &
        & bind(c, name = 'axyMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function axyMultiVecDeviceDoubleComplex
  end interface

  interface axybzMultiVecDevice
    function axybzMultiVecDeviceDoubleComplex(n, alpha, deviceVecA, deviceVecB, beta, deviceVecZ) result(res) &
        & bind(c, name = 'axybzMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: deviceVecA, deviceVecB, deviceVecZ
      integer(c_int)  :: res
    end function axybzMultiVecDeviceDoubleComplex
  end interface

  interface absMultiVecDevice
    function absMultiVecDeviceDoubleComplex(n, alpha, deviceVecA) result(res) &
        & bind(c, name = 'absMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA
      integer(c_int)  :: res
    end function absMultiVecDeviceDoubleComplex

    function absMultiVecDeviceDoubleComplex2(n, alpha, deviceVecA, deviceVecB) result(res) &
        & bind(c, name = 'absMultiVecDeviceDoubleComplex2')
      use iso_c_binding
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
      integer(c_int)  :: res
    end function absMultiVecDeviceDoubleComplex2
  end interface

  interface inner_register
    module procedure inner_registerDoubleComplex
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterDoubleComplex
  end interface

contains
  function inner_registerDoubleComplex(buffer, dval) result(res)
    complex(c_double_complex), allocatable, target :: buffer(:)
    type(c_ptr)                         :: dval
    integer(c_int)  :: res

    complex(c_double_complex)  :: dummy
    res = registerMapped(c_loc(buffer), dval, size(buffer), dummy)        
  end function inner_registerDoubleComplex

  subroutine inner_unregisterDoubleComplex(buffer)
    complex(c_double_complex), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterDoubleComplex
end module psb_z_vectordev_mod