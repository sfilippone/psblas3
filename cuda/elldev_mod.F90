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

module elldev_mod
  use iso_c_binding 
  use core_mod 
  
  type, bind(c) :: elldev_parms
    integer(c_int) :: element_type
    integer(c_int) :: pitch
    integer(c_int) :: rows
    integer(c_int) :: columns
    integer(c_int) :: maxRowSize
    integer(c_int) :: avgRowSize
    integer(c_int) :: firstIndex
  end type elldev_parms

  interface 
    function FgetEllDeviceParams(rows, maxRowSize, nnzeros, columns, elementType, firstIndex) result(res) &
        & bind(c, name = 'getEllDeviceParams')
      use iso_c_binding
      import :: elldev_parms
      integer(c_int), value :: rows, maxRowSize, nnzeros, columns, elementType, firstIndex
      type(elldev_parms)  :: res
    end function FgetEllDeviceParams
  end interface

  interface 
    function FallocEllDevice(deviceMat, rows, maxRowSize, nnzeros, columns, elementType, firstIndex) result(res) &
        & bind(c, name = 'FallocEllDevice')
      use iso_c_binding
      type(c_ptr)           :: deviceMat
      integer(c_int), value :: rows, maxRowSize, nnzeros, columns, elementType, firstIndex
      integer(c_int)  :: res
    end function FallocEllDevice
  end interface

  interface writeEllDevice
    function writeEllDeviceFloat(deviceMat, val, ja, ldj, irn, idiag) result(res) &
        & bind(c, name = 'writeEllDeviceFloat')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_float)         :: val(ldj, *)
      integer(c_int)        :: ja(ldj, *), irn(*), idiag(*)
      integer(c_int), value :: ldj
      integer(c_int)  :: res
    end function writeEllDeviceFloat

    function writeEllDeviceDouble(deviceMat, val, ja, ldj, irn, idiag) result(res) &
        & bind(c, name = 'writeEllDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_double)        :: val(ldj, *)
      integer(c_int)        :: ja(ldj, *), irn(*), idiag(*)
      integer(c_int), value :: ldj
      integer(c_int)  :: res
    end function writeEllDeviceDouble

    function writeEllDeviceFloatComplex(deviceMat, val, ja, ldj, irn, idiag) result(res) &
        & bind(c, name = 'writeEllDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value       :: deviceMat
      complex(c_float_complex) :: val(ldj, *)
      integer(c_int)           :: ja(ldj, *), irn(*), idiag(*)
      integer(c_int), value    :: ldj
      integer(c_int)  :: res
    end function writeEllDeviceFloatComplex

    function writeEllDeviceDoubleComplex(deviceMat, val, ja, ldj, irn, idiag) result(res) &
        & bind(c, name = 'writeEllDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_double_complex) :: val(ldj, *)
      integer(c_int)            :: ja(ldj, *), irn(*), idiag(*)
      integer(c_int), value     :: ldj
      integer(c_int)  :: res
    end function writeEllDeviceDoubleComplex
  end interface

  interface readEllDevice 
    function readEllDeviceFloat(deviceMat, val, ja, ldj, irn, idiag) result(res) &
        & bind(c, name = 'readEllDeviceFloat')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_float)         :: val(ldj, *)
      integer(c_int)        :: ja(ldj, *), irn(*), idiag(*)
      integer(c_int), value :: ldj
      integer(c_int)  :: res
    end function readEllDeviceFloat

    function readEllDeviceDouble(deviceMat, val, ja, ldj, irn, idiag) result(res) &
        & bind(c, name = 'readEllDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_double)        :: val(ldj, *)
      integer(c_int)        :: ja(ldj, *), irn(*), idiag(*)
      integer(c_int), value :: ldj
      integer(c_int)  :: res
    end function readEllDeviceDouble

    function readEllDeviceFloatComplex(deviceMat, val, ja, ldj, irn, idiag) result(res) &
        & bind(c, name = 'readEllDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value       :: deviceMat
      complex(c_float_complex) :: val(ldj, *)
      integer(c_int)           :: ja(ldj, *), irn(*), idiag(*)
      integer(c_int), value    :: ldj
      integer(c_int)  :: res
    end function readEllDeviceFloatComplex

    function readEllDeviceDoubleComplex(deviceMat, val, ja, ldj, irn, idiag) result(res) &
        & bind(c, name = 'readEllDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      integer(c_int), value     :: ldj
      complex(c_double_complex) :: val(ldj, *)
      integer(c_int)            :: ja(ldj, *), irn(*), idiag(*)
      integer(c_int)  :: res
    end function readEllDeviceDoubleComplex
  end interface

  interface 
    subroutine freeEllDevice(deviceMat) &
        & bind(c, name = 'freeEllDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
    end subroutine freeEllDevice
  end interface

  interface 
    subroutine zeroEllDevice(deviceMat) &
        & bind(c, name = 'zeroEllDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
    end subroutine zeroEllDevice
  end interface

  interface 
    subroutine resetEllTimer() &
        & bind(c, name = 'resetEllTimer')
      use iso_c_binding
    end subroutine resetEllTimer
  end interface

  interface 
    function getEllTimer() result(res) &
        & bind(c, name = 'getEllTimer')
      use iso_c_binding
      real(c_double)  :: res
    end function getEllTimer
  end interface

  interface 
    function getEllDevicePitch(deviceMat) result(res) &
        & bind(c, name = 'getEllDevicePitch')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      integer(c_int)  :: res
    end function getEllDevicePitch
  end interface

  interface 
    function getEllDeviceMaxRowSize(deviceMat) &
        & bind(c, name = 'getEllDeviceMaxRowSize') result(res)
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      integer(c_int)  :: res
    end function getEllDeviceMaxRowSize
  end interface

  interface psi_CopyCooToElg
    function psiCopyCooToElgFloat(nr, nc, nza, hacksz, ldv, nzm, irn, &
                               & idisp, ja, val, deviceMat) result(res) &
        & bind(c, name = 'psiCopyCooToElgFloat')
      use iso_c_binding
      integer(c_int), value :: nr, nc, nza, hacksz, ldv, nzm
      integer(c_int)        :: irn(*), idisp(*), ja(*)
      real(c_float)         :: val(*)
      type(c_ptr), value    :: deviceMat
      integer(c_int)  :: res
    end function psiCopyCooToElgFloat

    function psiCopyCooToElgDouble(nr, nc, nza, hacksz, ldv, nzm, irn, &
                                 & idisp, ja, val, deviceMat) result(res) &
        & bind(c, name = 'psiCopyCooToElgDouble')
      use iso_c_binding
      integer(c_int), value   :: nr, nc, nza, hacksz, ldv, nzm
      integer(c_int)          :: irn(*), idisp(*), ja(*)
      real(c_double)          :: val(*)
      type(c_ptr), value      :: deviceMat
      integer(c_int)  :: res
    end function psiCopyCooToElgDouble

    function psiCopyCooToElgFloatComplex(nr, nc, nza, hacksz, ldv, nzm, irn, &
                                       & idisp, ja, val, deviceMat) result(res) &
        & bind(c, name = 'psiCopyCooToElgFloatComplex')
      use iso_c_binding
      integer(c_int), value     :: nr, nc, nza, hacksz, ldv, nzm
      integer(c_int)            :: irn(*), idisp(*), ja(*)
      complex(c_float_complex)  :: val(*)
      type(c_ptr), value        :: deviceMat
      integer(c_int)  :: res
    end function psiCopyCooToElgFloatComplex

    function psiCopyCooToElgDoubleComplex(nr, nc, nza, hacksz, ldv, nzm, irn, &
                                       & idisp, ja, val, deviceMat) result(res) &
        & bind(c, name = 'psiCopyCooToElgDoubleComplex')
      use iso_c_binding
      integer(c_int), value     :: nr, nc, nza, hacksz, ldv, nzm
      integer(c_int)            :: irn(*), idisp(*), ja(*)
      complex(c_double_complex) :: val(*)
      type(c_ptr), value        :: deviceMat
      integer(c_int)  :: res
    end function psiCopyCooToElgDoubleComplex
  end interface

  interface csputEllDeviceFloat
    function dev_csputEllDeviceFloat(deviceMat, nnz, ia, ja, val) result(res) &
        & bind(c, name = 'dev_csputEllDeviceFloat')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, ia, ja, val
      integer(c_int), value :: nnz
      integer(c_int)  :: res
    end function dev_csputEllDeviceFloat
  end interface

  interface csputEllDeviceDouble
    function dev_csputEllDeviceDouble(deviceMat, nnz, ia, ja, val) result(res) &
        & bind(c, name = 'dev_csputEllDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat , ia, ja, val
      integer(c_int), value :: nnz
      integer(c_int)  :: res
    end function dev_csputEllDeviceDouble
  end interface

  interface csputEllDeviceFloatComplex
    function dev_csputEllDeviceFloatComplex(deviceMat, nnz, ia, ja, val) result(res) &
        & bind(c, name = 'dev_csputEllDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat , ia, ja, val
      integer(c_int), value :: nnz
      integer(c_int)  :: res
    end function dev_csputEllDeviceFloatComplex
  end interface

  interface csputEllDeviceDoubleComplex
    function dev_csputEllDeviceDoubleComplex(deviceMat, nnz, ia, ja, val) result(res) &
        & bind(c, name = 'dev_csputEllDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat , ia, ja, val
      integer(c_int), value :: nnz
      integer(c_int)  :: res
    end function dev_csputEllDeviceDoubleComplex
  end interface

  interface spmvEllDevice
    function spmvEllDeviceFloatF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvEllDeviceFloatF')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int)  :: res
    end function spmvEllDeviceFloatF

    function spmvEllDeviceFloatMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvEllDeviceFloatMV')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int), value :: idx_y
      integer(c_int)  :: res
    end function spmvEllDeviceFloatMV

    function spmvEllDeviceFloatVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvEllDeviceFloatVM')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int), value :: idx_x
      integer(c_int)  :: res
    end function spmvEllDeviceFloatVM

    function spmvEllDeviceFloatMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvEllDeviceFloatMMi')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value :: alpha, beta
      integer(c_int), value :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvEllDeviceFloatMMi

    function spmvEllDeviceDoubleF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvEllDeviceDoubleF')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int)  :: res
    end function spmvEllDeviceDoubleF

    function spmvEllDeviceDoubleMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvEllDeviceDoubleMV')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_y
      integer(c_int)  :: res
    end function spmvEllDeviceDoubleMV

    function spmvEllDeviceDoubleVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvEllDeviceDoubleVM')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_x
      integer(c_int)  :: res
    end function spmvEllDeviceDoubleVM

    function spmvEllDeviceDoubleMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvEllDeviceDoubleMMi')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvEllDeviceDoubleMMi

    function spmvEllDeviceFloatComplexF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvEllDeviceFloatComplexF')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y 
      complex(c_float_complex), value :: alpha, beta
      integer(c_int)  :: res
    end function spmvEllDeviceFloatComplexF

    function spmvEllDeviceFloatComplexMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvEllDeviceFloatComplexMV')
      use iso_c_binding
      type(c_ptr), value               :: deviceMat, x, y
      complex(c_float_complex), value  :: alpha, beta
      integer(c_int), value            :: idx_y
      integer(c_int)  :: res
    end function spmvEllDeviceFloatComplexMV

    function spmvEllDeviceFloatComplexVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvEllDeviceFloatComplexVM')
      use iso_c_binding
      type(c_ptr), value               :: deviceMat, x, y
      complex(c_float_complex), value  :: alpha, beta
      integer(c_int), value            :: idx_x
      integer(c_int)  :: res
    end function spmvEllDeviceFloatComplexVM

    function spmvEllDeviceFloatComplexMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvEllDeviceFloatComplexMMi')
      use iso_c_binding
      type(c_ptr), value               :: deviceMat, x, y
      complex(c_float_complex), value  :: alpha, beta
      integer(c_int), value            :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvEllDeviceFloatComplexMMi

    function spmvEllDeviceDoubleComplexF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvEllDeviceDoubleComplexF')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y 
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int)  :: res
    end function spmvEllDeviceDoubleComplexF

    function spmvEllDeviceDoubleComplexMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvEllDeviceDoubleComplexMV')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_y
      integer(c_int)  :: res
    end function spmvEllDeviceDoubleComplexMV

    function spmvEllDeviceDoubleComplexVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvEllDeviceDoubleComplexVM')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_x
      integer(c_int)  :: res
    end function spmvEllDeviceDoubleComplexVM

    function spmvEllDeviceDoubleComplexMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvEllDeviceDoubleComplexMMi')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvEllDeviceDoubleComplexMMi
  end interface
end module elldev_mod
