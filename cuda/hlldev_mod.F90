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
  

module hlldev_mod
  use iso_c_binding 
  use core_mod 

  type, bind(c) :: hlldev_parms
    integer(c_int) :: element_type
    integer(c_int) :: hackSize
    integer(c_int) :: rows
    integer(c_int) :: avgNzr
    integer(c_int) :: allocsize
    integer(c_int) :: firstIndex
  end type hlldev_parms

  interface 
    function bldHllDeviceParams(hksize, rows, nzeros, allocsize, elementType, firstIndex) result(res) &
        & bind(c, name = 'bldHllDeviceParams')
      use iso_c_binding
      import :: hlldev_parms
      integer(c_int), value :: hksize, rows, nzeros, allocsize, elementType, firstIndex
      type(hlldev_parms)  :: res
    end function BldHllDeviceParams
  end interface

  interface 
    function getHllDeviceParams(deviceMat, hksize, rows, nzeros, allocsize, &
                              & hackOffsLength, firstIndex, avgnzr) result(res) &
        & bind(c, name = 'getHllDeviceParams')
      use iso_c_binding
      import :: hlldev_parms
      type(c_ptr), value  :: deviceMat
      integer(c_int)      :: hksize, rows, nzeros, allocsize, hackOffsLength, firstIndex, avgnzr
      integer(c_int)  :: res
    end function GetHllDeviceParams
  end interface

  interface 
    function FallocHllDevice(deviceMat, hksize, rows, nzeros, allocsize,  &
                          & elementType, firstIndex) result(res) &
        & bind(c, name = 'FallocHllDevice')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: hksize, rows, nzeros, allocsize, elementType, firstIndex
      type(c_ptr)           :: deviceMat
    end function FallocHllDevice
  end interface

  interface writeHllDevice 
    function writeHllDeviceFloat(deviceMat, val, ja, hkoffs, irn, idiag) result(res) &
        & bind(c, name = 'writeHllDeviceFloat')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      real(c_float)       :: val(*)
      integer(c_int)      :: ja(*), irn(*), hkoffs(*), idiag(*)
      integer(c_int)  :: res
    end function writeHllDeviceFloat

    function writeHllDeviceDouble(deviceMat, val, ja, hkoffs, irn, idiag) result(res) &
        & bind(c, name = 'writeHllDeviceDouble')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      real(c_double)      :: val(*)
      integer(c_int)      :: ja(*), irn(*), hkoffs(*), idiag(*)
      integer(c_int)  :: res
    end function writeHllDeviceDouble

    function writeHllDeviceFloatComplex(deviceMat, val, ja, hkoffs, irn, idiag) result(res) &
        & bind(c, name = 'writeHllDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_float_complex)  :: val(*)
      integer(c_int)            :: ja(*), irn(*), hkoffs(*), idiag(*)
      integer(c_int)  :: res
    end function writeHllDeviceFloatComplex

    function writeHllDeviceDoubleComplex(deviceMat, val, ja, hkoffs, irn, idiag) result(res) &
        & bind(c, name = 'writeHllDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_double_complex) :: val(*)
      integer(c_int)            :: ja(*), irn(*), hkoffs(*), idiag(*)
      integer(c_int)  :: res
    end function writeHllDeviceDoubleComplex
  end interface

  interface readHllDevice 
    function readHllDeviceFloat(deviceMat, val, ja, hkoffs, irn, idiag) result(res) &
        & bind(c, name = 'readHllDeviceFloat')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      real(c_float)       :: val(*)
      integer(c_int)      :: ja(*), irn(*), hkoffs(*), idiag(*)
      integer(c_int)  :: res
    end function readHllDeviceFloat

    function readHllDeviceDouble(deviceMat, val, ja, hkoffs, irn, idiag) result(res) &
        & bind(c, name = 'readHllDeviceDouble')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      real(c_double)      :: val(*)
      integer(c_int)      :: ja(*), irn(*), hkoffs(*), idiag(*)
      integer(c_int)  :: res
    end function readHllDeviceDouble

    function readHllDeviceFloatComplex(deviceMat, val, ja, hkoffs, irn, idiag) result(res) &
        & bind(c, name = 'readHllDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value       :: deviceMat
      complex(c_float_complex) :: val(*)
      integer(c_int)           :: ja(*), irn(*), hkoffs(*), idiag(*)
      integer(c_int)  :: res
    end function readHllDeviceFloatComplex

    function readHllDeviceDoubleComplex(deviceMat, val, ja, hkoffs, irn, idiag) result(res) &
        & bind(c, name = 'readHllDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_double_complex) :: val(*)
      integer(c_int)            :: ja(*), irn(*), hkoffs(*), idiag(*)
      integer(c_int)  :: res
    end function readHllDeviceDoubleComplex
  end interface

  interface 
    subroutine freeHllDevice(deviceMat) &
        & bind(c, name = 'freeHllDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
    end subroutine freeHllDevice
  end interface

  interface psi_CopyCooToHlg
    function psiCopyCooToHlgFloat(nr, nc, nza, hacksz, noffs, isz, irn, &
                                & hoffs, idisp, ja, val, deviceMat) result(res) &
        & bind(c, name = 'psiCopyCooToHlgFloat')
      use iso_c_binding
      integer(c_int), value :: nr, nc, nza, hacksz, noffs, isz
      integer(c_int)        :: irn(*), hoffs(*), idisp(*), ja(*)
      real(c_float)         :: val(*)
      type(c_ptr), value    :: deviceMat
      integer(c_int)  :: res
    end function psiCopyCooToHlgFloat

    function psiCopyCooToHlgDouble(nr, nc, nza, hacksz, noffs, isz, irn, &
                                & hoffs, idisp, ja, val, deviceMat) result(res) &
        & bind(c, name = 'psiCopyCooToHlgDouble')
      use iso_c_binding
      integer(c_int), value :: nr, nc, nza, hacksz, noffs, isz
      integer(c_int)        :: irn(*), hoffs(*), idisp(*), ja(*)
      real(c_double)        :: val(*)
      type(c_ptr), value    :: deviceMat
      integer(c_int)  :: res
    end function psiCopyCooToHlgDouble

    function psiCopyCooToHlgFloatComplex(nr, nc, nza, hacksz, noffs, isz, irn, &
                                        & hoffs, idisp, ja, val, deviceMat) result(res) &
        & bind(c, name = 'psiCopyCooToHlgFloatComplex')
      use iso_c_binding
      integer(c_int), value     :: nr, nc, nza, hacksz, noffs, isz
      integer(c_int)            :: irn(*), idisp(*), ja(*), hoffs(*)
      complex(c_float_complex)  :: val(*)
      type(c_ptr), value        :: deviceMat
      integer(c_int)  :: res
    end function psiCopyCooToHlgFloatComplex

    function psiCopyCooToHlgDoubleComplex(nr, nc, nza, hacksz, noffs, isz, irn, &
                                        & hoffs, idisp, ja, val, deviceMat) result(res) &
        & bind(c, name = 'psiCopyCooToHlgDoubleComplex')
      use iso_c_binding
      integer(c_int), value     :: nr, nc, nza, hacksz, noffs, isz
      integer(c_int)	          :: irn(*), idisp(*), ja(*), hoffs(*)
      complex(c_double_complex) :: val(*)
      type(c_ptr), value        :: deviceMat
      integer(c_int)  :: res
    end function psiCopyCooToHlgDoubleComplex
  end interface

  !interface 
  ! function  getHllDevicePitch(deviceMat) &
  !      & bind(c, name = 'getHllDevicePitch') result(res)
  !   use iso_c_binding
  !   type(c_ptr), value  :: deviceMat
  !   integer(c_int)      :: res
  !  end function getHllDevicePitch
  !end interface

  !interface 
  !  function  getHllDeviceMaxRowSize(deviceMat) &
  !       & bind(c, name = 'getHllDeviceMaxRowSize') result(res)
  !    use iso_c_binding
  !    type(c_ptr), value  :: deviceMat
  !    integer(c_int)      :: res
  !  end function getHllDeviceMaxRowSize
  !end interface

  interface spmvHllDevice
    function spmvHllDeviceFloatF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvHllDeviceFloatF')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int)  :: res
    end function spmvHllDeviceFloatF

    function spmvHllDeviceFloatMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvHllDeviceFloatMV')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int), value :: idx_y
      integer(c_int)  :: res
    end function spmvHllDeviceFloatMV

    function spmvHllDeviceFloatVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleVM')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int), value :: idx_x
      integer(c_int)  :: res
    end function spmvHllDeviceFloatVM

    function spmvHllDeviceFloatMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleMMi')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int), value :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvHllDeviceFloatMMi

    function spmvHllDeviceDoubleF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleF')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int)  :: res
    end function spmvHllDeviceDoubleF

    function spmvHllDeviceDoubleMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleMV')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_y
      integer(c_int)  :: res
    end function spmvHllDeviceDoubleMV

    function spmvHllDeviceDoubleVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleVM')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_x
      integer(c_int)  :: res
    end function spmvHllDeviceDoubleVM

    function spmvHllDeviceDoubleMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleMMi')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvHllDeviceDoubleMMi

    function spmvHllDeviceFloatComplexF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvHllDeviceFloatComplexF')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y
      complex(c_float_complex), value :: alpha, beta
      integer(c_int)  :: res
    end function spmvHllDeviceFloatComplexF

    function spmvHllDeviceFloatComplexMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvHllDeviceFloatComplexMV')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y
      complex(c_float_complex), value :: alpha, beta
      integer(c_int), value           :: idx_y
      integer(c_int)  :: res
    end function spmvHllDeviceFloatComplexMV

    function spmvHllDeviceFloatComplexVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvHllDeviceFloatComplexVM')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y
      complex(c_float_complex), value :: alpha, beta
      integer(c_int), value           :: idx_x
      integer(c_int)  :: res
    end function spmvHllDeviceFloatComplexVM

    function spmvHllDeviceFloatComplexMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvHllDeviceFloatComplexMMi')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y
      complex(c_float_complex), value :: alpha, beta
      integer(c_int), value           :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvHllDeviceFloatComplexMMi

    function spmvHllDeviceDoubleComplexF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleComplexF')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int)  :: res
    end function spmvHllDeviceDoubleComplexF

    function spmvHllDeviceDoubleComplexMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleComplexMV')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_y
      integer(c_int)  :: res
    end function spmvHllDeviceDoubleComplexMV

    function spmvHllDeviceDoubleComplexVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleComplexVM')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_x
      integer(c_int)  :: res
    end function spmvHllDeviceDoubleComplexVM

    function spmvHllDeviceDoubleComplexMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvHllDeviceDoubleComplexMMi')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvHllDeviceDoubleComplexMMi
  end interface
end module hlldev_mod
