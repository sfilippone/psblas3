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
  

module diagdev_mod
  use iso_c_binding 
  use core_mod 

  type, bind(c) :: diagdev_parms
    integer(c_int) :: element_type
    integer(c_int) :: rows
    integer(c_int) :: columns
    integer(c_int) :: firstIndex
 end type diagdev_parms

  interface 
    function FgetDiagDeviceParams(rows, columns, elementType, firstIndex) result(res) &
        & bind(c, name = 'getDiagDeviceParams')
      use iso_c_binding
      import :: diagdev_parms
      integer(c_int), value :: rows, columns, elementType, firstIndex
      type(diagdev_parms) :: res
    end function FgetDiagDeviceParams
  end interface
  
  interface 
    function FallocDiagDevice(deviceMat, rows, columns, elementType, firstIndex) result(res) &
        & bind(c, name = 'FallocDiagDevice')
      use iso_c_binding
      type(c_ptr)           :: deviceMat
      integer(c_int), value :: rows, columns, elementType, firstIndex
      integer(c_int)  :: res
    end function FallocDiagDevice
  end interface

  interface writeDiagDevice
    function writeDiagDeviceFloat(deviceMat, a, off, n) result(res) &
        & bind(c, name = 'writeDiagDeviceFloat')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_float)         :: a(n, *)
      integer(c_int)        :: off(*)!, irn(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function writeDiagDeviceFloat

    function writeDiagDeviceDouble(deviceMat, a, off, n) result(res) &
        & bind(c, name = 'writeDiagDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_double)        :: a(n, *)
      integer(c_int)        :: off(*)
      integer(c_int), value :: n
      integer(c_int)  :: res
    end function writeDiagDeviceDouble

    function writeDiagDeviceFloatComplex(deviceMat, a, off, n) result(res) &
        & bind(c, name = 'writeDiagDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value       :: deviceMat
      complex(c_float_complex) :: a(n, *)
      integer(c_int)           :: off(*)!, irn(*)
      integer(c_int), value    :: n
      integer(c_int)  :: res
    end function writeDiagDeviceFloatComplex

    function writeDiagDeviceDoubleComplex(deviceMat, a, off, n) result(res) &
        & bind(c, name = 'writeDiagDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_double_complex) :: a(n, *)
      integer(c_int)            :: off(*)!, irn(*)
      integer(c_int), value     :: n
      integer(c_int)  :: res
    end function writeDiagDeviceDoubleComplex
  end interface

  interface readDiagDevice 
    function readDiagDeviceFloat(deviceMat, a, off, n) result(res) &
        & bind(c, name = 'readDiagDeviceFloat')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      real(c_float)       :: a(n, *)
      integer(c_int)      :: off(*)!, irn(*)
      integer(c_int)  :: res
    end function readDiagDeviceFloat

    function readDiagDeviceDouble(deviceMat, a, off, n) result(res) &
        & bind(c, name = 'readDiagDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_double)        :: a(n, *)
      integer(c_int)        :: off(*)
      integer(c_int), value :: n
      integer(c_int)   :: res
    end function readDiagDeviceDouble

    function readDiagDeviceFloatComplex(deviceMat, a, off, n) result(res) &
        & bind(c, name = 'readDiagDeviceFloatComplex')
      use iso_c_binding
      type(c_ptr), value       :: deviceMat
      complex(c_float_complex) :: a(n, *)
      integer(c_int)           :: off(*)!, irn(*)
      integer(c_int), value    :: n
      integer(c_int)  :: res
    end function readDiagDeviceFloatComplex

    function readDiagDeviceDoubleComplex(deviceMat, a, off, n) result(res) &
        & bind(c, name = 'readDiagDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_double_complex) :: a(n, *)
      integer(c_int)            :: off(*)!, irn(*)
      integer(c_int), value     :: n
      integer(c_int)  :: res
    end function readDiagDeviceDoubleComplex
  end interface

  interface 
    subroutine freeDiagDevice(deviceMat) &
        & bind(c, name = 'freeDiagDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
    end subroutine freeDiagDevice
  end interface

  interface 
    subroutine resetDiagTimer() bind(c, name = 'resetDiagTimer')
      use iso_c_binding
    end subroutine resetDiagTimer
  end interface

  interface 
    function getDiagTimer() result(res) &
        & bind(c, name = 'getDiagTimer')
      use iso_c_binding
      real(c_double)      :: res
    end function getDiagTimer
  end interface

  interface 
    function getDiagDevicePitch(deviceMat) result(res) &
        & bind(c, name = 'getDiagDevicePitch')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      integer(c_int)  :: res
    end function getDiagDevicePitch
  end interface

  interface 
    function getDiagDeviceMaxRowSize(deviceMat) result(res) &
        & bind(c, name = 'getDiagDeviceMaxRowSize')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      integer(c_int)  :: res
    end function getDiagDeviceMaxRowSize
  end interface

  interface spmvDiagDevice
    function spmvDiagDeviceFloatF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvDiagDeviceFloatF')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int)  :: res
    end function spmvDiagDeviceFloatF

    function spmvDiagDeviceFloatMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDiagDeviceFloatMV')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int), value :: idx_y
      integer(c_int)  :: res
    end function spmvDiagDeviceFloatMV

    function spmvDiagDeviceFloatVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvDiagDeviceFloatVM')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int), value :: idx_x
      integer(c_int)  :: res
    end function spmvDiagDeviceFloatVM

    function spmvDiagDeviceFloatMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDiagDeviceFloatMMi')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_float), value  :: alpha, beta
      integer(c_int), value :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvDiagDeviceFloatMMi

    function spmvDiagDeviceDoubleF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvDiagDeviceDoubleF')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int)  :: res
    end function spmvDiagDeviceDoubleF

    function spmvDiagDeviceDoubleMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDiagDeviceDoubleMV')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_y
      integer(c_int)  :: res
    end function spmvDiagDeviceDoubleMV

    function spmvDiagDeviceDoubleVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvDiagDeviceDoubleVM')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_x
      integer(c_int)  :: res
    end function spmvDiagDeviceDoubleVM
    
    function spmvDiagDeviceDoubleMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDiagDeviceDoubleMMi')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat, x, y
      real(c_double), value :: alpha, beta
      integer(c_int), value :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvDiagDeviceDoubleMMi
    
    function spmvDiagDeviceFloatComplexF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvDiagDeviceFloatComplexF')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y
      complex(c_float_complex), value :: alpha, beta
      integer(c_int)  :: res
    end function spmvDiagDeviceFloatComplexF

    function spmvDiagDeviceFloatComplexMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDiagDeviceFloatComplexMV')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y
      complex(c_float_complex), value :: alpha, beta
      integer(c_int), value           :: idx_y
      integer(c_int)  :: res
    end function spmvDiagDeviceFloatComplexMV

    function spmvDiagDeviceFloatComplexVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvDiagDeviceFloatComplexVM')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y
      complex(c_float_complex), value :: alpha, beta
      integer(c_int), value           :: idx_x
      integer(c_int)  :: res
    end function spmvDiagDeviceFloatComplexVM

    function spmvDiagDeviceFloatComplexMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDiagDeviceFloatComplexMMi')
      use iso_c_binding
      type(c_ptr), value              :: deviceMat, x, y
      complex(c_float_complex), value :: alpha, beta
      integer(c_int), value           :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvDiagDeviceFloatComplexMMi

    function spmvDiagDeviceDoubleComplexF(deviceMat, alpha, x, beta, y) result(res) &
        & bind(c, name = 'spmvDiagDeviceDoubleComplexF')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int)  :: res
    end function spmvDiagDeviceDoubleComplexF
    
    function spmvDiagDeviceDoubleComplexMV(deviceMat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDiagDeviceDoubleComplexMV')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_y
      integer(c_int)  :: res
    end function spmvDiagDeviceDoubleComplexMV

    function spmvDiagDeviceDoubleComplexVM(deviceMat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvDiagDeviceDoubleComplexVM')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_x
      integer(c_int)  :: res
    end function spmvDiagDeviceDoubleComplexVM

    function spmvDiagDeviceDoubleComplexMMi(deviceMat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDiagDeviceDoubleComplexMMi')
      use iso_c_binding
      type(c_ptr), value                :: deviceMat, x, y
      complex(c_double_complex), value  :: alpha, beta
      integer(c_int), value             :: idx_x, idx_y
      integer(c_int)  :: res
    end function spmvDiagDeviceDoubleComplexMMi
  end interface spmvDiagDevice
end module diagdev_mod