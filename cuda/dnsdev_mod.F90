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

module dnsdev_mod
  use iso_c_binding 
  use core_mod 
  
  type, bind(c) :: dnsdev_parms
    integer(c_int) :: element_type
    integer(c_int) :: pitch
    integer(c_int) :: rows
    integer(c_int) :: columns
    integer(c_int) :: maxRowSize
    integer(c_int) :: avgRowSize
    integer(c_int) :: firstIndex
  end type dnsdev_parms

  interface 
    function FgetDnsDeviceParams(rows, columns, elementType, firstIndex) result(res) &
        & bind(c, name = 'getDnsDeviceParams')
      use iso_c_binding
      import :: dnsdev_parms
      integer(c_int), value :: rows, columns, elementType, firstIndex
      type(dnsdev_parms)  :: res
    end function FgetDnsDeviceParams
  end interface

  interface 
    function FallocDnsDevice(deviceMat, rows, columns, elementType, firstIndex) result(res) &
        & bind(c, name = 'FallocDnsDevice')
      use iso_c_binding
      type(c_ptr)           :: deviceMat
      integer(c_int), value :: rows, columns, elementType, firstIndex
      integer(c_int)  :: res
    end function FallocDnsDevice
  end interface

  interface writeDnsDevice
    function writeDnsDeviceFloat(deviceMat, val, lda, nc) result(res) &
        & bind(c, name = 'writeDnsDeviceFloat')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_float)         :: val(lda, *)
      integer(c_int), value :: lda, nc
      integer(c_int)  :: res
    end function writeDnsDeviceFloat

    function writeDnsDeviceDouble(deviceMat, val, lda, nc) result(res) &
        & bind(c, name = 'writeDnsDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_double)        :: val(lda, *)
      integer(c_int), value :: lda, nc
      integer(c_int)  :: res
    end function writeDnsDeviceDouble

    function writeDnsDeviceFloatComplex(deviceMat, val, lda, nc) result(res) &
        & bind(c, name = 'writeDnsDeviceFloatComple')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_float_complex)  :: val(lda, *)
      integer(c_int), value     :: lda, nc
      integer(c_int)  :: res
    end function writeDnsDeviceFloatComplex

    function writeDnsDeviceDoubleComplex(deviceMat, val, lda, nc) result(res) &
        & bind(c, name = 'writeDnsDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_double_complex) :: val(lda, *)
      integer(c_int), value     :: lda, nc
      integer(c_int)  :: res
    end function writeDnsDeviceDoubleComplex
  end interface

  interface readDnsDevice
    function readDnsDeviceFloat(deviceMat, val, lda, nc) result(res) &
        & bind(c, name = 'readDnsDeviceFloat')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_float)         :: val(lda, *)
      integer(c_int), value :: lda, nc
      integer(c_int)  :: res
    end function readDnsDeviceFloat

    function readDnsDeviceDouble(deviceMat, val, lda, nc) result(res) &
        & bind(c, name = 'readDnsDeviceDouble')
      use iso_c_binding
      type(c_ptr), value    :: deviceMat
      real(c_double)        :: val(lda, *)
      integer(c_int), value :: lda, nc
      integer(c_int)  :: res
    end function readDnsDeviceDouble

    function readDnsDeviceFloatComplex(deviceMat, val, lda, nc) result(res) &
        & bind(c, name = 'readDnsDeviceFloatComple')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_float_complex)  :: val(lda, *)
      integer(c_int), value     :: lda, nc
      integer(c_int)  :: res
    end function readDnsDeviceFloatComplex

    function readDnsDeviceDoubleComplex(deviceMat, val, lda, nc) result(res) &
        & bind(c, name = 'readDnsDeviceDoubleComplex')
      use iso_c_binding
      type(c_ptr), value        :: deviceMat
      complex(c_double_complex) :: val(lda, *)
      integer(c_int), value     :: lda, nc
      integer(c_int)  :: res
    end function readDnsDeviceDoubleComplex
  end interface

  interface 
    subroutine freeDnsDevice(deviceMat) &
        & bind(c, name = 'freeDnsDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
    end subroutine freeDnsDevice
  end interface

  interface 
    subroutine resetDnsTimer() bind(c, name = 'resetDnsTimer')
      use iso_c_binding
    end subroutine resetDnsTimer
  end interface

  interface 
    function getDnsTimer() result(res) &
        & bind(c, name = 'getDnsTimer')
      use iso_c_binding
      real(c_double)  :: res
    end function getDnsTimer
  end interface

  interface 
    function getDnsDevicePitch(deviceMat) result(res) &
        & bind(c, name = 'getDnsDevicePitch')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      integer(c_int)  :: res
    end function getDnsDevicePitch
  end interface

  !!$ interface csputDnsDeviceFloat
  !!$   function dev_csputDnsDeviceFloat(deviceMat, nnz, ia, ja, val) &
  !!$        & result(res) bind(c, name = 'dev_csputDnsDeviceFloat')
  !!$     use iso_c_binding
  !!$     integer(c_int)		:: res
  !!$     type(c_ptr), value 	:: deviceMat , ia, ja, val
  !!$     integer(c_int), value     :: nnz 
  !!$   end function dev_csputDnsDeviceFloat
  !!$ end interface
  !!$
  !!$ interface csputDnsDeviceDouble
  !!$   function dev_csputDnsDeviceDouble(deviceMat, nnz, ia, ja, val) &
  !!$        & result(res) bind(c, name = 'dev_csputDnsDeviceDouble')
  !!$     use iso_c_binding
  !!$     integer(c_int)		:: res
  !!$     type(c_ptr), value 	:: deviceMat , ia, ja, val
  !!$     integer(c_int), value     :: nnz 
  !!$   end function dev_csputDnsDeviceDouble
  !!$ end interface
  !!$
  !!$ interface csputDnsDeviceFloatComplex
  !!$   function dev_csputDnsDeviceFloatComplex(deviceMat, nnz, ia, ja, val) &
  !!$        & result(res) bind(c, name = 'dev_csputDnsDeviceFloatComplex')
  !!$     use iso_c_binding
  !!$     integer(c_int)		:: res
  !!$     type(c_ptr), value 	:: deviceMat , ia, ja, val
  !!$     integer(c_int), value     :: nnz 
  !!$   end function dev_csputDnsDeviceFloatComplex
  !!$ end interface
  !!$
  !!$ interface csputDnsDeviceDoubleComplex
  !!$   function dev_csputDnsDeviceDoubleComplex(deviceMat, nnz, ia, ja, val) &
  !!$        & result(res) bind(c, name = 'dev_csputDnsDeviceDoubleComplex')
  !!$     use iso_c_binding
  !!$     integer(c_int)		:: res
  !!$     type(c_ptr), value 	:: deviceMat , ia, ja, val
  !!$     integer(c_int), value     :: nnz 
  !!$   end function dev_csputDnsDeviceDoubleComplex
  !!$ end interface

  interface spmvDnsDevice
    function spmvDnsDeviceFloatF(transa, m, n, k, alpha, deviceMat, x, beta, y) result(res) &
        & bind(c, name = 'spmvDnsDeviceFloatF')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k
      real(c_float)             :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceFloatF

    function spmvDnsDeviceFloatMV(transa, m, n, k, alpha, deviceMat, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDnsDeviceFloatMV')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_y
      real(c_float)             :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceFloatMV

    function spmvDnsDeviceFloatVM(transa, m, n, k, alpha, deviceMat, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvDnsDeviceFloatVM')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_x
      real(c_float)             :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceFloatVM

    function spmvDnsDeviceFloatMMi(transa, m, n, k, alpha, deviceMat, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDnsDeviceFloatMMi')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_x, idx_y
      real(c_float)             :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceFloatMMi

    function spmvDnsDeviceDoubleF(transa, m, n, k, alpha, deviceMat, x, beta, y) result(res) &
        & bind(c, name = 'spmvDnsDeviceDoubleF')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k
      real(c_double)            :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceDoubleF

    function spmvDnsDeviceDoubleMV(transa, m, n, k, alpha, deviceMat, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDnsDeviceDoubleMV')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_y
      real(c_double)            :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceDoubleMV

    function spmvDnsDeviceDoubleVM(transa, m, n, k, alpha, deviceMat, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvDnsDeviceDoubleVM')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_x
      real(c_double)            :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceDoubleVM

    function spmvDnsDeviceDoubleMMi(transa, m, n, k, alpha, deviceMat, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDnsDeviceDoubleMMi')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_x, idx_y
      real(c_double)            :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceDoubleMMi

    function spmvDnsDeviceFloatComplexF(transa, m, n, k, alpha, deviceMat, x, beta, y) result(res) &
        & bind(c, name = 'spmvDnsDeviceFloatComplexF')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k
      complex(c_float_complex)  :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceFloatComplexF

    function spmvDnsDeviceFloatComplexMV(transa, m, n, k, alpha, deviceMat, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDnsDeviceFloatComplexMV')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_y
      complex(c_float_complex)  :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceFloatComplexMV

    function spmvDnsDeviceFloatComplexVM(transa, m, n, k, alpha, deviceMat, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvDnsDeviceFloatComplexVM')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_x
      complex(c_float_complex)  :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceFloatComplexVM

    function spmvDnsDeviceFloatComplexMMi(transa, m, n, k, alpha, deviceMat, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDnsDeviceFloatComplexMMi')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_x, idx_y
      complex(c_float_complex)  :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceFloatComplexMMi

    function spmvDnsDeviceDoubleComplexF(transa, m, n, k, alpha, deviceMat, x, beta, y) result(res) &
        & bind(c, name = 'spmvDnsDeviceDoubleComplexF')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k
      complex(c_double_complex) :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceDoubleComplexF

    function spmvDnsDeviceDoubleComplexMV(transa, m, n, k, alpha, deviceMat, x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDnsDeviceDoubleComplexMV')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_y
      complex(c_double_complex) :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceDoubleComplexMV

    function spmvDnsDeviceDoubleComplexVM(transa, m, n, k, alpha, deviceMat, x, idx_x, beta, y) result(res) &
        & bind(c, name = 'spmvDnsDeviceDoubleComplexVM')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_x
      complex(c_double_complex) :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceDoubleComplexVM

    function spmvDnsDeviceDoubleComplexMMi(transa, m, n, k, alpha, deviceMat, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = 'spmvDnsDeviceDoubleComplexMMi')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k, idx_x, idx_y
      complex(c_double_complex) :: alpha, beta
      type(c_ptr), value        :: deviceMat, x, y
      integer(c_int)  :: res
    end function spmvDnsDeviceDoubleComplexMMi
  end interface
end module dnsdev_mod
