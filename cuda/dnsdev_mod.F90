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
!         software without specific written permission.
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
    function FgetDnsDeviceParams(rows, columns, elementType, firstIndex) &
         & result(res) bind(c,name='getDnsDeviceParams')
      use iso_c_binding
      import :: dnsdev_parms
      type(dnsdev_parms)    :: res
      integer(c_int), value :: rows,columns,elementType,firstIndex
    end function FgetDnsDeviceParams
  end interface


  interface 
    function FallocDnsDevice(deviceMat,rows,columns,&
         & elementType,firstIndex) &
         & result(res) bind(c,name='FallocDnsDevice')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: rows,columns,elementType,firstIndex
      type(c_ptr)           :: deviceMat
    end function FallocDnsDevice
  end interface


  interface writeDnsDevice

    function writeDnsDeviceFloat(deviceMat,val,lda,nc) &
         & result(res) bind(c,name='writeDnsDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: res
      type(c_ptr), value    :: deviceMat
      integer(c_int), value :: lda,nc
      real(c_float)        :: val(lda,*)
    end function writeDnsDeviceFloat


    function writeDnsDeviceDouble(deviceMat,val,lda,nc) &
         & result(res) bind(c,name='writeDnsDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: res
      type(c_ptr), value    :: deviceMat
      integer(c_int), value :: lda,nc
      real(c_double)        :: val(lda,*)
    end function writeDnsDeviceDouble


    function writeDnsDeviceFloatComplex(deviceMat,val,lda,nc) &
         & result(res) bind(c,name='writeDnsDeviceFloatComple')
      use iso_c_binding
      integer(c_int)        :: res
      type(c_ptr), value    :: deviceMat
      integer(c_int), value :: lda,nc
      complex(c_float_complex)        :: val(lda,*)
    end function writeDnsDeviceFloatComplex


    function writeDnsDeviceDoubleComplex(deviceMat,val,lda,nc) &
         & result(res) bind(c,name='writeDnsDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: res
      type(c_ptr), value    :: deviceMat
      integer(c_int), value :: lda,nc
      complex(c_double_complex)        :: val(lda,*)
    end function writeDnsDeviceDoubleComplex

  end interface

  interface readDnsDevice

    function readDnsDeviceFloat(deviceMat,val,lda,nc) &
         & result(res) bind(c,name='readDnsDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: res
      type(c_ptr), value    :: deviceMat
      integer(c_int), value :: lda,nc
      real(c_float)        :: val(lda,*)
    end function readDnsDeviceFloat


    function readDnsDeviceDouble(deviceMat,val,lda,nc) &
         & result(res) bind(c,name='readDnsDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: res
      type(c_ptr), value    :: deviceMat
      integer(c_int), value :: lda,nc
      real(c_double)        :: val(lda,*)
    end function readDnsDeviceDouble


    function readDnsDeviceFloatComplex(deviceMat,val,lda,nc) &
         & result(res) bind(c,name='readDnsDeviceFloatComple')
      use iso_c_binding
      integer(c_int)        :: res
      type(c_ptr), value    :: deviceMat
      integer(c_int), value :: lda,nc
      complex(c_float_complex)        :: val(lda,*)
    end function readDnsDeviceFloatComplex


    function readDnsDeviceDoubleComplex(deviceMat,val,lda,nc) &
         & result(res) bind(c,name='readDnsDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: res
      type(c_ptr), value    :: deviceMat
      integer(c_int), value :: lda,nc
      complex(c_double_complex)        :: val(lda,*)
    end function readDnsDeviceDoubleComplex

  end interface

  interface 
    subroutine  freeDnsDevice(deviceMat) &
         & bind(c,name='freeDnsDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
    end subroutine freeDnsDevice
  end interface

  interface 
    subroutine resetDnsTimer() bind(c,name='resetDnsTimer')
      use iso_c_binding
    end subroutine resetDnsTimer
  end interface
  interface 
    function  getDnsTimer() &
         & bind(c,name='getDnsTimer') result(res)
      use iso_c_binding
      real(c_double)      :: res
    end function getDnsTimer
  end interface


  interface 
    function  getDnsDevicePitch(deviceMat) &
         & bind(c,name='getDnsDevicePitch') result(res)
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
      integer(c_int)      :: res
    end function getDnsDevicePitch
  end interface

!!$  interface csputDnsDeviceFloat
!!$    function dev_csputDnsDeviceFloat(deviceMat, nnz, ia, ja, val) &
!!$         & result(res) bind(c,name='dev_csputDnsDeviceFloat')
!!$      use iso_c_binding
!!$      integer(c_int)		:: res
!!$      type(c_ptr), value 	:: deviceMat , ia, ja, val
!!$      integer(c_int), value     :: nnz 
!!$    end function dev_csputDnsDeviceFloat
!!$  end interface
!!$
!!$  interface csputDnsDeviceDouble
!!$    function dev_csputDnsDeviceDouble(deviceMat, nnz, ia, ja, val) &
!!$         & result(res) bind(c,name='dev_csputDnsDeviceDouble')
!!$      use iso_c_binding
!!$      integer(c_int)		:: res
!!$      type(c_ptr), value 	:: deviceMat , ia, ja, val
!!$      integer(c_int), value     :: nnz 
!!$    end function dev_csputDnsDeviceDouble
!!$  end interface
!!$
!!$  interface csputDnsDeviceFloatComplex
!!$    function dev_csputDnsDeviceFloatComplex(deviceMat, nnz, ia, ja, val) &
!!$         & result(res) bind(c,name='dev_csputDnsDeviceFloatComplex')
!!$      use iso_c_binding
!!$      integer(c_int)		:: res
!!$      type(c_ptr), value 	:: deviceMat , ia, ja, val
!!$      integer(c_int), value     :: nnz 
!!$    end function dev_csputDnsDeviceFloatComplex
!!$  end interface
!!$
!!$  interface csputDnsDeviceDoubleComplex
!!$    function dev_csputDnsDeviceDoubleComplex(deviceMat, nnz, ia, ja, val) &
!!$         & result(res) bind(c,name='dev_csputDnsDeviceDoubleComplex')
!!$      use iso_c_binding
!!$      integer(c_int)		:: res
!!$      type(c_ptr), value 	:: deviceMat , ia, ja, val
!!$      integer(c_int), value     :: nnz 
!!$    end function dev_csputDnsDeviceDoubleComplex
!!$  end interface

  interface spmvDnsDevice
    function spmvDnsDeviceFloat(transa,m,n,k,alpha,deviceMat,x,beta,y) &
         & result(res) bind(c,name='spmvDnsDeviceFloat')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value    :: m, n, k
      integer(c_int)		:: res
      type(c_ptr), value	:: deviceMat, x, y 
      real(c_float)             :: alpha,  beta
    end function spmvDnsDeviceFloat

    function spmvDnsDeviceDouble(transa,m,n,k,alpha,deviceMat,x,beta,y) &
         & result(res) bind(c,name='spmvDnsDeviceDouble')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k
      integer(c_int)		:: res
      type(c_ptr), value	:: deviceMat, x, y 
      real(c_double)            :: alpha,  beta
    end function spmvDnsDeviceDouble

    function spmvDnsDeviceFloatComplex(transa,m,n,k,alpha,deviceMat,x,beta,y) &
         & result(res) bind(c,name='spmvDnsDeviceFloatComplex')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k
      integer(c_int)	        :: res
      type(c_ptr), value        :: deviceMat, x, y 
      complex(c_float_complex)  :: alpha,  beta
    end function spmvDnsDeviceFloatComplex

    function spmvDnsDeviceDoubleComplex(transa,m,n,k,alpha,deviceMat,x,beta,y) &
         & result(res) bind(c,name='spmvDnsDeviceDoubleComplex')
      use iso_c_binding
      character(c_char), value  :: transa
      integer(c_int), value     :: m, n, k
      integer(c_int)	        :: res
      type(c_ptr), value        :: deviceMat, x, y 
      complex(c_double_complex) :: alpha,  beta
    end function spmvDnsDeviceDoubleComplex

  end interface

end module dnsdev_mod
