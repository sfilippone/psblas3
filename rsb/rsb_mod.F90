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
  

module rsb_mod
  use rsb
  use iso_c_binding 

#ifdef HAVE_RSB

 interface Rsb_from_coo
    function Rsb_double_from_coo(rsbMat,va,ia,ja,nnz,nr,nc,br,bc) &
         & result(res) bind(c,name='Rsb_double_from_coo')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr) :: rsbMat
      real(c_double) :: va(*)
      integer(c_int) :: ia(*),ja(*)
      integer(c_int),value :: nnz,nr,nc,br,bc
    end function Rsb_double_from_coo
 end interface Rsb_from_coo

 interface
    function Rsb_get_nzeros(rsbMat) &
         & result(res) bind(c,name='Rsb_getNZeros')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr),value :: rsbMat
    end function Rsb_get_nzeros
 end interface

 interface Rsb_spmv
    function Rsb_double_spmv(rsbMat,x,alfa,y,beta,trans) &
         & result(res) bind(c,name='Rsb_double_spmv')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr),value :: rsbMat
      real(c_double) :: x(*),y(*)
      real(c_double),value :: alfa,beta
      character(c_char),value :: trans
    end function Rsb_double_spmv
 end interface Rsb_spmv

  interface 
    subroutine  freeRsbMat(rsbMat) &
         & bind(c,name='freeRsbMat')
      use iso_c_binding
      type(c_ptr), value  :: rsbMat
    end subroutine freeRsbMat
  end interface

  ! interface writeEllDevice
 
  !   function writeEllDeviceFloat(deviceMat,val,ja,ldj,irn) &
  !        & result(res) bind(c,name='writeEllDeviceFloat')
  !     use iso_c_binding
  !     integer(c_int)      :: res
  !     type(c_ptr), value  :: deviceMat
  !     integer(c_int), value :: ldj
  !     real(c_float)       :: val(ldj,*)
  !     integer(c_int)      :: ja(ldj,*),irn(*)
  !   end function writeEllDeviceFloat

  !   function writeEllDeviceDouble(deviceMat,val,ja,ldj,irn) &
  !        & result(res) bind(c,name='writeEllDeviceDouble')
  !     use iso_c_binding
  !     integer(c_int)      :: res
  !     type(c_ptr), value  :: deviceMat
  !     integer(c_int), value :: ldj
  !     real(c_double)      :: val(ldj,*)
  !     integer(c_int)      :: ja(ldj,*),irn(*)
  !   end function writeEllDeviceDouble

  !   function writeEllDeviceFloatComplex(deviceMat,val,ja,ldj,irn) &
  !        & result(res) bind(c,name='writeEllDeviceFloatComplex')
  !     use iso_c_binding
  !     integer(c_int)           :: res
  !     type(c_ptr), value       :: deviceMat
  !     integer(c_int), value    :: ldj
  !     complex(c_float_complex) :: val(ldj,*)
  !     integer(c_int)           :: ja(ldj,*),irn(*)
  !   end function writeEllDeviceFloatComplex

  !   function writeEllDeviceDoubleComplex(deviceMat,val,ja,ldj,irn) &
  !        & result(res) bind(c,name='writeEllDeviceDoubleComplex')
  !     use iso_c_binding
  !     integer(c_int)            :: res
  !     type(c_ptr), value        :: deviceMat
  !     integer(c_int), value     :: ldj
  !     complex(c_double_complex) :: val(ldj,*)
  !     integer(c_int)            :: ja(ldj,*),irn(*)
  !   end function writeEllDeviceDoubleComplex

  ! end interface writeEllDevice

  ! interface readEllDevice 

  !   function readEllDeviceFloat(deviceMat,val,ja,ldj,irn) &
  !        & result(res) bind(c,name='readEllDeviceFloat')
  !     use iso_c_binding
  !     integer(c_int)      :: res
  !     type(c_ptr), value  :: deviceMat
  !     integer(c_int), value :: ldj
  !     real(c_float)       :: val(ldj,*)
  !     integer(c_int)      :: ja(ldj,*),irn(*)
  !   end function readEllDeviceFloat

  !   function readEllDeviceDouble(deviceMat,val,ja,ldj,irn) &
  !        & result(res) bind(c,name='readEllDeviceDouble')
  !     use iso_c_binding
  !     integer(c_int)      :: res
  !     type(c_ptr), value  :: deviceMat
  !     integer(c_int), value :: ldj
  !     real(c_double)      :: val(ldj,*)
  !     integer(c_int)      :: ja(ldj,*),irn(*)
  !   end function readEllDeviceDouble

  !   function readEllDeviceFloatComplex(deviceMat,val,ja,ldj,irn) &
  !        & result(res) bind(c,name='readEllDeviceFloatComplex')
  !     use iso_c_binding
  !     integer(c_int)           :: res
  !     type(c_ptr), value       :: deviceMat
  !     integer(c_int), value    :: ldj
  !     complex(c_float_complex) :: val(ldj,*)
  !     integer(c_int)           :: ja(ldj,*),irn(*)
  !   end function readEllDeviceFloatComplex

  !   function readEllDeviceDoubleComplex(deviceMat,val,ja,ldj,irn) &
  !        & result(res) bind(c,name='readEllDeviceDoubleComplex')
  !     use iso_c_binding
  !     integer(c_int)           :: res
  !     type(c_ptr), value       :: deviceMat
  !     integer(c_int), value    :: ldj
  !     complex(c_double_complex) :: val(ldj,*)
  !     integer(c_int)           :: ja(ldj,*),irn(*)
  !   end function readEllDeviceDoubleComplex

  ! end interface readEllDevice

  ! interface 
  !   subroutine resetEllTimer() bind(c,name='resetEllTimer')
  !     use iso_c_binding
  !   end subroutine resetEllTimer
  ! end interface
  ! interface 
  !   function  getEllTimer() &
  !        & bind(c,name='getEllTimer') result(res)
  !     use iso_c_binding
  !     real(c_double)      :: res
  !   end function getEllTimer
  ! end interface


  ! interface 
  !   function  getEllDevicePitch(deviceMat) &
  !        & bind(c,name='getEllDevicePitch') result(res)
  !     use iso_c_binding
  !     type(c_ptr), value  :: deviceMat
  !     integer(c_int)      :: res
  !   end function getEllDevicePitch
  ! end interface

  ! interface 
  !   function  getEllDeviceMaxRowSize(deviceMat) &
  !        & bind(c,name='getEllDeviceMaxRowSize') result(res)
  !     use iso_c_binding
  !     type(c_ptr), value  :: deviceMat
  !     integer(c_int)      :: res
  !   end function getEllDeviceMaxRowSize
  ! end interface


  ! interface spmvEllDevice
  !   function spmvEllDeviceFloat(deviceMat,alpha,x,beta,y) &
  !        & result(res) bind(c,name='spmvEllDeviceFloat')
  !     use iso_c_binding
  !     integer(c_int)		:: res
  !     type(c_ptr), value 	:: deviceMat, x, y
  !     real(c_float),value     	:: alpha, beta
  !   end function spmvEllDeviceFloat
  !   function spmvEllDeviceDouble(deviceMat,alpha,x,beta,y) &
  !        & result(res) bind(c,name='spmvEllDeviceDouble')
  !     use iso_c_binding
  !     integer(c_int)		:: res
  !     type(c_ptr), value	:: deviceMat, x, y 
  !     real(c_double),value     	:: alpha,  beta
  !   end function spmvEllDeviceDouble
  !   function spmvEllDeviceFloatComplex(deviceMat,alpha,x,beta,y) &
  !        & result(res) bind(c,name='spmvEllDeviceFloatComplex')
  !     use iso_c_binding
  !     integer(c_int)		     :: res
  !     type(c_ptr), value	     :: deviceMat, x, y 
  !     complex(c_float_complex),value :: alpha,  beta
  !   end function spmvEllDeviceFloatComplex
  !   function spmvEllDeviceDoubleComplex(deviceMat,alpha,x,beta,y) &
  !        & result(res) bind(c,name='spmvEllDeviceDoubleComplex')
  !     use iso_c_binding
  !     integer(c_int)		      :: res
  !     type(c_ptr), value	      :: deviceMat, x, y 
  !     complex(c_double_complex),value :: alpha,  beta
  !   end function spmvEllDeviceDoubleComplex
  ! end interface spmvEllDevice
    
#endif


end module rsb_mod
