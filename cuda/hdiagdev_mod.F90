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
  

module hdiagdev_mod
  use iso_c_binding 
  use core_mod 

  type, bind(c) :: hdiagdev_parms
    integer(c_int) :: element_type
    integer(c_int) :: rows
    integer(c_int) :: columns
    integer(c_int) :: hackSize
    integer(c_int) :: hackCount
    integer(c_int) :: allocationHeight    
 end type hdiagdev_parms

#ifdef HAVE_SPGPU  

 ! interface computeHdiaHacksCount
 !    function computeHdiaHacksCountDouble(allocationHeight,hackOffsets,hackSize, &
 !         & diaValues,diaValuesPitch,diags,rows)&
 !         & result(res) bind(c,name='computeHdiaHackOffsetsDouble')
 !      use iso_c_binding
 !      integer(c_int)        :: res
 !      integer(c_int), value :: rows,diags,diaValuesPitch,hackSize,elementType
 !      real(c_double)        :: diaValues(rows,:)
 !      integer(c_int)        :: hackOffsets,allocationHeight
 !    end function computeHdiaHacksCountDouble
 ! end interface computeHdiaHacksCount

  interface 
    function FgetHdiagDeviceParams(rows, columns, allocationHeight,hackSize, &
         & hackCount, elementType) &
         & result(res) bind(c,name='getHdiagDeviceParams')
      use iso_c_binding
      import :: hdiagdev_parms
      type(hdiagdev_parms)    :: res
      integer(c_int), value :: rows,columns,allocationHeight,&
           & elementType,hackSize,hackCount
    end function FgetHdiagDeviceParams
  end interface
  

  interface 
    function FallocHdiagDevice(deviceMat,rows,columns,allocationHeight,&
         & hackSize,hackCount,elementType) &
         & result(res) bind(c,name='FallocHdiagDevice')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: rows,columns,allocationHeight,hackSize,&
           & hackCount,elementType
      type(c_ptr)           :: deviceMat
    end function FallocHdiagDevice
  end interface


  interface 
    function sizeofHdiagDeviceDouble(deviceMat) &
         & result(res) bind(c,name='sizeofHdiagDeviceDouble') 
      use iso_c_binding
      integer(c_long_long) :: res
      type(c_ptr), value  :: deviceMat
    end function sizeofHdiagDeviceDouble
  end interface

  interface writeHdiagDevice

    function writeHdiagDeviceFloat(deviceMat,val,hdiaOffsets, hackOffsets) &
         & result(res) bind(c,name='writeHdiagDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceMat
      real(c_float)      :: val(*)
      integer(c_int)      :: hdiaOffsets(*), hackOffsets(*)
    end function writeHdiagDeviceFloat

    function writeHdiagDeviceDouble(deviceMat,val,hdiaOffsets, hackOffsets) &
         & result(res) bind(c,name='writeHdiagDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceMat
      real(c_double)      :: val(*)
      integer(c_int)      :: hdiaOffsets(*), hackOffsets(*)
    end function writeHdiagDeviceDouble

  end interface writeHdiagDevice

!!$  interface readHdiagDevice 
!!$
!!$    function readHdiagDeviceFloat(deviceMat,val,ja,ldj,irn) &
!!$         & result(res) bind(c,name='readHdiagDeviceFloat')
!!$      use iso_c_binding
!!$      integer(c_int)      :: res
!!$      type(c_ptr), value  :: deviceMat
!!$      integer(c_int), value :: ldj
!!$      real(c_float)       :: val(ldj,*)
!!$      integer(c_int)      :: ja(ldj,*),irn(*)
!!$    end function readHdiagDeviceFloat
!!$
!!$    function readHdiagDeviceDouble(deviceMat,a,off,n) &
!!$         & result(res) bind(c,name='readHdiagDeviceDouble')
!!$      use iso_c_binding
!!$      integer(c_int)      :: res
!!$      type(c_ptr), value  :: deviceMat
!!$      integer(c_int),value :: n
!!$      real(c_double)      :: a(n,*)
!!$      integer(c_int)      :: off(*)
!!$    end function readHdiagDeviceDouble
!!$
!!$    function readHdiagDeviceFloatComplex(deviceMat,val,ja,ldj,irn) &
!!$         & result(res) bind(c,name='readHdiagDeviceFloatComplex')
!!$      use iso_c_binding
!!$      integer(c_int)           :: res
!!$      type(c_ptr), value       :: deviceMat
!!$      integer(c_int), value    :: ldj
!!$      complex(c_float_complex) :: val(ldj,*)
!!$      integer(c_int)           :: ja(ldj,*),irn(*)
!!$    end function readHdiagDeviceFloatComplex
!!$
!!$    function readHdiagDeviceDoubleComplex(deviceMat,val,ja,ldj,irn) &
!!$         & result(res) bind(c,name='readHdiagDeviceDoubleComplex')
!!$      use iso_c_binding
!!$      integer(c_int)           :: res
!!$      type(c_ptr), value       :: deviceMat
!!$      integer(c_int), value    :: ldj
!!$      complex(c_double_complex) :: val(ldj,*)
!!$      integer(c_int)           :: ja(ldj,*),irn(*)
!!$    end function readHdiagDeviceDoubleComplex
!!$
!!$  end interface readHdiagDevice
!!$
  interface 
    subroutine  freeHdiagDevice(deviceMat) &
         & bind(c,name='freeHdiagDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceMat
    end subroutine freeHdiagDevice
  end interface


  interface spmvHdiagDevice
    function spmvHdiagDeviceFloat(deviceMat,alpha,x,beta,y) &
         & result(res) bind(c,name='spmvHdiagDeviceFloat')
      use iso_c_binding
      integer(c_int)		:: res
      type(c_ptr), value 	:: deviceMat, x, y
      real(c_float),value     	:: alpha, beta
    end function spmvHdiagDeviceFloat
    function spmvHdiagDeviceDouble(deviceMat,alpha,x,beta,y) &
         & result(res) bind(c,name='spmvHdiagDeviceDouble')
      use iso_c_binding
      integer(c_int)		:: res
      type(c_ptr), value	:: deviceMat, x, y 
      real(c_double),value     	:: alpha,  beta
    end function spmvHdiagDeviceDouble
!!$    function spmvHdiagDeviceFloatComplex(deviceMat,alpha,x,beta,y) &
!!$         & result(res) bind(c,name='spmvHdiagDeviceFloatComplex')
!!$      use iso_c_binding
!!$      integer(c_int)		     :: res
!!$      type(c_ptr), value	     :: deviceMat, x, y 
!!$      complex(c_float_complex),value :: alpha,  beta
!!$    end function spmvHdiagDeviceFloatComplex
!!$    function spmvHdiagDeviceDoubleComplex(deviceMat,alpha,x,beta,y) &
!!$         & result(res) bind(c,name='spmvHdiagDeviceDoubleComplex')
!!$      use iso_c_binding
!!$      integer(c_int)		      :: res
!!$      type(c_ptr), value	      :: deviceMat, x, y 
!!$      complex(c_double_complex),value :: alpha,  beta
!!$    end function spmvHdiagDeviceDoubleComplex
  end interface spmvHdiagDevice
    
#endif  

end module hdiagdev_mod
