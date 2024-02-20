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
module psb_d_vectordev_mod

  use psb_base_vectordev_mod
 
  interface registerMapped
    function registerMappedDouble(buf,d_p,n,dummy) &
         & result(res) bind(c,name='registerMappedDouble')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: buf
      type(c_ptr) :: d_p
      integer(c_int),value :: n
      real(c_double), value :: dummy
    end function registerMappedDouble
  end interface

  interface writeMultiVecDevice 
    function writeMultiVecDeviceDouble(deviceVec,hostVec) &
         & result(res) bind(c,name='writeMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)             :: res
      type(c_ptr), value         :: deviceVec
      real(c_double)   :: hostVec(*)
    end function writeMultiVecDeviceDouble
    function writeMultiVecDeviceDoubleR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='writeMultiVecDeviceDoubleR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      real(c_double)      :: hostVec(ld,*)
    end function writeMultiVecDeviceDoubleR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceDouble(deviceVec,hostVec) &
         & result(res) bind(c,name='readMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)           :: res
      type(c_ptr), value       :: deviceVec
      real(c_double) :: hostVec(*)
    end function readMultiVecDeviceDouble
    function readMultiVecDeviceDoubleR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='readMultiVecDeviceDoubleR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      real(c_double)      :: hostVec(ld,*)
    end function readMultiVecDeviceDoubleR2
  end interface 

  interface allocateDouble
    function allocateDouble(didx,n) &
         & result(res) bind(c,name='allocateDouble') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: n
      integer(c_int)  :: res
    end function allocateDouble
    function allocateMultiDouble(didx,m,n) &
         & result(res) bind(c,name='allocateMultiDouble') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: m,n
      integer(c_int)  :: res
    end function allocateMultiDouble
  end interface

  interface writeDouble
    function writeDouble(didx,hidx,n) &
         & result(res) bind(c,name='writeDouble')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      real(c_double)       :: hidx(*)
      integer(c_int),value :: n
    end function writeDouble
    function writeDoubleFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='writeDoubleFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      real(c_double)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function writeDoubleFirst
    function writeMultiDouble(didx,hidx,m,n) &
         & result(res) bind(c,name='writeMultiDouble')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      real(c_double)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function writeMultiDouble
  end interface
  
  interface readDouble
    function readDouble(didx,hidx,n) &
         & result(res) bind(c,name='readDouble')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      real(c_double)       :: hidx(*)
      integer(c_int),value :: n
    end function readDouble
    function readDoubleFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='readDoubleFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      real(c_double)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function readDoubleFirst
    function readMultiDouble(didx,hidx,m,n) &
         & result(res) bind(c,name='readMultiDouble')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      real(c_double)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function readMultiDouble
  end interface
  
  interface
    subroutine  freeDouble(didx) &
         & bind(c,name='freeDouble')
      use iso_c_binding
      type(c_ptr), value :: didx
    end subroutine freeDouble
  end interface
  

  interface setScalDevice
    function setScalMultiVecDeviceDouble(val, first, last, &
         & indexBase, deviceVecX) result(res) &
         & bind(c,name='setscalMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: first,last,indexbase
      real(c_double), value :: val
      type(c_ptr),   value  :: deviceVecX
    end function setScalMultiVecDeviceDouble
  end interface

  interface 
    function geinsMultiVecDeviceDouble(n,deviceVecIrl,deviceVecVal,&
         & dupl,indexbase,deviceVecX) &
         & result(res) bind(c,name='geinsMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n, dupl,indexbase
      type(c_ptr), value  :: deviceVecIrl, deviceVecVal, deviceVecX
    end function geinsMultiVecDeviceDouble
  end interface

  ! New gather functions

  interface 
    function igathMultiVecDeviceDouble(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceDouble
  end interface

  interface 
    function igathMultiVecDeviceDoubleVecIdx(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceDoubleVecIdx')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceDoubleVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceDouble(deviceVec, vectorId, & 
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      real(c_double),value :: beta
    end function iscatMultiVecDeviceDouble
  end interface

  interface 
    function iscatMultiVecDeviceDoubleVecIdx(deviceVec, vectorId, &
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceDoubleVecIdx')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      real(c_double),value :: beta
    end function iscatMultiVecDeviceDoubleVecIdx
  end interface


  interface scalMultiVecDevice
    function scalMultiVecDeviceDouble(alpha,deviceVecA) &
         & result(val) bind(c,name='scalMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: res
      real(c_double), value :: alpha
      type(c_ptr), value    :: deviceVecA
    end function scalMultiVecDeviceDouble
  end interface

  interface dotMultiVecDevice
    function dotMultiVecDeviceDouble(res, n,deviceVecA,deviceVecB) &
         & result(val) bind(c,name='dotMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_double) :: res
      type(c_ptr), value    :: deviceVecA, deviceVecB
    end function dotMultiVecDeviceDouble
  end interface
    
  interface nrm2MultiVecDevice
    function nrm2MultiVecDeviceDouble(res,n,deviceVecA) &
         & result(val) bind(c,name='nrm2MultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_double)         :: res
      type(c_ptr), value    :: deviceVecA
    end function nrm2MultiVecDeviceDouble
  end interface

  interface amaxMultiVecDevice
    function amaxMultiVecDeviceDouble(res,n,deviceVecA) &
         & result(val) bind(c,name='amaxMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_double) :: res
      type(c_ptr), value    :: deviceVecA
    end function amaxMultiVecDeviceDouble
  end interface

  interface asumMultiVecDevice
    function asumMultiVecDeviceDouble(res,n,deviceVecA) &
         & result(val) bind(c,name='asumMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_double) :: res
      type(c_ptr), value    :: deviceVecA
    end function asumMultiVecDeviceDouble
  end interface

  interface axpbyMultiVecDevice
    function axpbyMultiVecDeviceDouble(n,alpha,deviceVecA,beta,deviceVecB) &
         & result(res) bind(c,name='axpbyMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      real(c_double), value :: alpha, beta
      type(c_ptr), value  :: deviceVecA, deviceVecB
    end function axpbyMultiVecDeviceDouble
  end interface

  interface abgdxyzMultiVecDevice
    function abgdxyzMultiVecDeviceDouble(n,alpha,beta,gamma,delta,deviceVecX,&
         & deviceVecY,deviceVecZ) &
         & result(res) bind(c,name='abgdxyzMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      real(c_double), value :: alpha, beta,gamma,delta
      type(c_ptr), value  :: deviceVecX, deviceVecY, deviceVecZ
    end function abgdxyzMultiVecDeviceDouble
  end interface

  interface axyMultiVecDevice
    function axyMultiVecDeviceDouble(n,alpha,deviceVecA,deviceVecB) &
         & result(res) bind(c,name='axyMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value       :: deviceVecA, deviceVecB
    end function axyMultiVecDeviceDouble
  end interface

  interface axybzMultiVecDevice
    function axybzMultiVecDeviceDouble(n,alpha,deviceVecA,deviceVecB,beta,deviceVecZ) &
         & result(res) bind(c,name='axybzMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)              :: res
      integer(c_int), value       :: n
      real(c_double), value     :: alpha, beta
      type(c_ptr), value          :: deviceVecA, deviceVecB,deviceVecZ
    end function axybzMultiVecDeviceDouble
  end interface


  interface absMultiVecDevice
    function absMultiVecDeviceDouble(n,alpha,deviceVecA) &
         & result(res) bind(c,name='absMultiVecDeviceDouble')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value    :: deviceVecA
    end function absMultiVecDeviceDouble
    function absMultiVecDeviceDouble2(n,alpha,deviceVecA,deviceVecB) &
         & result(res) bind(c,name='absMultiVecDeviceDouble2')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      real(c_double), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
    end function absMultiVecDeviceDouble2
  end interface

  interface inner_register
    module procedure inner_registerDouble
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterDouble
  end interface

contains


  function inner_registerDouble(buffer,dval) result(res)
    real(c_double), allocatable, target :: buffer(:)
    type(c_ptr)            :: dval
    integer(c_int)         :: res
    real(c_double)         :: dummy
    res = registerMapped(c_loc(buffer),dval,size(buffer), dummy)        
  end function inner_registerDouble

  subroutine inner_unregisterDouble(buffer)
    real(c_double), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterDouble

end module psb_d_vectordev_mod
