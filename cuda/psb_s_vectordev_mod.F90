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
module psb_s_vectordev_mod

  use psb_base_vectordev_mod
 
  interface registerMapped
    function registerMappedFloat(buf,d_p,n,dummy) &
         & result(res) bind(c,name='registerMappedFloat')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: buf
      type(c_ptr) :: d_p
      integer(c_int),value :: n
      real(c_float), value :: dummy
    end function registerMappedFloat
  end interface

  interface writeMultiVecDevice 
    function writeMultiVecDeviceFloat(deviceVec,hostVec) &
         & result(res) bind(c,name='writeMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)             :: res
      type(c_ptr), value         :: deviceVec
      real(c_float)   :: hostVec(*)
    end function writeMultiVecDeviceFloat
    function writeMultiVecDeviceFloatR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='writeMultiVecDeviceFloatR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      real(c_float)      :: hostVec(ld,*)
    end function writeMultiVecDeviceFloatR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceFloat(deviceVec,hostVec) &
         & result(res) bind(c,name='readMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)           :: res
      type(c_ptr), value       :: deviceVec
      real(c_float) :: hostVec(*)
    end function readMultiVecDeviceFloat
    function readMultiVecDeviceFloatR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='readMultiVecDeviceFloatR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      real(c_float)      :: hostVec(ld,*)
    end function readMultiVecDeviceFloatR2
  end interface 

  interface allocateFloat
    function allocateFloat(didx,n) &
         & result(res) bind(c,name='allocateFloat') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: n
      integer(c_int)  :: res
    end function allocateFloat
    function allocateMultiFloat(didx,m,n) &
         & result(res) bind(c,name='allocateMultiFloat') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: m,n
      integer(c_int)  :: res
    end function allocateMultiFloat
  end interface

  interface writeFloat
    function writeFloat(didx,hidx,n) &
         & result(res) bind(c,name='writeFloat')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      real(c_float)       :: hidx(*)
      integer(c_int),value :: n
    end function writeFloat
    function writeFloatFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='writeFloatFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      real(c_float)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function writeFloatFirst
    function writeMultiFloat(didx,hidx,m,n) &
         & result(res) bind(c,name='writeMultiFloat')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      real(c_float)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function writeMultiFloat
  end interface
  
  interface readFloat
    function readFloat(didx,hidx,n) &
         & result(res) bind(c,name='readFloat')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      real(c_float)       :: hidx(*)
      integer(c_int),value :: n
    end function readFloat
    function readFloatFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='readFloatFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      real(c_float)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function readFloatFirst
    function readMultiFloat(didx,hidx,m,n) &
         & result(res) bind(c,name='readMultiFloat')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      real(c_float)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function readMultiFloat
  end interface
  
  interface
    subroutine  freeFloat(didx) &
         & bind(c,name='freeFloat')
      use iso_c_binding
      type(c_ptr), value :: didx
    end subroutine freeFloat
  end interface
  

  interface setScalDevice
    function setScalMultiVecDeviceFloat(val, first, last, &
         & indexBase, deviceVecX) result(res) &
         & bind(c,name='setscalMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: first,last,indexbase
      real(c_float), value :: val
      type(c_ptr),   value  :: deviceVecX
    end function setScalMultiVecDeviceFloat
  end interface

  interface 
    function geinsMultiVecDeviceFloat(n,deviceVecIrl,deviceVecVal,&
         & dupl,indexbase,deviceVecX) &
         & result(res) bind(c,name='geinsMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n, dupl,indexbase
      type(c_ptr), value  :: deviceVecIrl, deviceVecVal, deviceVecX
    end function geinsMultiVecDeviceFloat
  end interface

  ! New gather functions

  interface 
    function igathMultiVecDeviceFloat(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceFloat
  end interface

  interface 
    function igathMultiVecDeviceFloatVecIdx(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceFloatVecIdx')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceFloatVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceFloat(deviceVec, vectorId, & 
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      real(c_float),value :: beta
    end function iscatMultiVecDeviceFloat
  end interface

  interface 
    function iscatMultiVecDeviceFloatVecIdx(deviceVec, vectorId, &
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceFloatVecIdx')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      real(c_float),value :: beta
    end function iscatMultiVecDeviceFloatVecIdx
  end interface


  interface scalMultiVecDevice
    function scalMultiVecDeviceFloat(alpha,deviceVecA) &
         & result(val) bind(c,name='scalMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: res
      real(c_float), value :: alpha
      type(c_ptr), value    :: deviceVecA
    end function scalMultiVecDeviceFloat
  end interface

  interface dotMultiVecDevice
    function dotMultiVecDeviceFloat(res, n,deviceVecA,deviceVecB) &
         & result(val) bind(c,name='dotMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_float) :: res
      type(c_ptr), value    :: deviceVecA, deviceVecB
    end function dotMultiVecDeviceFloat
  end interface
    
  interface nrm2MultiVecDevice
    function nrm2MultiVecDeviceFloat(res,n,deviceVecA) &
         & result(val) bind(c,name='nrm2MultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_float)         :: res
      type(c_ptr), value    :: deviceVecA
    end function nrm2MultiVecDeviceFloat
  end interface

  interface amaxMultiVecDevice
    function amaxMultiVecDeviceFloat(res,n,deviceVecA) &
         & result(val) bind(c,name='amaxMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_float) :: res
      type(c_ptr), value    :: deviceVecA
    end function amaxMultiVecDeviceFloat
  end interface

  interface asumMultiVecDevice
    function asumMultiVecDeviceFloat(res,n,deviceVecA) &
         & result(val) bind(c,name='asumMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_float) :: res
      type(c_ptr), value    :: deviceVecA
    end function asumMultiVecDeviceFloat
  end interface

  interface axpbyMultiVecDevice
    function axpbyMultiVecDeviceFloat(n,alpha,deviceVecA,beta,deviceVecB) &
         & result(res) bind(c,name='axpbyMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      real(c_float), value :: alpha, beta
      type(c_ptr), value  :: deviceVecA, deviceVecB
    end function axpbyMultiVecDeviceFloat
  end interface

  interface abgdxyzMultiVecDevice
    function abgdxyzMultiVecDeviceFloat(n,alpha,beta,gamma,delta,deviceVecX,&
         & deviceVecY,deviceVecZ) &
         & result(res) bind(c,name='abgdxyzMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      real(c_float), value :: alpha, beta,gamma,delta
      type(c_ptr), value  :: deviceVecX, deviceVecY, deviceVecZ
    end function abgdxyzMultiVecDeviceFloat
  end interface

  interface axyMultiVecDevice
    function axyMultiVecDeviceFloat(n,alpha,deviceVecA,deviceVecB) &
         & result(res) bind(c,name='axyMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      real(c_float), value :: alpha
      type(c_ptr), value       :: deviceVecA, deviceVecB
    end function axyMultiVecDeviceFloat
  end interface

  interface axybzMultiVecDevice
    function axybzMultiVecDeviceFloat(n,alpha,deviceVecA,deviceVecB,beta,deviceVecZ) &
         & result(res) bind(c,name='axybzMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)              :: res
      integer(c_int), value       :: n
      real(c_float), value     :: alpha, beta
      type(c_ptr), value          :: deviceVecA, deviceVecB,deviceVecZ
    end function axybzMultiVecDeviceFloat
  end interface


  interface absMultiVecDevice
    function absMultiVecDeviceFloat(n,alpha,deviceVecA) &
         & result(res) bind(c,name='absMultiVecDeviceFloat')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      real(c_float), value :: alpha
      type(c_ptr), value    :: deviceVecA
    end function absMultiVecDeviceFloat
    function absMultiVecDeviceFloat2(n,alpha,deviceVecA,deviceVecB) &
         & result(res) bind(c,name='absMultiVecDeviceFloat2')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      real(c_float), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
    end function absMultiVecDeviceFloat2
  end interface

  interface inner_register
    module procedure inner_registerFloat
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterFloat
  end interface

contains


  function inner_registerFloat(buffer,dval) result(res)
    real(c_float), allocatable, target :: buffer(:)
    type(c_ptr)            :: dval
    integer(c_int)         :: res
    real(c_float)         :: dummy
    res = registerMapped(c_loc(buffer),dval,size(buffer), dummy)        
  end function inner_registerFloat

  subroutine inner_unregisterFloat(buffer)
    real(c_float), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterFloat

end module psb_s_vectordev_mod
