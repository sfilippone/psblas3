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
  

module psb_i_vectordev_mod

  use psb_base_vectordev_mod
 
#ifdef HAVE_SPGPU  
  
  interface registerMapped
    function registerMappedInt(buf,d_p,n,dummy) &
         & result(res) bind(c,name='registerMappedInt')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: buf
      type(c_ptr) :: d_p
      integer(c_int),value :: n
      integer(c_int), value :: dummy
    end function registerMappedInt
  end interface

  interface writeMultiVecDevice 
    function writeMultiVecDeviceInt(deviceVec,hostVec) &
         & result(res) bind(c,name='writeMultiVecDeviceInt')
      use iso_c_binding
      integer(c_int)             :: res
      type(c_ptr), value         :: deviceVec
      integer(c_int)   :: hostVec(*)
    end function writeMultiVecDeviceInt
    function writeMultiVecDeviceIntR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='writeMultiVecDeviceIntR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      integer(c_int)      :: hostVec(ld,*)
    end function writeMultiVecDeviceIntR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceInt(deviceVec,hostVec) &
         & result(res) bind(c,name='readMultiVecDeviceInt')
      use iso_c_binding
      integer(c_int)           :: res
      type(c_ptr), value       :: deviceVec
      integer(c_int) :: hostVec(*)
    end function readMultiVecDeviceInt
    function readMultiVecDeviceIntR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='readMultiVecDeviceIntR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      integer(c_int)      :: hostVec(ld,*)
    end function readMultiVecDeviceIntR2
  end interface 

  interface allocateInt
    function allocateInt(didx,n) &
         & result(res) bind(c,name='allocateInt') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: n
      integer(c_int)  :: res
    end function allocateInt
    function allocateMultiInt(didx,m,n) &
         & result(res) bind(c,name='allocateMultiInt') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: m,n
      integer(c_int)  :: res
    end function allocateMultiInt
  end interface

  interface writeInt
    function writeInt(didx,hidx,n) &
         & result(res) bind(c,name='writeInt')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      integer(c_int)       :: hidx(*)
      integer(c_int),value :: n
    end function writeInt
    function writeIntFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='writeIntFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      integer(c_int)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function writeIntFirst
    function writeMultiInt(didx,hidx,m,n) &
         & result(res) bind(c,name='writeMultiInt')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      integer(c_int)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function writeMultiInt
  end interface
  
  interface readInt
    function readInt(didx,hidx,n) &
         & result(res) bind(c,name='readInt')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      integer(c_int)       :: hidx(*)
      integer(c_int),value :: n
    end function readInt
    function readIntFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='readIntFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      integer(c_int)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function readIntFirst
    function readMultiInt(didx,hidx,m,n) &
         & result(res) bind(c,name='readMultiInt')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      integer(c_int)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function readMultiInt
  end interface
  
  interface
    subroutine  freeInt(didx) &
         & bind(c,name='freeInt')
      use iso_c_binding
      type(c_ptr), value :: didx
    end subroutine freeInt
  end interface
  

  interface setScalDevice
    function setScalMultiVecDeviceInt(val, first, last, &
         & indexBase, deviceVecX) result(res) &
         & bind(c,name='setscalMultiVecDeviceInt')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: first,last,indexbase
      integer(c_int), value :: val
      type(c_ptr),   value  :: deviceVecX
    end function setScalMultiVecDeviceInt
  end interface

  interface 
    function geinsMultiVecDeviceInt(n,deviceVecIrl,deviceVecVal,&
         & dupl,indexbase,deviceVecX) &
         & result(res) bind(c,name='geinsMultiVecDeviceInt')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n, dupl,indexbase
      type(c_ptr), value  :: deviceVecIrl, deviceVecVal, deviceVecX
    end function geinsMultiVecDeviceInt
  end interface

  ! New gather functions

  interface 
    function igathMultiVecDeviceInt(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceInt')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceInt
  end interface

  interface 
    function igathMultiVecDeviceIntVecIdx(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceIntVecIdx')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceIntVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceInt(deviceVec, vectorId, & 
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceInt')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      integer(c_int),value :: beta
    end function iscatMultiVecDeviceInt
  end interface

  interface 
    function iscatMultiVecDeviceIntVecIdx(deviceVec, vectorId, &
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceIntVecIdx')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      integer(c_int),value :: beta
    end function iscatMultiVecDeviceIntVecIdx
  end interface



  interface inner_register
    module procedure inner_registerInt
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterInt
  end interface

contains


  function inner_registerInt(buffer,dval) result(res)
    integer(c_int), allocatable, target :: buffer(:)
    type(c_ptr)            :: dval
    integer(c_int)         :: res
    integer(c_int)         :: dummy
    res = registerMapped(c_loc(buffer),dval,size(buffer), dummy)        
  end function inner_registerInt

  subroutine inner_unregisterInt(buffer)
    integer(c_int), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterInt

#endif  

end module psb_i_vectordev_mod
