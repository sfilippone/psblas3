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
  

module psb_base_vectordev_mod
  use iso_c_binding 
  use core_mod

  type, bind(c) :: multivec_dev_parms
    integer(c_int) :: count
    integer(c_int) :: element_type
    integer(c_int) :: pitch
    integer(c_int) :: size
  end type multivec_dev_parms
 
#ifdef HAVE_SPGPU  


  interface 
    function FallocMultiVecDevice(deviceVec,count,Size,elementType) &
         & result(res) bind(c,name='FallocMultiVecDevice')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: count,Size,elementType
      type(c_ptr)           :: deviceVec
    end function FallocMultiVecDevice
  end interface


  interface 
    subroutine  unregisterMapped(buf) &
         & bind(c,name='unregisterMapped')
      use iso_c_binding
      type(c_ptr), value :: buf
    end subroutine unregisterMapped
  end interface

  interface 
    subroutine  freeMultiVecDevice(deviceVec) &
         & bind(c,name='freeMultiVecDevice')
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
    end subroutine freeMultiVecDevice
  end interface

  interface 
    function  getMultiVecDeviceSize(deviceVec) &
         & bind(c,name='getMultiVecDeviceSize') result(res)
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      integer(c_int)      :: res
    end function getMultiVecDeviceSize
  end interface

  interface 
    function  getMultiVecDeviceCount(deviceVec) &
         & bind(c,name='getMultiVecDeviceCount') result(res)
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      integer(c_int)      :: res
    end function getMultiVecDeviceCount
  end interface

  interface 
    function  getMultiVecDevicePitch(deviceVec) &
         & bind(c,name='getMultiVecDevicePitch') result(res)
      use iso_c_binding
      type(c_ptr), value  :: deviceVec
      integer(c_int)      :: res
    end function getMultiVecDevicePitch
  end interface

#endif  


end module psb_base_vectordev_mod
