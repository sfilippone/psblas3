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
  

module psb_c_vectordev_mod

  use psb_base_vectordev_mod
 
  interface registerMapped
    function registerMappedFloatComplex(buf,d_p,n,dummy) &
         & result(res) bind(c,name='registerMappedFloatComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: buf
      type(c_ptr) :: d_p
      integer(c_int),value :: n
      complex(c_float_complex), value :: dummy
    end function registerMappedFloatComplex
  end interface

  interface writeMultiVecDevice 
    function writeMultiVecDeviceFloatComplex(deviceVec,hostVec) &
         & result(res) bind(c,name='writeMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)             :: res
      type(c_ptr), value         :: deviceVec
      complex(c_float_complex)   :: hostVec(*)
    end function writeMultiVecDeviceFloatComplex
    function writeMultiVecDeviceFloatComplexR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='writeMultiVecDeviceFloatComplexR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      complex(c_float_complex)      :: hostVec(ld,*)
    end function writeMultiVecDeviceFloatComplexR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceFloatComplex(deviceVec,hostVec) &
         & result(res) bind(c,name='readMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)           :: res
      type(c_ptr), value       :: deviceVec
      complex(c_float_complex) :: hostVec(*)
    end function readMultiVecDeviceFloatComplex
    function readMultiVecDeviceFloatComplexR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='readMultiVecDeviceFloatComplexR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      complex(c_float_complex)      :: hostVec(ld,*)
    end function readMultiVecDeviceFloatComplexR2
  end interface 

  interface allocateFloatComplex
    function allocateFloatComplex(didx,n) &
         & result(res) bind(c,name='allocateFloatComplex') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: n
      integer(c_int)  :: res
    end function allocateFloatComplex
    function allocateMultiFloatComplex(didx,m,n) &
         & result(res) bind(c,name='allocateMultiFloatComplex') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: m,n
      integer(c_int)  :: res
    end function allocateMultiFloatComplex
  end interface

  interface writeFloatComplex
    function writeFloatComplex(didx,hidx,n) &
         & result(res) bind(c,name='writeFloatComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      complex(c_float_complex)       :: hidx(*)
      integer(c_int),value :: n
    end function writeFloatComplex
    function writeFloatComplexFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='writeFloatComplexFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      complex(c_float_complex)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function writeFloatComplexFirst
    function writeMultiFloatComplex(didx,hidx,m,n) &
         & result(res) bind(c,name='writeMultiFloatComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      complex(c_float_complex)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function writeMultiFloatComplex
  end interface
  
  interface readFloatComplex
    function readFloatComplex(didx,hidx,n) &
         & result(res) bind(c,name='readFloatComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      complex(c_float_complex)       :: hidx(*)
      integer(c_int),value :: n
    end function readFloatComplex
    function readFloatComplexFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='readFloatComplexFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      complex(c_float_complex)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function readFloatComplexFirst
    function readMultiFloatComplex(didx,hidx,m,n) &
         & result(res) bind(c,name='readMultiFloatComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      complex(c_float_complex)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function readMultiFloatComplex
  end interface
  
  interface
    subroutine  freeFloatComplex(didx) &
         & bind(c,name='freeFloatComplex')
      use iso_c_binding
      type(c_ptr), value :: didx
    end subroutine freeFloatComplex
  end interface
  

  interface setScalDevice
    function setScalMultiVecDeviceFloatComplex(val, first, last, &
         & indexBase, deviceVecX) result(res) &
         & bind(c,name='setscalMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: first,last,indexbase
      complex(c_float_complex), value :: val
      type(c_ptr),   value  :: deviceVecX
    end function setScalMultiVecDeviceFloatComplex
  end interface

  interface 
    function geinsMultiVecDeviceFloatComplex(n,deviceVecIrl,deviceVecVal,&
         & dupl,indexbase,deviceVecX) &
         & result(res) bind(c,name='geinsMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n, dupl,indexbase
      type(c_ptr), value  :: deviceVecIrl, deviceVecVal, deviceVecX
    end function geinsMultiVecDeviceFloatComplex
  end interface

  ! New gather functions

  interface 
    function igathMultiVecDeviceFloatComplex(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceFloatComplex
  end interface

  interface 
    function igathMultiVecDeviceFloatComplexVecIdx(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceFloatComplexVecIdx')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceFloatComplexVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceFloatComplex(deviceVec, vectorId, & 
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      complex(c_float_complex),value :: beta
    end function iscatMultiVecDeviceFloatComplex
  end interface

  interface 
    function iscatMultiVecDeviceFloatComplexVecIdx(deviceVec, vectorId, &
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceFloatComplexVecIdx')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      complex(c_float_complex),value :: beta
    end function iscatMultiVecDeviceFloatComplexVecIdx
  end interface


  interface scalMultiVecDevice
    function scalMultiVecDeviceFloatComplex(alpha,deviceVecA) &
         & result(val) bind(c,name='scalMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)        :: res
      complex(c_float_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA
    end function scalMultiVecDeviceFloatComplex
  end interface

  interface dotMultiVecDevice
    function dotMultiVecDeviceFloatComplex(res, n,deviceVecA,deviceVecB) &
         & result(val) bind(c,name='dotMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      complex(c_float_complex) :: res
      type(c_ptr), value    :: deviceVecA, deviceVecB
    end function dotMultiVecDeviceFloatComplex
  end interface
    
  interface nrm2MultiVecDeviceComplex
    function nrm2MultiVecDeviceFloatComplex(res,n,deviceVecA) &
         & result(val) bind(c,name='nrm2MultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_float)         :: res
      type(c_ptr), value    :: deviceVecA
    end function nrm2MultiVecDeviceFloatComplex
  end interface

  interface amaxMultiVecDeviceComplex
    function amaxMultiVecDeviceFloatComplex(res,n,deviceVecA) &
         & result(val) bind(c,name='amaxMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_float) :: res
      type(c_ptr), value    :: deviceVecA
    end function amaxMultiVecDeviceFloatComplex
  end interface

  interface asumMultiVecDeviceComplex
    function asumMultiVecDeviceFloatComplex(res,n,deviceVecA) &
         & result(val) bind(c,name='asumMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_float) :: res
      type(c_ptr), value    :: deviceVecA
    end function asumMultiVecDeviceFloatComplex
  end interface


  interface axpbyMultiVecDevice
    function axpbyMultiVecDeviceFloatComplex(n,alpha,deviceVecA,beta,deviceVecB) &
         & result(res) bind(c,name='axpbyMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha, beta
      type(c_ptr), value  :: deviceVecA, deviceVecB
    end function axpbyMultiVecDeviceFloatComplex
  end interface

  interface axyMultiVecDevice
    function axyMultiVecDeviceFloatComplex(n,alpha,deviceVecA,deviceVecB) &
         & result(res) bind(c,name='axyMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha
      type(c_ptr), value       :: deviceVecA, deviceVecB
    end function axyMultiVecDeviceFloatComplex
  end interface

  interface axybzMultiVecDevice
    function axybzMultiVecDeviceFloatComplex(n,alpha,deviceVecA,deviceVecB,beta,deviceVecZ) &
         & result(res) bind(c,name='axybzMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)              :: res
      integer(c_int), value       :: n
      complex(c_float_complex), value     :: alpha, beta
      type(c_ptr), value          :: deviceVecA, deviceVecB,deviceVecZ
    end function axybzMultiVecDeviceFloatComplex
  end interface


  interface absMultiVecDevice
    function absMultiVecDeviceFloatComplex(n,alpha,deviceVecA) &
         & result(res) bind(c,name='absMultiVecDeviceFloatComplex')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA
    end function absMultiVecDeviceFloatComplex
    function absMultiVecDeviceFloatComplex2(n,alpha,deviceVecA,deviceVecB) &
         & result(res) bind(c,name='absMultiVecDeviceFloatComplex2')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      complex(c_float_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
    end function absMultiVecDeviceFloatComplex2
  end interface

  interface inner_register
    module procedure inner_registerFloatComplex
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterFloatComplex
  end interface

contains


  function inner_registerFloatComplex(buffer,dval) result(res)
    complex(c_float_complex), allocatable, target :: buffer(:)
    type(c_ptr)            :: dval
    integer(c_int)         :: res
    complex(c_float_complex)         :: dummy
    res = registerMapped(c_loc(buffer),dval,size(buffer), dummy)        
  end function inner_registerFloatComplex

  subroutine inner_unregisterFloatComplex(buffer)
    complex(c_float_complex), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterFloatComplex

end module psb_c_vectordev_mod
