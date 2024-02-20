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
module psb_z_vectordev_mod

  use psb_base_vectordev_mod
 
  interface registerMapped
    function registerMappedDoubleComplex(buf,d_p,n,dummy) &
         & result(res) bind(c,name='registerMappedDoubleComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: buf
      type(c_ptr) :: d_p
      integer(c_int),value :: n
      complex(c_double_complex), value :: dummy
    end function registerMappedDoubleComplex
  end interface

  interface writeMultiVecDevice 
    function writeMultiVecDeviceDoubleComplex(deviceVec,hostVec) &
         & result(res) bind(c,name='writeMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)             :: res
      type(c_ptr), value         :: deviceVec
      complex(c_double_complex)   :: hostVec(*)
    end function writeMultiVecDeviceDoubleComplex
    function writeMultiVecDeviceDoubleComplexR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='writeMultiVecDeviceDoubleComplexR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      complex(c_double_complex)      :: hostVec(ld,*)
    end function writeMultiVecDeviceDoubleComplexR2
  end interface 

  interface readMultiVecDevice
    function readMultiVecDeviceDoubleComplex(deviceVec,hostVec) &
         & result(res) bind(c,name='readMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)           :: res
      type(c_ptr), value       :: deviceVec
      complex(c_double_complex) :: hostVec(*)
    end function readMultiVecDeviceDoubleComplex
    function readMultiVecDeviceDoubleComplexR2(deviceVec,hostVec,ld) &
         & result(res) bind(c,name='readMultiVecDeviceDoubleComplexR2')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int), value :: ld
      complex(c_double_complex)      :: hostVec(ld,*)
    end function readMultiVecDeviceDoubleComplexR2
  end interface 

  interface allocateDoubleComplex
    function allocateDoubleComplex(didx,n) &
         & result(res) bind(c,name='allocateDoubleComplex') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: n
      integer(c_int)  :: res
    end function allocateDoubleComplex
    function allocateMultiDoubleComplex(didx,m,n) &
         & result(res) bind(c,name='allocateMultiDoubleComplex') 
      use iso_c_binding
      type(c_ptr) :: didx
      integer(c_int),value :: m,n
      integer(c_int)  :: res
    end function allocateMultiDoubleComplex
  end interface

  interface writeDoubleComplex
    function writeDoubleComplex(didx,hidx,n) &
         & result(res) bind(c,name='writeDoubleComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      complex(c_double_complex)       :: hidx(*)
      integer(c_int),value :: n
    end function writeDoubleComplex
    function writeDoubleComplexFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='writeDoubleComplexFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      complex(c_double_complex)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function writeDoubleComplexFirst
    function writeMultiDoubleComplex(didx,hidx,m,n) &
         & result(res) bind(c,name='writeMultiDoubleComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      complex(c_double_complex)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function writeMultiDoubleComplex
  end interface
  
  interface readDoubleComplex
    function readDoubleComplex(didx,hidx,n) &
         & result(res) bind(c,name='readDoubleComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      complex(c_double_complex)       :: hidx(*)
      integer(c_int),value :: n
    end function readDoubleComplex
    function readDoubleComplexFirst(first,didx,hidx,n,IndexBase) &
         & result(res) bind(c,name='readDoubleComplexFirst')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value   :: didx
      complex(c_double_complex)     :: hidx(*)
      integer(c_int),value :: n, first, IndexBase
    end function readDoubleComplexFirst
    function readMultiDoubleComplex(didx,hidx,m,n) &
         & result(res) bind(c,name='readMultiDoubleComplex')
      use iso_c_binding
      integer(c_int) :: res
      type(c_ptr), value :: didx
      complex(c_double_complex)       :: hidx(m,*)
      integer(c_int),value :: m,n
    end function readMultiDoubleComplex
  end interface
  
  interface
    subroutine  freeDoubleComplex(didx) &
         & bind(c,name='freeDoubleComplex')
      use iso_c_binding
      type(c_ptr), value :: didx
    end subroutine freeDoubleComplex
  end interface
  

  interface setScalDevice
    function setScalMultiVecDeviceDoubleComplex(val, first, last, &
         & indexBase, deviceVecX) result(res) &
         & bind(c,name='setscalMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: first,last,indexbase
      complex(c_double_complex), value :: val
      type(c_ptr),   value  :: deviceVecX
    end function setScalMultiVecDeviceDoubleComplex
  end interface

  interface 
    function geinsMultiVecDeviceDoubleComplex(n,deviceVecIrl,deviceVecVal,&
         & dupl,indexbase,deviceVecX) &
         & result(res) bind(c,name='geinsMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n, dupl,indexbase
      type(c_ptr), value  :: deviceVecIrl, deviceVecVal, deviceVecX
    end function geinsMultiVecDeviceDoubleComplex
  end interface

  ! New gather functions

  interface 
    function igathMultiVecDeviceDoubleComplex(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceDoubleComplex
  end interface

  interface 
    function igathMultiVecDeviceDoubleComplexVecIdx(deviceVec, vectorId, n, first, idx, &
         & hfirst, hostVec, indexBase) &
         & result(res) bind(c,name='igathMultiVecDeviceDoubleComplexVecIdx')
      use iso_c_binding
      integer(c_int)      :: res
      type(c_ptr), value  :: deviceVec
      integer(c_int),value:: vectorId
      integer(c_int),value:: first, n, hfirst
      type(c_ptr),value	  :: idx
      type(c_ptr),value   :: hostVec
      integer(c_int),value:: indexBase
    end function igathMultiVecDeviceDoubleComplexVecIdx
  end interface

  interface 
    function iscatMultiVecDeviceDoubleComplex(deviceVec, vectorId, & 
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      complex(c_double_complex),value :: beta
    end function iscatMultiVecDeviceDoubleComplex
  end interface

  interface 
    function iscatMultiVecDeviceDoubleComplexVecIdx(deviceVec, vectorId, &
         & first, n, idx, hfirst, hostVec, indexBase, beta) &
         & result(res) bind(c,name='iscatMultiVecDeviceDoubleComplexVecIdx')
      use iso_c_binding
      integer(c_int)         :: res
      type(c_ptr), value     :: deviceVec
      integer(c_int),value   :: vectorId
      integer(c_int),value   :: first, n, hfirst
      type(c_ptr), value     :: idx
      type(c_ptr), value     :: hostVec
      integer(c_int),value   :: indexBase
      complex(c_double_complex),value :: beta
    end function iscatMultiVecDeviceDoubleComplexVecIdx
  end interface


  interface scalMultiVecDevice
    function scalMultiVecDeviceDoubleComplex(alpha,deviceVecA) &
         & result(val) bind(c,name='scalMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: res
      complex(c_double_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA
    end function scalMultiVecDeviceDoubleComplex
  end interface

  interface dotMultiVecDevice
    function dotMultiVecDeviceDoubleComplex(res, n,deviceVecA,deviceVecB) &
         & result(val) bind(c,name='dotMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      complex(c_double_complex) :: res
      type(c_ptr), value    :: deviceVecA, deviceVecB
    end function dotMultiVecDeviceDoubleComplex
  end interface
    
  interface nrm2MultiVecDeviceComplex
    function nrm2MultiVecDeviceDoubleComplex(res,n,deviceVecA) &
         & result(val) bind(c,name='nrm2MultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_double)         :: res
      type(c_ptr), value    :: deviceVecA
    end function nrm2MultiVecDeviceDoubleComplex
  end interface

  interface amaxMultiVecDeviceComplex
    function amaxMultiVecDeviceDoubleComplex(res,n,deviceVecA) &
         & result(val) bind(c,name='amaxMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_double) :: res
      type(c_ptr), value    :: deviceVecA
    end function amaxMultiVecDeviceDoubleComplex
  end interface

  interface asumMultiVecDeviceComplex
    function asumMultiVecDeviceDoubleComplex(res,n,deviceVecA) &
         & result(val) bind(c,name='asumMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: val
      integer(c_int), value :: n
      real(c_double) :: res
      type(c_ptr), value    :: deviceVecA
    end function asumMultiVecDeviceDoubleComplex
  end interface

  interface axpbyMultiVecDevice
    function axpbyMultiVecDeviceDoubleComplex(n,alpha,deviceVecA,beta,deviceVecB) &
         & result(res) bind(c,name='axpbyMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value  :: deviceVecA, deviceVecB
    end function axpbyMultiVecDeviceDoubleComplex
  end interface

  interface abgdxyzMultiVecDevice
    function abgdxyzMultiVecDeviceDoubleComplex(n,alpha,beta,gamma,delta,deviceVecX,&
         & deviceVecY,deviceVecZ) &
         & result(res) bind(c,name='abgdxyzMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)      :: res
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha, beta,gamma,delta
      type(c_ptr), value  :: deviceVecX, deviceVecY, deviceVecZ
    end function abgdxyzMultiVecDeviceDoubleComplex
  end interface

  interface axyMultiVecDevice
    function axyMultiVecDeviceDoubleComplex(n,alpha,deviceVecA,deviceVecB) &
         & result(res) bind(c,name='axyMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha
      type(c_ptr), value       :: deviceVecA, deviceVecB
    end function axyMultiVecDeviceDoubleComplex
  end interface

  interface axybzMultiVecDevice
    function axybzMultiVecDeviceDoubleComplex(n,alpha,deviceVecA,deviceVecB,beta,deviceVecZ) &
         & result(res) bind(c,name='axybzMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)              :: res
      integer(c_int), value       :: n
      complex(c_double_complex), value     :: alpha, beta
      type(c_ptr), value          :: deviceVecA, deviceVecB,deviceVecZ
    end function axybzMultiVecDeviceDoubleComplex
  end interface


  interface absMultiVecDevice
    function absMultiVecDeviceDoubleComplex(n,alpha,deviceVecA) &
         & result(res) bind(c,name='absMultiVecDeviceDoubleComplex')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA
    end function absMultiVecDeviceDoubleComplex
    function absMultiVecDeviceDoubleComplex2(n,alpha,deviceVecA,deviceVecB) &
         & result(res) bind(c,name='absMultiVecDeviceDoubleComplex2')
      use iso_c_binding
      integer(c_int)        :: res
      integer(c_int), value :: n
      complex(c_double_complex), value :: alpha
      type(c_ptr), value    :: deviceVecA, deviceVecB
    end function absMultiVecDeviceDoubleComplex2
  end interface

  interface inner_register
    module procedure inner_registerDoubleComplex
  end interface
  
  interface inner_unregister
    module procedure inner_unregisterDoubleComplex
  end interface

contains


  function inner_registerDoubleComplex(buffer,dval) result(res)
    complex(c_double_complex), allocatable, target :: buffer(:)
    type(c_ptr)            :: dval
    integer(c_int)         :: res
    complex(c_double_complex)         :: dummy
    res = registerMapped(c_loc(buffer),dval,size(buffer), dummy)        
  end function inner_registerDoubleComplex

  subroutine inner_unregisterDoubleComplex(buffer)
    complex(c_double_complex), allocatable, target :: buffer(:)

    call  unregisterMapped(c_loc(buffer))
  end subroutine inner_unregisterDoubleComplex

end module psb_z_vectordev_mod
