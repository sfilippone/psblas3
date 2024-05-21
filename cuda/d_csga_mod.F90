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
module d_csga_mod
  use d_cusparse_mod
  integer, parameter :: MAX_NNZ_PER_WG = 4096
  integer, parameter :: MAX_GRID_SIZE  = 65536
  
  type, bind(c) :: d_CAmat
    type(c_ptr) :: Mat = c_null_ptr
  end type d_CAmat

  interface CSGADeviceFree
    function d_CSGADeviceFree(Mat) &
         & bind(c,name="d_CSGADeviceFree") result(res)
      use iso_c_binding
      import  d_CAmat
      type(d_CAmat)  :: Mat
      integer(c_int) :: res
    end function d_CSGADeviceFree
  end interface

  interface CSGADeviceAlloc
    function d_CSGADeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="d_CSGADeviceAlloc") result(res)
      use iso_c_binding
      import  d_CAmat
      type(d_CAmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function d_CSGADeviceAlloc
  end interface
  
  interface CSGADeviceSetMatDiagType
    function d_CSGADeviceSetMatDiagType(Mat,type) &
         & bind(c,name="d_CSGADeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  d_CAmat
      type(d_CAmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_CSGADeviceSetMatDiagType
  end interface
    
  interface CSGADeviceSetMatFillMode
    function d_CSGADeviceSetMatFillMode(Mat,type) &
         & bind(c,name="d_CSGADeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  d_CAmat
      type(d_CAmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_CSGADeviceSetMatFillMode
  end interface
  
  interface CSGAHost2Device
    function d_CSGAHost2Device(Mat,m,n,nz,irp,ja,val,rowBlocks) &
         &  bind(c,name="d_CSGAHost2Device") result(res)
      use iso_c_binding
      import  d_CAmat
      type(d_CAmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*), rowBlocks(*)
      real(c_double)         :: val(*)
      integer(c_int)        :: res
    end function d_CSGAHost2Device
  end interface
  
  interface CSGADevice2Host
    function d_CSGADevice2Host(Mat,m,n,nz,irp,ja,val,rowBlocks) &
         &  bind(c,name="d_CSGADevice2Host") result(res)
      use iso_c_binding
      import  d_CAmat
      type(d_CAmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*), rowBlocks(*)
      real(c_double)         :: val(*)
      integer(c_int)        :: res
    end function d_CSGADevice2Host
  end interface
  

  

end module d_csga_mod
