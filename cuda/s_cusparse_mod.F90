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
  

module s_cusparse_mod
  use base_cusparse_mod
  
  type, bind(c) :: s_Cmat
    type(c_ptr) :: Mat = c_null_ptr
  end type s_Cmat
  
#if CUDA_SHORT_VERSION <= 10 
  type, bind(c) :: s_Hmat
    type(c_ptr) :: Mat = c_null_ptr
  end type s_Hmat
#endif
  
  interface CSRGDeviceFree
    function s_CSRGDeviceFree(Mat) &
         & bind(c,name="s_CSRGDeviceFree") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)   :: Mat
      integer(c_int) :: res
    end function s_CSRGDeviceFree
  end interface
  
  interface CSRGDeviceSetMatType
    function s_CSRGDeviceSetMatType(Mat,type) &
         & bind(c,name="s_CSRGDeviceSetMatType") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function s_CSRGDeviceSetMatType
  end interface
  
  interface CSRGDeviceSetMatFillMode
    function s_CSRGDeviceSetMatFillMode(Mat,type) &
         & bind(c,name="s_CSRGDeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function s_CSRGDeviceSetMatFillMode
  end interface
  
  interface CSRGDeviceSetMatDiagType
    function s_CSRGDeviceSetMatDiagType(Mat,type) &
         & bind(c,name="s_CSRGDeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function s_CSRGDeviceSetMatDiagType
  end interface
  
  interface CSRGDeviceSetMatIndexBase
    function s_CSRGDeviceSetMatIndexBase(Mat,type) &
         & bind(c,name="s_CSRGDeviceSetMatIndexBase") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function s_CSRGDeviceSetMatIndexBase
  end interface
  
  interface CSRGDeviceCsrsmAnalysis
    function s_CSRGDeviceCsrsmAnalysis(Mat) &
         & bind(c,name="s_CSRGDeviceCsrsmAnalysis") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int)        :: res
    end function s_CSRGDeviceCsrsmAnalysis
  end interface
  
  interface CSRGDeviceAlloc
    function s_CSRGDeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="s_CSRGDeviceAlloc") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function s_CSRGDeviceAlloc
  end interface

  interface CSRGDeviceGetParms
    function s_CSRGDeviceGetParms(Mat,nr,nc,nz) &
         & bind(c,name="s_CSRGDeviceGetParms") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int)        :: nr, nc, nz
      integer(c_int)        :: res
    end function s_CSRGDeviceGetParms
  end interface
  
  interface spsvCSRGDevice
    function s_spsvCSRGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="s_spsvCSRGDevice") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_float), value :: alpha,beta
      integer(c_int)        :: res
    end function s_spsvCSRGDevice
  end interface
  
  interface spmvCSRGDevice
    function s_spmvCSRGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="s_spmvCSRGDevice") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_float), value :: alpha,beta
      integer(c_int)        :: res
    end function s_spmvCSRGDevice
  end interface
  
  interface CSRGHost2Device
    function s_CSRGHost2Device(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="s_CSRGHost2Device") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_float)         :: val(*)
      integer(c_int)        :: res
    end function s_CSRGHost2Device
  end interface
  
  interface CSRGDevice2Host
    function s_CSRGDevice2Host(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="s_CSRGDevice2Host") result(res)
      use iso_c_binding
      import  s_Cmat
      type(s_Cmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_float)         :: val(*)
      integer(c_int)        :: res
    end function s_CSRGDevice2Host
  end interface
  
#if CUDA_SHORT_VERSION <= 10
  interface HYBGDeviceAlloc
    function s_HYBGDeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="s_HYBGDeviceAlloc") result(res)
      use iso_c_binding
      import  s_hmat
      type(s_Hmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function s_HYBGDeviceAlloc
  end interface
  
  interface HYBGDeviceFree
    function s_HYBGDeviceFree(Mat) &
         & bind(c,name="s_HYBGDeviceFree") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)   :: Mat
      integer(c_int) :: res
    end function s_HYBGDeviceFree
  end interface
  
  interface HYBGDeviceSetMatType 
    function s_HYBGDeviceSetMatType(Mat,type) &
         & bind(c,name="s_HYBGDeviceSetMatType") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function s_HYBGDeviceSetMatType
  end interface
  
  interface HYBGDeviceSetMatFillMode
    function s_HYBGDeviceSetMatFillMode(Mat,type) &
         & bind(c,name="s_HYBGDeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function s_HYBGDeviceSetMatFillMode
  end interface
  
  interface HYBGDeviceSetMatDiagType
    function s_HYBGDeviceSetMatDiagType(Mat,type) &
         & bind(c,name="s_HYBGDeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function s_HYBGDeviceSetMatDiagType
  end interface
  
  interface HYBGDeviceSetMatIndexBase
    function s_HYBGDeviceSetMatIndexBase(Mat,type) &
         & bind(c,name="s_HYBGDeviceSetMatIndexBase") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function s_HYBGDeviceSetMatIndexBase
  end interface
  
  interface HYBGDeviceHybsmAnalysis
    function s_HYBGDeviceHybsmAnalysis(Mat) &
         & bind(c,name="s_HYBGDeviceHybsmAnalysis") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)          :: Mat
      integer(c_int)        :: res
    end function s_HYBGDeviceHybsmAnalysis
  end interface
  
  interface spsvHYBGDevice
    function s_spsvHYBGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="s_spsvHYBGDevice") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_float), value  :: alpha,beta
      integer(c_int)        :: res
    end function s_spsvHYBGDevice
  end interface
  
  interface spmvHYBGDevice
    function s_spmvHYBGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="s_spmvHYBGDevice") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_float), value  :: alpha,beta
      integer(c_int)        :: res
    end function s_spmvHYBGDevice
  end interface
  
  interface HYBGHost2Device
    function s_HYBGHost2Device(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="s_HYBGHost2Device") result(res)
      use iso_c_binding
      import  s_Hmat
      type(s_Hmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_float)         :: val(*)
      integer(c_int)        :: res
    end function s_HYBGHost2Device
  end interface
#endif
  
end module s_cusparse_mod
