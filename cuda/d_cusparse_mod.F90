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
  

module d_cusparse_mod
  use base_cusparse_mod
  
  type, bind(c) :: d_Cmat
    type(c_ptr) :: Mat = c_null_ptr
  end type d_Cmat
  
#if CUDA_SHORT_VERSION <= 10 
  type, bind(c) :: d_Hmat
    type(c_ptr) :: Mat = c_null_ptr
  end type d_Hmat
#endif
  
  
#if defined(HAVE_CUDA) && defined(HAVE_SPGPU)
  
  interface CSRGDeviceFree
    function d_CSRGDeviceFree(Mat) &
         & bind(c,name="d_CSRGDeviceFree") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)   :: Mat
      integer(c_int) :: res
    end function d_CSRGDeviceFree
  end interface
  
  interface CSRGDeviceSetMatType
    function d_CSRGDeviceSetMatType(Mat,type) &
         & bind(c,name="d_CSRGDeviceSetMatType") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_CSRGDeviceSetMatType
  end interface
  
  interface CSRGDeviceSetMatFillMode
    function d_CSRGDeviceSetMatFillMode(Mat,type) &
         & bind(c,name="d_CSRGDeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_CSRGDeviceSetMatFillMode
  end interface
  
  interface CSRGDeviceSetMatDiagType
    function d_CSRGDeviceSetMatDiagType(Mat,type) &
         & bind(c,name="d_CSRGDeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_CSRGDeviceSetMatDiagType
  end interface
  
  interface CSRGDeviceSetMatIndexBase
    function d_CSRGDeviceSetMatIndexBase(Mat,type) &
         & bind(c,name="d_CSRGDeviceSetMatIndexBase") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_CSRGDeviceSetMatIndexBase
  end interface
  
  interface CSRGDeviceCsrsmAnalysis
    function d_CSRGDeviceCsrsmAnalysis(Mat) &
         & bind(c,name="d_CSRGDeviceCsrsmAnalysis") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int)        :: res
    end function d_CSRGDeviceCsrsmAnalysis
  end interface
  
  interface CSRGDeviceAlloc
    function d_CSRGDeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="d_CSRGDeviceAlloc") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function d_CSRGDeviceAlloc
  end interface

  interface CSRGDeviceGetParms
    function d_CSRGDeviceGetParms(Mat,nr,nc,nz) &
         & bind(c,name="d_CSRGDeviceGetParms") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int)        :: nr, nc, nz
      integer(c_int)        :: res
    end function d_CSRGDeviceGetParms
  end interface
  
  interface spsvCSRGDevice
    function d_spsvCSRGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="d_spsvCSRGDevice") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_double), value :: alpha,beta
      integer(c_int)        :: res
    end function d_spsvCSRGDevice
  end interface
  
  interface spmvCSRGDevice
    function d_spmvCSRGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="d_spmvCSRGDevice") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_double), value :: alpha,beta
      integer(c_int)        :: res
    end function d_spmvCSRGDevice
  end interface
  
  interface CSRGHost2Device
    function d_CSRGHost2Device(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="d_CSRGHost2Device") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_double)         :: val(*)
      integer(c_int)        :: res
    end function d_CSRGHost2Device
  end interface
  
  interface CSRGDevice2Host
    function d_CSRGDevice2Host(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="d_CSRGDevice2Host") result(res)
      use iso_c_binding
      import  d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_double)         :: val(*)
      integer(c_int)        :: res
    end function d_CSRGDevice2Host
  end interface
  
#if CUDA_SHORT_VERSION <= 10
  interface HYBGDeviceAlloc
    function d_HYBGDeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="d_HYBGDeviceAlloc") result(res)
      use iso_c_binding
      import  d_hmat
      type(d_Hmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function d_HYBGDeviceAlloc
  end interface
  
  interface HYBGDeviceFree
    function d_HYBGDeviceFree(Mat) &
         & bind(c,name="d_HYBGDeviceFree") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)   :: Mat
      integer(c_int) :: res
    end function d_HYBGDeviceFree
  end interface
  
  interface HYBGDeviceSetMatType 
    function d_HYBGDeviceSetMatType(Mat,type) &
         & bind(c,name="d_HYBGDeviceSetMatType") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_HYBGDeviceSetMatType
  end interface
  
  interface HYBGDeviceSetMatFillMode
    function d_HYBGDeviceSetMatFillMode(Mat,type) &
         & bind(c,name="d_HYBGDeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_HYBGDeviceSetMatFillMode
  end interface
  
  interface HYBGDeviceSetMatDiagType
    function d_HYBGDeviceSetMatDiagType(Mat,type) &
         & bind(c,name="d_HYBGDeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_HYBGDeviceSetMatDiagType
  end interface
  
  interface HYBGDeviceSetMatIndexBase
    function d_HYBGDeviceSetMatIndexBase(Mat,type) &
         & bind(c,name="d_HYBGDeviceSetMatIndexBase") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function d_HYBGDeviceSetMatIndexBase
  end interface
  
  interface HYBGDeviceHybsmAnalysis
    function d_HYBGDeviceHybsmAnalysis(Mat) &
         & bind(c,name="d_HYBGDeviceHybsmAnalysis") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int)        :: res
    end function d_HYBGDeviceHybsmAnalysis
  end interface
  
  interface spsvHYBGDevice
    function d_spsvHYBGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="d_spsvHYBGDevice") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_double), value  :: alpha,beta
      integer(c_int)        :: res
    end function d_spsvHYBGDevice
  end interface
  
  interface spmvHYBGDevice
    function d_spmvHYBGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="d_spmvHYBGDevice") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_double), value  :: alpha,beta
      integer(c_int)        :: res
    end function d_spmvHYBGDevice
  end interface
  
  interface HYBGHost2Device
    function d_HYBGHost2Device(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="d_HYBGHost2Device") result(res)
      use iso_c_binding
      import  d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_double)         :: val(*)
      integer(c_int)        :: res
    end function d_HYBGHost2Device
  end interface
#endif
  
#endif
  
end module d_cusparse_mod
