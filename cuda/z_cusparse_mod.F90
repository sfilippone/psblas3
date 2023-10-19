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
  

module z_cusparse_mod
  use base_cusparse_mod
  
  type, bind(c) :: z_Cmat
    type(c_ptr) :: Mat = c_null_ptr
  end type z_Cmat
  
#if CUDA_SHORT_VERSION <= 10 
  type, bind(c) :: z_Hmat
    type(c_ptr) :: Mat = c_null_ptr
  end type z_Hmat
#endif
  
  
#if defined(HAVE_CUDA) && defined(HAVE_SPGPU)
  
  interface CSRGDeviceFree
    function z_CSRGDeviceFree(Mat) &
         & bind(c,name="z_CSRGDeviceFree") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)   :: Mat
      integer(c_int) :: res
    end function z_CSRGDeviceFree
  end interface
  
  interface CSRGDeviceSetMatType
    function z_CSRGDeviceSetMatType(Mat,type) &
         & bind(c,name="z_CSRGDeviceSetMatType") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function z_CSRGDeviceSetMatType
  end interface
  
  interface CSRGDeviceSetMatFillMode
    function z_CSRGDeviceSetMatFillMode(Mat,type) &
         & bind(c,name="z_CSRGDeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function z_CSRGDeviceSetMatFillMode
  end interface
  
  interface CSRGDeviceSetMatDiagType
    function z_CSRGDeviceSetMatDiagType(Mat,type) &
         & bind(c,name="z_CSRGDeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function z_CSRGDeviceSetMatDiagType
  end interface
  
  interface CSRGDeviceSetMatIndexBase
    function z_CSRGDeviceSetMatIndexBase(Mat,type) &
         & bind(c,name="z_CSRGDeviceSetMatIndexBase") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function z_CSRGDeviceSetMatIndexBase
  end interface
  
  interface CSRGDeviceCsrsmAnalysis
    function z_CSRGDeviceCsrsmAnalysis(Mat) &
         & bind(c,name="z_CSRGDeviceCsrsmAnalysis") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int)        :: res
    end function z_CSRGDeviceCsrsmAnalysis
  end interface
  
  interface CSRGDeviceAlloc
    function z_CSRGDeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="z_CSRGDeviceAlloc") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function z_CSRGDeviceAlloc
  end interface

  interface CSRGDeviceGetParms
    function z_CSRGDeviceGetParms(Mat,nr,nc,nz) &
         & bind(c,name="z_CSRGDeviceGetParms") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int)        :: nr, nc, nz
      integer(c_int)        :: res
    end function z_CSRGDeviceGetParms
  end interface 
  
  interface spsvCSRGDevice
    function z_spsvCSRGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="z_spsvCSRGDevice") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_double_complex), value :: alpha,beta
      integer(c_int)        :: res
    end function z_spsvCSRGDevice
  end interface
  
  interface spmvCSRGDevice
    function z_spmvCSRGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="z_spmvCSRGDevice") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_double_complex), value :: alpha,beta
      integer(c_int)        :: res
    end function z_spmvCSRGDevice
  end interface
  
  interface CSRGHost2Device
    function z_CSRGHost2Device(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="z_CSRGHost2Device") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_double_complex)         :: val(*)
      integer(c_int)        :: res
    end function z_CSRGHost2Device
  end interface
  
  interface CSRGDevice2Host
    function z_CSRGDevice2Host(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="z_CSRGDevice2Host") result(res)
      use iso_c_binding
      import  z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_double_complex)         :: val(*)
      integer(c_int)        :: res
    end function z_CSRGDevice2Host
  end interface
  
#if CUDA_SHORT_VERSION <= 10
  interface HYBGDeviceAlloc
    function z_HYBGDeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="z_HYBGDeviceAlloc") result(res)
      use iso_c_binding
      import  z_hmat
      type(z_Hmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function z_HYBGDeviceAlloc
  end interface
  
  interface HYBGDeviceFree
    function z_HYBGDeviceFree(Mat) &
         & bind(c,name="z_HYBGDeviceFree") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)   :: Mat
      integer(c_int) :: res
    end function z_HYBGDeviceFree
  end interface
  
  interface HYBGDeviceSetMatType 
    function z_HYBGDeviceSetMatType(Mat,type) &
         & bind(c,name="z_HYBGDeviceSetMatType") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function z_HYBGDeviceSetMatType
  end interface
  
  interface HYBGDeviceSetMatFillMode
    function z_HYBGDeviceSetMatFillMode(Mat,type) &
         & bind(c,name="z_HYBGDeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function z_HYBGDeviceSetMatFillMode
  end interface
  
  interface HYBGDeviceSetMatDiagType
    function z_HYBGDeviceSetMatDiagType(Mat,type) &
         & bind(c,name="z_HYBGDeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function z_HYBGDeviceSetMatDiagType
  end interface
  
  interface HYBGDeviceSetMatIndexBase
    function z_HYBGDeviceSetMatIndexBase(Mat,type) &
         & bind(c,name="z_HYBGDeviceSetMatIndexBase") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function z_HYBGDeviceSetMatIndexBase
  end interface
  
  interface HYBGDeviceHybsmAnalysis
    function z_HYBGDeviceHybsmAnalysis(Mat) &
         & bind(c,name="z_HYBGDeviceHybsmAnalysis") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int)        :: res
    end function z_HYBGDeviceHybsmAnalysis
  end interface
  
  interface spsvHYBGDevice
    function z_spsvHYBGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="z_spsvHYBGDevice") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_double_complex), value  :: alpha,beta
      integer(c_int)        :: res
    end function z_spsvHYBGDevice
  end interface
  
  interface spmvHYBGDevice
    function z_spmvHYBGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="z_spmvHYBGDevice") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_double_complex), value  :: alpha,beta
      integer(c_int)        :: res
    end function z_spmvHYBGDevice
  end interface
  
  interface HYBGHost2Device
    function z_HYBGHost2Device(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="z_HYBGHost2Device") result(res)
      use iso_c_binding
      import  z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_double_complex)         :: val(*)
      integer(c_int)        :: res
    end function z_HYBGHost2Device
  end interface
#endif
  
#endif
  
end module z_cusparse_mod
