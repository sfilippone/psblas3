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
  

module c_cusparse_mod
  use base_cusparse_mod
  
  type, bind(c) :: c_Cmat
    type(c_ptr) :: Mat = c_null_ptr
  end type c_Cmat
  
#if CUDA_SHORT_VERSION <= 10 
  type, bind(c) :: c_Hmat
    type(c_ptr) :: Mat = c_null_ptr
  end type c_Hmat
#endif
  
  interface CSRGDeviceFree
    function c_CSRGDeviceFree(Mat) &
         & bind(c,name="c_CSRGDeviceFree") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)   :: Mat
      integer(c_int) :: res
    end function c_CSRGDeviceFree
  end interface
  
  interface CSRGDeviceSetMatType
    function c_CSRGDeviceSetMatType(Mat,type) &
         & bind(c,name="c_CSRGDeviceSetMatType") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function c_CSRGDeviceSetMatType
  end interface
  
  interface CSRGDeviceSetMatFillMode
    function c_CSRGDeviceSetMatFillMode(Mat,type) &
         & bind(c,name="c_CSRGDeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function c_CSRGDeviceSetMatFillMode
  end interface
  
  interface CSRGDeviceSetMatDiagType
    function c_CSRGDeviceSetMatDiagType(Mat,type) &
         & bind(c,name="c_CSRGDeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function c_CSRGDeviceSetMatDiagType
  end interface
  
  interface CSRGDeviceSetMatIndexBase
    function c_CSRGDeviceSetMatIndexBase(Mat,type) &
         & bind(c,name="c_CSRGDeviceSetMatIndexBase") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function c_CSRGDeviceSetMatIndexBase
  end interface
  
#if CUDA_SHORT_VERSION <= 10  
  interface CSRGDeviceCsrsmAnalysis
    function c_CSRGDeviceCsrsmAnalysis(Mat) &
         & bind(c,name="c_CSRGDeviceCsrsmAnalysis") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int)        :: res
    end function c_CSRGDeviceCsrsmAnalysis
  end interface
#else
  interface CSRGIsNullSvBuffer
    function c_CSRGIsNullSvBuffer(Mat) &
         & bind(c,name="c_CSRGIsNullSvBuffer") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int)        :: res
    end function c_CSRGIsNullSvBuffer
  end interface
#endif

  interface CSRGDeviceAlloc
    function c_CSRGDeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="c_CSRGDeviceAlloc") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function c_CSRGDeviceAlloc
  end interface

  interface CSRGDeviceGetParms
    function c_CSRGDeviceGetParms(Mat,nr,nc,nz) &
         & bind(c,name="c_CSRGDeviceGetParms") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int)        :: nr, nc, nz
      integer(c_int)        :: res
    end function c_CSRGDeviceGetParms
  end interface  
  
  interface spsvCSRGDevice
    function c_spsvCSRGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="c_spsvCSRGDevice") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_float_complex), value :: alpha,beta
      integer(c_int)        :: res
    end function c_spsvCSRGDevice
  end interface
  
  interface spmvCSRGDevice
    function c_spmvCSRGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="c_spmvCSRGDevice") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_float_complex), value :: alpha,beta
      integer(c_int)        :: res
    end function c_spmvCSRGDevice
  end interface
  
  interface CSRGHost2Device
    function c_CSRGHost2Device(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="c_CSRGHost2Device") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_float_complex)         :: val(*)
      integer(c_int)        :: res
    end function c_CSRGHost2Device
  end interface
  
  interface CSRGDevice2Host
    function c_CSRGDevice2Host(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="c_CSRGDevice2Host") result(res)
      use iso_c_binding
      import  c_Cmat
      type(c_Cmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_float_complex)         :: val(*)
      integer(c_int)        :: res
    end function c_CSRGDevice2Host
  end interface

#if CUDA_SHORT_VERSION <=10  
  interface HYBGDeviceAlloc
    function c_HYBGDeviceAlloc(Mat,nr,nc,nz) &
         & bind(c,name="c_HYBGDeviceAlloc") result(res)
      use iso_c_binding
      import  c_hmat
      type(c_Hmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)        :: res
    end function c_HYBGDeviceAlloc
  end interface
  
  interface HYBGDeviceFree
    function c_HYBGDeviceFree(Mat) &
         & bind(c,name="c_HYBGDeviceFree") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)   :: Mat
      integer(c_int) :: res
    end function c_HYBGDeviceFree
  end interface
  
  interface HYBGDeviceSetMatType 
    function c_HYBGDeviceSetMatType(Mat,type) &
         & bind(c,name="c_HYBGDeviceSetMatType") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function c_HYBGDeviceSetMatType
  end interface
  
  interface HYBGDeviceSetMatFillMode
    function c_HYBGDeviceSetMatFillMode(Mat,type) &
         & bind(c,name="c_HYBGDeviceSetMatFillMode") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function c_HYBGDeviceSetMatFillMode
  end interface
  
  interface HYBGDeviceSetMatDiagType
    function c_HYBGDeviceSetMatDiagType(Mat,type) &
         & bind(c,name="c_HYBGDeviceSetMatDiagType") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function c_HYBGDeviceSetMatDiagType
  end interface
  
  interface HYBGDeviceSetMatIndexBase
    function c_HYBGDeviceSetMatIndexBase(Mat,type) &
         & bind(c,name="c_HYBGDeviceSetMatIndexBase") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)          :: Mat
      integer(c_int),value  :: type
      integer(c_int)        :: res
    end function c_HYBGDeviceSetMatIndexBase
  end interface
  
  interface HYBGDeviceHybsmAnalysis
    function c_HYBGDeviceHybsmAnalysis(Mat) &
         & bind(c,name="c_HYBGDeviceHybsmAnalysis") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)          :: Mat
      integer(c_int)        :: res
    end function c_HYBGDeviceHybsmAnalysis
  end interface
  
  interface spsvHYBGDevice
    function c_spsvHYBGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="c_spsvHYBGDevice") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_float_complex), value  :: alpha,beta
      integer(c_int)        :: res
    end function c_spsvHYBGDevice
  end interface
  
  interface spmvHYBGDevice
    function c_spmvHYBGDevice(Mat,alpha,x,beta,y) &
         &  bind(c,name="c_spmvHYBGDevice") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_float_complex), value  :: alpha,beta
      integer(c_int)        :: res
    end function c_spmvHYBGDevice
  end interface
  
  interface HYBGHost2Device
    function c_HYBGHost2Device(Mat,m,n,nz,irp,ja,val) &
         &  bind(c,name="c_HYBGHost2Device") result(res)
      use iso_c_binding
      import  c_Hmat
      type(c_Hmat)          :: Mat
      integer(c_int), value :: m,n,nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_float_complex)         :: val(*)
      integer(c_int)        :: res
    end function c_HYBGHost2Device
  end interface
#endif
  
end module c_cusparse_mod
