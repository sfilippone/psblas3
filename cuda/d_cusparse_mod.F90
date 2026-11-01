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
!         software without specific prior written permission.
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
  
#if PSB_CUDA_SHORT_VERSION <= 10 
  type, bind(c) :: d_Hmat
    type(c_ptr) :: Mat = c_null_ptr
  end type d_Hmat
#endif
  
  interface CSRGDeviceFree
    function d_CSRGDeviceFree(Mat) result(res) &
        & bind(c, name = "d_CSRGDeviceFree")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)   :: Mat
      integer(c_int) :: res
    end function d_CSRGDeviceFree
  end interface
  
  interface CSRGDeviceSetMatType
    function d_CSRGDeviceSetMatType(Mat, type) result(res) &
        & bind(c, name = "d_CSRGDeviceSetMatType")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function d_CSRGDeviceSetMatType
  end interface
  
  interface CSRGDeviceSetMatFillMode
    function d_CSRGDeviceSetMatFillMode(Mat, type) result(res) &
        & bind(c, name = "d_CSRGDeviceSetMatFillMode")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function d_CSRGDeviceSetMatFillMode
  end interface
  
  interface CSRGDeviceSetMatDiagType
    function d_CSRGDeviceSetMatDiagType(Mat, type) result(res) &
        & bind(c, name = "d_CSRGDeviceSetMatDiagType")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function d_CSRGDeviceSetMatDiagType
  end interface
  
  interface CSRGDeviceSetMatIndexBase
    function d_CSRGDeviceSetMatIndexBase(Mat, type) result(res) &
        & bind(c, name = "d_CSRGDeviceSetMatIndexBase")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function d_CSRGDeviceSetMatIndexBase
  end interface
  
#if PSB_CUDA_SHORT_VERSION <= 10  
  interface CSRGDeviceCsrsmAnalysis
    function d_CSRGDeviceCsrsmAnalysis(Mat) result(res) &
        & bind(c, name = "d_CSRGDeviceCsrsmAnalysis")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)    :: Mat
      integer(c_int)  :: res
    end function d_CSRGDeviceCsrsmAnalysis
  end interface
#else
  interface CSRGIsNullSvBuffer
    function d_CSRGIsNullSvBuffer(Mat) result(res) &
        & bind(c, name = "d_CSRGIsNullSvBuffer")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)    :: Mat
      integer(c_int)  :: res
    end function d_CSRGIsNullSvBuffer
  end interface
#endif

  interface CSRGDeviceAlloc
    function d_CSRGDeviceAlloc(Mat, nr, nc, nz) result(res) &
        & bind(c, name = "d_CSRGDeviceAlloc")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)  :: res
    end function d_CSRGDeviceAlloc
  end interface

  interface CSRGDeviceGetParms
    function d_CSRGDeviceGetParms(Mat, nr, nc, nz) result(res) &
        & bind(c, name = "d_CSRGDeviceGetParms")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)    :: Mat
      integer(c_int)  :: nr, nc, nz
      integer(c_int)  :: res
    end function d_CSRGDeviceGetParms
  end interface
  
  interface spsvCSRGDevice
    function d_spsvCSRGDevice(Mat, alpha, x, beta, y) result(res) &
        & bind(c, name = "d_spsvCSRGDevice")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      real(c_double), value :: alpha, beta
      integer(c_int)  :: res
    end function d_spsvCSRGDevice
  end interface
  
  interface spmvCSRGDevice
    function d_spmvCSRGDeviceF(mat, alpha, x, beta, y) result(res) &
        & bind(c, name = "d_spmvCSRGDeviceF")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: mat
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int)  :: res
    end function d_spmvCSRGDeviceF

    function d_spmvCSRGDeviceMV(mat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = "d_spmvCSRGDeviceMV")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: mat
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int), value :: idx_y
      integer(c_int)  :: res
    end function d_spmvCSRGDeviceMV

    function d_spmvCSRGDeviceVM(mat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = "d_spmvCSRGDeviceVM")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: mat
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int), value :: idx_x
      integer(c_int)  :: res
    end function d_spmvCSRGDeviceVM
    
    function d_spmvCSRGDeviceMMi(mat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = "d_spmvCSRGDeviceMMi")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: mat
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int), value :: idx_x, idx_y
      integer(c_int)  :: res
    end function d_spmvCSRGDeviceMMi
  end interface
  
  interface CSRGHost2Device
    function d_CSRGHost2Device(Mat, m, n, nz, irp, ja, val) result(res) &
        & bind(c, name = "d_CSRGHost2Device")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: m, n, nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_double)        :: val(*)
      integer(c_int)  :: res
    end function d_CSRGHost2Device
  end interface
  
  interface CSRGDevice2Host
    function d_CSRGDevice2Host(Mat, m, n, nz, irp, ja, val) result(res) &
        & bind(c, name = "d_CSRGDevice2Host")
      use iso_c_binding
      import d_Cmat
      type(d_Cmat)          :: Mat
      integer(c_int), value :: m, n, nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_double)        :: val(*)
      integer(c_int)  :: res
    end function d_CSRGDevice2Host
  end interface
  
#if PSB_CUDA_SHORT_VERSION <= 10
  interface HYBGDeviceAlloc
    function d_HYBGDeviceAlloc(Mat, nr, nc, nz) result(res) &
        & bind(c, name = "d_HYBGDeviceAlloc")
      use iso_c_binding
      import d_hmat
      type(d_Hmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)  :: res
    end function d_HYBGDeviceAlloc
  end interface
  
  interface HYBGDeviceFree
    function d_HYBGDeviceFree(Mat) result(res) &
        & bind(c, name = "d_HYBGDeviceFree")
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)   :: Mat
      integer(c_int) :: res
    end function d_HYBGDeviceFree
  end interface
  
  interface HYBGDeviceSetMatType 
    function d_HYBGDeviceSetMatType(Mat, type) result(res) &
        & bind(c, name = "d_HYBGDeviceSetMatType")
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function d_HYBGDeviceSetMatType
  end interface
  
  interface HYBGDeviceSetMatFillMode
    function d_HYBGDeviceSetMatFillMode(Mat, type) result(res) &
        & bind(c, name = "d_HYBGDeviceSetMatFillMode")
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function d_HYBGDeviceSetMatFillMode
  end interface
  
  interface HYBGDeviceSetMatDiagType
    function d_HYBGDeviceSetMatDiagType(Mat, type) result(res) &
        & bind(c, name = "d_HYBGDeviceSetMatDiagType")
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function d_HYBGDeviceSetMatDiagType
  end interface
  
  interface HYBGDeviceSetMatIndexBase
    function d_HYBGDeviceSetMatIndexBase(Mat, type) result(res) &
        & bind(c, name = "d_HYBGDeviceSetMatIndexBase")
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function d_HYBGDeviceSetMatIndexBase
  end interface
  
  interface HYBGDeviceHybsmAnalysis
    function d_HYBGDeviceHybsmAnalysis(Mat) result(res) &
        & bind(c, name = "d_HYBGDeviceHybsmAnalysis")
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)  :: Mat
      integer(c_int)  :: res
    end function d_HYBGDeviceHybsmAnalysis
  end interface
  
  interface spsvHYBGDevice
    function d_spsvHYBGDevice(Mat, alpha, x, beta, y) result(res) &
        & bind(c, name = "d_spsvHYBGDevice")
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)          :: Mat
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int)  :: res
    end function d_spsvHYBGDevice
  end interface
  
  interface spmvHYBGDevice
    function d_spmvHYBGDevice(Mat, alpha, x, beta, y) result(res) &
        & bind(c, name = "d_spmvHYBGDevice")
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)          :: Mat
      real(c_double), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int)  :: res
    end function d_spmvHYBGDevice
  end interface
  
  interface HYBGHost2Device
    function d_HYBGHost2Device(Mat, m, n, nz, irp, ja, val) &
        & bind(c, name = "d_HYBGHost2Device") result(res)
      use iso_c_binding
      import d_Hmat
      type(d_Hmat)          :: Mat
      integer(c_int), value :: m, n, nz
      integer(c_int)        :: irp(*), ja(*)
      real(c_double)        :: val(*)
      integer(c_int)  :: res
    end function d_HYBGHost2Device
  end interface
#endif
end module d_cusparse_mod
