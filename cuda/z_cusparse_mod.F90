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

module z_cusparse_mod
  use base_cusparse_mod
  
  type, bind(c) :: z_Cmat
    type(c_ptr) :: Mat = c_null_ptr
  end type z_Cmat
  
#if PSB_CUDA_SHORT_VERSION <= 10 
  type, bind(c) :: z_Hmat
    type(c_ptr) :: Mat = c_null_ptr
  end type z_Hmat
#endif
  
  interface CSRGDeviceFree
    function z_CSRGDeviceFree(Mat) result(res) &
        & bind(c, name = "z_CSRGDeviceFree")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)   :: Mat
      integer(c_int) :: res
    end function z_CSRGDeviceFree
  end interface
  
  interface CSRGDeviceSetMatType
    function z_CSRGDeviceSetMatType(Mat, type) result(res) &
        & bind(c, name = "z_CSRGDeviceSetMatType")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function z_CSRGDeviceSetMatType
  end interface
  
  interface CSRGDeviceSetMatFillMode
    function z_CSRGDeviceSetMatFillMode(Mat, type) result(res) &
        & bind(c, name = "z_CSRGDeviceSetMatFillMode")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function z_CSRGDeviceSetMatFillMode
  end interface
  
  interface CSRGDeviceSetMatDiagType
    function z_CSRGDeviceSetMatDiagType(Mat, type) result(res) &
        & bind(c, name = "z_CSRGDeviceSetMatDiagType")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function z_CSRGDeviceSetMatDiagType
  end interface
  
  interface CSRGDeviceSetMatIndexBase
    function z_CSRGDeviceSetMatIndexBase(Mat, type) result(res) &
        & bind(c, name = "z_CSRGDeviceSetMatIndexBase")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function z_CSRGDeviceSetMatIndexBase
  end interface
  
#if PSB_CUDA_SHORT_VERSION <= 10  
  interface CSRGDeviceCsrsmAnalysis
    function z_CSRGDeviceCsrsmAnalysis(Mat) result(res) &
        & bind(c, name = "z_CSRGDeviceCsrsmAnalysis")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)    :: Mat
      integer(c_int)  :: res
    end function z_CSRGDeviceCsrsmAnalysis
  end interface
#else
  interface CSRGIsNullSvBuffer
    function z_CSRGIsNullSvBuffer(Mat) result(res) &
        & bind(c, name = "z_CSRGIsNullSvBuffer")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)    :: Mat
      integer(c_int)  :: res
    end function z_CSRGIsNullSvBuffer
  end interface
#endif

  interface CSRGDeviceAlloc
    function z_CSRGDeviceAlloc(Mat, nr, nc, nz) result(res) &
        & bind(c, name = "z_CSRGDeviceAlloc")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)  :: res
    end function z_CSRGDeviceAlloc
  end interface

  interface CSRGDeviceGetParms
    function z_CSRGDeviceGetParms(Mat, nr, nc, nz) result(res) &
        & bind(c, name = "z_CSRGDeviceGetParms")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)    :: Mat
      integer(c_int)  :: nr, nc, nz
      integer(c_int)  :: res
    end function z_CSRGDeviceGetParms
  end interface
  
  interface spsvCSRGDevice
    function z_spsvCSRGDevice(Mat, alpha, x, beta, y) result(res) &
        & bind(c, name = "z_spsvCSRGDevice")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: Mat
      type(c_ptr), value    :: x
      type(c_ptr), value    :: y
      complex(c_double_complex), value :: alpha, beta
      integer(c_int)  :: res
    end function z_spsvCSRGDevice
  end interface
  
  interface spmvCSRGDevice
    function z_spmvCSRGDeviceF(mat, alpha, x, beta, y) result(res) &
        & bind(c, name = "z_spmvCSRGDeviceF")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: mat
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int)  :: res
    end function z_spmvCSRGDeviceF

    function z_spmvCSRGDeviceMV(mat, alpha, x, beta, y, idx_y) result(res) &
        & bind(c, name = "z_spmvCSRGDeviceMV")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: mat
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int), value :: idx_y
      integer(c_int)  :: res
    end function z_spmvCSRGDeviceMV

    function z_spmvCSRGDeviceVM(mat, alpha, x, idx_x, beta, y) result(res) &
        & bind(c, name = "z_spmvCSRGDeviceVM")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: mat
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int), value :: idx_x
      integer(c_int)  :: res
    end function z_spmvCSRGDeviceVM
    
    function z_spmvCSRGDeviceMMi(mat, alpha, x, idx_x, beta, y, idx_y) result(res) &
        & bind(c, name = "z_spmvCSRGDeviceMMi")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: mat
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int), value :: idx_x, idx_y
      integer(c_int)  :: res
    end function z_spmvCSRGDeviceMMi
  end interface
  
  interface CSRGHost2Device
    function z_CSRGHost2Device(Mat, m, n, nz, irp, ja, val) result(res) &
        & bind(c, name = "z_CSRGHost2Device")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: m, n, nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_double_complex)        :: val(*)
      integer(c_int)  :: res
    end function z_CSRGHost2Device
  end interface
  
  interface CSRGDevice2Host
    function z_CSRGDevice2Host(Mat, m, n, nz, irp, ja, val) result(res) &
        & bind(c, name = "z_CSRGDevice2Host")
      use iso_c_binding
      import z_Cmat
      type(z_Cmat)          :: Mat
      integer(c_int), value :: m, n, nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_double_complex)        :: val(*)
      integer(c_int)  :: res
    end function z_CSRGDevice2Host
  end interface
  
#if PSB_CUDA_SHORT_VERSION <= 10
  interface HYBGDeviceAlloc
    function z_HYBGDeviceAlloc(Mat, nr, nc, nz) result(res) &
        & bind(c, name = "z_HYBGDeviceAlloc")
      use iso_c_binding
      import z_hmat
      type(z_Hmat)          :: Mat
      integer(c_int), value :: nr, nc, nz
      integer(c_int)  :: res
    end function z_HYBGDeviceAlloc
  end interface
  
  interface HYBGDeviceFree
    function z_HYBGDeviceFree(Mat) result(res) &
        & bind(c, name = "z_HYBGDeviceFree")
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)   :: Mat
      integer(c_int) :: res
    end function z_HYBGDeviceFree
  end interface
  
  interface HYBGDeviceSetMatType 
    function z_HYBGDeviceSetMatType(Mat, type) result(res) &
        & bind(c, name = "z_HYBGDeviceSetMatType")
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function z_HYBGDeviceSetMatType
  end interface
  
  interface HYBGDeviceSetMatFillMode
    function z_HYBGDeviceSetMatFillMode(Mat, type) result(res) &
        & bind(c, name = "z_HYBGDeviceSetMatFillMode")
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function z_HYBGDeviceSetMatFillMode
  end interface
  
  interface HYBGDeviceSetMatDiagType
    function z_HYBGDeviceSetMatDiagType(Mat, type) result(res) &
        & bind(c, name = "z_HYBGDeviceSetMatDiagType")
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function z_HYBGDeviceSetMatDiagType
  end interface
  
  interface HYBGDeviceSetMatIndexBase
    function z_HYBGDeviceSetMatIndexBase(Mat, type) result(res) &
        & bind(c, name = "z_HYBGDeviceSetMatIndexBase")
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int), value :: type
      integer(c_int)  :: res
    end function z_HYBGDeviceSetMatIndexBase
  end interface
  
  interface HYBGDeviceHybsmAnalysis
    function z_HYBGDeviceHybsmAnalysis(Mat) result(res) &
        & bind(c, name = "z_HYBGDeviceHybsmAnalysis")
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)  :: Mat
      integer(c_int)  :: res
    end function z_HYBGDeviceHybsmAnalysis
  end interface
  
  interface spsvHYBGDevice
    function z_spsvHYBGDevice(Mat, alpha, x, beta, y) result(res) &
        & bind(c, name = "z_spsvHYBGDevice")
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)          :: Mat
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int)  :: res
    end function z_spsvHYBGDevice
  end interface
  
  interface spmvHYBGDevice
    function z_spmvHYBGDevice(Mat, alpha, x, beta, y) result(res) &
        & bind(c, name = "z_spmvHYBGDevice")
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)          :: Mat
      complex(c_double_complex), value :: alpha, beta
      type(c_ptr), value    :: x, y
      integer(c_int)  :: res
    end function z_spmvHYBGDevice
  end interface
  
  interface HYBGHost2Device
    function z_HYBGHost2Device(Mat, m, n, nz, irp, ja, val) &
        & bind(c, name = "z_HYBGHost2Device") result(res)
      use iso_c_binding
      import z_Hmat
      type(z_Hmat)          :: Mat
      integer(c_int), value :: m, n, nz
      integer(c_int)        :: irp(*), ja(*)
      complex(c_double_complex)        :: val(*)
      integer(c_int)  :: res
    end function z_HYBGHost2Device
  end interface
#endif
end module z_cusparse_mod
