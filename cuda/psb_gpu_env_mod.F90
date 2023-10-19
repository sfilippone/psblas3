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
  

module psb_gpu_env_mod
  use psb_const_mod
  use iso_c_binding
  use base_cusparse_mod
!  interface psb_gpu_init
!    module procedure  psb_gpu_init
!  end interface
#if defined(HAVE_CUDA)
  use core_mod
  
  interface 
    function psb_gpuGetHandle() &
         & result(res) bind(c,name='psb_gpuGetHandle')
      use iso_c_binding   
      type(c_ptr)		:: res
    end function psb_gpuGetHandle
  end interface

  interface 
    function psb_gpuGetStream() &
         & result(res) bind(c,name='psb_gpuGetStream')
      use iso_c_binding   
      type(c_ptr)		:: res
    end function psb_gpuGetStream
  end interface

  interface 
    function psb_C_gpu_init(dev) &
         & result(res) bind(c,name='gpuInit')
      use iso_c_binding   
      integer(c_int),value	:: dev
      integer(c_int)		:: res
    end function psb_C_gpu_init
  end interface

  interface 
    function psb_cuda_getDeviceCount() &
         & result(res) bind(c,name='getDeviceCount')
      use iso_c_binding   
      integer(c_int)		:: res
    end function psb_cuda_getDeviceCount
  end interface

  interface 
    function psb_cuda_getDevice() &
         & result(res) bind(c,name='getDevice')
      use iso_c_binding   
      integer(c_int)		:: res
    end function psb_cuda_getDevice
  end interface

  interface 
    function psb_cuda_setDevice(dev) &
         & result(res) bind(c,name='setDevice')
      use iso_c_binding   
      integer(c_int), value	:: dev
      integer(c_int)		:: res
    end function psb_cuda_setDevice
  end interface
  

  interface 
    subroutine psb_gpuCreateHandle() &
         & bind(c,name='psb_gpuCreateHandle')
      use iso_c_binding   
    end subroutine psb_gpuCreateHandle
  end interface

  interface 
    subroutine psb_gpuSetStream(handle,stream) &
         & bind(c,name='psb_gpuSetStream')
      use iso_c_binding   
      type(c_ptr), value :: handle, stream
    end subroutine psb_gpuSetStream
  end interface

  interface 
    subroutine psb_gpuDestroyHandle() &
         & bind(c,name='psb_gpuDestroyHandle')
      use iso_c_binding   
    end subroutine psb_gpuDestroyHandle
  end interface

  interface 
    subroutine psb_cudaReset() &
         & bind(c,name='cudaReset')
      use iso_c_binding   
    end subroutine psb_cudaReset
  end interface

  interface 
    subroutine psb_gpuClose() &
         & bind(c,name='gpuClose')
      use iso_c_binding   
    end subroutine psb_gpuClose
  end interface
#endif

  interface 
    function psb_C_DeviceHasUVA() &
         & result(res) bind(c,name='DeviceHasUVA')
      use iso_c_binding   
      integer(c_int)		:: res
    end function psb_C_DeviceHasUVA
  end interface

  interface 
    function psb_C_get_MultiProcessors() &
         & result(res) bind(c,name='getGPUMultiProcessors')
      use iso_c_binding
      integer(c_int) :: res
    end function psb_C_get_MultiProcessors
    function psb_C_get_MemoryBusWidth() &
         & result(res) bind(c,name='getGPUMemoryBusWidth')
      use iso_c_binding
      integer(c_int) :: res
    end function psb_C_get_MemoryBusWidth
    function psb_C_get_MemoryClockRate() &
         & result(res) bind(c,name='getGPUMemoryClockRate')
      use iso_c_binding
      integer(c_int) :: res
    end function psb_C_get_MemoryClockRate
    function psb_C_get_WarpSize() &
         & result(res) bind(c,name='getGPUWarpSize')
      use iso_c_binding
      integer(c_int) :: res
    end function psb_C_get_WarpSize
    function psb_C_get_MaxThreadsPerMP() &
         & result(res) bind(c,name='getGPUMaxThreadsPerMP')
      use iso_c_binding
      integer(c_int) :: res
    end function psb_C_get_MaxThreadsPerMP
    function psb_C_get_MaxRegistersPerBlock() &
         & result(res) bind(c,name='getGPUMaxRegistersPerBlock')
      use iso_c_binding
      integer(c_int) :: res
    end function psb_C_get_MaxRegistersPerBlock
  end interface
  interface 
    subroutine psb_C_cpy_NameString(cstring) &
         & bind(c,name='cpyGPUNameString')
      use iso_c_binding
      character(c_char) :: cstring(*) 
    end subroutine psb_C_cpy_NameString
  end interface

  logical, private :: gpu_do_maybe_free_buffer = .false.

Contains
  
  function psb_gpu_get_maybe_free_buffer() result(res)
    logical :: res
    res = gpu_do_maybe_free_buffer
  end function psb_gpu_get_maybe_free_buffer

  subroutine psb_gpu_set_maybe_free_buffer(val)
    logical, intent(in) :: val
    gpu_do_maybe_free_buffer = val
  end subroutine psb_gpu_set_maybe_free_buffer
  
  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Environment handling 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!


  subroutine psb_gpu_init(ctxt,dev)
    use psb_penv_mod
    use psb_const_mod
    use psb_error_mod
    type(psb_ctxt_type), intent(in) :: ctxt
    integer, intent(in), optional   :: dev

    integer :: np, npavail, iam, info, count, dev_
    Integer(Psb_ipk_)  :: err_act

    info = psb_success_
    call psb_erractionsave(err_act)
#if defined (HAVE_CUDA)
#if defined(SERIAL_MPI) 
    iam = 0
#else
    call psb_info(ctxt,iam,np)
#endif

    count = psb_cuda_getDeviceCount()

    if (present(dev)) then 
      info = psb_C_gpu_init(dev)
    else
      if (count >0) then 
        dev_ = mod(iam,count)
      else
        dev_ = 0
      end if
      info = psb_C_gpu_init(dev_)
    end if
    if (info == 0) info = initFcusparse()
    if (info /= 0) then 
      call psb_errpush(psb_err_internal_error_,'psb_gpu_init')
      goto 9999
    end if
    call psb_gpuCreateHandle()
#endif
    call psb_erractionrestore(err_act)
    return
9999 call psb_error_handler(ctxt,err_act)

    return

  end subroutine psb_gpu_init


  subroutine psb_gpu_DeviceSync()
#if defined(HAVE_CUDA)
    call psb_cudaSync()
#endif
  end subroutine psb_gpu_DeviceSync

  function psb_gpu_getDeviceCount() result(res)
    integer :: res
#if defined(HAVE_CUDA)
    res = psb_cuda_getDeviceCount()
#else 
    res = 0
#endif
  end function psb_gpu_getDeviceCount

  subroutine psb_gpu_exit()
    integer :: res
    res =  closeFcusparse()
    call psb_gpuClose()
    call psb_cudaReset()
  end subroutine psb_gpu_exit

  function psb_gpu_DeviceHasUVA() result(res)
    logical :: res
    res =  (psb_C_DeviceHasUVA() == 1)
  end function psb_gpu_DeviceHasUVA

  function psb_gpu_MultiProcessors() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MultiProcessors()
  end function psb_gpu_MultiProcessors

  function psb_gpu_MaxRegistersPerBlock() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MaxRegistersPerBlock()
  end function psb_gpu_MaxRegistersPerBlock

  function psb_gpu_MaxThreadsPerMP() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MaxThreadsPerMP()
  end function psb_gpu_MaxThreadsPerMP

  function psb_gpu_WarpSize() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_WarpSize()
  end function psb_gpu_WarpSize

  function psb_gpu_MemoryClockRate() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MemoryClockRate()
  end function psb_gpu_MemoryClockRate

  function psb_gpu_MemoryBusWidth() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MemoryBusWidth()
  end function psb_gpu_MemoryBusWidth

  function psb_gpu_MemoryPeakBandwidth() result(res)     
    real(psb_dpk_) :: res
    ! Formula here: 2*ClockRate(KHz)*BusWidth(bit)
    ! normalization: bit/byte, KHz/MHz
    ! output: MBytes/s
    res =  2.d0*0.125d0*1.d-3*psb_C_get_MemoryBusWidth()*psb_C_get_MemoryClockRate()
  end function psb_gpu_MemoryPeakBandwidth

  function psb_gpu_DeviceName() result(res)     
    character(len=256) :: res
    character :: cstring(256)
    call psb_C_cpy_NameString(cstring)
    call stringc2f(cstring,res)
  end function psb_gpu_DeviceName


  subroutine stringc2f(cstring,fstring)
    character(c_char)        :: cstring(*)
    character(len=*)         :: fstring
    integer :: i
    
    i = 1
    do 
      if (cstring(i) == c_null_char) exit
      if (i > len(fstring)) exit
      fstring(i:i) = cstring(i)
      i = i + 1 
    end do
    do 
      if (i > len(fstring)) exit
      fstring(i:i) = " "
      i = i + 1 
    end do
    return
  end subroutine stringc2f

end module psb_gpu_env_mod
