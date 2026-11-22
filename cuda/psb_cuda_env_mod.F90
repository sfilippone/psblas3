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
  

module psb_cuda_env_mod
  use psb_const_mod
  use iso_c_binding
  use base_cusparse_mod
!  interface psb_cuda_init
!    module procedure  psb_cuda_init
!  end interface
  use core_mod
  
  interface 
    function psb_cudaGetHandle() &
         & result(res) bind(c,name='psb_cudaGetHandle')
      use iso_c_binding   
      type(c_ptr)		:: res
    end function psb_cudaGetHandle
  end interface

  interface 
    function psb_cudaGetStream() &
         & result(res) bind(c,name='psb_cudaGetStream')
      use iso_c_binding   
      type(c_ptr)		:: res
    end function psb_cudaGetStream
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
    function psb_cuda_inner_getDeviceCount() &
         & result(res) bind(c,name='getDeviceCount')
      use iso_c_binding   
      integer(c_int)		:: res
    end function psb_cuda_inner_getDeviceCount
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
    subroutine psb_cudaCreateHandle() &
         & bind(c,name='psb_cudaCreateHandle')
      use iso_c_binding   
    end subroutine psb_cudaCreateHandle
  end interface

  interface 
    subroutine psb_cudaSetStream(handle,stream) &
         & bind(c,name='psb_cudaSetStream')
      use iso_c_binding   
      type(c_ptr), value :: handle, stream
    end subroutine psb_cudaSetStream
  end interface

  interface 
    subroutine psb_cudaDestroyHandle() &
         & bind(c,name='psb_cudaDestroyHandle')
      use iso_c_binding   
    end subroutine psb_cudaDestroyHandle
  end interface

  interface 
    subroutine psb_cuda_innerReset() &
         & bind(c,name='cudaReset')
      use iso_c_binding   
    end subroutine psb_cuda_innerReset
  end interface

  interface 
    subroutine psb_cuda_innerClose() &
         & bind(c,name='gpuClose')
      use iso_c_binding   
    end subroutine psb_cuda_innerClose
  end interface

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
#if (PSB_CUDA_VERSION >= 13000)    
    function psb_C_get_MemoryClockRate(dev) &
         & result(res) bind(c,name='getGPUMemoryClockRate')
      use iso_c_binding
      integer(c_int), value :: dev
      integer(c_int) :: res
    end function psb_C_get_MemoryClockRate
#else
    function psb_C_get_MemoryClockRate() &
         & result(res) bind(c,name='getGPUMemoryClockRate')
      use iso_c_binding
      integer(c_int) :: res
    end function psb_C_get_MemoryClockRate
#endif
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

  logical, private, save :: gpu_do_maybe_free_buffer = .false.
  integer(c_int), private, save :: myDevice
  !
  ! Halo exchange on device buffers with CUDA-aware MPI. Set at psb_cuda_init
  ! from PSB_CUDA_AWARE (0/1/auto) and what the MPI library reports.
  !
  logical, private, save :: gpu_dev_comm = .false.

  interface
    function psi_c_cuda_aware_mpi() result(res) bind(c,name='psi_c_cuda_aware_mpi')
      use iso_c_binding
      integer(c_int) :: res
    end function psi_c_cuda_aware_mpi
    function psi_c_node_local_rank() result(res) bind(c,name='psi_c_node_local_rank')
      use iso_c_binding
      integer(c_int) :: res
    end function psi_c_node_local_rank
    function psi_c_cuda_aware_selftest(want) result(res) bind(c,name='psi_c_cuda_aware_selftest')
      use iso_c_binding
      integer(c_int), value :: want
      integer(c_int) :: res
    end function psi_c_cuda_aware_selftest
    function psi_c_cuda_comm_buf_acquire(nbytes,buf) result(res) &
         & bind(c,name='psi_c_cuda_comm_buf_acquire')
      use iso_c_binding
      integer(c_size_t), value :: nbytes
      type(c_ptr) :: buf
      integer(c_int) :: res
    end function psi_c_cuda_comm_buf_acquire
    subroutine psi_c_cuda_comm_buf_release() bind(c,name='psi_c_cuda_comm_buf_release')
    end subroutine psi_c_cuda_comm_buf_release
    subroutine psi_c_cuda_comm_buf_free() bind(c,name='psi_c_cuda_comm_buf_free')
    end subroutine psi_c_cuda_comm_buf_free
    function psi_c_cuda_ipc_warmup() result(res) bind(c,name='psi_c_cuda_ipc_warmup')
      use iso_c_binding
      real(c_double) :: res
    end function psi_c_cuda_ipc_warmup
  end interface

Contains

  function psb_cuda_get_dev_comm() result(res)
    logical :: res
    res = gpu_dev_comm
  end function psb_cuda_get_dev_comm

  subroutine psb_cuda_set_dev_comm(val)
    logical, intent(in) :: val
    gpu_dev_comm = val
  end subroutine psb_cuda_set_dev_comm

  !
  ! Shared device communication buffer (one per process) for the halo
  ! exchange on device buffers: .true. and BUF set if taken, .false. if
  ! CUDA-aware exchange is off, the buffer is held by another exchange or
  ! it cannot be allocated (the caller then uses the host buffers).
  !
  function psb_cuda_comm_buf_acquire(nbytes,buf) result(res)
    integer(c_size_t), intent(in) :: nbytes
    type(c_ptr), intent(out) :: buf
    logical :: res
    buf = c_null_ptr
    res = .false.
    if (.not.gpu_dev_comm) return
    res = (psi_c_cuda_comm_buf_acquire(nbytes,buf) == 0)
  end function psb_cuda_comm_buf_acquire

  subroutine psb_cuda_comm_buf_release()
    call psi_c_cuda_comm_buf_release()
  end subroutine psb_cuda_comm_buf_release
  
  function psb_cuda_get_maybe_free_buffer() result(res)
    logical :: res
    res = gpu_do_maybe_free_buffer
  end function psb_cuda_get_maybe_free_buffer

  subroutine psb_cuda_set_maybe_free_buffer(val)
    logical, intent(in) :: val
    gpu_do_maybe_free_buffer = val
  end subroutine psb_cuda_set_maybe_free_buffer
  
  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Environment handling 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!


  subroutine psb_cuda_init(ctxt,dev)
    use psb_penv_mod
    use psb_const_mod
    use psb_error_mod
    type(psb_ctxt_type), intent(in) :: ctxt
    integer, intent(in), optional   :: dev

    integer :: np, npavail, iam, info, count, dev_, mode
    integer(c_int) :: want
    real(c_double) :: twarm
    Integer(Psb_ipk_)  :: err_act

    info = psb_success_
    call psb_erractionsave(err_act)
#if defined(PSB_SERIAL_MPI) 
    iam = 0
#else
    call psb_info(ctxt,iam,np)
#endif

    count = psb_cuda_getDeviceCount()

    if (present(dev)) then 
      info = psb_C_gpu_init(dev)
      myDevice = dev
    else
      if (count >0) then
        !
        ! Device chosen by the rank inside the node, so that it is correct
        ! whatever the placement of the ranks; with one visible device per
        ! process (e.g. bsc_affinity) this is always device 0. Making all the
        ! devices of the node visible lets CUDA-aware MPI use CUDA IPC/NVLink
        ! between the GPUs of the same node.
        !
        dev_ = mod(psi_c_node_local_rank(),count)
      else
        dev_ = 0
      end if
      info = psb_C_gpu_init(dev_)
      myDevice = dev_
    end if
    if (info == 0) info = initFcusparse()
    if (info /= 0) then 
      call psb_errpush(psb_err_internal_error_,'psb_cuda_init')
      goto 9999
    end if
    call psb_cudaCreateHandle()
    !
    ! Halo exchange on device buffers only if requested/allowed (mode) AND a
    ! self test with device buffers succeeds on every process; otherwise the
    ! library silently keeps the host-buffer path (rollback).
    !
    mode = psi_c_cuda_aware_mpi()
    if ((mode == 1).or.(mode == 3)) then
      want = 1
    else
      want = 0
    end if
    ! Collective on all processes, decision consistent everywhere
    gpu_dev_comm = (psi_c_cuda_aware_selftest(want) == 1)
    ! Open the device-to-device paths inside the node now (cheap on an
    ! almost empty device) instead of at the first halo exchange
    twarm = -1.0_c_double
    if (gpu_dev_comm) twarm = psi_c_cuda_ipc_warmup()
    if (iam == 0) then
      if (gpu_dev_comm) then
        write(psb_out_unit,'(a)') ' PSBLAS CUDA: halo exchange on device buffers (CUDA-aware MPI)'
        if (twarm >= 0.0_c_double) write(psb_out_unit,'(a,f10.1,a)') &
             & ' PSBLAS CUDA: IPC warm-up inside the node:',twarm,' ms'
      else
        write(psb_out_unit,'(a)') ' PSBLAS CUDA: halo exchange through host buffers'
        if (mode == 2) write(psb_out_unit,'(a)') &
             & ' PSBLAS CUDA: PSB_CUDA_AWARE=1 ignored, MPI reports no CUDA support'
        if ((mode == 1).or.(mode == 3)) write(psb_out_unit,'(a)') &
             & ' PSBLAS CUDA: CUDA-aware MPI self test failed, using host buffers'
      end if
    end if
    call psb_erractionrestore(err_act)
    return
9999 call psb_error_handler(ctxt,err_act)

    return

  end subroutine psb_cuda_init


  subroutine psb_cuda_DeviceSync()
    call psb_cudaSync()
  end subroutine psb_cuda_DeviceSync

  function psb_cuda_getDeviceCount() result(res)
    integer :: res
    res = psb_cuda_inner_getDeviceCount()
  end function psb_cuda_getDeviceCount

  subroutine psb_cuda_exit()
    integer :: res
    call psi_c_cuda_comm_buf_free()
    res =  closeFcusparse()
    call psb_cuda_innerClose()
    call psb_cuda_innerReset()
  end subroutine psb_cuda_exit

  function psb_cuda_DeviceHasUVA() result(res)
    logical :: res
    res =  (psb_C_DeviceHasUVA() == 1)
  end function psb_cuda_DeviceHasUVA

  function psb_cuda_MultiProcessors() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MultiProcessors()
  end function psb_cuda_MultiProcessors

  function psb_cuda_MaxRegistersPerBlock() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MaxRegistersPerBlock()
  end function psb_cuda_MaxRegistersPerBlock

  function psb_cuda_MaxThreadsPerMP() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MaxThreadsPerMP()
  end function psb_cuda_MaxThreadsPerMP

  function psb_cuda_WarpSize() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_WarpSize()
  end function psb_cuda_WarpSize

  function psb_cuda_MemoryClockRate() result(res)     
    integer(psb_ipk_) :: res
#if (PSB_CUDA_VERSION >= 13000)    
    res =  psb_C_get_MemoryClockRate(myDevice)
#else
    res =  psb_C_get_MemoryClockRate()
#endif
  end function psb_cuda_MemoryClockRate

  function psb_cuda_MemoryBusWidth() result(res)     
    integer(psb_ipk_) :: res
    res =  psb_C_get_MemoryBusWidth()
  end function psb_cuda_MemoryBusWidth

  function psb_cuda_MemoryPeakBandwidth() result(res)     
    real(psb_dpk_) :: res
    ! Formula here: 2*ClockRate(KHz)*BusWidth(bit)
    ! normalization: bit/byte, KHz/MHz
    ! output: MBytes/s
#if PSB_CUDA_VERSION >= 13000
    res =  2.d0*0.125d0*1.d-3*psb_C_get_MemoryBusWidth()*psb_C_get_MemoryClockRate(myDevice)
#else
    res =  2.d0*0.125d0*1.d-3*psb_C_get_MemoryBusWidth()*psb_C_get_MemoryClockRate()
#endif
  end function psb_cuda_MemoryPeakBandwidth

  function psb_cuda_DeviceName() result(res)     
    character(len=256) :: res
    character :: cstring(256)
    call psb_C_cpy_NameString(cstring)
    call stringc2f(cstring,res)
  end function psb_cuda_DeviceName


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

end module psb_cuda_env_mod
