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
  

module psb_s_cuda_vect_mod
  use iso_c_binding
  use psb_const_mod
  use psb_error_mod
  use psb_s_vect_mod
  use psb_cuda_env_mod
  use psb_i_vect_mod
  use psb_i_cuda_vect_mod
  use psb_i_vectordev_mod
  use psb_s_vectordev_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 
  
  type, extends(psb_s_base_vect_type) ::  psb_s_vect_cuda
    integer     :: state      = is_host
    type(c_ptr) :: deviceVect = c_null_ptr
    real(c_float), allocatable :: pinned_buffer(:)
    type(c_ptr) :: dt_p_buf = c_null_ptr
    real(c_float), allocatable :: buffer(:)
    type(c_ptr) :: dt_buf = c_null_ptr
    integer :: dt_buf_sz = 0
    type(c_ptr) :: i_buf = c_null_ptr
    integer :: i_buf_sz = 0
  contains
    procedure, pass(x) :: get_nrows => s_cuda_get_nrows
    procedure, nopass  :: get_fmt   => s_cuda_get_fmt

    procedure, pass(x) :: all      => s_cuda_all
    procedure, pass(x) :: zero     => s_cuda_zero
    procedure, pass(x) :: asb_m    => s_cuda_asb_m
    procedure, pass(x) :: sync     => s_cuda_sync
    procedure, pass(x) :: sync_space => s_cuda_sync_space
    procedure, pass(x) :: bld_x    => s_cuda_bld_x
    procedure, pass(x) :: bld_mn   => s_cuda_bld_mn
    procedure, pass(x) :: free     => s_cuda_free
    procedure, pass(x) :: ins_a    => s_cuda_ins_a
    procedure, pass(x) :: ins_v    => s_cuda_ins_v
    procedure, pass(x) :: is_host  => s_cuda_is_host
    procedure, pass(x) :: is_dev   => s_cuda_is_dev
    procedure, pass(x) :: is_sync  => s_cuda_is_sync
    procedure, pass(x) :: set_host => s_cuda_set_host
    procedure, pass(x) :: set_dev  => s_cuda_set_dev
    procedure, pass(x) :: set_sync => s_cuda_set_sync
    procedure, pass(x) :: set_scal => s_cuda_set_scal
!!$    procedure, pass(x) :: set_vect => s_cuda_set_vect
    procedure, pass(x) :: gthzv_x  => s_cuda_gthzv_x
    procedure, pass(y) :: sctb     => s_cuda_sctb
    procedure, pass(y) :: sctb_x   => s_cuda_sctb_x
    procedure, pass(x) :: gthzbuf  => s_cuda_gthzbuf
    procedure, pass(y) :: sctb_buf => s_cuda_sctb_buf
    procedure, pass(x) :: new_buffer   => s_cuda_new_buffer
    procedure, nopass  :: device_wait  => s_cuda_device_wait
    procedure, pass(x) :: free_buffer  => s_cuda_free_buffer
    procedure, pass(x) :: maybe_free_buffer  => s_cuda_maybe_free_buffer
    procedure, pass(x) :: dot_v    => s_cuda_dot_v
    procedure, pass(x) :: dot_a    => s_cuda_dot_a
    procedure, pass(y) :: axpby_v  => s_cuda_axpby_v
    procedure, pass(y) :: axpby_a  => s_cuda_axpby_a
    procedure, pass(z) :: abgdxyz  => s_cuda_abgdxyz
    procedure, pass(y) :: mlt_v    => s_cuda_mlt_v
    procedure, pass(y) :: mlt_a    => s_cuda_mlt_a
    procedure, pass(z) :: mlt_a_2  => s_cuda_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => s_cuda_mlt_v_2
    procedure, pass(x) :: scal     => s_cuda_scal
    procedure, pass(x) :: nrm2     => s_cuda_nrm2
    procedure, pass(x) :: amax     => s_cuda_amax
    procedure, pass(x) :: asum     => s_cuda_asum
    procedure, pass(x) :: absval1  => s_cuda_absval1
    procedure, pass(x) :: absval2  => s_cuda_absval2

    final              :: s_cuda_vect_finalize
  end type psb_s_vect_cuda

  public  :: psb_s_vect_cuda_
  private :: constructor
  interface psb_s_vect_cuda_
    module procedure constructor
  end interface psb_s_vect_cuda_

contains
  
  function constructor(x) result(this)
    real(psb_spk_)       :: x(:)
    type(psb_s_vect_cuda) :: this
    integer(psb_ipk_) :: info

    this%v = x
    call this%asb(size(x),info)

  end function constructor
    
  subroutine s_cuda_device_wait()
    call psb_cudaSync()
  end subroutine s_cuda_device_wait

  subroutine s_cuda_new_buffer(n,x,info)
    use psb_realloc_mod
    use psb_cuda_env_mod
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info


    if (psb_cuda_DeviceHasUVA()) then       
      if (allocated(x%combuf)) then 
        if (size(x%combuf)<n) then 
          call inner_unregister(x%combuf)
          deallocate(x%combuf,stat=info)
        end if
      end if
      if (.not.allocated(x%combuf)) then 
        call psb_realloc(n,x%combuf,info)
        if (info == 0) info = inner_register(x%combuf,x%dt_p_buf)        
      end if
    else
      call psb_realloc(n,x%combuf,info)
    end if
    if (c_associated(x%dt_buf)) then 
      if (x%dt_buf_sz < n) then
        call freeFloat(x%dt_buf)
        x%dt_buf = c_null_ptr
        x%dt_buf_sz = 0
      end if
    end if
    if (.not.c_associated(x%dt_buf)) then
      info       =  allocateFloat(x%dt_buf,n)
      x%dt_buf_sz = n
    end if
    if (c_associated(x%i_buf)) then 
      if (x%i_buf_sz < n) then
        call freeInt(x%i_buf)
        x%i_buf = c_null_ptr
        x%i_buf_sz = 0
      end if
    end if
    if (.not.c_associated(x%i_buf)) then
      info       =  allocateInt(x%i_buf,n)
      x%i_buf_sz = n
    end if

  end subroutine s_cuda_new_buffer

  subroutine s_cuda_maybe_free_buffer(x,info)
    use psb_realloc_mod
    use psb_cuda_env_mod
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    info = 0 
    if (psb_cuda_get_maybe_free_buffer()) &
         & call x%free_buffer(info)
  end subroutine s_cuda_maybe_free_buffer

  subroutine s_cuda_free_buffer(x,info)
    use psb_realloc_mod
    use psb_cuda_env_mod
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info
    
    if (allocated(x%pinned_buffer)) then 
      call inner_unregister(x%pinned_buffer)
      deallocate(x%pinned_buffer, stat=info)
    end if
    if (allocated(x%combuf)) then 
      if (psb_cuda_DeviceHasUVA()) &
           & call inner_unregister(x%combuf)
      deallocate(x%combuf, stat=info)
    end if
    if (allocated(x%buffer)) then 
      deallocate(x%buffer, stat=info)
    end if
    if (c_associated(x%dt_buf)) then 
      call freeFloat(x%dt_buf)
      x%dt_buf  = c_null_ptr
    end if
    if (c_associated(x%i_buf)) then 
      call freeInt(x%i_buf)
      x%i_buf = c_null_ptr
    end if
    x%dt_buf_sz=0
    x%i_buf_sz=0

  end subroutine s_cuda_free_buffer

  subroutine s_cuda_gthzv_x(i,n,idx,x,y)
    use psb_cuda_env_mod
    use psi_serial_mod
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    real(psb_spk_) ::  y(:)
    class(psb_s_vect_cuda) :: x
    integer ::  info, ni

    info = 0 

    select type(ii=> idx) 
    class is (psb_i_vect_cuda) 
      if (ii%is_host()) call ii%sync()
      if (x%is_host())  call x%sync()

      if (psb_cuda_DeviceHasUVA()) then 
        !
        ! Only need a sync in this branch; in the others
        ! cudamemCpy acts as a sync point.
        !
        if (allocated(x%pinned_buffer)) then  
          if (size(x%pinned_buffer) < n) then 
            call inner_unregister(x%pinned_buffer)
            deallocate(x%pinned_buffer, stat=info)
          end if
        end if

        if (.not.allocated(x%pinned_buffer)) then
          allocate(x%pinned_buffer(n),stat=info)
          if (info == 0) info = inner_register(x%pinned_buffer,x%dt_p_buf)        
          if (info /= 0) &
               & write(0,*) 'Error from inner_register ',info
        endif
        info = igathMultiVecDeviceFloatVecIdx(x%deviceVect,&
             & 0, n, i, ii%deviceVect, 1, x%dt_p_buf, 1)
        call psb_cudaSync()
        y(1:n) = x%pinned_buffer(1:n)

      else
        if (allocated(x%buffer)) then 
          if (size(x%buffer) < n) then 
            deallocate(x%buffer, stat=info)
          end if
        end if

        if (.not.allocated(x%buffer)) then
          allocate(x%buffer(n),stat=info)
        end if

        if (x%dt_buf_sz < n) then 
          if (c_associated(x%dt_buf)) then 
            call freeFloat(x%dt_buf)
            x%dt_buf = c_null_ptr
          end if
          info =  allocateFloat(x%dt_buf,n)
          x%dt_buf_sz=n
        end if
        if (info == 0) &
             & info = igathMultiVecDeviceFloatVecIdx(x%deviceVect,&
             & 0, n, i, ii%deviceVect, 1, x%dt_buf, 1)
        if (info == 0) &
             &  info = readFloat(x%dt_buf,y,n)

      endif

    class default
      ! Do not go for brute force, but move the index vector
      ni = size(ii%v)

      if (x%i_buf_sz < ni) then 
        if (c_associated(x%i_buf)) then 
          call freeInt(x%i_buf)
          x%i_buf = c_null_ptr
        end if
        info =  allocateInt(x%i_buf,ni)
        x%i_buf_sz=ni
      end if
      if (allocated(x%buffer)) then 
        if (size(x%buffer) < n) then 
          deallocate(x%buffer, stat=info)
        end if
      end if

      if (.not.allocated(x%buffer)) then
        allocate(x%buffer(n),stat=info)
      end if

      if (x%dt_buf_sz < n) then 
        if (c_associated(x%dt_buf)) then 
          call freeFloat(x%dt_buf)
            x%dt_buf = c_null_ptr
          end if
        info =  allocateFloat(x%dt_buf,n)
        x%dt_buf_sz=n
      end if

      if (info == 0) &
           & info = writeInt(x%i_buf,ii%v,ni)
      if (info == 0) &
           & info = igathMultiVecDeviceFloat(x%deviceVect,&
           & 0, n, i, x%i_buf, 1, x%dt_buf, 1)
      if (info == 0) &
           &  info = readFloat(x%dt_buf,y,n)

    end select
    
  end subroutine s_cuda_gthzv_x

  subroutine s_cuda_gthzbuf(i,n,idx,x)
    use psb_cuda_env_mod
    use psi_serial_mod
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    class(psb_s_vect_cuda)       :: x
    integer ::  info, ni

    info = 0 
!!$    write(0,*) 'Starting gth_zbuf'
    if (.not.allocated(x%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'gthzbuf')
      return
    end if

    select type(ii=> idx) 
    class is (psb_i_vect_cuda) 
      if (ii%is_host()) call ii%sync()
      if (x%is_host())  call x%sync()

      if (psb_cuda_DeviceHasUVA()) then 
        info = igathMultiVecDeviceFloatVecIdx(x%deviceVect,&
             & 0, n, i, ii%deviceVect, i,x%dt_p_buf, 1)

      else
        info = igathMultiVecDeviceFloatVecIdx(x%deviceVect,&
             & 0, n, i, ii%deviceVect, i,x%dt_buf, 1)
        if (info == 0) &
             &  info = readFloat(i,x%dt_buf,x%combuf(i:),n,1)
      endif

    class default
      ! Do not go for brute force, but move the index vector
      ni = size(ii%v)
      info = 0 
      if (.not.c_associated(x%i_buf)) then 
        info =  allocateInt(x%i_buf,ni)
        x%i_buf_sz=ni
      end if
      if (info == 0) &
           & info = writeInt(i,x%i_buf,ii%v(i:),n,1)

      if (info == 0) &
           & info = igathMultiVecDeviceFloat(x%deviceVect,&
           & 0, n, i, x%i_buf, i,x%dt_buf, 1)
      
      if (info == 0) &
           &  info = readFloat(i,x%dt_buf,x%combuf(i:),n,1)

    end select

  end subroutine s_cuda_gthzbuf

  subroutine s_cuda_sctb(n,idx,x,beta,y)
    implicit none
    !use psb_const_mod
    integer(psb_ipk_)     :: n, idx(:)
    real(psb_spk_)        :: beta, x(:)
    class(psb_s_vect_cuda) :: y
    integer(psb_ipk_)     :: info

    if (n == 0) return
    
    if (y%is_dev())  call y%sync()
          
    call y%psb_s_base_vect_type%sctb(n,idx,x,beta)
    call y%set_host()

  end subroutine s_cuda_sctb

  subroutine s_cuda_sctb_x(i,n,idx,x,beta,y)
    use psb_cuda_env_mod
    use psi_serial_mod
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    real(psb_spk_) :: beta, x(:)
    class(psb_s_vect_cuda) :: y
    integer :: info, ni

    select type(ii=> idx) 
    class is (psb_i_vect_cuda) 
      if (ii%is_host()) call ii%sync()
      if (y%is_host())  call y%sync()
      
      ! 
      if (psb_cuda_DeviceHasUVA()) then 
        if (allocated(y%pinned_buffer)) then  
          if (size(y%pinned_buffer) < n) then 
            call inner_unregister(y%pinned_buffer)
            deallocate(y%pinned_buffer, stat=info)
          end if
        end if

        if (.not.allocated(y%pinned_buffer)) then
          allocate(y%pinned_buffer(n),stat=info)
          if (info == 0) info = inner_register(y%pinned_buffer,y%dt_p_buf)        
          if (info /= 0) &
               & write(0,*) 'Error from inner_register ',info
        endif
        y%pinned_buffer(1:n) = x(1:n) 
        info = iscatMultiVecDeviceFloatVecIdx(y%deviceVect,&
             & 0, n, i, ii%deviceVect, 1, y%dt_p_buf, 1,beta)
      else
        
        if (allocated(y%buffer)) then 
          if (size(y%buffer) < n) then 
            deallocate(y%buffer, stat=info)
          end if
        end if
        
        if (.not.allocated(y%buffer)) then
          allocate(y%buffer(n),stat=info)
        end if

        if (y%dt_buf_sz < n) then 
          if (c_associated(y%dt_buf)) then 
            call freeFloat(y%dt_buf)
            y%dt_buf = c_null_ptr
          end if
          info =  allocateFloat(y%dt_buf,n)
          y%dt_buf_sz=n
        end if
        info = writeFloat(y%dt_buf,x,n)
        info = iscatMultiVecDeviceFloatVecIdx(y%deviceVect,&
             & 0, n, i, ii%deviceVect, 1, y%dt_buf, 1,beta)

      end if
      
    class default
      ni = size(ii%v)
      
      if (y%i_buf_sz < ni) then 
        if (c_associated(y%i_buf)) then 
          call freeInt(y%i_buf)
          y%i_buf = c_null_ptr
        end if
        info =  allocateInt(y%i_buf,ni)
        y%i_buf_sz=ni
      end if
      if (allocated(y%buffer)) then 
        if (size(y%buffer) < n) then 
          deallocate(y%buffer, stat=info)
        end if
      end if

      if (.not.allocated(y%buffer)) then
        allocate(y%buffer(n),stat=info)
      end if

      if (y%dt_buf_sz < n) then 
        if (c_associated(y%dt_buf)) then 
          call freeFloat(y%dt_buf)
          y%dt_buf = c_null_ptr
        end if
        info =  allocateFloat(y%dt_buf,n)
        y%dt_buf_sz=n
      end if

      if (info == 0) &
           & info = writeInt(y%i_buf,ii%v(i:i+n-1),n)
      info = writeFloat(y%dt_buf,x,n)
      info = iscatMultiVecDeviceFloat(y%deviceVect,&
           & 0, n, 1, y%i_buf, 1, y%dt_buf, 1,beta)


    end select
    !
    !  Need a sync here to make sure we are not reallocating
    !  the buffers before iscatMulti has finished.
    !
    call psb_cudaSync()       
    call y%set_dev()

  end subroutine s_cuda_sctb_x

  subroutine s_cuda_sctb_buf(i,n,idx,beta,y)
    use psi_serial_mod
    use psb_cuda_env_mod
    implicit none 
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    real(psb_spk_) :: beta
    class(psb_s_vect_cuda) :: y
    integer(psb_ipk_) :: info, ni
    
!!$    write(0,*) 'Starting sctb_buf'
    if (.not.allocated(y%combuf)) then 
      call psb_errpush(psb_err_alloc_dealloc_,'sctb_buf')
      return
    end if
    

    select type(ii=> idx) 
    class is (psb_i_vect_cuda) 
              
      if (ii%is_host()) call ii%sync()
      if (y%is_host())  call y%sync()
      if (psb_cuda_DeviceHasUVA()) then 
        info = iscatMultiVecDeviceFloatVecIdx(y%deviceVect,&
             & 0, n, i, ii%deviceVect, i, y%dt_p_buf, 1,beta)
      else 
        info = writeFloat(i,y%dt_buf,y%combuf(i:),n,1)
        info = iscatMultiVecDeviceFloatVecIdx(y%deviceVect,&
             & 0, n, i, ii%deviceVect, i, y%dt_buf, 1,beta)

      end if

    class default
      !call y%sct(n,ii%v(i:),x,beta)
      ni = size(ii%v)
      info = 0 
      if (.not.c_associated(y%i_buf)) then 
        info =  allocateInt(y%i_buf,ni)
        y%i_buf_sz=ni
      end if
      if (info == 0) &
           & info = writeInt(i,y%i_buf,ii%v(i:),n,1)
      if (info == 0) &
           & info = writeFloat(i,y%dt_buf,y%combuf(i:),n,1)
      if (info == 0) info = iscatMultiVecDeviceFloat(y%deviceVect,&
           & 0, n, i, y%i_buf, i, y%dt_buf, 1,beta)
    end select
!!$    write(0,*) 'Done sctb_buf'

  end subroutine s_cuda_sctb_buf


  subroutine s_cuda_bld_x(x,this)
    use psb_base_mod
    real(psb_spk_), intent(in)           :: this(:)
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this),x%v,info)
    if (info /= 0) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,'s_cuda_bld_x',&
           & i_err=(/size(this),izero,izero,izero,izero/))
    end if
    x%v(:)  = this(:) 
    call x%set_host()
    call x%sync()

  end subroutine s_cuda_bld_x

  subroutine s_cuda_bld_mn(x,n)
    integer(psb_mpk_), intent(in) :: n
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%all(n,info)
    if (info /= 0) then 
      call psb_errpush(info,'s_cuda_bld_n',i_err=(/n,n,n,n,n/))
    end if
    
  end subroutine s_cuda_bld_mn

  subroutine s_cuda_set_host(x)
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    
    x%state = is_host
  end subroutine s_cuda_set_host

  subroutine s_cuda_set_dev(x)
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    
    x%state = is_dev
  end subroutine s_cuda_set_dev

  subroutine s_cuda_set_sync(x)
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    
    x%state = is_sync
  end subroutine s_cuda_set_sync

  function s_cuda_is_dev(x) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(in) :: x
    logical  :: res
  
    res = (x%state == is_dev)
  end function s_cuda_is_dev
  
  function s_cuda_is_host(x) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(in) :: x
    logical  :: res

    res = (x%state == is_host)
  end function s_cuda_is_host

  function s_cuda_is_sync(x) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(in) :: x
    logical  :: res

    res = (x%state == is_sync)
  end function s_cuda_is_sync

  
  function s_cuda_get_nrows(x) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v)
  end function s_cuda_get_nrows

  function s_cuda_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'sGPU'
  end function s_cuda_get_fmt
  
  subroutine s_cuda_all(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)      :: n
    class(psb_s_vect_cuda), intent(out) :: x
    integer(psb_ipk_), intent(out)     :: info
    
    call psb_realloc(n,x%v,info)
    if (info == 0) call x%set_host()
    if (info == 0) call x%sync_space(info)
    if (info /= 0) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,'s_cuda_all',&
           & i_err=(/n,n,n,n,n/))
    end if
  end subroutine s_cuda_all

  subroutine s_cuda_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    ! Since we are overwriting, make sure to do it
    ! on the GPU side
    call x%set_dev()
    call x%set_scal(szero)
  end subroutine s_cuda_zero

  subroutine s_cuda_asb_m(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_mpk_), intent(in)        :: n
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_mpk_) :: nd
    
    if (x%is_dev()) then 
      nd  = getMultiVecDeviceSize(x%deviceVect)
      if (nd < n) then 
        call x%sync()
        call x%psb_s_base_vect_type%asb(n,info)      
        if (info == psb_success_) call x%sync_space(info)
        call x%set_host()
      end if
    else   !
      if (x%get_nrows()<n) then 
        call x%psb_s_base_vect_type%asb(n,info)      
        if (info == psb_success_) call x%sync_space(info)
        call x%set_host()      
      end if
    end if

  end subroutine s_cuda_asb_m

  subroutine s_cuda_sync_space(x,info)
    use psb_base_mod, only : psb_realloc
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(out)       :: info 
    integer(psb_ipk_) :: nh, nd
    
    info = 0
    if (x%is_dev()) then 
      ! 
      if (.not.allocated(x%v)) then 
        nh = 0
      else
        nh    = size(x%v)
      end if
      nd  = getMultiVecDeviceSize(x%deviceVect)
      if (nh < nd ) then 
        call psb_realloc(nd,x%v,info)
      end if
    else  !    if (x%is_host()) then 
      if (.not.allocated(x%v)) then 
        nh = 0
      else
        nh    = size(x%v)
      end if
      if (c_associated(x%deviceVect)) then 
        nd  = getMultiVecDeviceSize(x%deviceVect)
        if (nd < nh ) then 
          call freeMultiVecDevice(x%deviceVect)
          x%deviceVect=c_null_ptr
        end if
      end if
      if (.not.c_associated(x%deviceVect)) then 
        info = FallocMultiVecDevice(x%deviceVect,1,nh,spgpu_type_float)
        if  (info /= 0) then 
          if (info == spgpu_outofmem) then 
            info = psb_err_alloc_request_
          end if
        end if
      end if
    end if
    
  end subroutine s_cuda_sync_space

  subroutine s_cuda_sync(x)
    use psb_base_mod, only : psb_realloc
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_) :: n,info
    
    info = 0
    if (x%is_host()) then 
      if (.not.c_associated(x%deviceVect)) then 
        n    = size(x%v)
        info = FallocMultiVecDevice(x%deviceVect,1,n,spgpu_type_float)
      end if
      if (info == 0) &
           & info = writeMultiVecDevice(x%deviceVect,x%v)
    else if (x%is_dev()) then 
      n    = getMultiVecDeviceSize(x%deviceVect)
      if (.not.allocated(x%v)) then 
!!$        write(0,*) 'Incoherent situation : x%v not allocated'
        call psb_realloc(n,x%v,info)
      end if
      if ((n > size(x%v)).or.(n > x%get_nrows())) then 
!!$        write(0,*) 'Incoherent situation : sizes',n,size(x%v),x%get_nrows()
        call psb_realloc(n,x%v,info)
      end if
      info = readMultiVecDevice(x%deviceVect,x%v)
    end if
    if (info == 0)  call x%set_sync()
    if (info /= 0) then
      info=psb_err_internal_error_
      call psb_errpush(info,'s_cuda_sync')
    end if
    
  end subroutine s_cuda_sync

  subroutine s_cuda_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_vect_cuda), intent(inout)  :: x
    integer(psb_ipk_), intent(out)        :: info
    
    info = 0  
    if (allocated(x%v)) deallocate(x%v, stat=info)
    if (c_associated(x%deviceVect)) then
!!$      write(0,*)'d_cuda_free Calling freeMultiVecDevice'
      call freeMultiVecDevice(x%deviceVect)
      x%deviceVect=c_null_ptr
    end if
    call x%free_buffer(info)
    call x%set_sync()
  end subroutine s_cuda_free

  subroutine s_cuda_set_scal(x,val,first,last)
    class(psb_s_vect_cuda), intent(inout) :: x
    real(psb_spk_), intent(in)           :: val
    integer(psb_ipk_), optional :: first, last
    
    integer(psb_ipk_) :: info, first_, last_
    
    first_ = 1
    last_  = x%get_nrows()
    if (present(first)) first_ = max(1,first)
    if (present(last))  last_  = min(last,last_)
    
    info = setScalDevice(val,first_,last_,1,x%deviceVect)
    call x%set_dev()
    
  end subroutine s_cuda_set_scal
!!$
!!$  subroutine s_cuda_set_vect(x,val)
!!$    class(psb_s_vect_cuda), intent(inout) :: x
!!$    real(psb_spk_), intent(in)           :: val(:)
!!$    integer(psb_ipk_) :: nr
!!$    integer(psb_ipk_) :: info
!!$
!!$    if (x%is_dev()) call x%sync()
!!$    call x%psb_s_base_vect_type%set_vect(val)
!!$    call x%set_host()
!!$
!!$  end subroutine s_cuda_set_vect



  function s_cuda_dot_v(n,x,y) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(inout)       :: x
    class(psb_s_base_vect_type), intent(inout) :: y
    integer(psb_ipk_), intent(in)              :: n
    real(psb_spk_)                :: res
    real(psb_spk_), external      :: ddot
    integer(psb_ipk_) :: info
    
    res = szero
    !
    ! Note: this is the gpu implementation.
    !  When we get here, we are sure that X is of
    !  TYPE psb_s_vect
    !
    select type(yy => y)
    type is (psb_s_base_vect_type)
      if (x%is_dev()) call x%sync()
      res = ddot(n,x%v,1,yy%v,1)
    type is (psb_s_vect_cuda)
      if (x%is_host()) call x%sync()
      if (yy%is_host()) call yy%sync()
      info = dotMultiVecDevice(res,n,x%deviceVect,yy%deviceVect)
      if (info /= 0) then 
        info = psb_err_internal_error_
        call psb_errpush(info,'s_cuda_dot_v')
      end if

    class default
      ! y%sync is done in dot_a
      call x%sync()      
      res = y%dot(n,x%v)
    end select

  end function s_cuda_dot_v

  function s_cuda_dot_a(n,x,y) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    real(psb_spk_), intent(in)           :: y(:)
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                :: res
    real(psb_spk_), external      :: ddot
    
    if (x%is_dev()) call x%sync()
    res = ddot(n,y,1,x%v,1)

  end function s_cuda_dot_a
    
  subroutine s_cuda_axpby_v(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)              :: m
    class(psb_s_base_vect_type), intent(inout) :: x
    class(psb_s_vect_cuda), intent(inout)       :: y
    real(psb_spk_), intent (in)                :: alpha, beta
    integer(psb_ipk_), intent(out)             :: info
    integer(psb_ipk_) :: nx, ny

    info = psb_success_

    select type(xx => x)
    type is (psb_s_vect_cuda)
      ! Do something different here 
      if ((beta /= szero).and.y%is_host())&
           &  call y%sync()
      if (xx%is_host()) call xx%sync()
      nx = getMultiVecDeviceSize(xx%deviceVect)
      ny = getMultiVecDeviceSize(y%deviceVect)
      if ((nx<m).or.(ny<m)) then
        info = psb_err_internal_error_
      else
        info = axpbyMultiVecDevice(m,alpha,xx%deviceVect,beta,y%deviceVect)
      end if
      call y%set_dev()
    class default
      ! Do it on the host side
      if ((alpha /= szero).and.(x%is_dev()))&
         & call x%sync()
      call y%axpby(m,alpha,x%v,beta,info)
    end select

  end subroutine s_cuda_axpby_v

  subroutine s_cuda_abgdxyz(m,alpha, beta, gamma,delta,x, y, z, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_s_base_vect_type), intent(inout)  :: x
    class(psb_s_base_vect_type), intent(inout)  :: y
    class(psb_s_vect_cuda), intent(inout)  :: z
    real(psb_spk_), intent (in)       :: alpha, beta, gamma, delta
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: nx, ny, nz
    logical :: gpu_done

    info = psb_success_
    
    if (.true.) then
      gpu_done = .false.
      select type(xx => x)
      class is (psb_s_vect_cuda)
        select type(yy => y)
        class is (psb_s_vect_cuda)
          select type(zz => z)
          class is (psb_s_vect_cuda)
            ! Do something different here 
            if ((beta /= szero).and.yy%is_host())&
                 &  call yy%sync()
            if ((delta /= szero).and.zz%is_host())&
                 &  call zz%sync()
            if (xx%is_host()) call xx%sync()
            nx = getMultiVecDeviceSize(xx%deviceVect)
            ny = getMultiVecDeviceSize(yy%deviceVect)
            nz = getMultiVecDeviceSize(zz%deviceVect)
            if ((nx<m).or.(ny<m).or.(nz<m)) then
              info = psb_err_internal_error_
            else
              info = abgdxyzMultiVecDevice(m,alpha,beta,gamma,delta,&
                   & xx%deviceVect,yy%deviceVect,zz%deviceVect)
            end if
            call yy%set_dev()
            call zz%set_dev()
            gpu_done = .true.
          end select
        end select
      end select
      
      if (.not.gpu_done) then 
        if (x%is_host()) call x%sync()
        if (y%is_host()) call y%sync()
        if (z%is_host()) call z%sync()
        call y%axpby(m,alpha,x,beta,info)
        call z%axpby(m,gamma,y,delta,info)
      end if
    else

      if (x%is_host()) call x%sync()
      if (y%is_host()) call y%sync()
      if (z%is_host()) call z%sync()
      call y%axpby(m,alpha,x,beta,info)
      call z%axpby(m,gamma,y,delta,info)
    end if

  end subroutine s_cuda_abgdxyz

  subroutine s_cuda_xyzw(m,a,b,c,d,e,f,x, y, z,w, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_s_base_vect_type), intent(inout)  :: x
    class(psb_s_base_vect_type), intent(inout)  :: y
    class(psb_s_base_vect_type), intent(inout)  :: z
    class(psb_s_vect_cuda), intent(inout)  :: w
    real(psb_spk_), intent (in)       :: a,b,c,d,e,f
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: nx, ny, nz, nw
    logical :: gpu_done

    info = psb_success_

    gpu_done = .false.
    if ((a==szero).or.(b==szero).or. &
         & (c==szero).or.(d==szero).or.&
         & (e==szero).or.(f==szero)) then
      write(0,*) 'XYZW assumes  a,b,c,d,e,f are all nonzero'
    else
      select type(xx => x)
      class is (psb_s_vect_cuda)
        select type(yy => y)
        class is (psb_s_vect_cuda)
          select type(zz => z)
          class is (psb_s_vect_cuda)
            ! Do something different here 
            if (xx%is_host()) call xx%sync()
            if (yy%is_host()) call yy%sync()
            if (zz%is_host()) call zz%sync()
            if (w%is_host())  call w%sync()
            nx = getMultiVecDeviceSize(xx%deviceVect)
            ny = getMultiVecDeviceSize(yy%deviceVect)
            nz = getMultiVecDeviceSize(zz%deviceVect)
            nw = getMultiVecDeviceSize(w%deviceVect)
            if ((nx<m).or.(ny<m).or.(nz<m).or.(nw<m)) then
              info = psb_err_internal_error_
            else
              info = xyzwMultiVecDevice(m,a,b,c,d,e,f,&
                   & xx%deviceVect,yy%deviceVect,zz%deviceVect,w%deviceVect)
            end if
            call yy%set_dev()
            call zz%set_dev()
            call w%set_dev()
            gpu_done = .true.
          end select
        end select
      end select
      
      if (.not.gpu_done) then 
        if (x%is_host()) call x%sync()
        if (y%is_host()) call y%sync()
        if (z%is_host()) call z%sync()
        if (w%is_host()) call w%sync()
        call y%axpby(m,a,x,b,info)
        call z%axpby(m,c,y,d,info)
        call w%axpby(m,e,z,f,info)
      end if
    end if

  end subroutine s_cuda_xyzw
  
  subroutine s_cuda_axpby_a(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)        :: m
    real(psb_spk_), intent(in)           :: x(:)
    class(psb_s_vect_cuda), intent(inout) :: y
    real(psb_spk_), intent (in)          :: alpha, beta
    integer(psb_ipk_), intent(out)       :: info
    
    if ((beta /= szero).and.(y%is_dev()))&
         & call y%sync()
    call psb_geaxpby(m,alpha,x,beta,y%v,info)
    call y%set_host()
  end subroutine s_cuda_axpby_a

  subroutine s_cuda_mlt_v(x, y, info)
    use psi_serial_mod
    implicit none 
    class(psb_s_base_vect_type), intent(inout) :: x
    class(psb_s_vect_cuda), intent(inout)       :: y
    integer(psb_ipk_), intent(out)             :: info

    integer(psb_ipk_) :: i, n
    
    info = 0    
    n = min(x%get_nrows(),y%get_nrows())
    select type(xx => x)
    type is (psb_s_base_vect_type)
      if (y%is_dev()) call y%sync()
      do i=1, n
        y%v(i) = y%v(i) * xx%v(i)
      end do
      call y%set_host()
    type is (psb_s_vect_cuda)
      ! Do something different here 
      if (y%is_host())  call y%sync()
      if (xx%is_host()) call xx%sync()
      info = axyMultiVecDevice(n,sone,xx%deviceVect,y%deviceVect)
      call y%set_dev()
    class default
      if (xx%is_dev()) call xx%sync()
      if (y%is_dev())  call y%sync()
      call y%mlt(xx%v,info)
      call y%set_host()
    end select

  end subroutine s_cuda_mlt_v

  subroutine s_cuda_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none 
    real(psb_spk_), intent(in)           :: x(:)
    class(psb_s_vect_cuda), intent(inout) :: y
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: i, n
    
    info = 0    
    if (y%is_dev()) call y%sync()
    call y%psb_s_base_vect_type%mlt(x,info)
    ! set_host() is invoked in the base method
  end subroutine s_cuda_mlt_a

  subroutine s_cuda_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    real(psb_spk_), intent(in)           :: alpha,beta
    real(psb_spk_), intent(in)           :: x(:)
    real(psb_spk_), intent(in)           :: y(:)
    class(psb_s_vect_cuda), intent(inout) :: z
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: i, n
    
    info = 0    
    if (z%is_dev()) call z%sync()
    call z%psb_s_base_vect_type%mlt(alpha,x,y,beta,info)
    ! set_host() is invoked in the base method
  end subroutine s_cuda_mlt_a_2

  subroutine s_cuda_mlt_v_2(alpha,x,y, beta,z,info,conjgx,conjgy)
    use psi_serial_mod
    use psb_string_mod
    implicit none 
    real(psb_spk_), intent(in)                 :: alpha,beta
    class(psb_s_base_vect_type), intent(inout) :: x
    class(psb_s_base_vect_type), intent(inout) :: y
    class(psb_s_vect_cuda), intent(inout)       :: z
    integer(psb_ipk_), intent(out)             :: info
    character(len=1), intent(in), optional     :: conjgx, conjgy
    integer(psb_ipk_) :: i, n
    logical :: conjgx_, conjgy_

    if (.false.) then 
      ! These are present just for coherence with the
      ! complex versions; they do nothing here. 
      conjgx_=.false.
      if (present(conjgx)) conjgx_ = (psb_toupper(conjgx)=='C')
      conjgy_=.false.
      if (present(conjgy)) conjgy_ = (psb_toupper(conjgy)=='C')
    end if
    
    n = min(x%get_nrows(),y%get_nrows(),z%get_nrows())
    
    !
    ! Need to reconsider BETA in the GPU side
    !  of things.
    !
    info = 0    
    select type(xx => x) 
    type is (psb_s_vect_cuda)
      select type (yy => y) 
      type is (psb_s_vect_cuda)
        if (xx%is_host()) call xx%sync()
        if (yy%is_host()) call yy%sync()
        if ((beta /= szero).and.(z%is_host())) call z%sync()
        info = axybzMultiVecDevice(n,alpha,xx%deviceVect,&
             & yy%deviceVect,beta,z%deviceVect)
        call z%set_dev()
      class default
        if (xx%is_dev()) call xx%sync()
        if (yy%is_dev()) call yy%sync()
        if ((beta /= szero).and.(z%is_dev())) call z%sync()
        call z%psb_s_base_vect_type%mlt(alpha,xx,yy,beta,info)
        call z%set_host()
      end select
      
    class default
      if (x%is_dev()) call x%sync()
      if (y%is_dev()) call y%sync()
      if ((beta /= szero).and.(z%is_dev())) call z%sync()
      call z%psb_s_base_vect_type%mlt(alpha,x,y,beta,info)
      call z%set_host()
    end select
  end subroutine s_cuda_mlt_v_2

  subroutine s_cuda_scal(alpha, x)
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    real(psb_spk_), intent (in)          :: alpha
    integer(psb_ipk_) :: info
    
    if (x%is_host()) call x%sync()
    info = scalMultiVecDevice(alpha,x%deviceVect)
    call x%set_dev()
  end subroutine s_cuda_scal


  function s_cuda_nrm2(n,x) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                       :: res
    integer(psb_ipk_) :: info
    ! WARNING: this should be changed. 
    if (x%is_host()) call x%sync()
    info = nrm2MultiVecDevice(res,n,x%deviceVect)
    
  end function s_cuda_nrm2
  
  function s_cuda_amax(n,x) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                :: res
    integer(psb_ipk_) :: info

    if (x%is_host()) call x%sync()
    info = amaxMultiVecDevice(res,n,x%deviceVect)

  end function s_cuda_amax

  function s_cuda_asum(n,x) result(res)
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n
    real(psb_spk_)                :: res
    integer(psb_ipk_) :: info

    if (x%is_host()) call x%sync()
    info = asumMultiVecDevice(res,n,x%deviceVect)

  end function s_cuda_asum

  subroutine  s_cuda_absval1(x) 
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: info

    if (x%is_host()) call x%sync()
    n=x%get_nrows()
    info = absMultiVecDevice(n,sone,x%deviceVect)

  end subroutine s_cuda_absval1

  subroutine  s_cuda_absval2(x,y) 
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    class(psb_s_base_vect_type), intent(inout) :: y
    integer(psb_ipk_) :: n
    integer(psb_ipk_) :: info

    n=min(x%get_nrows(),y%get_nrows())
    select type (yy=> y) 
    class is (psb_s_vect_cuda)
      if (x%is_host()) call x%sync()    
      if (yy%is_host()) call yy%sync()
      info = absMultiVecDevice(n,sone,x%deviceVect,yy%deviceVect)
    class default
      if (x%is_dev()) call x%sync()    
      if (y%is_dev()) call y%sync()
      call x%psb_s_base_vect_type%absval(y)
    end select
  end subroutine s_cuda_absval2


  subroutine s_cuda_vect_finalize(x)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    type(psb_s_vect_cuda), intent(inout)  :: x
    integer(psb_ipk_)        :: info
    
    info = 0
    call x%free(info)
  end subroutine s_cuda_vect_finalize

  subroutine s_cuda_ins_v(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_s_vect_cuda), intent(inout)        :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    class(psb_i_base_vect_type), intent(inout)  :: irl
    class(psb_s_base_vect_type), intent(inout)  :: val
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, isz
    logical :: done_cuda

    info = 0
    if (psb_errstatus_fatal()) return 

    done_cuda = .false. 
    select type(virl => irl)
    class is (psb_i_vect_cuda) 
      select type(vval => val)
      class is (psb_s_vect_cuda) 
        if (vval%is_host()) call vval%sync()
        if (virl%is_host()) call virl%sync()
        if (x%is_host())    call x%sync()
        info = geinsMultiVecDeviceFloat(n,virl%deviceVect,&
             & vval%deviceVect,dupl,1,x%deviceVect)
        call x%set_dev()
        done_cuda=.true.
      end select
    end select

    if (.not.done_cuda) then 
      if (irl%is_dev()) call irl%sync()
      if (val%is_dev()) call val%sync()
      call x%ins(n,irl%v,val%v,dupl,info)
    end if

    if (info /= 0) then 
      call psb_errpush(info,'cuda_vect_ins')
      return
    end if

  end subroutine s_cuda_ins_v
  
  subroutine s_cuda_ins_a(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_s_vect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n, dupl
    integer(psb_ipk_), intent(in)        :: irl(:)
    real(psb_spk_), intent(in)           :: val(:)
    integer(psb_ipk_), intent(out)       :: info

    integer(psb_ipk_) :: i

    info = 0
    if (x%is_dev()) call x%sync()
    call x%psb_s_base_vect_type%ins(n,irl,val,dupl,info)
    call x%set_host()

  end subroutine s_cuda_ins_a

end module psb_s_cuda_vect_mod


!
! Multivectors
! 



module psb_s_cuda_multivect_mod
  use iso_c_binding
  use psb_const_mod
  use psb_error_mod
  use psb_s_multivect_mod
  use psb_s_base_multivect_mod
  use psb_cuda_env_mod
  use psb_i_multivect_mod
  use psb_i_cuda_multivect_mod
  use psb_s_vectordev_mod

  integer(psb_ipk_), parameter, private :: is_host = -1
  integer(psb_ipk_), parameter, private :: is_sync = 0 
  integer(psb_ipk_), parameter, private :: is_dev  = 1 
  
  type, extends(psb_s_base_multivect_type) ::  psb_s_multivect_cuda
    integer(psb_ipk_)  :: state      = is_host, m_nrows=0, m_ncols=0
    type(c_ptr) :: deviceVect = c_null_ptr
    real(c_double), allocatable :: buffer(:,:)
    type(c_ptr) :: dt_buf = c_null_ptr
  contains
    procedure, pass(x) :: get_nrows => s_cuda_multi_get_nrows
    procedure, pass(x) :: get_ncols => s_cuda_multi_get_ncols
    procedure, nopass  :: get_fmt   => s_cuda_multi_get_fmt
!!$    procedure, pass(x) :: dot_v    => s_cuda_multi_dot_v
!!$    procedure, pass(x) :: dot_a    => s_cuda_multi_dot_a
!!$    procedure, pass(y) :: axpby_v  => s_cuda_multi_axpby_v
!!$    procedure, pass(y) :: axpby_a  => s_cuda_multi_axpby_a
!!$    procedure, pass(y) :: mlt_v    => s_cuda_multi_mlt_v
!!$    procedure, pass(y) :: mlt_a    => s_cuda_multi_mlt_a
!!$    procedure, pass(z) :: mlt_a_2  => s_cuda_multi_mlt_a_2
!!$    procedure, pass(z) :: mlt_v_2  => s_cuda_multi_mlt_v_2
!!$    procedure, pass(x) :: scal     => s_cuda_multi_scal
!!$    procedure, pass(x) :: nrm2     => s_cuda_multi_nrm2
!!$    procedure, pass(x) :: amax     => s_cuda_multi_amax
!!$    procedure, pass(x) :: asum     => s_cuda_multi_asum
    procedure, pass(x) :: all      => s_cuda_multi_all
    procedure, pass(x) :: zero     => s_cuda_multi_zero
    procedure, pass(x) :: asb      => s_cuda_multi_asb
    procedure, pass(x) :: sync     => s_cuda_multi_sync
    procedure, pass(x) :: sync_space => s_cuda_multi_sync_space
    procedure, pass(x) :: bld_x    => s_cuda_multi_bld_x
    procedure, pass(x) :: bld_n    => s_cuda_multi_bld_n
    procedure, pass(x) :: free     => s_cuda_multi_free
    procedure, pass(x) :: ins      => s_cuda_multi_ins
    procedure, pass(x) :: is_host  => s_cuda_multi_is_host
    procedure, pass(x) :: is_dev   => s_cuda_multi_is_dev
    procedure, pass(x) :: is_sync  => s_cuda_multi_is_sync
    procedure, pass(x) :: set_host => s_cuda_multi_set_host
    procedure, pass(x) :: set_dev  => s_cuda_multi_set_dev
    procedure, pass(x) :: set_sync => s_cuda_multi_set_sync
    procedure, pass(x) :: set_scal => s_cuda_multi_set_scal
    procedure, pass(x) :: set_vect => s_cuda_multi_set_vect
!!$    procedure, pass(x) :: gthzv_x  => s_cuda_multi_gthzv_x
!!$    procedure, pass(y) :: sctb     => s_cuda_multi_sctb
!!$    procedure, pass(y) :: sctb_x   => s_cuda_multi_sctb_x
    final              :: s_cuda_multi_vect_finalize
  end type psb_s_multivect_cuda

  public  :: psb_s_multivect_cuda
  private :: constructor
  interface psb_s_multivect_cuda
    module procedure constructor
  end interface

contains
  
  function constructor(x) result(this)
    real(psb_spk_)       :: x(:,:)
    type(psb_s_multivect_cuda) :: this
    integer(psb_ipk_) :: info

    this%v = x
    call this%asb(size(x,1),size(x,2),info)

  end function constructor
    

!!$  subroutine s_cuda_multi_gthzv_x(i,n,idx,x,y)
!!$    use psi_serial_mod
!!$    integer(psb_ipk_) :: i,n
!!$    class(psb_i_base_multivect_type) :: idx
!!$    real(psb_spk_) ::  y(:)
!!$    class(psb_s_multivect_cuda) :: x
!!$
!!$    select type(ii=> idx) 
!!$    class is (psb_i_vect_cuda) 
!!$      if (ii%is_host()) call ii%sync()
!!$      if (x%is_host())  call x%sync()
!!$
!!$      if (allocated(x%buffer)) then 
!!$        if (size(x%buffer) < n) then 
!!$          call inner_unregister(x%buffer)
!!$          deallocate(x%buffer, stat=info)
!!$        end if
!!$      end if
!!$      
!!$      if (.not.allocated(x%buffer)) then
!!$        allocate(x%buffer(n),stat=info)
!!$        if (info == 0) info = inner_register(x%buffer,x%dt_buf)        
!!$      endif
!!$      info = igathMultiVecDeviceDouble(x%deviceVect,&
!!$           & 0, i, n, ii%deviceVect, x%dt_buf, 1)
!!$      call psb_cudaSync()
!!$      y(1:n) = x%buffer(1:n)
!!$      
!!$    class default
!!$      call x%gth(n,ii%v(i:),y)
!!$    end select
!!$
!!$
!!$  end subroutine s_cuda_multi_gthzv_x
!!$
!!$
!!$
!!$  subroutine s_cuda_multi_sctb(n,idx,x,beta,y)
!!$    implicit none
!!$    !use psb_const_mod
!!$    integer(psb_ipk_)     :: n, idx(:)
!!$    real(psb_spk_)        :: beta, x(:)
!!$    class(psb_s_multivect_cuda) :: y
!!$    integer(psb_ipk_)     :: info
!!$
!!$    if (n == 0) return
!!$    
!!$    if (y%is_dev())  call y%sync()
!!$          
!!$    call y%psb_s_base_multivect_type%sctb(n,idx,x,beta)
!!$    call y%set_host()
!!$
!!$  end subroutine s_cuda_multi_sctb
!!$
!!$  subroutine s_cuda_multi_sctb_x(i,n,idx,x,beta,y)
!!$    use psi_serial_mod
!!$    integer(psb_ipk_) :: i, n
!!$    class(psb_i_base_multivect_type) :: idx
!!$    real(psb_spk_) :: beta, x(:)
!!$    class(psb_s_multivect_cuda) :: y
!!$
!!$    select type(ii=> idx) 
!!$    class is (psb_i_vect_cuda) 
!!$      if (ii%is_host()) call ii%sync()
!!$      if (y%is_host())  call y%sync()
!!$
!!$      if (allocated(y%buffer)) then 
!!$        if (size(y%buffer) < n) then 
!!$          call inner_unregister(y%buffer)
!!$          deallocate(y%buffer, stat=info)
!!$        end if
!!$      end if
!!$      
!!$      if (.not.allocated(y%buffer)) then
!!$        allocate(y%buffer(n),stat=info)
!!$        if (info == 0) info = inner_register(y%buffer,y%dt_buf)        
!!$      endif
!!$      y%buffer(1:n) = x(1:n) 
!!$      info = iscatMultiVecDeviceDouble(y%deviceVect,&
!!$           & 0, i, n, ii%deviceVect, y%dt_buf, 1,beta)
!!$
!!$      call y%set_dev()
!!$      call psb_cudaSync()   
!!$      
!!$    class default
!!$      call y%sct(n,ii%v(i:),x,beta)
!!$    end select
!!$
!!$  end subroutine s_cuda_multi_sctb_x


  subroutine s_cuda_multi_bld_x(x,this)
    use psb_base_mod
    real(psb_spk_), intent(in)           :: this(:,:)
    class(psb_s_multivect_cuda), intent(inout) :: x
    integer(psb_ipk_) :: info, m, n
    
    m=size(this,1)
    n=size(this,2)
    x%m_nrows = m
    x%m_ncols = n
    call psb_realloc(m,n,x%v,info)
    if (info /= 0) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,'s_cuda_multi_bld_x',&
           & i_err=(/size(this,1),size(this,2),izero,izero,izero,izero/))
    end if
    x%v(1:m,1:n)  = this(1:m,1:n) 
    call x%set_host()
    call x%sync()

  end subroutine s_cuda_multi_bld_x

  subroutine s_cuda_multi_bld_n(x,m,n)
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_s_multivect_cuda), intent(inout) :: x
    integer(psb_ipk_) :: info

    call x%all(m,n,info)
    if (info /= 0) then 
      call psb_errpush(info,'s_cuda_multi_bld_n',i_err=(/m,n,n,n,n/))
    end if
    
  end subroutine s_cuda_multi_bld_n


  subroutine s_cuda_multi_set_host(x)
    implicit none 
    class(psb_s_multivect_cuda), intent(inout) :: x
    
    x%state = is_host
  end subroutine s_cuda_multi_set_host

  subroutine s_cuda_multi_set_dev(x)
    implicit none 
    class(psb_s_multivect_cuda), intent(inout) :: x
    
    x%state = is_dev
  end subroutine s_cuda_multi_set_dev

  subroutine s_cuda_multi_set_sync(x)
    implicit none 
    class(psb_s_multivect_cuda), intent(inout) :: x
    
    x%state = is_sync
  end subroutine s_cuda_multi_set_sync

  function s_cuda_multi_is_dev(x) result(res)
    implicit none 
    class(psb_s_multivect_cuda), intent(in) :: x
    logical  :: res
  
    res = (x%state == is_dev)
  end function s_cuda_multi_is_dev
  
  function s_cuda_multi_is_host(x) result(res)
    implicit none 
    class(psb_s_multivect_cuda), intent(in) :: x
    logical  :: res

    res = (x%state == is_host)
  end function s_cuda_multi_is_host

  function s_cuda_multi_is_sync(x) result(res)
    implicit none 
    class(psb_s_multivect_cuda), intent(in) :: x
    logical  :: res

    res = (x%state == is_sync)
  end function s_cuda_multi_is_sync

  
  function s_cuda_multi_get_nrows(x) result(res)
    implicit none 
    class(psb_s_multivect_cuda), intent(in) :: x
    integer(psb_ipk_) :: res

    res = x%m_nrows

  end function s_cuda_multi_get_nrows
  
  function s_cuda_multi_get_ncols(x) result(res)
    implicit none 
    class(psb_s_multivect_cuda), intent(in) :: x
    integer(psb_ipk_) :: res

    res = x%m_ncols

  end function s_cuda_multi_get_ncols

  function s_cuda_multi_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'sGPU'
  end function s_cuda_multi_get_fmt

!!$  function s_cuda_multi_dot_v(n,x,y) result(res)
!!$    implicit none 
!!$    class(psb_s_multivect_cuda), intent(inout)       :: x
!!$    class(psb_s_base_multivect_type), intent(inout) :: y
!!$    integer(psb_ipk_), intent(in)              :: n
!!$    real(psb_spk_)                :: res
!!$    real(psb_spk_), external      :: ddot
!!$    integer(psb_ipk_) :: info
!!$    
!!$    res = szero
!!$    !
!!$    ! Note: this is the gpu implementation.
!!$    !  When we get here, we are sure that X is of
!!$    !  TYPE psb_s_vect
!!$    !
!!$    select type(yy => y)
!!$    type is (psb_s_base_multivect_type)
!!$      if (x%is_dev()) call x%sync()
!!$      res = ddot(n,x%v,1,yy%v,1)
!!$    type is (psb_s_multivect_cuda)
!!$      if (x%is_host()) call x%sync()
!!$      if (yy%is_host()) call yy%sync()
!!$      info = dotMultiVecDevice(res,n,x%deviceVect,yy%deviceVect)
!!$      if (info /= 0) then 
!!$        info = psb_err_internal_error_
!!$        call psb_errpush(info,'s_cuda_multi_dot_v')
!!$      end if
!!$
!!$    class default
!!$      ! y%sync is done in dot_a
!!$      call x%sync()      
!!$      res = y%dot(n,x%v)
!!$    end select
!!$
!!$  end function s_cuda_multi_dot_v
!!$
!!$  function s_cuda_multi_dot_a(n,x,y) result(res)
!!$    implicit none 
!!$    class(psb_s_multivect_cuda), intent(inout) :: x
!!$    real(psb_spk_), intent(in)           :: y(:)
!!$    integer(psb_ipk_), intent(in)        :: n
!!$    real(psb_spk_)                :: res
!!$    real(psb_spk_), external      :: ddot
!!$    
!!$    if (x%is_dev()) call x%sync()
!!$    res = ddot(n,y,1,x%v,1)
!!$
!!$  end function s_cuda_multi_dot_a
!!$    
!!$  subroutine s_cuda_multi_axpby_v(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    integer(psb_ipk_), intent(in)              :: m
!!$    class(psb_s_base_multivect_type), intent(inout) :: x
!!$    class(psb_s_multivect_cuda), intent(inout)       :: y
!!$    real(psb_spk_), intent (in)                :: alpha, beta
!!$    integer(psb_ipk_), intent(out)             :: info
!!$    integer(psb_ipk_) :: nx, ny
!!$
!!$    info = psb_success_
!!$
!!$    select type(xx => x)
!!$    type is (psb_s_base_multivect_type)
!!$      if ((beta /= szero).and.(y%is_dev()))&
!!$           & call y%sync()
!!$      call psb_geaxpby(m,alpha,xx%v,beta,y%v,info)
!!$      call y%set_host()
!!$    type is (psb_s_multivect_cuda)
!!$      ! Do something different here 
!!$      if ((beta /= szero).and.y%is_host())&
!!$           &  call y%sync()
!!$      if (xx%is_host()) call xx%sync()
!!$      nx = getMultiVecDeviceSize(xx%deviceVect)
!!$      ny = getMultiVecDeviceSize(y%deviceVect)
!!$      if ((nx<m).or.(ny<m)) then
!!$        info = psb_err_internal_error_
!!$        info = psb_err_internal_error_
!!$      else
!!$        info = axpbyMultiVecDevice(m,alpha,xx%deviceVect,beta,y%deviceVect)
!!$      end if
!!$      call y%set_dev()
!!$    class default
!!$      call x%sync()
!!$      call y%axpby(m,alpha,x%v,beta,info)
!!$    end select
!!$
!!$  end subroutine s_cuda_multi_axpby_v
!!$
!!$  subroutine s_cuda_multi_axpby_a(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    integer(psb_ipk_), intent(in)        :: m
!!$    real(psb_spk_), intent(in)           :: x(:)
!!$    class(psb_s_multivect_cuda), intent(inout) :: y
!!$    real(psb_spk_), intent (in)          :: alpha, beta
!!$    integer(psb_ipk_), intent(out)       :: info
!!$
!!$    if (y%is_dev()) call y%sync()
!!$    call psb_geaxpby(m,alpha,x,beta,y%v,info)
!!$    call y%set_host()
!!$  end subroutine s_cuda_multi_axpby_a
!!$
!!$  subroutine s_cuda_multi_mlt_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    class(psb_s_base_multivect_type), intent(inout) :: x
!!$    class(psb_s_multivect_cuda), intent(inout)       :: y
!!$    integer(psb_ipk_), intent(out)             :: info
!!$
!!$    integer(psb_ipk_) :: i, n
!!$    
!!$    info = 0    
!!$    n = min(x%get_nrows(),y%get_nrows())
!!$    select type(xx => x)
!!$    type is (psb_s_base_multivect_type)
!!$      if (y%is_dev()) call y%sync()
!!$      do i=1, n
!!$        y%v(i) = y%v(i) * xx%v(i)
!!$      end do
!!$      call y%set_host()
!!$    type is (psb_s_multivect_cuda)
!!$      ! Do something different here 
!!$      if (y%is_host())  call y%sync()
!!$      if (xx%is_host()) call xx%sync()
!!$      info = axyMultiVecDevice(n,done,xx%deviceVect,y%deviceVect)
!!$      call y%set_dev()
!!$    class default
!!$      call xx%sync()
!!$      call y%mlt(xx%v,info)
!!$      call y%set_host()
!!$    end select
!!$
!!$  end subroutine s_cuda_multi_mlt_v
!!$
!!$  subroutine s_cuda_multi_mlt_a(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    real(psb_spk_), intent(in)           :: x(:)
!!$    class(psb_s_multivect_cuda), intent(inout) :: y
!!$    integer(psb_ipk_), intent(out)       :: info
!!$    integer(psb_ipk_) :: i, n
!!$    
!!$    info = 0    
!!$    call y%sync()
!!$    call y%psb_s_base_multivect_type%mlt(x,info)
!!$    call y%set_host()
!!$  end subroutine s_cuda_multi_mlt_a
!!$
!!$  subroutine s_cuda_multi_mlt_a_2(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none 
!!$    real(psb_spk_), intent(in)           :: alpha,beta
!!$    real(psb_spk_), intent(in)           :: x(:)
!!$    real(psb_spk_), intent(in)           :: y(:)
!!$    class(psb_s_multivect_cuda), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)       :: info
!!$    integer(psb_ipk_) :: i, n
!!$    
!!$    info = 0    
!!$    if (z%is_dev()) call z%sync()
!!$    call z%psb_s_base_multivect_type%mlt(alpha,x,y,beta,info)
!!$    call z%set_host()
!!$  end subroutine s_cuda_multi_mlt_a_2
!!$
!!$  subroutine s_cuda_multi_mlt_v_2(alpha,x,y, beta,z,info,conjgx,conjgy)
!!$    use psi_serial_mod
!!$    use psb_string_mod
!!$    implicit none 
!!$    real(psb_spk_), intent(in)                 :: alpha,beta
!!$    class(psb_s_base_multivect_type), intent(inout) :: x
!!$    class(psb_s_base_multivect_type), intent(inout) :: y
!!$    class(psb_s_multivect_cuda), intent(inout)       :: z
!!$    integer(psb_ipk_), intent(out)             :: info
!!$    character(len=1), intent(in), optional     :: conjgx, conjgy
!!$    integer(psb_ipk_) :: i, n
!!$    logical :: conjgx_, conjgy_
!!$
!!$    if (.false.) then 
!!$      ! These are present just for coherence with the
!!$      ! complex versions; they do nothing here. 
!!$      conjgx_=.false.
!!$      if (present(conjgx)) conjgx_ = (psb_toupper(conjgx)=='C')
!!$      conjgy_=.false.
!!$      if (present(conjgy)) conjgy_ = (psb_toupper(conjgy)=='C')
!!$    end if
!!$    
!!$    n = min(x%get_nrows(),y%get_nrows(),z%get_nrows())
!!$    
!!$    !
!!$    ! Need to reconsider BETA in the GPU side
!!$    !  of things.
!!$    !
!!$    info = 0    
!!$    select type(xx => x) 
!!$    type is (psb_s_multivect_cuda)
!!$      select type (yy => y) 
!!$      type is (psb_s_multivect_cuda)
!!$        if (xx%is_host()) call xx%sync()
!!$        if (yy%is_host()) call yy%sync()
!!$        ! Z state is irrelevant: it will be done on the GPU. 
!!$        info = axybzMultiVecDevice(n,alpha,xx%deviceVect,&
!!$             & yy%deviceVect,beta,z%deviceVect)
!!$        call z%set_dev()
!!$      class default
!!$        call xx%sync()
!!$        call yy%sync()
!!$        call z%psb_s_base_multivect_type%mlt(alpha,xx,yy,beta,info)
!!$        call z%set_host()
!!$      end select
!!$      
!!$    class default
!!$      call x%sync()
!!$      call y%sync()
!!$      call z%psb_s_base_multivect_type%mlt(alpha,x,y,beta,info)
!!$      call z%set_host()
!!$    end select
!!$  end subroutine s_cuda_multi_mlt_v_2


  subroutine s_cuda_multi_set_scal(x,val)
    class(psb_s_multivect_cuda), intent(inout) :: x
    real(psb_spk_), intent(in)           :: val
        
    integer(psb_ipk_) :: info

    if (x%is_dev()) call x%sync()
    call x%psb_s_base_multivect_type%set_scal(val)
    call x%set_host()
  end subroutine s_cuda_multi_set_scal

  subroutine s_cuda_multi_set_vect(x,val)
    class(psb_s_multivect_cuda), intent(inout) :: x
    real(psb_spk_), intent(in)           :: val(:,:)
    integer(psb_ipk_) :: nr
    integer(psb_ipk_) :: info

    if (x%is_dev()) call x%sync()
    call x%psb_s_base_multivect_type%set_vect(val)
    call x%set_host()

  end subroutine s_cuda_multi_set_vect



!!$  subroutine s_cuda_multi_scal(alpha, x)
!!$    implicit none 
!!$    class(psb_s_multivect_cuda), intent(inout) :: x
!!$    real(psb_spk_), intent (in)          :: alpha
!!$    
!!$    if (x%is_dev()) call x%sync()
!!$    call x%psb_s_base_multivect_type%scal(alpha)
!!$    call x%set_host()
!!$  end subroutine s_cuda_multi_scal
!!$
!!$
!!$  function s_cuda_multi_nrm2(n,x) result(res)
!!$    implicit none 
!!$    class(psb_s_multivect_cuda), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)        :: n
!!$    real(psb_spk_)                       :: res
!!$    integer(psb_ipk_) :: info
!!$    ! WARNING: this should be changed. 
!!$    if (x%is_host()) call x%sync()
!!$    info = nrm2MultiVecDevice(res,n,x%deviceVect)
!!$    
!!$  end function s_cuda_multi_nrm2
!!$  
!!$  function s_cuda_multi_amax(n,x) result(res)
!!$    implicit none 
!!$    class(psb_s_multivect_cuda), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)        :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (x%is_dev()) call x%sync()
!!$    res =  maxval(abs(x%v(1:n)))
!!$
!!$  end function s_cuda_multi_amax
!!$
!!$  function s_cuda_multi_asum(n,x) result(res)
!!$    implicit none 
!!$    class(psb_s_multivect_cuda), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)        :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (x%is_dev()) call x%sync()
!!$    res =  sum(abs(x%v(1:n)))
!!$
!!$  end function s_cuda_multi_asum
  
  subroutine s_cuda_multi_all(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)      :: m,n
    class(psb_s_multivect_cuda), intent(out) :: x
    integer(psb_ipk_), intent(out)     :: info
    
    call psb_realloc(m,n,x%v,info,pad=szero)
    x%m_nrows = m
    x%m_ncols = n
    if (info == 0) call x%set_host()
    if (info == 0) call x%sync_space(info)
    if (info /= 0) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,'s_cuda_multi_all',&
           & i_err=(/m,n,n,n,n/))
    end if
  end subroutine s_cuda_multi_all

  subroutine s_cuda_multi_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_s_multivect_cuda), intent(inout) :: x
    
    if (allocated(x%v)) x%v=szero
    call x%set_host()
  end subroutine s_cuda_multi_zero

  subroutine s_cuda_multi_asb(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)        :: m,n
    class(psb_s_multivect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(out)       :: info
    integer(psb_ipk_) :: nd, nc


    x%m_nrows = m
    x%m_ncols = n
    if (x%is_host()) then 
      call x%psb_s_base_multivect_type%asb(m,n,info)
      if (info == psb_success_) call x%sync_space(info)
    else if (x%is_dev()) then 
      nd  = getMultiVecDevicePitch(x%deviceVect)
      nc  = getMultiVecDeviceCount(x%deviceVect)
      if ((nd < m).or.(nc<n)) then 
        call x%sync()
        call x%psb_s_base_multivect_type%asb(m,n,info)      
        if (info == psb_success_) call x%sync_space(info)
        call x%set_host()
      end if
    end if
  end subroutine s_cuda_multi_asb

  subroutine s_cuda_multi_sync_space(x,info)
    use psb_realloc_mod
    implicit none 
    class(psb_s_multivect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(out)       :: info 
    integer(psb_ipk_) :: mh,nh,md,nd
    
    info = 0
    if (x%is_host()) then 
      if (allocated(x%v)) then 
        mh = size(x%v,1)
        nh = size(x%v,2)
      else
        mh=0
        nh=0
      end if
      if (c_associated(x%deviceVect)) then 
        md  = getMultiVecDevicePitch(x%deviceVect)
        nd  = getMultiVecDeviceCount(x%deviceVect)
        if ((md < mh).or.(nd<nh)) then 
          call freeMultiVecDevice(x%deviceVect)
          x%deviceVect=c_null_ptr
        end if
      end if

      if (.not.c_associated(x%deviceVect)) then 
        info = FallocMultiVecDevice(x%deviceVect,nh,mh,spgpu_type_float)
        if (info == 0) &
             & call psb_realloc(getMultiVecDevicePitch(x%deviceVect),&
             & getMultiVecDeviceCount(x%deviceVect),x%v,info,pad=szero)
        if  (info /= 0) then 
!!$          write(0,*) 'Error from FallocMultiVecDevice',info,n
          if (info == spgpu_outofmem) then 
            info = psb_err_alloc_request_
          end if
        end if
        
      end if
    else if (x%is_dev()) then 
      ! 
      if (allocated(x%v)) then 
        mh = size(x%v,1)
        nh = size(x%v,2)
      else
        mh=0
        nh=0
      end if
      md  = getMultiVecDevicePitch(x%deviceVect)
      nd  = getMultiVecDeviceCount(x%deviceVect)
      if ((mh /= md).or.(nh /= nd)) then 
        call psb_realloc(getMultiVecDevicePitch(x%deviceVect),&
             & getMultiVecDeviceCount(x%deviceVect),x%v,info,pad=szero)
      end if
      
    end if
    
  end subroutine s_cuda_multi_sync_space

  subroutine s_cuda_multi_sync(x)
    implicit none 
    class(psb_s_multivect_cuda), intent(inout) :: x
    integer(psb_ipk_) :: n,info
    
    info = 0
    if (x%is_host()) then 
      if (.not.c_associated(x%deviceVect)) then 
        call x%sync_space(info)
      end if
      if (info == 0) &
           & info = writeMultiVecDevice(x%deviceVect,x%v,size(x%v,1))
    else if (x%is_dev()) then 
      info = readMultiVecDevice(x%deviceVect,x%v,size(x%v,1))
    end if
    if (info == 0)  call x%set_sync()
    if (info /= 0) then
      info=psb_err_internal_error_
      call psb_errpush(info,'s_cuda_multi_sync')
    end if
    
  end subroutine s_cuda_multi_sync

  subroutine s_cuda_multi_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_s_multivect_cuda), intent(inout)  :: x
    integer(psb_ipk_), intent(out)        :: info
    
    info = 0
    if (c_associated(x%deviceVect)) then 
      call freeMultiVecDevice(x%deviceVect)
      x%deviceVect=c_null_ptr
    end if
    if (allocated(x%buffer)) then 
!!$      call inner_unregister(x%buffer)
      deallocate(x%buffer, stat=info)
    end if

    if (allocated(x%v)) deallocate(x%v, stat=info)
    call x%set_sync()
  end subroutine s_cuda_multi_free

  subroutine s_cuda_multi_vect_finalize(x)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    type(psb_s_multivect_cuda), intent(inout)  :: x
    integer(psb_ipk_)        :: info
    
    info = 0
    if (c_associated(x%deviceVect)) then 
      call freeMultiVecDevice(x%deviceVect)
      x%deviceVect=c_null_ptr
    end if
    if (allocated(x%buffer)) then 
!!$      call inner_unregister(x%buffer)
      deallocate(x%buffer, stat=info)
    end if

    if (allocated(x%v)) deallocate(x%v, stat=info)
    call x%set_sync()
  end subroutine s_cuda_multi_vect_finalize

  subroutine s_cuda_multi_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_s_multivect_cuda), intent(inout) :: x
    integer(psb_ipk_), intent(in)        :: n, dupl
    integer(psb_ipk_), intent(in)        :: irl(:)
    real(psb_spk_), intent(in)           :: val(:,:)
    integer(psb_ipk_), intent(out)       :: info

    integer(psb_ipk_) :: i

    info = 0
    if (x%is_dev()) call x%sync()
    call x%psb_s_base_multivect_type%ins(n,irl,val,dupl,info)
    call x%set_host()

  end subroutine s_cuda_multi_ins

end module psb_s_cuda_multivect_mod



