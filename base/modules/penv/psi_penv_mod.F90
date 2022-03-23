!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
!    
#if defined(SERIAL_MPI)
! Provide a fake mpi module just to keep the compiler(s) happy.
module mpi
  use psb_const_mod
  integer(psb_mpk_), parameter :: mpi_success          = 0
  integer(psb_mpk_), parameter :: mpi_request_null     = 0
  integer(psb_mpk_), parameter :: mpi_status_size      = 1
  integer(psb_mpk_), parameter :: mpi_integer          = 1
  integer(psb_mpk_), parameter :: mpi_integer8         = 2
  integer(psb_mpk_), parameter :: mpi_real             = 3
  integer(psb_mpk_), parameter :: mpi_double_precision = 4
  integer(psb_mpk_), parameter :: mpi_complex          = 5   
  integer(psb_mpk_), parameter :: mpi_double_complex   = 6 
  integer(psb_mpk_), parameter :: mpi_character        = 7
  integer(psb_mpk_), parameter :: mpi_logical          = 8
  integer(psb_mpk_), parameter :: mpi_integer2         = 9
  integer(psb_mpk_), parameter :: mpi_integer4         = 10
  integer(psb_mpk_), parameter :: mpi_comm_null        = -1
  integer(psb_mpk_), parameter :: mpi_comm_world       = 1
  
  real(psb_dpk_), external :: mpi_wtime
end module mpi
#endif    


module psi_penv_mod
  use psb_const_mod

  integer(psb_mpk_), parameter:: psb_int_tag      = 543987
  integer(psb_mpk_), parameter:: psb_real_tag     = psb_int_tag      + 1
  integer(psb_mpk_), parameter:: psb_double_tag   = psb_real_tag     + 1
  integer(psb_mpk_), parameter:: psb_complex_tag  = psb_double_tag   + 1
  integer(psb_mpk_), parameter:: psb_dcomplex_tag = psb_complex_tag  + 1
  integer(psb_mpk_), parameter:: psb_logical_tag  = psb_dcomplex_tag + 1
  integer(psb_mpk_), parameter:: psb_char_tag     = psb_logical_tag  + 1
  integer(psb_mpk_), parameter:: psb_int8_tag     = psb_char_tag     + 1
  integer(psb_mpk_), parameter:: psb_int2_tag     = psb_int8_tag     + 1
  integer(psb_mpk_), parameter:: psb_int4_tag     = psb_int2_tag     + 1
  integer(psb_mpk_), parameter:: psb_long_tag     = psb_int4_tag     + 1

  integer(psb_mpk_), parameter:: psb_int_swap_tag      = psb_int_tag      + psb_int_tag
  integer(psb_mpk_), parameter:: psb_real_swap_tag     = psb_real_tag     + psb_int_tag
  integer(psb_mpk_), parameter:: psb_double_swap_tag   = psb_double_tag   + psb_int_tag
  integer(psb_mpk_), parameter:: psb_complex_swap_tag  = psb_complex_tag  + psb_int_tag
  integer(psb_mpk_), parameter:: psb_dcomplex_swap_tag = psb_dcomplex_tag + psb_int_tag
  integer(psb_mpk_), parameter:: psb_logical_swap_tag  = psb_logical_tag  + psb_int_tag
  integer(psb_mpk_), parameter:: psb_char_swap_tag     = psb_char_tag     + psb_int_tag
  integer(psb_mpk_), parameter:: psb_int8_swap_tag     = psb_int8_tag     + psb_int_tag
  integer(psb_mpk_), parameter:: psb_int2_swap_tag     = psb_int2_tag     + psb_int_tag
  integer(psb_mpk_), parameter:: psb_int4_swap_tag     = psb_int4_tag     + psb_int_tag
  integer(psb_mpk_), parameter:: psb_long_swap_tag     = psb_long_tag     + psb_int_tag


  
  integer(psb_mpk_), private, parameter:: psb_int_type      = 987543
  integer(psb_mpk_), private, parameter:: psb_real_type     = psb_int_type      + 1
  integer(psb_mpk_), private, parameter:: psb_double_type   = psb_real_type     + 1
  integer(psb_mpk_), private, parameter:: psb_complex_type  = psb_double_type   + 1
  integer(psb_mpk_), private, parameter:: psb_dcomplex_type = psb_complex_type  + 1
  integer(psb_mpk_), private, parameter:: psb_logical_type  = psb_dcomplex_type + 1
  integer(psb_mpk_), private, parameter:: psb_char_type     = psb_logical_type  + 1
  integer(psb_mpk_), private, parameter:: psb_int8_type     = psb_char_type     + 1
  integer(psb_mpk_), private, parameter:: psb_int2_type     = psb_int8_type     + 1
  integer(psb_mpk_), private, parameter:: psb_int4_type     = psb_int2_type     + 1
  integer(psb_mpk_), private, parameter:: psb_long_type     = psb_int4_type     + 1

  type psb_buffer_node
    integer(psb_mpk_)   :: request
    type(psb_ctxt_type) :: ctxt 
    integer(psb_mpk_) :: buffer_type
    integer(psb_epk_), allocatable     :: int8buf(:)
    integer(psb_i2pk_), allocatable    :: int2buf(:)
    integer(psb_mpk_), allocatable     :: int4buf(:)
    real(psb_spk_), allocatable        :: realbuf(:)
    real(psb_dpk_), allocatable        :: doublebuf(:)
    complex(psb_spk_), allocatable     :: complexbuf(:)
    complex(psb_dpk_), allocatable     :: dcomplbuf(:)
    logical, allocatable               :: logbuf(:)
    character(len=1), allocatable      :: charbuf(:)
    type(psb_buffer_node), pointer :: prev=>null(), next=>null()
  end type psb_buffer_node

  type psb_buffer_queue
    type(psb_buffer_node), pointer :: head=>null(), tail=>null()
  end type psb_buffer_queue

  interface psi_snd
    module procedure&
         & psi_msnd, psi_esnd,&
         & psi_ssnd, psi_dsnd,&
         & psi_csnd, psi_zsnd,&
         & psi_logsnd, psi_hsnd,&
         & psi_i2snd
  end interface


  interface psb_init
    module procedure  psb_init_mpik
  end interface

  interface psb_exit
    module procedure  psb_exit_mpik
  end interface

  interface psb_abort
    module procedure  psb_abort_mpik
  end interface

  interface psb_info
    module procedure psb_info_mpik
  end interface
#if defined(IPK4) && defined(LPK8)
  interface psb_info
    module procedure psb_info_epk
  end interface
#endif
  
  interface psb_barrier
    module procedure  psb_barrier_mpik
  end interface
  
  interface psb_wtime
    module procedure  psb_wtime
  end interface psb_wtime

  interface psb_get_mpi_comm
    module procedure psb_m_get_mpi_comm !, psb_e_get_mpi_comm
  end interface psb_get_mpi_comm
  
  interface psb_get_mpi_rank
    module procedure psb_m_get_mpi_rank!, psb_e_get_mpi_rank
  end interface psb_get_mpi_rank

#if defined(SERIAL_MPI)
  integer(psb_mpk_), private, save :: nctxt=0

#else 

  integer(psb_mpk_), save :: mpi_iamx_op, mpi_iamn_op
  integer(psb_mpk_), save :: mpi_mamx_op, mpi_mamn_op
  integer(psb_mpk_), save :: mpi_eamx_op, mpi_eamn_op
  integer(psb_mpk_), save :: mpi_samx_op, mpi_samn_op
  integer(psb_mpk_), save :: mpi_damx_op, mpi_damn_op
  integer(psb_mpk_), save :: mpi_camx_op, mpi_camn_op
  integer(psb_mpk_), save :: mpi_zamx_op, mpi_zamn_op
  integer(psb_mpk_), save :: mpi_snrm2_op, mpi_dnrm2_op

  type(psb_buffer_queue), save :: psb_mesg_queue 

#endif

  private :: psi_get_sizes,  psi_register_mpi_extras
  private :: psi_iamx_op, psi_iamn_op 
  private :: psi_mamx_op, psi_mamn_op 
  private :: psi_eamx_op, psi_eamn_op 
  private :: psi_samx_op, psi_samn_op 
  private :: psi_damx_op, psi_damn_op 
  private :: psi_camx_op, psi_camn_op 
  private :: psi_zamx_op, psi_zamn_op 
  private :: psi_snrm2_op, psi_dnrm2_op 


contains
  subroutine psb_init_queue(mesg_queue,info)
    implicit none 
    type(psb_buffer_queue), intent(inout) :: mesg_queue
    integer(psb_ipk_), intent(out)                  :: info

    info = 0
    if ((.not.associated(mesg_queue%head)).and.&
         & (.not.associated(mesg_queue%tail))) then 
      ! Nothing to do
      return
    end if

    if ((.not.associated(mesg_queue%head)).or.&
         & (.not.associated(mesg_queue%tail))) then 
      ! If we are here one is associated, the other is not.
      ! This is impossible. 
      info = -1
      write(psb_err_unit,*) 'Wrong status on init '
      return
    end if

  end subroutine psb_init_queue

  subroutine psb_wait_buffer(node, info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_buffer_node), intent(inout) :: node
    integer(psb_ipk_), intent(out) :: info 
    integer(psb_mpk_) :: status(mpi_status_size),minfo
    minfo = mpi_success
    call mpi_wait(node%request,status,minfo)
    info=minfo
  end subroutine psb_wait_buffer

  subroutine psb_test_buffer(node, flag, info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_buffer_node), intent(inout) :: node
    logical, intent(out) :: flag
    integer(psb_ipk_), intent(out) :: info 
    integer(psb_mpk_) :: status(mpi_status_size), minfo
    minfo = mpi_success
#if defined(SERIAL_MPI)
    flag  = .true.
#else
    call mpi_test(node%request,flag,status,minfo)
#endif
    info=minfo
  end subroutine psb_test_buffer
  

  subroutine psb_close_context(mesg_queue,ctxt)
    type(psb_buffer_queue), intent(inout) :: mesg_queue
    type(psb_ctxt_type), intent(in) :: ctxt
    integer(psb_ipk_) :: info
    type(psb_buffer_node), pointer :: node, nextnode

    node => mesg_queue%head
    do 
      if (.not.associated(node)) exit
      nextnode => node%next
      if (psb_cmp_ctxt(node%ctxt,ctxt)) then 
        call psb_wait_buffer(node,info)
        call psb_delete_node(mesg_queue,node)
      end if
      node => nextnode
    end do
  end subroutine psb_close_context

  subroutine psb_close_all_context(mesg_queue)
    type(psb_buffer_queue), intent(inout) :: mesg_queue
    type(psb_buffer_node), pointer :: node, nextnode
    integer(psb_ipk_) :: info
    
    node => mesg_queue%head
    do 
      if (.not.associated(node)) exit
      nextnode => node%next
      call psb_wait_buffer(node,info)
      call psb_delete_node(mesg_queue,node)
      node => nextnode
    end do
  end subroutine psb_close_all_context


  subroutine psb_delete_node(mesg_queue,node)
    type(psb_buffer_queue), intent(inout) :: mesg_queue
    type(psb_buffer_node), pointer   :: node
    type(psb_buffer_node), pointer  :: prevnode
    
    if (.not.associated(node)) then 
      return
    end if
    prevnode => node%prev
    if (associated(mesg_queue%head,node)) mesg_queue%head => node%next
    if (associated(mesg_queue%tail,node)) mesg_queue%tail => prevnode
    if (associated(prevnode)) prevnode%next => node%next
    if (associated(node%next)) node%next%prev => prevnode
    deallocate(node)
    
  end subroutine psb_delete_node

  subroutine psb_insert_node(mesg_queue,node)
    type(psb_buffer_queue), intent(inout) :: mesg_queue
    type(psb_buffer_node), pointer   :: node

    node%next => null()
    node%prev => null()
    if ((.not.associated(mesg_queue%head)).and.&
         & (.not.associated(mesg_queue%tail))) then 
      mesg_Queue%head => node
      mesg_queue%tail => node
      return
    end if
    mesg_queue%tail%next => node
    node%prev => mesg_queue%tail
    mesg_queue%tail => node

  end subroutine psb_insert_node

  subroutine psb_test_nodes(mesg_queue)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node, nextnode
    integer(psb_ipk_) :: info
    logical :: flag
    
    node => mesg_queue%head
    do 
      if (.not.associated(node)) exit
      nextnode => node%next
      call psb_test_buffer(node,flag,info)
      if (flag) then 
        call psb_delete_node(mesg_queue,node)
      end if
      node => nextnode
    end do
  end subroutine psb_test_nodes

  ! !!!!!!!!!!!!!!!!!
  !
  ! Inner send. Basic idea:
  !  the input buffer is MOVE_ALLOCed
  !  to a node in the mesg queue, then it is sent.
  !  Thus the calling process should guarantee that
  !  the buffer is dispensable, i.e. the user data
  !  has already been copied. 
  !
  ! !!!!!!!!!!!!!!!!!
  subroutine psi_msnd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: tag, dest
    integer(psb_mpk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_int_type
    call move_alloc(buffer,node%int4buf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%int4buf,size(node%int4buf),psb_mpi_mpk_,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)

  end subroutine psi_msnd


  subroutine psi_esnd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: tag, dest
    integer(psb_epk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_int8_type
    call move_alloc(buffer,node%int8buf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%int8buf,size(node%int8buf),psb_mpi_epk_,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo 
    call psb_insert_node(mesg_queue,node)
    call psb_test_nodes(mesg_queue)

  end subroutine psi_esnd

  subroutine psi_i2snd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: tag, dest
    integer(psb_i2pk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_int2_type
    call move_alloc(buffer,node%int2buf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%int2buf,size(node%int2buf),psb_mpi_i2pk_,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    call psb_test_nodes(mesg_queue)

  end subroutine psi_i2snd

  subroutine psi_ssnd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: tag, dest
    real(psb_spk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm

    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_real_type
    call move_alloc(buffer,node%realbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%realbuf,size(node%realbuf),psb_mpi_r_spk_,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)  
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_ssnd

  subroutine psi_dsnd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_)   :: tag, dest
    real(psb_dpk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm

    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_double_type
    call move_alloc(buffer,node%doublebuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%doublebuf,size(node%doublebuf),psb_mpi_r_dpk_,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_dsnd
    
  subroutine psi_csnd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: tag, dest
    complex(psb_spk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm

    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_complex_type
    call move_alloc(buffer,node%complexbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%complexbuf,size(node%complexbuf),psb_mpi_c_spk_,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo 
    call psb_insert_node(mesg_queue,node)
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_csnd

  subroutine psi_zsnd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: tag, dest
    complex(psb_dpk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_dcomplex_type
    call move_alloc(buffer,node%dcomplbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%dcomplbuf,size(node%dcomplbuf),psb_mpi_c_dpk_,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_zsnd


  subroutine psi_logsnd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: tag, dest
    logical, allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_logical_type
    call move_alloc(buffer,node%logbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%logbuf,size(node%logbuf),mpi_logical,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_logsnd


  subroutine psi_hsnd(ctxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: tag, dest
    character(len=1), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpk_) :: minfo, icomm
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%ctxt     = ctxt
    icomm            = psb_get_mpi_comm(ctxt)
    node%buffer_type = psb_char_type
    call move_alloc(buffer,node%charbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%charbuf,size(node%charbuf),mpi_character,&
         & dest,tag,icomm,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_hsnd

  
  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Environment handling 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!

  subroutine psi_get_sizes()
    use psb_const_mod
    use iso_c_binding
    
    real(psb_dpk_), target     :: dv(2) 
    real(psb_spk_), target     :: sv(2) 
    integer(psb_i2pk_), target :: i2v(2)
    integer(psb_mpk_), target  :: mv(2)
    integer(psb_ipk_), target  :: iv(2)
    integer(psb_lpk_), target  :: lv(2)
    integer(psb_epk_), target  :: ev(2)
    interface
      subroutine psi_c_diffadd(p1, p2, val) &
           & bind(c,name="psi_c_diffadd")
        use iso_c_binding
        import :: psb_mpk_
        type(c_ptr), value :: p1, p2
        integer(psb_mpk_) :: val
      end subroutine psi_c_diffadd
    end interface
    
    call psi_c_diffadd(c_loc(sv(1)),c_loc(sv(2)),psb_sizeof_sp)
    call psi_c_diffadd(c_loc(dv(1)),c_loc(dv(2)),psb_sizeof_dp)
    call psi_c_diffadd(c_loc(i2v(1)),c_loc(i2v(2)),psb_sizeof_i2p)
    call psi_c_diffadd(c_loc(mv(1)),c_loc(mv(2)),psb_sizeof_mp)
    call psi_c_diffadd(c_loc(iv(1)),c_loc(iv(2)),psb_sizeof_ip)
    call psi_c_diffadd(c_loc(lv(1)),c_loc(lv(2)),psb_sizeof_lp)
    call psi_c_diffadd(c_loc(ev(1)),c_loc(ev(2)),psb_sizeof_ep)

  end subroutine psi_get_sizes

  subroutine  psi_register_mpi_extras(info)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpk_) :: info
    
    info = 0
#if 0
    if (info == 0) call mpi_type_create_f90_integer(psb_ipk_, psb_mpi_ipk_ ,info)
    if (info == 0) call mpi_type_create_f90_integer(psb_lpk_, psb_mpi_lpk_ ,info)
    if (info == 0) call mpi_type_create_f90_integer(psb_mpk_, psb_mpi_mpk_ ,info)
    if (info == 0) call mpi_type_create_f90_integer(psb_epk_, psb_mpi_lpk_ ,info)
    if (info == 0) call mpi_type_create_f90_real(psb_spk_p_,psb_spk_r_, psb_mpi_r_spk_,info)
    if (info == 0) call mpi_type_create_f90_real(psb_dpk_p_,psb_dpk_r_, psb_mpi_r_dpk_,info)
    if (info == 0) call mpi_type_create_f90_complex(psb_spk_p_,psb_spk_r_, psb_mpi_c_spk_,info)
    if (info == 0) call mpi_type_create_f90_complex(psb_dpk_p_,psb_dpk_r_, psb_mpi_c_dpk_,info)
#else
#if defined(IPK4) && defined(LPK4)
    psb_mpi_ipk_ = mpi_integer4
    psb_mpi_lpk_ = mpi_integer4
#elif defined(IPK4) && defined(LPK8)
    psb_mpi_ipk_ = mpi_integer4
    psb_mpi_lpk_ = mpi_integer8
#elif defined(IPK8) && defined(LPK8)
    psb_mpi_ipk_ = mpi_integer8
    psb_mpi_lpk_ = mpi_integer8
#else
    ! This should never happen
    write(psb_err_unit,*) 'Warning: an impossible IPK/LPK combination.'
    write(psb_err_unit,*) 'Something went wrong at configuration time.'
    psb_mpi_ipk_ = -1
    psb_mpi_lpk_ = -1
#endif
    psb_mpi_i2pk_ = mpi_integer2
    psb_mpi_mpk_  = mpi_integer4
    psb_mpi_epk_  = mpi_integer8
    psb_mpi_r_spk_  = mpi_real
    psb_mpi_r_dpk_  = mpi_double_precision
    psb_mpi_c_spk_  = mpi_complex
    psb_mpi_c_dpk_  = mpi_double_complex
#endif

#if defined(SERIAL_MPI)
#else 
    if (info == 0) call mpi_op_create(psi_mamx_op,.true.,mpi_mamx_op,info)
    if (info == 0) call mpi_op_create(psi_mamn_op,.true.,mpi_mamn_op,info)
    if (info == 0) call mpi_op_create(psi_eamx_op,.true.,mpi_eamx_op,info)
    if (info == 0) call mpi_op_create(psi_eamn_op,.true.,mpi_eamn_op,info)
    if (info == 0) call mpi_op_create(psi_samx_op,.true.,mpi_samx_op,info)
    if (info == 0) call mpi_op_create(psi_samn_op,.true.,mpi_samn_op,info)
    if (info == 0) call mpi_op_create(psi_damx_op,.true.,mpi_damx_op,info)
    if (info == 0) call mpi_op_create(psi_damn_op,.true.,mpi_damn_op,info)
    if (info == 0) call mpi_op_create(psi_camx_op,.true.,mpi_camx_op,info)
    if (info == 0) call mpi_op_create(psi_camn_op,.true.,mpi_camn_op,info)
    if (info == 0) call mpi_op_create(psi_zamx_op,.true.,mpi_zamx_op,info)
    if (info == 0) call mpi_op_create(psi_zamn_op,.true.,mpi_zamn_op,info)
    if (info == 0) call mpi_op_create(psi_snrm2_op,.true.,mpi_snrm2_op,info)
    if (info == 0) call mpi_op_create(psi_dnrm2_op,.true.,mpi_dnrm2_op,info)
#endif

  end subroutine psi_register_mpi_extras

#if defined(IPK4) && defined(LPK8)
  subroutine psb_info_epk(ctxt,iam,np)

    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_epk_), intent(out) :: iam, np

    !
    ! Simple caching scheme, keep track
    ! of the last CTXT encountered. 
    !
    integer(psb_mpk_), save :: lam, lnp
    call psb_info(ctxt,lam,lnp)
    iam = lam
    np  = lnp
  end subroutine psb_info_epk
#endif
  
  subroutine psb_init_mpik(ctxt,np,basectxt,ids)
    use psb_const_mod
    use psb_error_mod
    use psb_mat_mod
    use psb_vect_mod
! !$    use psb_rsb_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(out) :: ctxt
    type(psb_ctxt_type), intent(in), optional :: basectxt
    integer(psb_mpk_), intent(in), optional :: np, ids(:)

    integer(psb_mpk_) :: i, isnullcomm, icomm
    integer(psb_mpk_), allocatable :: iids(:) 
    logical :: initialized    
    integer(psb_mpk_) :: np_, npavail, iam, info, basecomm, basegroup, newgroup
    character(len=20), parameter :: name='psb_init'
    integer(psb_ipk_) :: iinfo
    !    
    call psb_set_debug_unit(psb_err_unit)

#if defined(SERIAL_MPI) 
    ctxt%ctxt = nctxt ! allocate on assignment
    nctxt = nctxt + 1

    call psi_register_mpi_extras(info)
    call psi_get_sizes()

#else    
    call mpi_initialized(initialized,info)
    if ((.not.initialized).or.(info /= mpi_success)) then 
      if (info == mpi_success) call mpi_init(info) 
      if (info /= mpi_success) then
        write(psb_err_unit,*) 'Error in initalizing MPI, bailing out',info 
        stop 
      end if
    end if

    if (present(basectxt)) then
      if (allocated(basectxt%ctxt)) then 
        basecomm = basectxt%ctxt
      else
        basecomm = mpi_comm_world
      end if
    else
      basecomm = mpi_comm_world
    end if

    if (present(np)) then 
      if (np < 1) then 
        iinfo=psb_err_initerror_neugh_procs_
        call psb_errpush(iinfo,name)
        call psb_error()
        !ctxt = mpi_comm_null
        return
      endif
      call mpi_comm_size(basecomm,np_,info)
      if (np_ < np) then 
        iinfo=psb_err_initerror_neugh_procs_
        call psb_errpush(iinfo,name)
        call psb_error()
        !ctxt = mpi_comm_null
        return
      endif
      call mpi_comm_group(basecomm,basegroup,info)
      if (present(ids)) then 
        if (size(ids)<np) then 
          write(psb_err_unit,*) 'Error in init: too few ids in input'
          !ctxt%ctxt = mpi_comm_null
          return
        end if
        do i=1, np 
          if ((ids(i)<0).or.(ids(i)>np_)) then 
            write(psb_err_unit,*) 'Error in init: invalid rank in input'
            !ctxt%ctxt = mpi_comm_null
            return
          end if
        end do
        call mpi_group_incl(basegroup,np,ids,newgroup,info)
        if (info /= mpi_success) then 
          !ctxt%ctxt = mpi_comm_null 
          return
        endif
      else
        allocate(iids(np),stat=info)
        if (info /= 0) then 
          !ctxt%ctxt = mpi_comm_null
          return
        endif
        do i=1, np
          iids(i) = i-1
        end do
        call mpi_group_incl(basegroup,np,iids,newgroup,info)
        if (info /= mpi_success) then 
          !ctxt = mpi_comm_null 
          return
        endif
        deallocate(iids)
      end if
      
      call mpi_comm_create(basecomm,newgroup,icomm,info)
   
    else
      if (basecomm /= mpi_comm_null) then 
        call mpi_comm_dup(basecomm,icomm,info)
      else 
        ! ctxt = mpi_comm_null
      end if
    endif
    if (info == 0) then
      ctxt%ctxt = icomm ! allocate on assignment
    end if
    call psi_register_mpi_extras(info)
    call psi_get_sizes()
    !if (ctxt == mpi_comm_null) return
    if (.not.allocated(ctxt%ctxt)) return 
#endif
    call psb_init_vect_defaults()
    call psb_init_mat_defaults()
    ! !$    call psb_rsb_init(info)
    ! !$    if (info.ne.psb_rsb_const_success) then 
    ! !$      if (info.eq.psb_rsb_const_not_available) then 
    ! !$        info=psb_success_ ! rsb is not present
    ! !$      else
    ! !$        ! rsb failed to initialize, and we issue an internal error.
    ! !$        ! or shall we tolerate this ?
    ! !$        info=psb_err_internal_error_
    ! !$        call psb_errpush(info,name)
    ! !$        call psb_error(ctxt)
    ! !$      endif
    ! !$    endif

  end subroutine psb_init_mpik

  subroutine psb_exit_mpik(ctxt,close)
    use psb_mat_mod
    use psb_vect_mod
! !$    use psb_rsb_mod
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(inout) :: ctxt
    logical, intent(in), optional :: close
    logical  :: close_
    integer(psb_mpk_) :: info
    character(len=20), parameter :: name='psb_exit'

    info = 0
    if (present(close)) then 
      close_ = close
    else
      close_ = .true.
    end if
! !$    if (close_) call psb_rsb_exit(info)
! !$    if (info.ne.psb_rsb_const_success) then 
! !$      if (info.eq.psb_rsb_const_not_available) then 
! !$        info=psb_success_ ! rsb is not present
! !$      else
! !$        info=psb_err_internal_error_ ! rsb failed to exit, and we issue an internal error. or  shall we tolerate this ?
! !$        call psb_errpush(info,name)
! !$        call psb_error(ctxt)
! !$      endif
! !$    endif
#if defined(SERIAL_MPI)
    ! Under serial mode, CLOSE has no effect, but reclaim
    ! the used ctxt number. 
    nctxt = max(0, nctxt - 1)    
#else 
    if (close_) then 
      call psb_close_all_context(psb_mesg_queue)
    else
      call psb_close_context(psb_mesg_queue,ctxt)
    end if
    !if ((ctxt /= mpi_comm_null).and.(ctxt /= mpi_comm_world)) then
    if (allocated(ctxt%ctxt)) then
      !write(0,*) ctxt%ctxt,mpi_comm_world,mpi_comm_null
      if ((ctxt%ctxt /= mpi_comm_world).and.(ctxt%ctxt /= mpi_comm_null)) &
           & call mpi_comm_Free(ctxt%ctxt,info)
    end if
    if (close_) then 
      if (info == 0) call mpi_op_free(mpi_mamx_op,info)
      if (info == 0) call mpi_op_free(mpi_mamn_op,info)
      if (info == 0) call mpi_op_free(mpi_eamx_op,info)
      if (info == 0) call mpi_op_free(mpi_eamn_op,info)
      if (info == 0) call mpi_op_free(mpi_samx_op,info)
      if (info == 0) call mpi_op_free(mpi_samn_op,info)
      if (info == 0) call mpi_op_free(mpi_damx_op,info)
      if (info == 0) call mpi_op_free(mpi_damn_op,info)
      if (info == 0) call mpi_op_free(mpi_camx_op,info)
      if (info == 0) call mpi_op_free(mpi_camn_op,info)
      if (info == 0) call mpi_op_free(mpi_zamx_op,info)
      if (info == 0) call mpi_op_free(mpi_zamn_op,info)
      if (info == 0) call mpi_op_free(mpi_snrm2_op,info)
      if (info == 0) call mpi_op_free(mpi_dnrm2_op,info)
    end if
    
    if (close_) call mpi_finalize(info)

#endif
    if (close_) call psb_clear_vect_defaults()
    if (close_) call psb_clear_mat_defaults()

  end subroutine psb_exit_mpik


  subroutine psb_barrier_mpik(ctxt)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type), intent(in) :: ctxt

    integer(psb_mpk_) :: info
#if !defined(SERIAL_MPI)
    if (allocated(ctxt%ctxt)) then 
      if (ctxt%ctxt /= mpi_comm_null) call mpi_barrier(ctxt%ctxt, info)
    end if
#endif    

  end subroutine psb_barrier_mpik

  function psb_wtime()
    use psb_const_mod
!    use mpi_constants
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    real(psb_dpk_) :: psb_wtime

    psb_wtime = mpi_wtime()
  end function psb_wtime

  subroutine psb_abort_mpik(ctxt,errc)

    type(psb_ctxt_type), intent(in) :: ctxt
    integer(psb_mpk_), intent(in), optional :: errc
    
    integer(psb_mpk_) :: code, info 

#if defined(SERIAL_MPI) 
    stop 
#else    
    if (present(errc)) then 
      code = errc
    else
      code = -1 
    endif

    if (allocated(ctxt%ctxt)) call mpi_abort(ctxt%ctxt,code,info)
#endif    

  end subroutine psb_abort_mpik


  subroutine psb_info_mpik(ctxt,iam,np)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif

    type(psb_ctxt_type), intent(in)  :: ctxt
    integer(psb_mpk_), intent(out) :: iam, np
    integer(psb_mpk_) :: info
    !
    ! Simple caching scheme, keep track
    ! of the last CTXT encountered. 
    !  
    integer(psb_mpk_), save :: lctxt=-1, lam, lnp

    !
    ! Note. There is NO way to ask MPI to check if
    ! a communicator handle is valid or not. Any
    ! call with an invalid handle will result in
    ! an error being thrown, and what happend
    ! depends on whether or not the error handler
    ! has been reset, which is a rather heavy-handed
    ! approach.
    ! This is why we  transformed ICTXT
    ! into an opaque object containing an ALLOCATABLE
    ! (be it an integer or a TYPE(MPI_COMM) object)
    ! and use its allocation status to record whether
    ! it's valid or not. 
    !
    
#if defined(SERIAL_MPI) 
    iam = 0
    np  = 1
#else    
    iam = -1
    np  = -1
    if (allocated(ctxt%ctxt)) then 
      if (ctxt%ctxt == lctxt) then
        iam = lam
        np  = lnp
      else
        if (ctxt%ctxt /= mpi_comm_null) then 
          call mpi_comm_size(ctxt%ctxt,np,info) 
          if (info /= mpi_success) np = -1 
          if (info == mpi_success) call mpi_comm_rank(ctxt%ctxt,iam,info) 
          if (info /= mpi_success) iam = -1 
        end if
        lctxt = ctxt%ctxt
        lam   = iam
        lnp   = np
      end if
    end if
#endif    
  end subroutine psb_info_mpik


  function psb_m_get_mpi_comm(ctxt) result(comm)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: comm
    comm = mpi_comm_null
    if (allocated(ctxt%ctxt)) comm = ctxt%ctxt
  end function psb_m_get_mpi_comm

  function psb_m_get_mpi_rank(ctxt,id) result(rank)
    integer(psb_mpk_) :: rank
    integer(psb_mpk_) :: id
    type(psb_ctxt_type) :: ctxt

    rank = id
  end function psb_m_get_mpi_rank

  subroutine psb_get_mpicomm(ctxt,comm)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: comm
    comm = mpi_comm_null
    if (allocated(ctxt%ctxt)) comm = ctxt%ctxt
  end subroutine psb_get_mpicomm

  subroutine psb_get_rank(rank,ctxt,id)
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: rank,id

    rank = psb_get_mpi_rank(ctxt,id) 
  end subroutine psb_get_rank


  ! !!!!!!!!!!!!!!!!!!!!!!
  !
  ! Base binary  operations
  !
  ! Note: len & type are always default integer.
  !
  ! !!!!!!!!!!!!!!!!!!!!!!
  subroutine psi_mamx_op(inv, outv,len,type) 
    integer(psb_mpk_) :: inv(len), outv(len)
    integer(psb_mpk_) :: len,type
    integer(psb_mpk_) :: i

    do i=1, len
      if (abs(inv(i)) > abs(outv(i))) outv(i) = inv(i)
    end do
  end subroutine psi_mamx_op

  subroutine psi_mamn_op(inv, outv,len,type) 
    integer(psb_mpk_) :: inv(len), outv(len)
    integer(psb_mpk_) :: len,type
    integer(psb_mpk_) :: i

    do i=1, len
      if (abs(inv(i)) < abs(outv(i))) outv(i) = inv(i)
    end do
  end subroutine psi_mamn_op

  subroutine psi_eamx_op(inv, outv,len,type) 
    integer(psb_epk_) :: inv(len), outv(len)
    integer(psb_mpk_) :: len,type
    integer(psb_mpk_) :: i

    do i=1, len
      if (abs(inv(i)) > abs(outv(i))) outv(i) = inv(i)
    end do
  end subroutine psi_eamx_op

  subroutine psi_eamn_op(inv, outv,len,type) 
    integer(psb_epk_) :: inv(len), outv(len)
    integer(psb_mpk_) :: len,type
    integer(psb_mpk_) :: i

    do i=1, len
      if (abs(inv(i)) < abs(outv(i))) outv(i) = inv(i)
    end do
  end subroutine psi_eamn_op

  subroutine psi_samx_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_spk_), intent(in)    :: vin(len)
    real(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) < abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_samx_op

  subroutine psi_samn_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_spk_), intent(in)    :: vin(len)
    real(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) > abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_samn_op

  subroutine psi_damx_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_dpk_), intent(in)    :: vin(len)
    real(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) < abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_damx_op

  subroutine psi_damn_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_dpk_), intent(in)    :: vin(len)
    real(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) > abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_damn_op

  subroutine psi_camx_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    complex(psb_spk_), intent(in)    :: vin(len)
    complex(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) < abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_camx_op

  subroutine psi_camn_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    complex(psb_spk_), intent(in)    :: vin(len)
    complex(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) > abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_camn_op

  subroutine psi_zamx_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    complex(psb_dpk_), intent(in)    :: vin(len)
    complex(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) < abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_zamx_op

  subroutine psi_zamn_op(vin,vinout,len,itype)
    integer(psb_mpk_), intent(in)           :: len, itype
    complex(psb_dpk_), intent(in)    :: vin(len)
    complex(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    do i=1, len
      if (abs(vinout(i)) > abs(vin(i))) vinout(i) = vin(i)
    end do
  end subroutine psi_zamn_op

  subroutine psi_snrm2_op(vin,vinout,len,itype)
    implicit none 
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_spk_), intent(in)    :: vin(len)
    real(psb_spk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    real(psb_spk_) :: w, z
    do i=1, len
      w = max( vin(i), vinout(i) )
      z = min( vin(i), vinout(i) )
      if ( z == szero ) then
        vinout(i) = w
      else
        vinout(i) = w*sqrt( sone+( z / w )**2 )
      end if
    end do
  end subroutine psi_snrm2_op

  subroutine psi_dnrm2_op(vin,vinout,len,itype)
    implicit none 
    integer(psb_mpk_), intent(in)           :: len, itype
    real(psb_dpk_), intent(in)    :: vin(len)
    real(psb_dpk_), intent(inout) :: vinout(len)

    integer(psb_mpk_) :: i
    real(psb_dpk_) :: w, z
    do i=1, len
      w = max( vin(i), vinout(i) )
      z = min( vin(i), vinout(i) )
      if ( z == dzero ) then
        vinout(i) = w
      else
        vinout(i) = w*sqrt( done+( z / w )**2 )
      end if
    end do
  end subroutine psi_dnrm2_op

  
end module psi_penv_mod
