!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
#if defined(SERIAL_MPI)
! Provide a fake mpi module just to keep the compiler(s) happy.
module mpi
  use psb_const_mod
  integer(psb_mpik_), parameter :: mpi_success          = 0
  integer(psb_mpik_), parameter :: mpi_request_null     = 0
  integer(psb_mpik_), parameter :: mpi_status_size      = 1
  integer(psb_mpik_), parameter :: mpi_integer          = 1
  integer(psb_mpik_), parameter :: mpi_integer8         = 2
  integer(psb_mpik_), parameter :: mpi_real             = 3
  integer(psb_mpik_), parameter :: mpi_double_precision = 4
  integer(psb_mpik_), parameter :: mpi_complex          = 5   
  integer(psb_mpik_), parameter :: mpi_double_complex   = 6 
  integer(psb_mpik_), parameter :: mpi_character        = 7
  integer(psb_mpik_), parameter :: mpi_logical          = 8
  integer(psb_mpik_), parameter :: mpi_integer2         = 9
  integer(psb_mpik_), parameter :: mpi_comm_null        = -1
  integer(psb_mpik_), parameter :: mpi_comm_world       = 1
  
  real(psb_dpk_), external :: mpi_wtime
end module mpi
#endif    

module psi_comm_buffers_mod
  use psb_const_mod

  integer(psb_mpik_), private, parameter:: psb_int_type      = 987543
  integer(psb_mpik_), private, parameter:: psb_real_type     = psb_int_type      + 1
  integer(psb_mpik_), private, parameter:: psb_double_type   = psb_real_type     + 1
  integer(psb_mpik_), private, parameter:: psb_complex_type  = psb_double_type   + 1
  integer(psb_mpik_), private, parameter:: psb_dcomplex_type = psb_complex_type  + 1
  integer(psb_mpik_), private, parameter:: psb_logical_type  = psb_dcomplex_type + 1
  integer(psb_mpik_), private, parameter:: psb_char_type     = psb_logical_type  + 1
  integer(psb_mpik_), private, parameter:: psb_int8_type     = psb_char_type     + 1
  integer(psb_mpik_), private, parameter:: psb_int2_type     = psb_int8_type     + 1
  integer(psb_mpik_), private, parameter:: psb_int4_type     = psb_int2_type     + 1


  type psb_buffer_node
    integer(psb_mpik_) :: request
    integer(psb_mpik_) :: icontxt 
    integer(psb_mpik_) :: buffer_type
    integer(psb_ipk_), allocatable        :: intbuf(:)
    integer(psb_long_int_k_), allocatable :: int8buf(:)
    integer(2), allocatable               :: int2buf(:)
    integer(psb_mpik_), allocatable       :: int4buf(:)
    real(psb_spk_), allocatable           :: realbuf(:)
    real(psb_dpk_), allocatable           :: doublebuf(:)
    complex(psb_spk_), allocatable        :: complexbuf(:)
    complex(psb_dpk_), allocatable        :: dcomplbuf(:)
    logical, allocatable                  :: logbuf(:)
    character(len=1), allocatable         :: charbuf(:)
    type(psb_buffer_node), pointer :: prev=>null(), next=>null()
  end type psb_buffer_node

  type psb_buffer_queue
    type(psb_buffer_node), pointer :: head=>null(), tail=>null()
  end type psb_buffer_queue


  interface psi_snd
    module procedure psi_isnd,&
         & psi_ssnd, psi_dsnd,&
         & psi_csnd, psi_zsnd,&
         & psi_lsnd, psi_hsnd
  end interface
#if defined(LONG_INTEGERS)
  interface psi_snd
    module procedure psi_i4snd
  end interface
#endif

#if !defined(LONG_INTEGERS)
  interface psi_snd
    module procedure psi_i8snd
  end interface
#endif

#if defined(SHORT_INTEGERS)
  interface psi_snd
    module procedure psi_i2snd
  end interface
#endif
  
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
    integer(psb_mpik_) :: status(mpi_status_size),minfo

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
    integer(psb_mpik_) :: status(mpi_status_size), minfo

    call mpi_test(node%request,flag,status,minfo)
    info=minfo
  end subroutine psb_test_buffer
  

  subroutine psb_close_context(mesg_queue,icontxt)
    type(psb_buffer_queue), intent(inout) :: mesg_queue
    integer(psb_mpik_), intent(in) :: icontxt
    integer(psb_ipk_) :: info
    type(psb_buffer_node), pointer :: node, nextnode

    node => mesg_queue%head
    do 
      if (.not.associated(node)) exit
      nextnode => node%next
      if (node%icontxt == icontxt) then 
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
  subroutine psi_isnd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    integer(psb_ipk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_int_type
    call move_alloc(buffer,node%intbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%intbuf,size(node%intbuf),psb_mpi_ipk_integer,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)

  end subroutine psi_isnd

#if defined(LONG_INTEGERS)
  subroutine psi_i4snd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    integer(psb_mpik_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_mpik_) :: info
    integer(psb_mpik_) :: minfo
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_int4_type
    call move_alloc(buffer,node%int4buf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%int4buf,size(node%int4buf),psb_mpi_def_integer,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo 
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)

  end subroutine psi_i4snd
#endif

#if !defined(LONG_INTEGERS)
  subroutine psi_i8snd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    integer(psb_long_int_k_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_int8_type
    call move_alloc(buffer,node%int8buf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%int8buf,size(node%int8buf),psb_mpi_lng_integer,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo 
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)

  end subroutine psi_i8snd
#endif

#if defined(SHORT_INTEGERS)
  subroutine psi_i2snd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    integer(2), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_int2_type
    call move_alloc(buffer,node%int2buf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%int2buf,size(node%int2buf),psb_mpi_def_integer2,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)

  end subroutine psi_i2snd
#endif

  subroutine psi_ssnd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    real(psb_spk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo

    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_real_type
    call move_alloc(buffer,node%realbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%realbuf,size(node%realbuf),mpi_real,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_ssnd

  subroutine psi_dsnd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    real(psb_dpk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo

    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_double_type
    call move_alloc(buffer,node%doublebuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%doublebuf,size(node%doublebuf),mpi_double_precision,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_dsnd
    
  subroutine psi_csnd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    complex(psb_spk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo

    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_complex_type
    call move_alloc(buffer,node%complexbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%complexbuf,size(node%complexbuf),mpi_complex,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo 
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_csnd

  subroutine psi_zsnd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    complex(psb_dpk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_dcomplex_type
    call move_alloc(buffer,node%dcomplbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%dcomplbuf,size(node%dcomplbuf),mpi_double_complex,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_zsnd


  subroutine psi_lsnd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    logical, allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_logical_type
    call move_alloc(buffer,node%logbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%logbuf,size(node%logbuf),mpi_logical,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_lsnd


  subroutine psi_hsnd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer(psb_mpik_) :: icontxt, tag, dest
    character(len=1), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer(psb_ipk_) :: info
    integer(psb_mpik_) :: minfo
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_char_type
    call move_alloc(buffer,node%charbuf)
    if (info /= 0) then 
      write(psb_err_unit,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%charbuf,size(node%charbuf),mpi_character,&
         & dest,tag,icontxt,node%request,minfo)
    info = minfo
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_hsnd


end module psi_comm_buffers_mod

