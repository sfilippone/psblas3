module psi_comm_buffers_mod
  use psb_const_mod

  integer, private, parameter:: psb_int_type      = 987543
  integer, private, parameter:: psb_real_type     = psb_int_type      + 1
  integer, private, parameter:: psb_double_type   = psb_real_type     + 1
  integer, private, parameter:: psb_complex_type  = psb_double_type   + 1
  integer, private, parameter:: psb_dcomplex_type = psb_complex_type  + 1
  integer, private, parameter:: psb_logical_type  = psb_dcomplex_type + 1
  integer, private, parameter:: psb_char_type     = psb_logical_type  + 1
  integer, private, parameter:: psb_int8_type     = psb_char_type     + 1


  type psb_buffer_node
    integer :: request
    integer :: icontxt 
    integer :: buffer_type
    integer(psb_int_k_), allocatable      :: intbuf(:)
    integer(psb_long_int_k_), allocatable :: int8buf(:)
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
#if !defined(LONG_INTEGERS)
  interface psi_snd
    module procedure psi_i8snd
  end interface
#endif
  
contains

  subroutine psb_init_queue(mesg_queue,info)
    type(psb_buffer_queue), intent(inout) :: mesg_queue
    type(psb_buffer_node), pointer :: item

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
      write(0,*) 'Wrong status on init '
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
    integer, intent(out) :: info 
    integer :: status(mpi_status_size)

    call mpi_wait(node%request,status,info)
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
    integer, intent(out) :: info 
    integer :: status(mpi_status_size)

    call mpi_test(node%request,flag,status,info)
  end subroutine psb_test_buffer
  

  subroutine psb_close_context(mesg_queue,icontxt)
    type(psb_buffer_queue), intent(inout) :: mesg_queue
    integer, intent(in) :: icontxt
    integer :: info
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
    integer :: info
    
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
    integer :: info
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
    integer :: icontxt, tag, dest
    integer(psb_int_k_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer :: info
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_int_type
    call move_alloc(buffer,node%intbuf)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%intbuf,size(node%intbuf),psb_mpi_integer,&
         & dest,tag,icontxt,node%request,info)
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)

  end subroutine psi_isnd

#if !defined(LONG_INTEGERS)
  subroutine psi_i8snd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer :: icontxt, tag, dest
    integer(psb_long_int_k_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer :: info
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_int8_type
    call move_alloc(buffer,node%int8buf)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%int8buf,size(node%int8buf),mpi_integer8,&
         & dest,tag,icontxt,node%request,info)
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)

  end subroutine psi_i8snd
#endif


  subroutine psi_ssnd(icontxt,tag,dest,buffer,mesg_queue)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    integer :: icontxt, tag, dest
    real(psb_spk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer :: info
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_real_type
    call move_alloc(buffer,node%realbuf)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%realbuf,size(node%realbuf),mpi_real,&
         & dest,tag,icontxt,node%request,info)
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
    integer :: icontxt, tag, dest
    real(psb_dpk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer :: info
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_double_type
    call move_alloc(buffer,node%doublebuf)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%doublebuf,size(node%doublebuf),mpi_double_precision,&
         & dest,tag,icontxt,node%request,info)
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
    integer :: icontxt, tag, dest
    complex(psb_spk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer :: info
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_complex_type
    call move_alloc(buffer,node%complexbuf)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%complexbuf,size(node%complexbuf),mpi_complex,&
         & dest,tag,icontxt,node%request,info)
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
    integer :: icontxt, tag, dest
    complex(psb_dpk_), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer :: info
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_dcomplex_type
    call move_alloc(buffer,node%dcomplbuf)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%dcomplbuf,size(node%dcomplbuf),mpi_double_complex,&
         & dest,tag,icontxt,node%request,info)
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
    integer :: icontxt, tag, dest
    logical, allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer :: info
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_logical_type
    call move_alloc(buffer,node%logbuf)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%logbuf,size(node%logbuf),mpi_logical,&
         & dest,tag,icontxt,node%request,info)
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
    integer :: icontxt, tag, dest
    character(len=1), allocatable, intent(inout) :: buffer(:)
    type(psb_buffer_queue) :: mesg_queue
    type(psb_buffer_node), pointer :: node
    integer :: info
    
    allocate(node, stat=info)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    node%icontxt     = icontxt
    node%buffer_type = psb_char_type
    call move_alloc(buffer,node%charbuf)
    if (info /= 0) then 
      write(0,*) 'Fatal memory error inside communication subsystem'
      return
    end if
    call mpi_isend(node%charbuf,size(node%charbuf),mpi_character,&
         & dest,tag,icontxt,node%request,info)
    call psb_insert_node(mesg_queue,node)
    
    call psb_test_nodes(mesg_queue)
    
  end subroutine psi_hsnd


end module psi_comm_buffers_mod

