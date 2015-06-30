!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
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
!
! package: psb_i_base_vect_mod
!
! This module contains the definition of the psb_i_base_vect type which
! is a container for dense vectors.
!  This is encapsulated instead of being just a simple array to allow for 
!  more complicated situations, such as GPU programming, where the memory
!  area we are interested in is not easily accessible from the host/Fortran
!  side. It is also meant to be encapsulated in an outer type, to allow
!  runtime switching as per the STATE design pattern, similar to the
!  sparse matrix types.
!
!
module psb_i_base_vect_mod
  
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod

  !> \namespace  psb_base_mod  \class psb_i_base_vect_type
  !! The psb_i_base_vect_type 
  !! defines a middle level  integer(psb_ipk_) encapsulated dense vector.
  !! The encapsulation is needed, in place of a simple array, to allow  
  !! for complicated situations, such as GPU programming, where the memory
  !!  area we are interested in is not easily accessible from the host/Fortran
  !!  side. It is also meant to be encapsulated in an outer type, to allow
  !!  runtime switching as per the STATE design pattern, similar to the
  !!  sparse matrix types.
  !!
  type psb_i_base_vect_type
    !> Values. 
    integer(psb_ipk_), allocatable :: v(:)
    integer(psb_ipk_), allocatable :: combuf(:) 
    integer(psb_ipk_), allocatable :: comid(:,:)
  contains
    !
    !  Constructors/allocators
    !
    procedure, pass(x) :: bld_x    => i_base_bld_x
    procedure, pass(x) :: bld_n    => i_base_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: all      => i_base_all
    procedure, pass(x) :: mold     => i_base_mold
    !
    ! Insert/set. Assembly and free.
    ! Assembly does almost nothing here, but is important
    ! in derived classes. 
    !
    procedure, pass(x) :: ins_a    => i_base_ins_a
    procedure, pass(x) :: ins_v    => i_base_ins_v
    generic, public    :: ins      => ins_a, ins_v
    procedure, pass(x) :: zero     => i_base_zero
    procedure, pass(x) :: asb      => i_base_asb
    procedure, pass(x) :: free     => i_base_free
    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder. 
    !
    procedure, pass(x) :: sync     => i_base_sync
    procedure, pass(x) :: is_host  => i_base_is_host
    procedure, pass(x) :: is_dev   => i_base_is_dev
    procedure, pass(x) :: is_sync  => i_base_is_sync
    procedure, pass(x) :: set_host => i_base_set_host
    procedure, pass(x) :: set_dev  => i_base_set_dev
    procedure, pass(x) :: set_sync => i_base_set_sync

    !
    ! These are for handling gather/scatter in new
    ! comm internals implementation.
    !
    procedure, nopass  :: use_buffer   => i_base_use_buffer
    procedure, pass(x) :: new_buffer   => i_base_new_buffer
    procedure, nopass  :: device_wait  => i_base_device_wait
    procedure, pass(x) :: free_buffer  => i_base_free_buffer
    procedure, pass(x) :: new_comid    => i_base_new_comid
    procedure, pass(x) :: free_comid   => i_base_free_comid

    !
    ! Basic info
    procedure, pass(x) :: get_nrows => i_base_get_nrows
    procedure, pass(x) :: sizeof    => i_base_sizeof
    procedure, nopass  :: get_fmt   => i_base_get_fmt
    !
    ! Set/get data from/to an external array; also
    ! overload assignment.
    !
    procedure, pass(x) :: get_vect => i_base_get_vect
    procedure, pass(x) :: set_scal => i_base_set_scal
    procedure, pass(x) :: set_vect => i_base_set_vect
    generic, public    :: set      => set_vect, set_scal
    !
    ! Gather/scatter. These are needed for MPI interfacing.
    ! May have to be reworked. 
    !
    procedure, pass(x) :: gthab    => i_base_gthab
    procedure, pass(x) :: gthzv    => i_base_gthzv
    procedure, pass(x) :: gthzv_x  => i_base_gthzv_x
    procedure, pass(x) :: gthzbuf  => i_base_gthzbuf
    generic, public    :: gth      => gthab, gthzv, gthzv_x, gthzbuf
    procedure, pass(y) :: sctb     => i_base_sctb
    procedure, pass(y) :: sctb_x   => i_base_sctb_x
    procedure, pass(y) :: sctb_buf => i_base_sctb_buf
    generic, public    :: sct      => sctb, sctb_x, sctb_buf



  end type psb_i_base_vect_type

  public  :: psb_i_base_vect
  private :: constructor, size_const
  interface psb_i_base_vect
    module procedure constructor, size_const
  end interface psb_i_base_vect

contains
  
  !
  ! Constructors. 
  !
  
  !> Function  constructor:
  !! \brief     Constructor from an array
  !!  \param   x(:)  input array to be copied
  !!
  function constructor(x) result(this)
    integer(psb_ipk_)   :: x(:)
    type(psb_i_base_vect_type) :: this
    integer(psb_ipk_) :: info

    this%v = x
    call this%asb(size(x,kind=psb_ipk_),info)
  end function constructor
    
  
  !> Function  constructor:
  !! \brief     Constructor from size
  !!  \param    n   Size of vector to be built. 
  !!
  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_i_base_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%asb(n,info)

  end function size_const
  
  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  subroutine i_base_bld_x(x,this)
    use psb_realloc_mod
    integer(psb_ipk_), intent(in) :: this(:)
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this),x%v,info)
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,'base_vect_bld')
      return
    end if
    x%v(:)  = this(:)

  end subroutine i_base_bld_x
    
  !
  ! Create with size, but no initialization
  !

  !> Function  bld_n:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated. 
  !!
  subroutine i_base_bld_n(x,n)
    use psb_realloc_mod
    integer(psb_ipk_), intent(in) :: n
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(n,x%v,info)
    call x%asb(n,info)

  end subroutine i_base_bld_n
  
  !> Function  base_all:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated. 
  !!  \param info  return code
  !!
  subroutine i_base_all(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)               :: n
    class(psb_i_base_vect_type), intent(out)    :: x
    integer(psb_ipk_), intent(out)              :: info
    
    call psb_realloc(n,x%v,info)
    
  end subroutine i_base_all

  !> Function  base_mold:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  subroutine i_base_mold(x, y, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_i_base_vect_type), intent(in)   :: x
    class(psb_i_base_vect_type), intent(out), allocatable :: y
    integer(psb_ipk_), intent(out)              :: info
    
    allocate(psb_i_base_vect_type :: y, stat=info)

  end subroutine i_base_mold

  !
  ! Insert a bunch of values at specified positions.
  !
  !> Function  base_ins:
  !! \memberof  psb_i_base_vect_type
  !! \brief Insert coefficients. 
  !!
  !!
  !!         Given  a list of N pairs
  !!           (IRL(i),VAL(i))
  !!         record a new coefficient in X such that
  !!            X(IRL(1:N)) = VAL(1:N).
  !!            
  !!         - the update operation will perform either
  !!               X(IRL(1:n)) = VAL(1:N)
  !!           or
  !!               X(IRL(1:n)) = X(IRL(1:n))+VAL(1:N)
  !!           according to the value of DUPLICATE.
  !!           
  !!           
  !!  \param n     number of pairs in input
  !!  \param irl(:)  the input row indices
  !!  \param val(:)  the input coefficients
  !!  \param dupl    how to treat duplicate entries
  !!  \param info  return code
  !!
  !
  subroutine i_base_ins_a(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_ipk_), intent(in)        :: val(:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, isz

    info = 0
    if (psb_errstatus_fatal()) return 

    if (.not.allocated(x%v)) then 
      info = psb_err_invalid_vect_state_
    else if (n > min(size(irl),size(val))) then 
      info = psb_err_invalid_input_

    else 
      isz = size(x%v)
      select case(dupl) 
      case(psb_dupl_ovwrt_) 
        do i = 1, n
          !loop over all val's rows

          ! row actual block row 
          if ((1 <= irl(i)).and.(irl(i) <= isz)) then
            ! this row belongs to me
            ! copy i-th row of block val in x
            x%v(irl(i)) = val(i)
          end if
        enddo

      case(psb_dupl_add_) 

        do i = 1, n
          !loop over all val's rows
          if ((1 <= irl(i)).and.(irl(i) <= isz)) then
            ! this row belongs to me
            ! copy i-th row of block val in x
            x%v(irl(i)) = x%v(irl(i)) +  val(i)
          end if
        enddo

      case default
        info = 321
        ! !$      call psb_errpush(info,name)
        ! !$      goto 9999
      end select
    end if
    call x%set_host()
    if (info /= 0) then 
      call psb_errpush(info,'base_vect_ins')
      return
    end if

  end subroutine i_base_ins_a

  subroutine i_base_ins_v(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    class(psb_i_base_vect_type), intent(inout)  :: irl
    class(psb_i_base_vect_type), intent(inout)  :: val
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, isz

    info = 0
    if (psb_errstatus_fatal()) return 

    if (irl%is_dev()) call irl%sync()
    if (val%is_dev()) call val%sync()
    if (x%is_dev())   call x%sync()
    call x%ins(n,irl%v,val%v,dupl,info)

    if (info /= 0) then 
      call psb_errpush(info,'base_vect_ins')
      return
    end if

  end subroutine i_base_ins_v


  !
  !> Function  base_zero
  !! \memberof  psb_i_base_vect_type
  !! \brief Zero out contents
  !!
  !
  subroutine i_base_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout)    :: x
    
    if (allocated(x%v)) x%v=izero
    call x%set_host()
  end subroutine i_base_zero

  
  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_asb:
  !! \memberof  psb_i_base_vect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!           
  !!  \param n     final size
  !!  \param info  return code
  !!
  !
 
  subroutine i_base_asb(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)              :: n
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info
    
    info = 0
    if (x%get_nrows() < n) &
         & call psb_realloc(n,x%v,info)
    if (info /= 0) &
         & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')
    call x%sync()
  end subroutine i_base_asb


  !
  !> Function  base_free:
  !! \memberof  psb_i_base_vect_type
  !! \brief Free vector
  !!           
  !!  \param info  return code
  !!
  !
  subroutine i_base_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info
    
    info = 0
    if (allocated(x%v)) deallocate(x%v, stat=info)
    if (info == 0) call x%free_buffer(info)
    if (info == 0) call x%free_comid(info)
    if (info /= 0) call & 
         & psb_errpush(psb_err_alloc_dealloc_,'vect_free')
        
  end subroutine i_base_free

  

  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  ! 
  !
  !> Function  base_sync:
  !! \memberof  psb_i_base_vect_type
  !! \brief Sync: base version is a no-op.
  !!           
  !
  subroutine i_base_sync(x)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    
  end subroutine i_base_sync

  !
  !> Function  base_set_host:
  !! \memberof  psb_i_base_vect_type
  !! \brief Set_host: base version is a no-op.
  !!           
  !
  subroutine i_base_set_host(x)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    
  end subroutine i_base_set_host

  !
  !> Function  base_set_dev:
  !! \memberof  psb_i_base_vect_type
  !! \brief Set_dev: base version is a no-op.
  !!           
  !
  subroutine i_base_set_dev(x)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    
  end subroutine i_base_set_dev

  !
  !> Function  base_set_sync:
  !! \memberof  psb_i_base_vect_type
  !! \brief Set_sync: base version is a no-op.
  !!           
  !
  subroutine i_base_set_sync(x)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    
  end subroutine i_base_set_sync

  !
  !> Function  base_is_dev:
  !! \memberof  psb_i_base_vect_type
  !! \brief Is  vector on external device    .
  !!           
  !
  function i_base_is_dev(x) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(in) :: x
    logical  :: res
  
    res = .false.
  end function i_base_is_dev
  
  !
  !> Function  base_is_host
  !! \memberof  psb_i_base_vect_type
  !! \brief Is  vector on standard memory    .
  !!           
  !
  function i_base_is_host(x) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function i_base_is_host

  !
  !> Function  base_is_sync
  !! \memberof  psb_i_base_vect_type
  !! \brief Is  vector on sync               .
  !!           
  !
  function i_base_is_sync(x) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function i_base_is_sync


  !
  ! Size info. 
  !
  !
  !> Function  base_get_nrows
  !! \memberof  psb_i_base_vect_type
  !! \brief  Number of entries
  !!           
  !
  function i_base_get_nrows(x) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v)

  end function i_base_get_nrows

  !
  !> Function  base_get_sizeof
  !! \memberof  psb_i_base_vect_type
  !! \brief  Size in bytes
  !!           
  !
  function i_base_sizeof(x) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(in) :: x
    integer(psb_long_int_k_) :: res
    
    ! Force 8-byte integers.
    res = (1_psb_long_int_k_ * psb_sizeof_int) * x%get_nrows()

  end function i_base_sizeof

  !
  !> Function  base_get_fmt
  !! \memberof  psb_i_base_vect_type
  !! \brief  Format
  !!           
  !
  function i_base_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'BASE'
  end function i_base_get_fmt
  

  !
  !
  !
  !> Function  base_get_vect
  !! \memberof  psb_i_base_vect_type
  !! \brief  Extract a copy of the contents
  !!
  !    
  function  i_base_get_vect(x) result(res)
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), allocatable                 :: res(:)
    integer(psb_ipk_) :: info
    
    if (.not.allocated(x%v)) return 
    if (.not.x%is_host()) call x%sync()
    allocate(res(x%get_nrows()),stat=info) 
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,'base_get_vect')
      return
    end if
    res(:) = x%v(:)
  end function i_base_get_vect
    
  !
  ! Reset all values 
  !
  !
  !> Function  base_set_scal
  !! \memberof  psb_i_base_vect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  subroutine i_base_set_scal(x,val,first,last)
    class(psb_i_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in) :: val
    integer(psb_ipk_), optional :: first, last
        
    integer(psb_ipk_) :: info, first_, last_

    first_=1
    last_=size(x%v)
    if (present(first)) first_ = max(1,first)
    if (present(last))  last_  = min(last,last_)
    
    if (x%is_dev()) call x%sync()
    x%v(first_:last_) = val
    call x%set_host()

  end subroutine i_base_set_scal


  !
  !> Function  base_set_vect
  !! \memberof  psb_i_base_vect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in 
  !!
  subroutine i_base_set_vect(x,val,first,last)
    class(psb_i_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in) :: val(:)
    integer(psb_ipk_), optional :: first, last
        
    integer(psb_ipk_) :: info, first_, last_, nr

    first_=1
    last_=min(psb_size(x%v),size(val))
    if (present(first)) first_ = max(1,first)
    if (present(last))  last_  = min(last,last_)

    if (allocated(x%v)) then 
      if (x%is_dev()) call x%sync()
      x%v(first_:last_) = val(1:last_-first_+1)
    else
      x%v = val
    end if
    call x%set_host()

  end subroutine i_base_set_vect


  
  
  !
  ! Gather: Y = beta * Y + alpha * X(IDX(:))
  !
  !
  !> Function  base_gthab
  !! \memberof  psb_i_base_vect_type
  !! \brief gather into an array
  !!    Y = beta * Y + alpha * X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param alpha
  !! \param beta
  subroutine i_base_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_ipk_) :: alpha, beta, y(:)
    class(psb_i_base_vect_type) :: x
    
    if (x%is_dev()) call x%sync()
    call psi_gth(n,idx,alpha,x%v,beta,y)

  end subroutine i_base_gthab
  !
  ! shortcut alpha=1 beta=0
  ! 
  !> Function  base_gthzv
  !! \memberof  psb_i_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  subroutine i_base_gthzv_x(i,n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    integer(psb_ipk_) ::  y(:)
    class(psb_i_base_vect_type) :: x
    
    if (idx%is_dev()) call idx%sync()
    call x%gth(n,idx%v(i:),y)

  end subroutine i_base_gthzv_x

  !
  ! New comm internals impl. 
  !
  subroutine i_base_gthzbuf(i,n,idx,x)
    use psi_serial_mod
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    class(psb_i_base_vect_type) :: x
    
    if (.not.allocated(x%combuf)) then 
      call psb_errpush(psb_err_alloc_dealloc_,'gthzbuf')
      return
    end if
    if (idx%is_dev()) call idx%sync()
    if (x%is_dev()) call x%sync()
    call x%gth(n,idx%v(i:),x%combuf(i:))

  end subroutine i_base_gthzbuf

  subroutine i_base_sctb_buf(i,n,idx,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    integer(psb_ipk_) :: beta
    class(psb_i_base_vect_type) :: y
    
    
    if (.not.allocated(y%combuf)) then 
      call psb_errpush(psb_err_alloc_dealloc_,'sctb_buf')
      return
    end if
    if (y%is_dev()) call y%sync()
    if (idx%is_dev()) call idx%sync()
    call y%sct(n,idx%v(i:),y%combuf(i:),beta)
    call y%set_host()

  end subroutine i_base_sctb_buf

  !
  !> Function  base_device_wait:
  !! \memberof  psb_i_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!           
  !
  subroutine i_base_device_wait()
    implicit none 
    
  end subroutine i_base_device_wait

  function i_base_use_buffer() result(res)
    logical :: res
    
    res = .true.
  end function i_base_use_buffer

  subroutine i_base_new_buffer(n,x,info)
    use psb_realloc_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    call psb_realloc(n,x%combuf,info)
  end subroutine i_base_new_buffer

  subroutine i_base_new_comid(n,x,info)
    use psb_realloc_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    call psb_realloc(n,2,x%comid,info)
  end subroutine i_base_new_comid


  subroutine i_base_free_buffer(x,info)
    use psb_realloc_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%combuf)) &
         &  deallocate(x%combuf,stat=info)
  end subroutine i_base_free_buffer

  subroutine i_base_free_comid(x,info)
    use psb_realloc_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%comid)) &
         &  deallocate(x%comid,stat=info)
  end subroutine i_base_free_comid


  !
  ! shortcut alpha=1 beta=0
  ! 
  !> Function  base_gthzv
  !! \memberof  psb_i_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  subroutine i_base_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_ipk_) ::  y(:)
    class(psb_i_base_vect_type) :: x
    
    if (x%is_dev()) call x%sync()
    call psi_gth(n,idx,x%v,y)

  end subroutine i_base_gthzv

  !
  ! Scatter: 
  ! Y(IDX(:)) = beta*Y(IDX(:)) + X(:)
  ! 
  !
  !> Function  base_sctb
  !! \memberof  psb_i_base_vect_type
  !! \brief scatter into a class(base_vect)
  !!    Y(IDX(:)) = beta * Y(IDX(:)) +  X(:)
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param beta
  !! \param x(:) 
  subroutine i_base_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_ipk_) :: beta, x(:)
    class(psb_i_base_vect_type) :: y
    
    if (y%is_dev()) call y%sync()
    call psi_sct(n,idx,x,beta,y%v)
    call y%set_host()

  end subroutine i_base_sctb

  subroutine i_base_sctb_x(i,n,idx,x,beta,y)
    use psi_serial_mod
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    integer(psb_ipk_) :: beta, x(:)
    class(psb_i_base_vect_type) :: y
    
    if (idx%is_dev()) call idx%sync()
    call y%sct(n,idx%v(i:),x,beta)
    call y%set_host()

  end subroutine i_base_sctb_x

end module psb_i_base_vect_mod





module psb_i_base_multivect_mod
  
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod

  !> \namespace  psb_base_mod  \class psb_i_base_vect_type
  !! The psb_i_base_vect_type 
  !! defines a middle level  integer(psb_ipk_) encapsulated dense vector.
  !! The encapsulation is needed, in place of a simple array, to allow  
  !! for complicated situations, such as GPU programming, where the memory
  !!  area we are interested in is not easily accessible from the host/Fortran
  !!  side. It is also meant to be encapsulated in an outer type, to allow
  !!  runtime switching as per the STATE design pattern, similar to the
  !!  sparse matrix types.
  !!
  private 
  public  :: psb_i_base_multivect, psb_i_base_multivect_type

  type psb_i_base_multivect_type
    !> Values. 
    integer(psb_ipk_), allocatable :: v(:,:)
  contains
    !
    !  Constructors/allocators
    !
    procedure, pass(x) :: bld_x    => i_base_mlv_bld_x
    procedure, pass(x) :: bld_n    => i_base_mlv_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: all      => i_base_mlv_all
    procedure, pass(x) :: mold     => i_base_mlv_mold
    !
    ! Insert/set. Assembly and free.
    ! Assembly does almost nothing here, but is important
    ! in derived classes. 
    !
    procedure, pass(x) :: ins      => i_base_mlv_ins
    procedure, pass(x) :: zero     => i_base_mlv_zero
    procedure, pass(x) :: asb      => i_base_mlv_asb
    procedure, pass(x) :: free     => i_base_mlv_free
    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder. 
    !
    procedure, pass(x) :: sync     => i_base_mlv_sync
    procedure, pass(x) :: is_host  => i_base_mlv_is_host
    procedure, pass(x) :: is_dev   => i_base_mlv_is_dev
    procedure, pass(x) :: is_sync  => i_base_mlv_is_sync
    procedure, pass(x) :: set_host => i_base_mlv_set_host
    procedure, pass(x) :: set_dev  => i_base_mlv_set_dev
    procedure, pass(x) :: set_sync => i_base_mlv_set_sync

    !
    ! Basic info
    procedure, pass(x) :: get_nrows => i_base_mlv_get_nrows
    procedure, pass(x) :: get_ncols => i_base_mlv_get_ncols
    procedure, pass(x) :: sizeof    => i_base_mlv_sizeof
    procedure, nopass  :: get_fmt   => i_base_mlv_get_fmt
    !
    ! Set/get data from/to an external array; also
    ! overload assignment.
    !
    procedure, pass(x) :: get_vect => i_base_mlv_get_vect
    procedure, pass(x) :: set_scal => i_base_mlv_set_scal
    procedure, pass(x) :: set_vect => i_base_mlv_set_vect
    generic, public    :: set      => set_vect, set_scal

  end type psb_i_base_multivect_type

  interface psb_i_base_multivect
    module procedure constructor, size_const
  end interface

contains
  
  !
  ! Constructors. 
  !
  
  !> Function  constructor:
  !! \brief     Constructor from an array
  !!  \param   x(:)  input array to be copied
  !!
  function constructor(x) result(this)
    integer(psb_ipk_)   :: x(:,:)
    type(psb_i_base_multivect_type) :: this
    integer(psb_ipk_) :: info

    this%v = x
    call this%asb(size(x,dim=1,kind=psb_ipk_),size(x,dim=2,kind=psb_ipk_),info)
  end function constructor
    
  
  !> Function  constructor:
  !! \brief     Constructor from size
  !!  \param    n   Size of vector to be built. 
  !!
  function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_i_base_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%asb(m,n,info)

  end function size_const
  
  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_i_base_multivect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  subroutine i_base_mlv_bld_x(x,this)
    use psb_realloc_mod
    integer(psb_ipk_), intent(in) :: this(:,:)
    class(psb_i_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this,1),size(this,2),x%v,info)
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,'base_mlv_vect_bld')
      return
    end if
    x%v(:,:)  = this(:,:)

  end subroutine i_base_mlv_bld_x
    
  !
  ! Create with size, but no initialization
  !

  !> Function  bld_n:
  !! \memberof  psb_i_base_multivect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated. 
  !!
  subroutine i_base_mlv_bld_n(x,m,n)
    use psb_realloc_mod
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_i_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(m,n,x%v,info)
    call x%asb(m,n,info)

  end subroutine i_base_mlv_bld_n
  
  !> Function  base_mlv_all:
  !! \memberof  psb_i_base_multivect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated. 
  !!  \param info  return code
  !!
  subroutine i_base_mlv_all(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)               :: m,n
    class(psb_i_base_multivect_type), intent(out) :: x
    integer(psb_ipk_), intent(out)              :: info
    
    call psb_realloc(m,n,x%v,info)
    
  end subroutine i_base_mlv_all

  !> Function  base_mlv_mold:
  !! \memberof  psb_i_base_multivect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  subroutine i_base_mlv_mold(x, y, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_i_base_multivect_type), intent(in)   :: x
    class(psb_i_base_multivect_type), intent(out), allocatable :: y
    integer(psb_ipk_), intent(out)              :: info
    
    allocate(psb_i_base_multivect_type :: y, stat=info)

  end subroutine i_base_mlv_mold

  !
  ! Insert a bunch of values at specified positions.
  !
  !> Function  base_mlv_ins:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Insert coefficients. 
  !!
  !!
  !!         Given  a list of N pairs
  !!           (IRL(i),VAL(i))
  !!         record a new coefficient in X such that
  !!            X(IRL(1:N)) = VAL(1:N).
  !!            
  !!         - the update operation will perform either
  !!               X(IRL(1:n)) = VAL(1:N)
  !!           or
  !!               X(IRL(1:n)) = X(IRL(1:n))+VAL(1:N)
  !!           according to the value of DUPLICATE.
  !!           
  !!           
  !!  \param n     number of pairs in input
  !!  \param irl(:)  the input row indices
  !!  \param val(:)  the input coefficients
  !!  \param dupl    how to treat duplicate entries
  !!  \param info  return code
  !!
  !
  subroutine i_base_mlv_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none 
    class(psb_i_base_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    integer(psb_ipk_), intent(in)        :: val(:,:)
    integer(psb_ipk_), intent(out)              :: info

    integer(psb_ipk_) :: i, isz

    info = 0
    if (psb_errstatus_fatal()) return 

    if (.not.allocated(x%v)) then 
      info = psb_err_invalid_vect_state_
    else if (n > min(size(irl),size(val))) then 
      info = psb_err_invalid_input_

    else 
      isz = size(x%v,1)
      select case(dupl) 
      case(psb_dupl_ovwrt_) 
        do i = 1, n
          !loop over all val's rows

          ! row actual block row 
          if ((1 <= irl(i)).and.(irl(i) <= isz)) then
            ! this row belongs to me
            ! copy i-th row of block val in x
            x%v(irl(i),:) = val(i,:)
          end if
        enddo

      case(psb_dupl_add_) 

        do i = 1, n
          !loop over all val's rows
          if ((1 <= irl(i)).and.(irl(i) <= isz)) then
            ! this row belongs to me
            ! copy i-th row of block val in x
            x%v(irl(i),:) = x%v(irl(i),:) +  val(i,:)
          end if
        enddo

      case default
        info = 321
! !$      call psb_errpush(info,name)
! !$      goto 9999
      end select
    end if
    if (info /= 0) then 
      call psb_errpush(info,'base_mlv_vect_ins')
      return
    end if

  end subroutine i_base_mlv_ins

  !
  !> Function  base_mlv_zero
  !! \memberof  psb_i_base_multivect_type
  !! \brief Zero out contents
  !!
  !
  subroutine i_base_mlv_zero(x)
    use psi_serial_mod
    implicit none 
    class(psb_i_base_multivect_type), intent(inout)    :: x
    
    if (allocated(x%v)) x%v=izero

  end subroutine i_base_mlv_zero

  
  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_mlv_asb:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!           
  !!  \param n     final size
  !!  \param info  return code
  !!
  !
 
  subroutine i_base_mlv_asb(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    integer(psb_ipk_), intent(in)              :: m,n
    class(psb_i_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info
    
    if ((x%get_nrows() < m).or.(x%get_ncols()<n)) &
         & call psb_realloc(m,n,x%v,info)
    if (info /= 0) &
         & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')

  end subroutine i_base_mlv_asb


  !
  !> Function  base_mlv_free:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Free vector
  !!           
  !!  \param info  return code
  !!
  !
  subroutine i_base_mlv_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none 
    class(psb_i_base_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info
    
    info = 0
    if (allocated(x%v)) deallocate(x%v, stat=info)
    if (info /= 0) call & 
         & psb_errpush(psb_err_alloc_dealloc_,'vect_free')
        
  end subroutine i_base_mlv_free

  

  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  ! 
  !
  !> Function  base_mlv_sync:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Sync: base version is a no-op.
  !!           
  !
  subroutine i_base_mlv_sync(x)
    implicit none 
    class(psb_i_base_multivect_type), intent(inout) :: x
    
  end subroutine i_base_mlv_sync

  !
  !> Function  base_mlv_set_host:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Set_host: base version is a no-op.
  !!           
  !
  subroutine i_base_mlv_set_host(x)
    implicit none 
    class(psb_i_base_multivect_type), intent(inout) :: x
    
  end subroutine i_base_mlv_set_host

  !
  !> Function  base_mlv_set_dev:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Set_dev: base version is a no-op.
  !!           
  !
  subroutine i_base_mlv_set_dev(x)
    implicit none 
    class(psb_i_base_multivect_type), intent(inout) :: x
    
  end subroutine i_base_mlv_set_dev

  !
  !> Function  base_mlv_set_sync:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Set_sync: base version is a no-op.
  !!           
  !
  subroutine i_base_mlv_set_sync(x)
    implicit none 
    class(psb_i_base_multivect_type), intent(inout) :: x
    
  end subroutine i_base_mlv_set_sync

  !
  !> Function  base_mlv_is_dev:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Is  vector on external device    .
  !!           
  !
  function i_base_mlv_is_dev(x) result(res)
    implicit none 
    class(psb_i_base_multivect_type), intent(in) :: x
    logical  :: res
  
    res = .false.
  end function i_base_mlv_is_dev
  
  !
  !> Function  base_mlv_is_host
  !! \memberof  psb_i_base_multivect_type
  !! \brief Is  vector on standard memory    .
  !!           
  !
  function i_base_mlv_is_host(x) result(res)
    implicit none 
    class(psb_i_base_multivect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function i_base_mlv_is_host

  !
  !> Function  base_mlv_is_sync
  !! \memberof  psb_i_base_multivect_type
  !! \brief Is  vector on sync               .
  !!           
  !
  function i_base_mlv_is_sync(x) result(res)
    implicit none 
    class(psb_i_base_multivect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function i_base_mlv_is_sync


  !
  ! Size info. 
  !
  !
  !> Function  base_mlv_get_nrows
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Number of entries
  !!           
  !
  function i_base_mlv_get_nrows(x) result(res)
    implicit none 
    class(psb_i_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v,1)

  end function i_base_mlv_get_nrows

  function i_base_mlv_get_ncols(x) result(res)
    implicit none 
    class(psb_i_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v,2)

  end function i_base_mlv_get_ncols
  
  !
  !> Function  base_mlv_get_sizeof
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Size in bytesa
  !!           
  !
  function i_base_mlv_sizeof(x) result(res)
    implicit none 
    class(psb_i_base_multivect_type), intent(in) :: x
    integer(psb_long_int_k_) :: res
    
    ! Force 8-byte integers.
    res = (1_psb_long_int_k_ * psb_sizeof_int) * x%get_nrows() * x%get_ncols()

  end function i_base_mlv_sizeof

  !
  !> Function  base_mlv_get_fmt
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Format
  !!           
  !
  function i_base_mlv_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'BASE'
  end function i_base_mlv_get_fmt
  

  !
  !
  !
  !> Function  base_mlv_get_vect
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Extract a copy of the contents
  !!
  !    
  function  i_base_mlv_get_vect(x) result(res)
    class(psb_i_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), allocatable                 :: res(:,:)
    integer(psb_ipk_) :: info,m,n
    m = x%get_nrows()
    n = x%get_ncols()
    if (.not.allocated(x%v)) return 
    call x%sync()
    allocate(res(m,n),stat=info) 
    if (info /= 0) then 
      call psb_errpush(psb_err_alloc_dealloc_,'base_mlv_get_vect')
      return
    end if
    res(1:m,1:n) = x%v(1:m,1:n)
  end function i_base_mlv_get_vect
    
  !
  ! Reset all values 
  !
  !
  !> Function  base_mlv_set_scal
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  subroutine i_base_mlv_set_scal(x,val)
    class(psb_i_base_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in) :: val
        
    integer(psb_ipk_) :: info
    x%v = val
    
  end subroutine i_base_mlv_set_scal

  !
  !> Function  base_mlv_set_vect
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in 
  !!
  subroutine i_base_mlv_set_vect(x,val)
    class(psb_i_base_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in) :: val(:,:)
    integer(psb_ipk_) :: nr
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then 
      nr = min(size(x%v,1),size(val,1))
      nc = min(size(x%v,2),size(val,2))
      
      x%v(1:nr,1:nc) = val(1:nr,1:nc)
    else
      x%v = val
    end if

  end subroutine i_base_mlv_set_vect

end module psb_i_base_multivect_mod

