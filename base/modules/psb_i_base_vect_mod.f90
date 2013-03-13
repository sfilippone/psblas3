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
    procedure, pass(x) :: ins      => i_base_ins
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
    ! Basic info
    procedure, pass(x) :: get_nrows => i_base_get_nrows
    procedure, pass(x) :: sizeof   => i_base_sizeof
    !
    ! Set/get data from/to an external array; also
    ! overload assignment.
    !
    procedure, pass(x) :: get_vect => i_base_get_vect
    procedure, pass(x) :: set_scal => i_base_set_scal
    procedure, pass(x) :: set_vect => i_base_set_vect
    generic, public    :: set      => set_vect, set_scal

    !
    ! Dot product and AXPBY
    !
    procedure, pass(x) :: dot_v    => i_base_dot_v
    procedure, pass(x) :: dot_a    => i_base_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => i_base_axpby_v
    procedure, pass(y) :: axpby_a  => i_base_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    !
    ! Vector by vector multiplication. Need all variants
    ! to handle multiple requirements from preconditioners
    !
    procedure, pass(y) :: mlt_v    => i_base_mlt_v
    procedure, pass(y) :: mlt_a    => i_base_mlt_a
    procedure, pass(z) :: mlt_a_2  => i_base_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => i_base_mlt_v_2
    procedure, pass(z) :: mlt_va   => i_base_mlt_va
    procedure, pass(z) :: mlt_av   => i_base_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2, mlt_v_2, mlt_av, mlt_va
    !
    ! Scaling and norms
    !
    procedure, pass(x) :: scal     => i_base_scal
    procedure, pass(x) :: nrm2     => i_base_nrm2
    procedure, pass(x) :: amax     => i_base_amax
    procedure, pass(x) :: asum     => i_base_asum
    !
    ! Gather/scatter. These are needed for MPI interfacing.
    ! May have to be reworked. 
    !
    procedure, pass(x) :: gthab    => i_base_gthab
    procedure, pass(x) :: gthzv    => i_base_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => i_base_sctb
    generic, public    :: sct      => sctb
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
  subroutine i_base_ins(n,irl,val,dupl,x,info)
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
    if (info /= 0) then 
      call psb_errpush(info,'base_vect_ins')
      return
    end if

  end subroutine i_base_ins

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
    
    if (x%get_nrows() < n) &
         & call psb_realloc(n,x%v,info)
    if (info /= 0) &
         & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')

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
    call x%sync()
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
  subroutine i_base_set_scal(x,val)
    class(psb_i_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in) :: val
        
    integer(psb_ipk_) :: info
    x%v = val
    
  end subroutine i_base_set_scal

  !
  !> Function  base_set_vect
  !! \memberof  psb_i_base_vect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in 
  !!
  subroutine i_base_set_vect(x,val)
    class(psb_i_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in) :: val(:)
    integer(psb_ipk_) :: nr
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then 
      nr = min(size(x%v),size(val))
      x%v(1:nr) = val(1:nr)
    else
      x%v = val
    end if

  end subroutine i_base_set_vect

  !
  ! Dot products 
  ! 
  !
  !> Function  base_dot_v
  !! \memberof  psb_i_base_vect_type
  !! \brief  Dot product by another base_vector
  !! \param n    Number of entries to be considere
  !! \param y    The other (base_vect) to be multiplied by
  !!
  function i_base_dot_v(n,x,y) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(in)           :: n
    integer(psb_ipk_)                :: res
    integer(psb_ipk_), external      :: idot
    
    res = izero
    !
    ! Note: this is the base implementation.
    !  When we get here, we are sure that X is of
    !  TYPE psb_i_base_vect.
    !  If Y is not, throw the burden on it, implicitly
    !  calling dot_a
    !
    select type(yy => y)
    type is (psb_i_base_vect_type)
      res = idot(n,x%v,1,y%v,1)
    class default
      res = y%dot(n,x%v)
    end select

  end function i_base_dot_v

  !
  ! Base workhorse is good old BLAS1
  !
  !
  !> Function  base_dot_a
  !! \memberof  psb_i_base_vect_type
  !! \brief  Dot product by a normal array
  !! \param n    Number of entries to be considere
  !! \param y(:) The array to be multiplied by
  !!
  function i_base_dot_a(n,x,y) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)    :: y(:)
    integer(psb_ipk_), intent(in)           :: n
    integer(psb_ipk_)                :: res
    integer(psb_ipk_), external      :: idot
    
    res = idot(n,y,1,x%v,1)

  end function i_base_dot_a
    
  !
  ! AXPBY is invoked via Y, hence the structure below. 
  !
  !
  !
  !> Function  base_axpby_v
  !! \memberof  psb_i_base_vect_type
  !! \brief AXPBY  by a (base_vect) y=alpha*x+beta*y
  !! \param m    Number of entries to be considere
  !! \param alpha scalar alpha
  !! \param x     The class(base_vect) to be added
  !! \param beta scalar alpha
  !! \param info   return code
  !!
  subroutine i_base_axpby_v(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)               :: m
    class(psb_i_base_vect_type), intent(inout)  :: x
    class(psb_i_base_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info
    
    select type(xx => x)
    type is (psb_i_base_vect_type)
      call psb_geaxpby(m,alpha,x%v,beta,y%v,info)
    class default
      call y%axpby(m,alpha,x%v,beta,info)
    end select

  end subroutine i_base_axpby_v

  !
  ! AXPBY is invoked via Y, hence the structure below. 
  !
  !
  !> Function  base_axpby_a
  !! \memberof  psb_i_base_vect_type
  !! \brief AXPBY  by a normal array y=alpha*x+beta*y
  !! \param m    Number of entries to be considere
  !! \param alpha scalar alpha
  !! \param x(:) The array to be added
  !! \param beta scalar alpha
  !! \param info   return code
  !!
  subroutine i_base_axpby_a(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)               :: m
    integer(psb_ipk_), intent(in)        :: x(:)
    class(psb_i_base_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info
    
    call psb_geaxpby(m,alpha,x,beta,y%v,info)
    
  end subroutine i_base_axpby_a

  
  !
  !  Multiple variants of two operations:
  !  Simple multiplication  Y(:) = X(:)*Y(:)
  !  blas-like:   Z(:) = alpha*X(:)*Y(:)+beta*Z(:)
  !
  !  Variants expanded according to the dynamic type
  !  of the involved entities
  !
  !
  !> Function  base_mlt_a
  !! \memberof  psb_i_base_vect_type
  !! \brief Vector entry-by-entry multiply  by a base_vect array y=x*y
  !! \param x   The class(base_vect) to be multiplied by
  !! \param info   return code
  !!
  subroutine i_base_mlt_v(x, y, info)
    use psi_serial_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout)  :: x
    class(psb_i_base_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info    
    integer(psb_ipk_) :: i, n

    info = 0
    select type(xx => x)
    type is (psb_i_base_vect_type)
      n = min(size(y%v), size(xx%v))
      do i=1, n 
        y%v(i) = y%v(i)*xx%v(i)
      end do
    class default
      call y%mlt(x%v,info)
    end select

  end subroutine i_base_mlt_v

  !
  !> Function  base_mlt_a
  !! \memberof  psb_i_base_vect_type
  !! \brief Vector entry-by-entry multiply  by a normal array y=x*y
  !! \param x(:) The array to be multiplied by
  !! \param info   return code
  !!
  subroutine i_base_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)        :: x(:)
    class(psb_i_base_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    n = min(size(y%v), size(x))
    do i=1, n 
      y%v(i) = y%v(i)*x(i)
    end do
    
  end subroutine i_base_mlt_a


  !
  !> Function  base_mlt_a_2
  !! \memberof  psb_i_base_vect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by normal arrays
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x(:) The array to be multiplied b
  !! \param y(:) The array to be multiplied by
  !! \param info   return code
  !!
  subroutine i_base_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)        :: alpha,beta
    integer(psb_ipk_), intent(in)        :: y(:)
    integer(psb_ipk_), intent(in)        :: x(:)
    class(psb_i_base_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0    
    n = min(size(z%v), size(x), size(y))
!!$    write(0,*) 'Mlt_a_2: ',n
    if (alpha == izero) then 
      if (beta == ione) then 
        return 
      else
        do i=1, n
          z%v(i) = beta*z%v(i)
        end do
      end if
    else
      if (alpha == ione) then 
        if (beta == izero) then 
          do i=1, n 
            z%v(i) = y(i)*x(i)
          end do
        else if (beta == ione) then 
          do i=1, n 
            z%v(i) = z%v(i) + y(i)*x(i)
          end do
        else 
          do i=1, n 
            z%v(i) = beta*z%v(i) + y(i)*x(i)
          end do
        end if
      else if (alpha == -ione) then 
        if (beta == izero) then 
          do i=1, n 
            z%v(i) = -y(i)*x(i)
          end do
        else if (beta == ione) then 
          do i=1, n 
            z%v(i) = z%v(i) - y(i)*x(i)
          end do
        else 
          do i=1, n 
            z%v(i) = beta*z%v(i) - y(i)*x(i)
          end do
        end if
      else
        if (beta == izero) then 
          do i=1, n 
            z%v(i) = alpha*y(i)*x(i)
          end do
        else if (beta == ione) then 
          do i=1, n 
            z%v(i) = z%v(i) + alpha*y(i)*x(i)
          end do
        else 
          do i=1, n 
            z%v(i) = beta*z%v(i) + alpha*y(i)*x(i)
          end do
        end if
      end if
    end if
  end subroutine i_base_mlt_a_2

  !
  !> Function  base_mlt_v_2
  !! \memberof  psb_i_base_vect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by class(base_vect)
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x  The class(base_vect) to be multiplied b
  !! \param y  The class(base_vect) to be multiplied by
  !! \param info   return code
  !!
  subroutine i_base_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
    use psi_serial_mod
    use psb_string_mod
    implicit none 
    integer(psb_ipk_), intent(in)        :: alpha,beta
    class(psb_i_base_vect_type), intent(inout)  :: x
    class(psb_i_base_vect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info    
    character(len=1), intent(in), optional     :: conjgx, conjgy
    integer(psb_ipk_) :: i, n
    logical :: conjgx_, conjgy_

    info = 0
    if (.not.psb_i_is_complex_) then
      call z%mlt(alpha,x%v,y%v,beta,info)
    else 
      conjgx_=.false.
      if (present(conjgx)) conjgx_ = (psb_toupper(conjgx)=='C')
      conjgy_=.false.
      if (present(conjgy)) conjgy_ = (psb_toupper(conjgy)=='C')
      if (conjgx_) x%v=(x%v)
      if (conjgy_) y%v=(y%v)
      call z%mlt(alpha,x%v,y%v,beta,info)
      if (conjgx_) x%v=(x%v)
      if (conjgy_) y%v=(y%v)
    end if
  end subroutine i_base_mlt_v_2

  subroutine i_base_mlt_av(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)        :: alpha,beta
    integer(psb_ipk_), intent(in)        :: x(:)
    class(psb_i_base_vect_type), intent(inout)  :: y
    class(psb_i_base_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info    
    integer(psb_ipk_) :: i, n

    info = 0
    
    call z%mlt(alpha,x,y%v,beta,info)

  end subroutine i_base_mlt_av

  subroutine i_base_mlt_va(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none 
    integer(psb_ipk_), intent(in)        :: alpha,beta
    integer(psb_ipk_), intent(in)        :: y(:)
    class(psb_i_base_vect_type), intent(inout)  :: x
    class(psb_i_base_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info    
    integer(psb_ipk_) :: i, n

    info = 0
    
    call z%mlt(alpha,y,x,beta,info)

  end subroutine i_base_mlt_va


  !
  ! Simple scaling 
  !
  !> Function  base_scal
  !! \memberof  psb_i_base_vect_type
  !! \brief Scale all entries  x = alpha*x
  !! \param alpha   The multiplier
  !!
  subroutine i_base_scal(alpha, x)
    use psi_serial_mod
    implicit none 
    class(psb_i_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent (in)       :: alpha
    
    if (allocated(x%v)) x%v = alpha*x%v

  end subroutine i_base_scal
  
  !
  ! Norms 1, 2 and infinity
  !
  !> Function  base_nrm2
  !! \memberof  psb_i_base_vect_type
  !! \brief 2-norm |x(1:n)|_2
  !! \param n  how many entries to consider
  function i_base_nrm2(n,x) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    integer(psb_ipk_)                :: res
    integer(psb_ipk_), external      :: inrm2
    
    res =  inrm2(n,x%v,1)

  end function i_base_nrm2
  
  !
  !> Function  base_amax
  !! \memberof  psb_i_base_vect_type
  !! \brief infinity-norm |x(1:n)|_\infty
  !! \param n  how many entries to consider
  function i_base_amax(n,x) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    integer(psb_ipk_)                :: res
    
    res =  maxval(abs(x%v(1:n)))

  end function i_base_amax

  !
  !> Function  base_asum
  !! \memberof  psb_i_base_vect_type
  !! \brief 1-norm |x(1:n)|_1
  !! \param n  how many entries to consider
  function i_base_asum(n,x) result(res)
    implicit none 
    class(psb_i_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    integer(psb_ipk_)                :: res
    
    res =  sum(abs(x%v(1:n)))

  end function i_base_asum
  
  
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
    
    call x%sync()
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
  subroutine i_base_gthzv(n,idx,x,y)
    use psi_serial_mod
    integer(psb_ipk_) :: n, idx(:)
    integer(psb_ipk_) ::  y(:)
    class(psb_i_base_vect_type) :: x
    
    call x%sync()
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
    
    call y%sync()
    call psi_sct(n,idx,x,beta,y%v)
    call y%set_host()

  end subroutine i_base_sctb

end module psb_i_base_vect_mod
