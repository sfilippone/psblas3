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
!
! package: psb_d_base_vect_mod
!
! This module contains the definition of the psb_d_base_vect type which
! is a container for dense vectors.
!  This is encapsulated instead of being just a simple array to allow for
!  more complicated situations, such as GPU programming, where the memory
!  area we are interested in is not easily accessible from the host/Fortran
!  side. It is also meant to be encapsulated in an outer type, to allow
!  runtime switching as per the STATE design pattern, similar to the
!  sparse matrix types.
!
!
module psb_d_base_vect_mod

  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_i_base_vect_mod
  use psb_l_base_vect_mod

  !> \namespace  psb_base_mod  \class psb_d_base_vect_type
  !! The psb_d_base_vect_type
  !! defines a middle level  real(psb_dpk_) encapsulated dense vector.
  !! The encapsulation is needed, in place of a simple array, to allow
  !! for complicated situations, such as GPU programming, where the memory
  !!  area we are interested in is not easily accessible from the host/Fortran
  !!  side. It is also meant to be encapsulated in an outer type, to allow
  !!  runtime switching as per the STATE design pattern, similar to the
  !!  sparse matrix types.
  !!
  type psb_d_persis_vect_type
    integer, allocatable :: init_request
    real(psb_dpk_), allocatable :: sndbuf(:)
    real(psb_dpk_), allocatable :: rcvbuf(:)
    ! required for packing and unpacking when using persistent communication
    integer(psb_ipk_), allocatable :: snd_count
    integer(psb_ipk_), allocatable :: rcv_count
    integer(psb_ipk_), allocatable :: snd_counts(:)
    integer(psb_ipk_), allocatable :: rcv_counts(:)
    integer(psb_ipk_), allocatable :: snd_ws(:)
    integer(psb_ipk_), allocatable :: rcv_ws(:)
    integer(psb_ipk_), allocatable :: snd_to(:)
    integer(psb_ipk_), allocatable :: rcv_from(:)
    integer(psb_ipk_), allocatable :: snd_displs(:) ! snd array displacements
    integer(psb_ipk_), allocatable :: rcv_displs(:) ! rcv array displacements
 end type psb_d_persis_vect_type

  type psb_d_base_vect_type
    !> Values.
    real(psb_dpk_), allocatable :: v(:)
    real(psb_dpk_), allocatable :: combuf(:)
    integer(psb_mpk_), allocatable :: comid(:,:)
    type(psb_d_persis_vect_type), allocatable  :: p
  contains
    !
    !  Constructors/allocators
    !
    procedure, pass(x) :: bld_x    => d_base_bld_x
    procedure, pass(x) :: bld_mn   => d_base_bld_mn
    procedure, pass(x) :: bld_en   => d_base_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
    procedure, pass(x) :: all      => d_base_all
    procedure, pass(x) :: mold     => d_base_mold
    !
    ! Insert/set. Assembly and free.
    ! Assembly does almost nothing here, but is important
    ! in derived classes.
    !
    procedure, pass(x) :: ins_a    => d_base_ins_a
    procedure, pass(x) :: ins_v    => d_base_ins_v
    generic, public    :: ins      => ins_a, ins_v
    procedure, pass(x) :: zero     => d_base_zero
    procedure, pass(x) :: asb_m    => d_base_asb_m
    procedure, pass(x) :: asb_e    => d_base_asb_e
    generic, public    :: asb      => asb_m, asb_e
    procedure, pass(x) :: free     => d_base_free
    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder.
    !
    procedure, pass(x) :: sync     => d_base_sync
    procedure, pass(x) :: is_host  => d_base_is_host
    procedure, pass(x) :: is_dev   => d_base_is_dev
    procedure, pass(x) :: is_sync  => d_base_is_sync
    procedure, pass(x) :: set_host => d_base_set_host
    procedure, pass(x) :: set_dev  => d_base_set_dev
    procedure, pass(x) :: set_sync => d_base_set_sync

    !
    ! These are for handling gather/scatter in new
    ! comm internals implementation.
    !
    procedure, nopass  :: use_buffer   => d_base_use_buffer
    procedure, pass(x) :: new_buffer   => d_base_new_buffer
    procedure, nopass  :: device_wait  => d_base_device_wait
    procedure, pass(x) :: maybe_free_buffer  => d_base_maybe_free_buffer
    procedure, pass(x) :: free_buffer  => d_base_free_buffer
    procedure, pass(x) :: new_comid    => d_base_new_comid
    procedure, pass(x) :: free_comid   => d_base_free_comid

    !
    ! Basic info
    procedure, pass(x) :: get_nrows => d_base_get_nrows
    procedure, pass(x) :: sizeof    => d_base_sizeof
    procedure, nopass  :: get_fmt   => d_base_get_fmt
    !
    ! Set/get data from/to an external array; also
    ! overload assignment.
    !
    procedure, pass(x) :: get_vect => d_base_get_vect
    procedure, pass(x) :: set_scal => d_base_set_scal
    procedure, pass(x) :: set_vect => d_base_set_vect
    generic, public    :: set      => set_vect, set_scal
    !
    ! Gather/scatter. These are needed for MPI interfacing.
    ! May have to be reworked.
    !
    procedure, pass(x) :: gthab    => d_base_gthab
    procedure, pass(x) :: gthzv    => d_base_gthzv
    procedure, pass(x) :: gthzv_x  => d_base_gthzv_x
    procedure, pass(x) :: gthzbuf  => d_base_gthzbuf
    generic, public    :: gth      => gthab, gthzv, gthzv_x, gthzbuf
    procedure, pass(y) :: sctb     => d_base_sctb
    procedure, pass(y) :: sctb_x   => d_base_sctb_x
    procedure, pass(y) :: sctb_buf => d_base_sctb_buf
    generic, public    :: sct      => sctb, sctb_x, sctb_buf


    !
    ! Dot product and AXPBY
    !
    procedure, pass(x) :: dot_v    => d_base_dot_v
    procedure, pass(x) :: dot_a    => d_base_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => d_base_axpby_v
    procedure, pass(y) :: axpby_a  => d_base_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    !
    ! Vector by vector multiplication. Need all variants
    ! to handle multiple requirements from preconditioners
    !
    procedure, pass(y) :: mlt_v    => d_base_mlt_v
    procedure, pass(y) :: mlt_a    => d_base_mlt_a
    procedure, pass(z) :: mlt_a_2  => d_base_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => d_base_mlt_v_2
    procedure, pass(z) :: mlt_va   => d_base_mlt_va
    procedure, pass(z) :: mlt_av   => d_base_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2, mlt_v_2, mlt_av, mlt_va
    !
    ! Scaling and norms
    !
    procedure, pass(x) :: scal     => d_base_scal
    procedure, pass(x) :: absval1  => d_base_absval1
    procedure, pass(x) :: absval2  => d_base_absval2
    generic, public    :: absval   => absval1, absval2
    procedure, pass(x) :: nrm2     => d_base_nrm2
    procedure, pass(x) :: amax     => d_base_amax
    procedure, pass(x) :: asum     => d_base_asum

  end type psb_d_base_vect_type

  public  :: psb_d_base_vect
  private :: constructor, size_const
  interface psb_d_base_vect
    module procedure constructor, size_const
  end interface psb_d_base_vect

contains

  !
  ! Constructors.
  !

  !> Function  constructor:
  !! \brief     Constructor from an array
  !!  \param   x(:)  input array to be copied
  !!
  function constructor(x) result(this)
    real(psb_dpk_)   :: x(:)
    type(psb_d_base_vect_type) :: this
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
    type(psb_d_base_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%asb(n,info)

  end function size_const

  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_d_base_vect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  subroutine d_base_bld_x(x,this)
    use psb_realloc_mod
    implicit none
    real(psb_dpk_), intent(in) :: this(:)
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this),x%v,info)
    if (info /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_,'base_vect_bld')
      return
    end if
    x%v(:)  = this(:)

  end subroutine d_base_bld_x

  !
  ! Create with size, but no initialization
  !

  !> Function  bld_mn:
  !! \memberof  psb_d_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  subroutine d_base_bld_mn(x,n)
    use psb_realloc_mod
    implicit none
    integer(psb_mpk_), intent(in) :: n
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(n,x%v,info)
    call x%asb(n,info)

  end subroutine d_base_bld_mn

  !> Function  bld_en:
  !! \memberof  psb_d_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  subroutine d_base_bld_en(x,n)
    use psb_realloc_mod
    implicit none
    integer(psb_epk_), intent(in) :: n
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(n,x%v,info)
    call x%asb(n,info)

  end subroutine d_base_bld_en

  !> Function  base_all:
  !! \memberof  psb_d_base_vect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated.
  !!  \param info  return code
  !!
  subroutine d_base_all(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: n
    class(psb_d_base_vect_type), intent(out)    :: x
    integer(psb_ipk_), intent(out)              :: info

    call psb_realloc(n,x%v,info)

  end subroutine d_base_all

  !> Function  base_mold:
  !! \memberof  psb_d_base_vect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  subroutine d_base_mold(x, y, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_d_base_vect_type), intent(in)   :: x
    class(psb_d_base_vect_type), intent(out), allocatable :: y
    integer(psb_ipk_), intent(out)              :: info

    allocate(psb_d_base_vect_type :: y, stat=info)

  end subroutine d_base_mold

  !
  ! Insert a bunch of values at specified positions.
  !
  !> Function  base_ins:
  !! \memberof  psb_d_base_vect_type
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
  subroutine d_base_ins_a(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    real(psb_dpk_), intent(in)        :: val(:)
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

  end subroutine d_base_ins_a

  subroutine d_base_ins_v(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    class(psb_i_base_vect_type), intent(inout)  :: irl
    class(psb_d_base_vect_type), intent(inout)  :: val
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

  end subroutine d_base_ins_v


  !
  !> Function  base_zero
  !! \memberof  psb_d_base_vect_type
  !! \brief Zero out contents
  !!
  !
  subroutine d_base_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout)    :: x

    if (allocated(x%v)) x%v=dzero
    call x%set_host()
  end subroutine d_base_zero


  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_asb:
  !! \memberof  psb_d_base_vect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !

  subroutine d_base_asb_m(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_mpk_), intent(in)              :: n
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    info = 0
    if (x%get_nrows() < n) &
         & call psb_realloc(n,x%v,info)
    if (info /= 0) &
         & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')
    call x%sync()
  end subroutine d_base_asb_m

  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_asb:
  !! \memberof  psb_d_base_vect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !

  subroutine d_base_asb_e(n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_epk_), intent(in)              :: n
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    info = 0
    if (x%get_nrows() < n) &
         & call psb_realloc(n,x%v,info)
    if (info /= 0) &
         & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')
    call x%sync()
  end subroutine d_base_asb_e

  !
  !> Function  base_free:
  !! \memberof  psb_d_base_vect_type
  !! \brief Free vector
  !!
  !!  \param info  return code
  !!
  !
  subroutine d_base_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) deallocate(x%v, stat=info)
    if (info == 0) call x%free_buffer(info)
    if (info == 0) call x%free_comid(info)
    if (info /= 0) call &
         & psb_errpush(psb_err_alloc_dealloc_,'vect_free')

  end subroutine d_base_free



  !
  !> Function  base_free_buffer:
  !! \memberof  psb_d_base_vect_type
  !! \brief Free aux buffer
  !!
  !!  \param info  return code
  !!
  !
  subroutine d_base_free_buffer(x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%combuf)) &
         &  deallocate(x%combuf,stat=info)
  end subroutine d_base_free_buffer

  !
  !> Function  base_maybe_free_buffer:
  !! \memberof  psb_d_base_vect_type
  !! \brief Conditionally Free aux buffer.
  !!        In some derived classes, e.g. GPU,
  !!        does not really frees to avoid  runtime
  !!        costs
  !!
  !!  \param info  return code
  !!
  !
  subroutine d_base_maybe_free_buffer(x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    info = 0
    if (psb_get_maybe_free_buffer())&
         &  call x%free_buffer(info)

  end subroutine d_base_maybe_free_buffer

  !
  !> Function  base_free_comid:
  !! \memberof  psb_d_base_vect_type
  !! \brief Free aux MPI communication id buffer
  !!
  !!  \param info  return code
  !!
  !
  subroutine d_base_free_comid(x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%comid)) &
         &  deallocate(x%comid,stat=info)
  end subroutine d_base_free_comid


  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  !
  !
  !> Function  base_sync:
  !! \memberof  psb_d_base_vect_type
  !! \brief Sync: base version is a no-op.
  !!
  !
  subroutine d_base_sync(x)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x

  end subroutine d_base_sync

  !
  !> Function  base_set_host:
  !! \memberof  psb_d_base_vect_type
  !! \brief Set_host: base version is a no-op.
  !!
  !
  subroutine d_base_set_host(x)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x

  end subroutine d_base_set_host

  !
  !> Function  base_set_dev:
  !! \memberof  psb_d_base_vect_type
  !! \brief Set_dev: base version is a no-op.
  !!
  !
  subroutine d_base_set_dev(x)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x

  end subroutine d_base_set_dev

  !
  !> Function  base_set_sync:
  !! \memberof  psb_d_base_vect_type
  !! \brief Set_sync: base version is a no-op.
  !!
  !
  subroutine d_base_set_sync(x)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x

  end subroutine d_base_set_sync

  !
  !> Function  base_is_dev:
  !! \memberof  psb_d_base_vect_type
  !! \brief Is  vector on external device    .
  !!
  !
  function d_base_is_dev(x) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(in) :: x
    logical  :: res

    res = .false.
  end function d_base_is_dev

  !
  !> Function  base_is_host
  !! \memberof  psb_d_base_vect_type
  !! \brief Is  vector on standard memory    .
  !!
  !
  function d_base_is_host(x) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function d_base_is_host

  !
  !> Function  base_is_sync
  !! \memberof  psb_d_base_vect_type
  !! \brief Is  vector on sync               .
  !!
  !
  function d_base_is_sync(x) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function d_base_is_sync


  !
  ! Size info.
  !
  !
  !> Function  base_get_nrows
  !! \memberof  psb_d_base_vect_type
  !! \brief  Number of entries
  !!
  !
  function d_base_get_nrows(x) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v)

  end function d_base_get_nrows

  !
  !> Function  base_get_sizeof
  !! \memberof  psb_d_base_vect_type
  !! \brief  Size in bytes
  !!
  !
  function d_base_sizeof(x) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(in) :: x
    integer(psb_epk_) :: res

    ! Force 8-byte integers.
    res = (1_psb_epk_ * psb_sizeof_dp) * x%get_nrows()

  end function d_base_sizeof

  !
  !> Function  base_get_fmt
  !! \memberof  psb_d_base_vect_type
  !! \brief  Format
  !!
  !
  function d_base_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'BASE'
  end function d_base_get_fmt


  !
  !
  !
  !> Function  base_get_vect
  !! \memberof  psb_d_base_vect_type
  !! \brief  Extract a copy of the contents
  !!
  !
  function  d_base_get_vect(x,n) result(res)
    class(psb_d_base_vect_type), intent(inout) :: x
    real(psb_dpk_), allocatable                 :: res(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional :: n
    ! Local variables
    integer(psb_ipk_) :: isz

    if (.not.allocated(x%v)) return
    if (.not.x%is_host()) call x%sync()
    isz = x%get_nrows()
    if (present(n)) isz = max(0,min(isz,n))
    allocate(res(isz),stat=info)
    if (info /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_,'base_get_vect')
      return
    end if
    res(1:isz) = x%v(1:isz)
  end function d_base_get_vect

  !
  ! Reset all values
  !
  !
  !> Function  base_set_scal
  !! \memberof  psb_d_base_vect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  subroutine d_base_set_scal(x,val,first,last)
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x
    real(psb_dpk_), intent(in) :: val
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info, first_, last_

    first_=1
    last_=size(x%v)
    if (present(first)) first_ = max(1,first)
    if (present(last))  last_  = min(last,last_)

    if (x%is_dev()) call x%sync()
    x%v(first_:last_) = val
    call x%set_host()

  end subroutine d_base_set_scal


  !
  !> Function  base_set_vect
  !! \memberof  psb_d_base_vect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in
  !!
  subroutine d_base_set_vect(x,val,first,last)
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x
    real(psb_dpk_), intent(in) :: val(:)
    integer(psb_ipk_), optional :: first, last

    integer(psb_ipk_) :: info, first_, last_, nr

    first_                     = 1
    if (present(first)) first_ = max(1,first)
    last_                      = min(psb_size(x%v),first_+size(val)-1)
    if (present(last))  last_  = min(last,last_)

    if (allocated(x%v)) then
      if (x%is_dev()) call x%sync()
      x%v(first_:last_) = val(1:last_-first_+1)
    else
      x%v = val
    end if
    call x%set_host()

  end subroutine d_base_set_vect


  !
  ! Overwrite with absolute value
  !
  !
  !> Function  base_absval1
  !! \memberof  psb_d_base_vect_type
  !! \brief  Set all entries to their respective absolute values.
  !!
  subroutine d_base_absval1(x)
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x

    if (allocated(x%v)) then
      if (x%is_dev()) call x%sync()
      x%v =  abs(x%v)
      call x%set_host()
    end if

  end subroutine d_base_absval1

  subroutine d_base_absval2(x,y)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    class(psb_d_base_vect_type), intent(inout) :: y
    integer(psb_ipk_) :: info
    if (.not.x%is_host()) call x%sync()
    if (allocated(x%v)) then
      call y%axpby(ione*min(x%get_nrows(),y%get_nrows()),done,x,dzero,info)
      call y%absval()
    end if

  end subroutine d_base_absval2

  !
  ! Dot products
  !
  !
  !> Function  base_dot_v
  !! \memberof  psb_d_base_vect_type
  !! \brief  Dot product by another base_vector
  !! \param n    Number of entries to be considered
  !! \param y    The other (base_vect) to be multiplied by
  !!
  function d_base_dot_v(n,x,y) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_)                :: res
    real(psb_dpk_), external      :: ddot

    res = dzero
    !
    ! Note: this is the base implementation.
    !  When we get here, we are sure that X is of
    !  TYPE psb_d_base_vect.
    !  If Y is not, throw the burden on it, implicitly
    !  calling dot_a
    !
    select type(yy => y)
    type is (psb_d_base_vect_type)
      res = ddot(n,x%v,1,y%v,1)
    class default
      res = y%dot(n,x%v)
    end select

  end function d_base_dot_v

  !
  ! Base workhorse is good old BLAS1
  !
  !
  !> Function  base_dot_a
  !! \memberof  psb_d_base_vect_type
  !! \brief  Dot product by a normal array
  !! \param n    Number of entries to be considered
  !! \param y(:) The array to be multiplied by
  !!
  function d_base_dot_a(n,x,y) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    real(psb_dpk_), intent(in)    :: y(:)
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_)                :: res
    real(psb_dpk_), external      :: ddot

    res = ddot(n,y,1,x%v,1)

  end function d_base_dot_a

  !
  ! AXPBY is invoked via Y, hence the structure below.
  !
  !
  !
  !> Function  base_axpby_v
  !! \memberof  psb_d_base_vect_type
  !! \brief AXPBY  by a (base_vect) y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x     The class(base_vect) to be added
  !! \param beta scalar alpha
  !! \param info   return code
  !!
  subroutine d_base_axpby_v(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_d_base_vect_type), intent(inout)  :: x
    class(psb_d_base_vect_type), intent(inout)  :: y
    real(psb_dpk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info

    if (x%is_dev()) call x%sync()

    call y%axpby(m,alpha,x%v,beta,info)

  end subroutine d_base_axpby_v

  !
  ! AXPBY is invoked via Y, hence the structure below.
  !
  !
  !> Function  base_axpby_a
  !! \memberof  psb_d_base_vect_type
  !! \brief AXPBY  by a normal array y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x(:) The array to be added
  !! \param beta scalar alpha
  !! \param info   return code
  !!
  subroutine d_base_axpby_a(m,alpha, x, beta, y, info)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_vect_type), intent(inout)  :: y
    real(psb_dpk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info

    if (y%is_dev()) call y%sync()
    call psb_geaxpby(m,alpha,x,beta,y%v,info)
    call y%set_host()

  end subroutine d_base_axpby_a


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
  !! \memberof  psb_d_base_vect_type
  !! \brief Vector entry-by-entry multiply  by a base_vect array y=x*y
  !! \param x   The class(base_vect) to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlt_v(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x
    class(psb_d_base_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (x%is_dev()) call x%sync()
    call y%mlt(x%v,info)

  end subroutine d_base_mlt_v

  !
  !> Function  base_mlt_a
  !! \memberof  psb_d_base_vect_type
  !! \brief Vector entry-by-entry multiply  by a normal array y=x*y
  !! \param x(:) The array to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlt_a(x, y, info)
    use psi_serial_mod
    implicit none
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_vect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (y%is_dev()) call y%sync()
    n = min(size(y%v), size(x))
    do i=1, n
      y%v(i) = y%v(i)*x(i)
    end do
    call y%set_host()

  end subroutine d_base_mlt_a


  !
  !> Function  base_mlt_a_2
  !! \memberof  psb_d_base_vect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by normal arrays
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x(:) The array to be multiplied b
  !! \param y(:) The array to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: y(:)
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (z%is_dev()) call z%sync()

    n = min(size(z%v), size(x), size(y))
    if (alpha == dzero) then
      if (beta == done) then
        return
      else
        do i=1, n
          z%v(i) = beta*z%v(i)
        end do
      end if
    else
      if (alpha == done) then
        if (beta == dzero) then
          do i=1, n
            z%v(i) = y(i)*x(i)
          end do
        else if (beta == done) then
          do i=1, n
            z%v(i) = z%v(i) + y(i)*x(i)
          end do
        else
          do i=1, n
            z%v(i) = beta*z%v(i) + y(i)*x(i)
          end do
        end if
      else if (alpha == -done) then
        if (beta == dzero) then
          do i=1, n
            z%v(i) = -y(i)*x(i)
          end do
        else if (beta == done) then
          do i=1, n
            z%v(i) = z%v(i) - y(i)*x(i)
          end do
        else
          do i=1, n
            z%v(i) = beta*z%v(i) - y(i)*x(i)
          end do
        end if
      else
        if (beta == dzero) then
          do i=1, n
            z%v(i) = alpha*y(i)*x(i)
          end do
        else if (beta == done) then
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
    call z%set_host()

  end subroutine d_base_mlt_a_2

  !
  !> Function  base_mlt_v_2
  !! \memberof  psb_d_base_vect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by class(base_vect)
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x  The class(base_vect) to be multiplied b
  !! \param y  The class(base_vect) to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
    use psi_serial_mod
    use psb_string_mod
    implicit none
    real(psb_dpk_), intent(in)        :: alpha,beta
    class(psb_d_base_vect_type), intent(inout)  :: x
    class(psb_d_base_vect_type), intent(inout)  :: y
    class(psb_d_base_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional     :: conjgx, conjgy
    integer(psb_ipk_) :: i, n
    logical :: conjgx_, conjgy_

    info = 0
    if (y%is_dev()) call y%sync()
    if (x%is_dev()) call x%sync()
    if (.not.psb_d_is_complex_) then
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
  end subroutine d_base_mlt_v_2

  subroutine d_base_mlt_av(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_vect_type), intent(inout)  :: y
    class(psb_d_base_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (y%is_dev()) call y%sync()
    call z%mlt(alpha,x,y%v,beta,info)

  end subroutine d_base_mlt_av

  subroutine d_base_mlt_va(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: y(:)
    class(psb_d_base_vect_type), intent(inout)  :: x
    class(psb_d_base_vect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    if (x%is_dev()) call x%sync()
    call z%mlt(alpha,y,x,beta,info)

  end subroutine d_base_mlt_va


  !
  ! Simple scaling
  !
  !> Function  base_scal
  !! \memberof  psb_d_base_vect_type
  !! \brief Scale all entries  x = alpha*x
  !! \param alpha   The multiplier
  !!
  subroutine d_base_scal(alpha, x)
    use psi_serial_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x
    real(psb_dpk_), intent (in)       :: alpha

    if (allocated(x%v)) then
      x%v = alpha*x%v
      call x%set_host()
    end if

  end subroutine d_base_scal

  !
  ! Norms 1, 2 and infinity
  !
  !> Function  base_nrm2
  !! \memberof  psb_d_base_vect_type
  !! \brief 2-norm |x(1:n)|_2
  !! \param n  how many entries to consider
  function d_base_nrm2(n,x) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_)                :: res
    real(psb_dpk_), external      :: dnrm2

    if (x%is_dev()) call x%sync()
    res =  dnrm2(n,x%v,1)

  end function d_base_nrm2

  !
  !> Function  base_amax
  !! \memberof  psb_d_base_vect_type
  !! \brief infinity-norm |x(1:n)|_\infty
  !! \param n  how many entries to consider
  function d_base_amax(n,x) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_)                :: res

    if (x%is_dev()) call x%sync()
    res =  maxval(abs(x%v(1:n)))

  end function d_base_amax

  !
  !> Function  base_asum
  !! \memberof  psb_d_base_vect_type
  !! \brief 1-norm |x(1:n)|_1
  !! \param n  how many entries to consider
  function d_base_asum(n,x) result(res)
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_)                :: res

    if (x%is_dev()) call x%sync()
    res =  sum(abs(x%v(1:n)))

  end function d_base_asum


  !
  ! Gather: Y = beta * Y + alpha * X(IDX(:))
  !
  !
  !> Function  base_gthab
  !! \memberof  psb_d_base_vect_type
  !! \brief gather into an array
  !!    Y = beta * Y + alpha * X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param alpha
  !! \param beta
  subroutine d_base_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: n, idx(:)
    real(psb_dpk_) :: alpha, beta, y(:)
    class(psb_d_base_vect_type) :: x

    if (x%is_dev()) call x%sync()
    call psi_gth(n,idx,alpha,x%v,beta,y)

  end subroutine d_base_gthab
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_gthzv
  !! \memberof  psb_d_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  subroutine d_base_gthzv_x(i,n,idx,x,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    real(psb_dpk_) ::  y(:)
    class(psb_d_base_vect_type) :: x

    if (idx%is_dev()) call idx%sync()
    call x%gth(n,idx%v(i:),y)

  end subroutine d_base_gthzv_x

  !
  ! New comm internals impl.
  !
  subroutine d_base_gthzbuf(i,n,idx,x)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    class(psb_d_base_vect_type) :: x

    if (.not.allocated(x%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'gthzbuf')
      return
    end if
    if (idx%is_dev()) call idx%sync()
    if (x%is_dev()) call x%sync()
    call x%gth(n,idx%v(i:),x%combuf(i:))

  end subroutine d_base_gthzbuf
  !
  !> Function  base_device_wait:
  !! \memberof  psb_d_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!
  !
  subroutine d_base_device_wait()
    implicit none

  end subroutine d_base_device_wait

  function d_base_use_buffer() result(res)
    logical :: res

    res = .true.
  end function d_base_use_buffer

  subroutine d_base_new_buffer(n,x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    call psb_realloc(n,x%combuf,info)
  end subroutine d_base_new_buffer

  subroutine d_base_new_comid(n,x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    call psb_realloc(n,2_psb_ipk_,x%comid,info)
  end subroutine d_base_new_comid


  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_gthzv
  !! \memberof  psb_d_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  subroutine d_base_gthzv(n,idx,x,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: n, idx(:)
    real(psb_dpk_) ::  y(:)
    class(psb_d_base_vect_type) :: x

    if (x%is_dev()) call x%sync()
    call psi_gth(n,idx,x%v,y)

  end subroutine d_base_gthzv

  !
  ! Scatter:
  ! Y(IDX(:)) = beta*Y(IDX(:)) + X(:)
  !
  !
  !> Function  base_sctb
  !! \memberof  psb_d_base_vect_type
  !! \brief scatter into a class(base_vect)
  !!    Y(IDX(:)) = beta * Y(IDX(:)) +  X(:)
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param beta
  !! \param x(:)
  subroutine d_base_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: n, idx(:)
    real(psb_dpk_) :: beta, x(:)
    class(psb_d_base_vect_type) :: y

    if (y%is_dev()) call y%sync()
    call psi_sct(n,idx,x,beta,y%v)
    call y%set_host()

  end subroutine d_base_sctb

  subroutine d_base_sctb_x(i,n,idx,x,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    real(psb_dpk_) :: beta, x(:)
    class(psb_d_base_vect_type) :: y

    if (idx%is_dev()) call idx%sync()
    call y%sct(n,idx%v(i:),x,beta)
    call y%set_host()

  end subroutine d_base_sctb_x

  subroutine d_base_sctb_buf(i,n,idx,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    real(psb_dpk_) :: beta
    class(psb_d_base_vect_type) :: y


    if (.not.allocated(y%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'sctb_buf')
      return
    end if
    if (y%is_dev()) call y%sync()
    if (idx%is_dev()) call idx%sync()
    call y%sct(n,idx%v(i:),y%combuf(i:),beta)
    call y%set_host()

  end subroutine d_base_sctb_buf

end module psb_d_base_vect_mod





module psb_d_base_multivect_mod

  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_d_base_vect_mod

  !> \namespace  psb_base_mod  \class psb_d_base_vect_type
  !! The psb_d_base_vect_type
  !! defines a middle level  integer(psb_ipk_) encapsulated dense vector.
  !! The encapsulation is needed, in place of a simple array, to allow
  !! for complicated situations, such as GPU programming, where the memory
  !!  area we are interested in is not easily accessible from the host/Fortran
  !!  side. It is also meant to be encapsulated in an outer type, to allow
  !!  runtime switching as per the STATE design pattern, similar to the
  !!  sparse matrix types.
  !!
  private
  public  :: psb_d_base_multivect, psb_d_base_multivect_type

  type psb_d_base_multivect_type
    !> Values.
    real(psb_dpk_), allocatable :: v(:,:)
    real(psb_dpk_), allocatable :: combuf(:)
    integer(psb_mpk_), allocatable :: comid(:,:)
  contains
    !
    !  Constructors/allocators
    !
    procedure, pass(x) :: bld_x    => d_base_mlv_bld_x
    procedure, pass(x) :: bld_n    => d_base_mlv_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: all      => d_base_mlv_all
    procedure, pass(x) :: mold     => d_base_mlv_mold
    !
    ! Insert/set. Assembly and free.
    ! Assembly does almost nothing here, but is important
    ! in derived classes.
    !
    procedure, pass(x) :: ins      => d_base_mlv_ins
    procedure, pass(x) :: zero     => d_base_mlv_zero
    procedure, pass(x) :: asb      => d_base_mlv_asb
    procedure, pass(x) :: free     => d_base_mlv_free
    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder.
    !
    procedure, pass(x) :: sync     => d_base_mlv_sync
    procedure, pass(x) :: is_host  => d_base_mlv_is_host
    procedure, pass(x) :: is_dev   => d_base_mlv_is_dev
    procedure, pass(x) :: is_sync  => d_base_mlv_is_sync
    procedure, pass(x) :: set_host => d_base_mlv_set_host
    procedure, pass(x) :: set_dev  => d_base_mlv_set_dev
    procedure, pass(x) :: set_sync => d_base_mlv_set_sync

    !
    ! Basic info
    procedure, pass(x) :: get_nrows => d_base_mlv_get_nrows
    procedure, pass(x) :: get_ncols => d_base_mlv_get_ncols
    procedure, pass(x) :: sizeof    => d_base_mlv_sizeof
    procedure, nopass  :: get_fmt   => d_base_mlv_get_fmt
    !
    ! Set/get data from/to an external array; also
    ! overload assignment.
    !
    procedure, pass(x) :: get_vect => d_base_mlv_get_vect
    procedure, pass(x) :: set_scal => d_base_mlv_set_scal
    procedure, pass(x) :: set_vect => d_base_mlv_set_vect
    generic, public    :: set      => set_vect, set_scal

    !
    ! Dot product and AXPBY
    !
    procedure, pass(x) :: dot_v    => d_base_mlv_dot_v
    procedure, pass(x) :: dot_a    => d_base_mlv_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => d_base_mlv_axpby_v
    procedure, pass(y) :: axpby_a  => d_base_mlv_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    !
    ! MultiVector by vector/multivector multiplication. Need all variants
    ! to handle multiple requirements from preconditioners
    !
    procedure, pass(y) :: mlt_mv   => d_base_mlv_mlt_mv
    procedure, pass(y) :: mlt_mv_v => d_base_mlv_mlt_mv_v
    procedure, pass(y) :: mlt_ar1  => d_base_mlv_mlt_ar1
    procedure, pass(y) :: mlt_ar2  => d_base_mlv_mlt_ar2
    procedure, pass(z) :: mlt_a_2  => d_base_mlv_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => d_base_mlv_mlt_v_2
!!$    procedure, pass(z) :: mlt_va   => d_base_mlv_mlt_va
!!$    procedure, pass(z) :: mlt_av   => d_base_mlv_mlt_av
    generic, public    :: mlt      => mlt_mv, mlt_mv_v, mlt_ar1, mlt_ar2, &
         & mlt_a_2, mlt_v_2 !, mlt_av, mlt_va
    !
    ! Scaling and norms
    !
    procedure, pass(x) :: scal     => d_base_mlv_scal
    procedure, pass(x) :: nrm2     => d_base_mlv_nrm2
    procedure, pass(x) :: amax     => d_base_mlv_amax
    procedure, pass(x) :: asum     => d_base_mlv_asum
    procedure, pass(x) :: absval1  => d_base_mlv_absval1
    procedure, pass(x) :: absval2  => d_base_mlv_absval2
    generic, public    :: absval   => absval1, absval2

    !
    ! These are for handling gather/scatter in new
    ! comm internals implementation.
    !
    procedure, nopass  :: use_buffer   => d_base_mlv_use_buffer
    procedure, pass(x) :: new_buffer   => d_base_mlv_new_buffer
    procedure, nopass  :: device_wait  => d_base_mlv_device_wait
    procedure, pass(x) :: maybe_free_buffer  => d_base_mlv_maybe_free_buffer
    procedure, pass(x) :: free_buffer  => d_base_mlv_free_buffer
    procedure, pass(x) :: new_comid    => d_base_mlv_new_comid
    procedure, pass(x) :: free_comid   => d_base_mlv_free_comid

    !
    ! Gather/scatter. These are needed for MPI interfacing.
    ! May have to be reworked.
    !
    procedure, pass(x) :: gthab    => d_base_mlv_gthab
    procedure, pass(x) :: gthzv    => d_base_mlv_gthzv
    procedure, pass(x) :: gthzm    => d_base_mlv_gthzm
    procedure, pass(x) :: gthzv_x  => d_base_mlv_gthzv_x
    procedure, pass(x) :: gthzbuf  => d_base_mlv_gthzbuf
    generic, public    :: gth      => gthab, gthzv, gthzm, gthzv_x, gthzbuf
    procedure, pass(y) :: sctb     => d_base_mlv_sctb
    procedure, pass(y) :: sctbr2   => d_base_mlv_sctbr2
    procedure, pass(y) :: sctb_x   => d_base_mlv_sctb_x
    procedure, pass(y) :: sctb_buf => d_base_mlv_sctb_buf
    generic, public    :: sct      => sctb, sctbr2, sctb_x, sctb_buf
  end type psb_d_base_multivect_type

  interface psb_d_base_multivect
    module procedure constructor, size_const
  end interface psb_d_base_multivect

contains

  !
  ! Constructors.
  !

  !> Function  constructor:
  !! \brief     Constructor from an array
  !!  \param   x(:)  input array to be copied
  !!
  function constructor(x) result(this)
    real(psb_dpk_)   :: x(:,:)
    type(psb_d_base_multivect_type) :: this
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
    type(psb_d_base_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%asb(m,n,info)

  end function size_const

  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_d_base_multivect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  subroutine d_base_mlv_bld_x(x,this)
    use psb_realloc_mod
    real(psb_dpk_), intent(in) :: this(:,:)
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(size(this,1),size(this,2),x%v,info)
    if (info /= 0) then
      call psb_errpush(psb_err_alloc_dealloc_,'base_mlv_vect_bld')
      return
    end if
    x%v(:,:)  = this(:,:)

  end subroutine d_base_mlv_bld_x

  !
  ! Create with size, but no initialization
  !

  !> Function  bld_n:
  !! \memberof  psb_d_base_multivect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  subroutine d_base_mlv_bld_n(x,m,n)
    use psb_realloc_mod
    integer(psb_ipk_), intent(in) :: m,n
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_) :: info

    call psb_realloc(m,n,x%v,info)
    call x%asb(m,n,info)

  end subroutine d_base_mlv_bld_n

  !> Function  base_mlv_all:
  !! \memberof  psb_d_base_multivect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated.
  !!  \param info  return code
  !!
  subroutine d_base_mlv_all(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m,n
    class(psb_d_base_multivect_type), intent(out) :: x
    integer(psb_ipk_), intent(out)              :: info

    call psb_realloc(m,n,x%v,info)

  end subroutine d_base_mlv_all

  !> Function  base_mlv_mold:
  !! \memberof  psb_d_base_multivect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  subroutine d_base_mlv_mold(x, y, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_d_base_multivect_type), intent(in)   :: x
    class(psb_d_base_multivect_type), intent(out), allocatable :: y
    integer(psb_ipk_), intent(out)              :: info

    allocate(psb_d_base_multivect_type :: y, stat=info)

  end subroutine d_base_mlv_mold

  !
  ! Insert a bunch of values at specified positions.
  !
  !> Function  base_mlv_ins:
  !! \memberof  psb_d_base_multivect_type
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
  subroutine d_base_mlv_ins(n,irl,val,dupl,x,info)
    use psi_serial_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(in)               :: n, dupl
    integer(psb_ipk_), intent(in)               :: irl(:)
    real(psb_dpk_), intent(in)        :: val(:,:)
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

  end subroutine d_base_mlv_ins

  !
  !> Function  base_mlv_zero
  !! \memberof  psb_d_base_multivect_type
  !! \brief Zero out contents
  !!
  !
  subroutine d_base_mlv_zero(x)
    use psi_serial_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout)    :: x

    if (allocated(x%v)) x%v=dzero

  end subroutine d_base_mlv_zero


  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_mlv_asb:
  !! \memberof  psb_d_base_multivect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !

  subroutine d_base_mlv_asb(m,n, x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in)              :: m,n
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if ((x%get_nrows() < m).or.(x%get_ncols()<n)) &
         & call psb_realloc(m,n,x%v,info)
    if (info /= 0) &
         & call psb_errpush(psb_err_alloc_dealloc_,'vect_asb')

  end subroutine d_base_mlv_asb


  !
  !> Function  base_mlv_free:
  !! \memberof  psb_d_base_multivect_type
  !! \brief Free vector
  !!
  !!  \param info  return code
  !!
  !
  subroutine d_base_mlv_free(x, info)
    use psi_serial_mod
    use psb_realloc_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout)  :: x
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (allocated(x%v)) deallocate(x%v, stat=info)
    if (info /= 0) call &
         & psb_errpush(psb_err_alloc_dealloc_,'vect_free')

  end subroutine d_base_mlv_free



  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  !
  !
  !> Function  base_mlv_sync:
  !! \memberof  psb_d_base_multivect_type
  !! \brief Sync: base version is a no-op.
  !!
  !
  subroutine d_base_mlv_sync(x)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x

  end subroutine d_base_mlv_sync

  !
  !> Function  base_mlv_set_host:
  !! \memberof  psb_d_base_multivect_type
  !! \brief Set_host: base version is a no-op.
  !!
  !
  subroutine d_base_mlv_set_host(x)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x

  end subroutine d_base_mlv_set_host

  !
  !> Function  base_mlv_set_dev:
  !! \memberof  psb_d_base_multivect_type
  !! \brief Set_dev: base version is a no-op.
  !!
  !
  subroutine d_base_mlv_set_dev(x)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x

  end subroutine d_base_mlv_set_dev

  !
  !> Function  base_mlv_set_sync:
  !! \memberof  psb_d_base_multivect_type
  !! \brief Set_sync: base version is a no-op.
  !!
  !
  subroutine d_base_mlv_set_sync(x)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x

  end subroutine d_base_mlv_set_sync

  !
  !> Function  base_mlv_is_dev:
  !! \memberof  psb_d_base_multivect_type
  !! \brief Is  vector on external device    .
  !!
  !
  function d_base_mlv_is_dev(x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(in) :: x
    logical  :: res

    res = .false.
  end function d_base_mlv_is_dev

  !
  !> Function  base_mlv_is_host
  !! \memberof  psb_d_base_multivect_type
  !! \brief Is  vector on standard memory    .
  !!
  !
  function d_base_mlv_is_host(x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function d_base_mlv_is_host

  !
  !> Function  base_mlv_is_sync
  !! \memberof  psb_d_base_multivect_type
  !! \brief Is  vector on sync               .
  !!
  !
  function d_base_mlv_is_sync(x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(in) :: x
    logical  :: res

    res = .true.
  end function d_base_mlv_is_sync


  !
  ! Size info.
  !
  !
  !> Function  base_mlv_get_nrows
  !! \memberof  psb_d_base_multivect_type
  !! \brief  Number of entries
  !!
  !
  function d_base_mlv_get_nrows(x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v,1)

  end function d_base_mlv_get_nrows

  function d_base_mlv_get_ncols(x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(in) :: x
    integer(psb_ipk_) :: res

    res = 0
    if (allocated(x%v)) res = size(x%v,2)

  end function d_base_mlv_get_ncols

  !
  !> Function  base_mlv_get_sizeof
  !! \memberof  psb_d_base_multivect_type
  !! \brief  Size in bytesa
  !!
  !
  function d_base_mlv_sizeof(x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(in) :: x
    integer(psb_epk_) :: res

    ! Force 8-byte integers.
    res = (1_psb_epk_ * psb_sizeof_ip) * x%get_nrows() * x%get_ncols()

  end function d_base_mlv_sizeof

  !
  !> Function  base_mlv_get_fmt
  !! \memberof  psb_d_base_multivect_type
  !! \brief  Format
  !!
  !
  function d_base_mlv_get_fmt() result(res)
    implicit none
    character(len=5) :: res
    res = 'BASE'
  end function d_base_mlv_get_fmt


  !
  !
  !
  !> Function  base_mlv_get_vect
  !! \memberof  psb_d_base_multivect_type
  !! \brief  Extract a copy of the contents
  !!
  !
  function  d_base_mlv_get_vect(x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    real(psb_dpk_), allocatable                 :: res(:,:)
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
  end function d_base_mlv_get_vect

  !
  ! Reset all values
  !
  !
  !> Function  base_mlv_set_scal
  !! \memberof  psb_d_base_multivect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  subroutine d_base_mlv_set_scal(x,val)
    implicit none
    class(psb_d_base_multivect_type), intent(inout)  :: x
    real(psb_dpk_), intent(in) :: val

    integer(psb_ipk_) :: info
    x%v = val

  end subroutine d_base_mlv_set_scal

  !
  !> Function  base_mlv_set_vect
  !! \memberof  psb_d_base_multivect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in
  !!
  subroutine d_base_mlv_set_vect(x,val)
    implicit none
    class(psb_d_base_multivect_type), intent(inout)  :: x
    real(psb_dpk_), intent(in) :: val(:,:)
    integer(psb_ipk_) :: nr, nc
    integer(psb_ipk_) :: info

    if (allocated(x%v)) then
      nr = min(size(x%v,1),size(val,1))
      nc = min(size(x%v,2),size(val,2))

      x%v(1:nr,1:nc) = val(1:nr,1:nc)
    else
      x%v = val
    end if

  end subroutine d_base_mlv_set_vect

  !
  ! Dot products
  !
  !
  !> Function  base_mlv_dot_v
  !! \memberof  psb_d_base_multivect_type
  !! \brief  Dot product by another base_mlv_vector
  !! \param n    Number of entries to be considered
  !! \param y    The other (base_mlv_vect) to be multiplied by
  !!
  function d_base_mlv_dot_v(n,x,y) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_), allocatable   :: res(:)
    real(psb_dpk_), external      :: ddot
    integer(psb_ipk_) :: j,nc

    if (x%is_dev()) call x%sync()
    res = dzero
    !
    ! Note: this is the base implementation.
    !  When we get here, we are sure that X is of
    !  TYPE psb_d_base_mlv_vect (or its class does not care).
    !  If Y is not, throw the burden on it, implicitly
    !  calling dot_a
    !
    select type(yy => y)
    type is (psb_d_base_multivect_type)
      if (y%is_dev()) call y%sync()
      nc = min(psb_size(x%v,2_psb_ipk_),psb_size(y%v,2_psb_ipk_))
      allocate(res(nc))
      do j=1,nc
        res(j) = ddot(n,x%v(:,j),1,y%v(:,j),1)
      end do
      class default
      res = y%dot(n,x%v)
    end select

  end function d_base_mlv_dot_v

  !
  ! Base workhorse is good old BLAS1
  !
  !
  !> Function  base_mlv_dot_a
  !! \memberof  psb_d_base_multivect_type
  !! \brief  Dot product by a normal array
  !! \param n    Number of entries to be considered
  !! \param y(:) The array to be multiplied by
  !!
  function d_base_mlv_dot_a(n,x,y) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    real(psb_dpk_), intent(in)    :: y(:,:)
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_), allocatable     :: res(:)
    real(psb_dpk_), external      :: ddot
    integer(psb_ipk_) :: j,nc

    if (x%is_dev()) call x%sync()
    nc = min(psb_size(x%v,2_psb_ipk_),size(y,2_psb_ipk_))
    allocate(res(nc))
    do j=1,nc
      res(j) = ddot(n,x%v(:,j),1,y(:,j),1)
    end do

  end function d_base_mlv_dot_a

  !
  ! AXPBY is invoked via Y, hence the structure below.
  !
  !
  !
  !> Function  base_mlv_axpby_v
  !! \memberof  psb_d_base_multivect_type
  !! \brief AXPBY  by a (base_mlv_vect) y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x     The class(base_mlv_vect) to be added
  !! \param beta scalar alpha
  !! \param info   return code
  !!
  subroutine d_base_mlv_axpby_v(m,alpha, x, beta, y, info, n)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    class(psb_d_base_multivect_type), intent(inout)  :: x
    class(psb_d_base_multivect_type), intent(inout)  :: y
    real(psb_dpk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_), intent(in), optional     :: n
    integer(psb_ipk_)  :: nc

    if (present(n)) then
      nc = n
    else
      nc = min(psb_size(x%v,2_psb_ipk_),psb_size(y%v,2_psb_ipk_))
    end if
    select type(xx => x)
    type is (psb_d_base_multivect_type)
      call psb_geaxpby(m,nc,alpha,x%v,beta,y%v,info)
      class default
      call y%axpby(m,alpha,x%v,beta,info,n=n)
    end select

  end subroutine d_base_mlv_axpby_v

  !
  ! AXPBY is invoked via Y, hence the structure below.
  !
  !
  !> Function  base_mlv_axpby_a
  !! \memberof  psb_d_base_multivect_type
  !! \brief AXPBY  by a normal array y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x(:) The array to be added
  !! \param beta scalar alpha
  !! \param info   return code
  !!
  subroutine d_base_mlv_axpby_a(m,alpha, x, beta, y, info,n)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_), intent(in)               :: m
    real(psb_dpk_), intent(in)        :: x(:,:)
    class(psb_d_base_multivect_type), intent(inout)  :: y
    real(psb_dpk_), intent (in)       :: alpha, beta
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_), intent(in), optional     :: n
    integer(psb_ipk_)  :: nc
    if (present(n)) then
      nc = n
    else
      nc = min(size(x,2),psb_size(y%v,2_psb_ipk_))
    end if

    call psb_geaxpby(m,nc,alpha,x,beta,y%v,info)

  end subroutine d_base_mlv_axpby_a


  !
  !  Multiple variants of two operations:
  !  Simple multiplication  Y(:.:) = X(:,:)*Y(:,:)
  !  blas-like:   Z(:) = alpha*X(:)*Y(:)+beta*Z(:)
  !
  !  Variants expanded according to the dynamic type
  !  of the involved entities
  !
  !
  !> Function  base_mlv_mlt_mv
  !! \memberof  psb_d_base_multivect_type
  !! \brief Multivector entry-by-entry multiply  by a base_mlv_multivect  y=x*y
  !! \param x   The class(base_mlv_vect) to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlv_mlt_mv(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout)  :: x
    class(psb_d_base_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (x%is_dev()) call x%sync()
    call y%mlt(x%v,info)

  end subroutine d_base_mlv_mlt_mv

  subroutine d_base_mlv_mlt_mv_v(x, y, info)
    use psi_serial_mod
    implicit none
    class(psb_d_base_vect_type), intent(inout)  :: x
    class(psb_d_base_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info

    info = 0
    if (x%is_dev()) call x%sync()
    call y%mlt(x%v,info)

  end subroutine d_base_mlv_mlt_mv_v

  !
  !> Function  base_mlv_mlt_ar1
  !! \memberof  psb_d_base_multivect_type
  !! \brief MultiVector entry-by-entry multiply  by a rank 1 array y=x*y
  !! \param x(:) The array to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlv_mlt_ar1(x, y, info)
    use psi_serial_mod
    implicit none
    real(psb_dpk_), intent(in)        :: x(:)
    class(psb_d_base_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, n

    info = 0
    n = min(psb_size(y%v,1_psb_ipk_), size(x))
    do i=1, n
      y%v(i,:) = y%v(i,:)*x(i)
    end do

  end subroutine d_base_mlv_mlt_ar1

  !
  !> Function  base_mlv_mlt_ar2
  !! \memberof  psb_d_base_multivect_type
  !! \brief MultiVector entry-by-entry multiply  by a rank 2 array y=x*y
  !! \param x(:,:) The array to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlv_mlt_ar2(x, y, info)
    use psi_serial_mod
    implicit none
    real(psb_dpk_), intent(in)        :: x(:,:)
    class(psb_d_base_multivect_type), intent(inout)  :: y
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, nr,nc

    info = 0
    nr   = min(psb_size(y%v,1_psb_ipk_), size(x,1))
    nc   = min(psb_size(y%v,2_psb_ipk_), size(x,2))
    y%v(1:nr,1:nc) = y%v(1:nr,1:nc)*x(1:nr,1:nc)

  end subroutine d_base_mlv_mlt_ar2


  !
  !> Function  base_mlv_mlt_a_2
  !! \memberof  psb_d_base_multivect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by normal arrays
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x(:) The array to be multiplied b
  !! \param y(:) The array to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlv_mlt_a_2(alpha,x,y,beta,z,info)
    use psi_serial_mod
    implicit none
    real(psb_dpk_), intent(in)        :: alpha,beta
    real(psb_dpk_), intent(in)        :: y(:,:)
    real(psb_dpk_), intent(in)        :: x(:,:)
    class(psb_d_base_multivect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    integer(psb_ipk_) :: i, nr, nc

    info = 0
    nr = min(psb_size(z%v,1_psb_ipk_), size(x,1), size(y,1))
    nc = min(psb_size(z%v,2_psb_ipk_), size(x,2), size(y,2))
    if (alpha == dzero) then
      if (beta == done) then
        return
      else
        z%v(1:nr,1:nc) = beta*z%v(1:nr,1:nc)
      end if
    else
      if (alpha == done) then
        if (beta == dzero) then
          z%v(1:nr,1:nc) = y(1:nr,1:nc)*x(1:nr,1:nc)
        else if (beta == done) then
          z%v(1:nr,1:nc) = z%v(1:nr,1:nc) + y(1:nr,1:nc)*x(1:nr,1:nc)
        else
          z%v(1:nr,1:nc) = beta*z%v(1:nr,1:nc) + y(1:nr,1:nc)*x(1:nr,1:nc)
        end if
      else if (alpha == -done) then
        if (beta == dzero) then
          z%v(1:nr,1:nc) = -y(1:nr,1:nc)*x(1:nr,1:nc)
        else if (beta == done) then
          z%v(1:nr,1:nc) = z%v(1:nr,1:nc) - y(1:nr,1:nc)*x(1:nr,1:nc)
        else
          z%v(1:nr,1:nc) = beta*z%v(1:nr,1:nc) - y(1:nr,1:nc)*x(1:nr,1:nc)
        end if
      else
        if (beta == dzero) then
          z%v(1:nr,1:nc) = alpha*y(1:nr,1:nc)*x(1:nr,1:nc)
        else if (beta == done) then
          z%v(1:nr,1:nc) = z%v(1:nr,1:nc) + alpha*y(1:nr,1:nc)*x(1:nr,1:nc)
        else
          z%v(1:nr,1:nc) = beta*z%v(1:nr,1:nc) + alpha*y(1:nr,1:nc)*x(1:nr,1:nc)
        end if
      end if
    end if
  end subroutine d_base_mlv_mlt_a_2

  !
  !> Function  base_mlv_mlt_v_2
  !! \memberof  psb_d_base_multivect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by class(base_mlv_vect)
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x  The class(base_mlv_vect) to be multiplied b
  !! \param y  The class(base_mlv_vect) to be multiplied by
  !! \param info   return code
  !!
  subroutine d_base_mlv_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
    use psi_serial_mod
    use psb_string_mod
    implicit none
    real(psb_dpk_), intent(in)        :: alpha,beta
    class(psb_d_base_multivect_type), intent(inout)  :: x
    class(psb_d_base_multivect_type), intent(inout)  :: y
    class(psb_d_base_multivect_type), intent(inout)  :: z
    integer(psb_ipk_), intent(out)              :: info
    character(len=1), intent(in), optional     :: conjgx, conjgy
    integer(psb_ipk_) :: i, n
    logical :: conjgx_, conjgy_

    info = 0
    if (x%is_dev()) call x%sync()
    if (y%is_dev()) call y%sync()
    if (z%is_dev()) call z%sync()
    if (.not.psb_d_is_complex_) then
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
  end subroutine d_base_mlv_mlt_v_2
!!$
!!$  subroutine d_base_mlv_mlt_av(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_dpk_), intent(in)        :: alpha,beta
!!$    real(psb_dpk_), intent(in)        :: x(:)
!!$    class(psb_d_base_multivect_type), intent(inout)  :: y
!!$    class(psb_d_base_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    call z%mlt(alpha,x,y%v,beta,info)
!!$
!!$  end subroutine d_base_mlv_mlt_av
!!$
!!$  subroutine d_base_mlv_mlt_va(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_dpk_), intent(in)        :: alpha,beta
!!$    real(psb_dpk_), intent(in)        :: y(:)
!!$    class(psb_d_base_multivect_type), intent(inout)  :: x
!!$    class(psb_d_base_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    call z%mlt(alpha,y,x,beta,info)
!!$
!!$  end subroutine d_base_mlv_mlt_va
!!$
!!$
  !
  ! Simple scaling
  !
  !> Function  base_mlv_scal
  !! \memberof  psb_d_base_multivect_type
  !! \brief Scale all entries  x = alpha*x
  !! \param alpha   The multiplier
  !!
  subroutine d_base_mlv_scal(alpha, x)
    use psi_serial_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout)  :: x
    real(psb_dpk_), intent (in)       :: alpha

    if (x%is_dev()) call x%sync()
    if (allocated(x%v)) x%v = alpha*x%v

  end subroutine d_base_mlv_scal

  !
  ! Norms 1, 2 and infinity
  !
  !> Function  base_mlv_nrm2
  !! \memberof  psb_d_base_multivect_type
  !! \brief 2-norm |x(1:n)|_2
  !! \param n  how many entries to consider
  function d_base_mlv_nrm2(n,x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_), allocatable    :: res(:)
    real(psb_dpk_), external      :: dnrm2
    integer(psb_ipk_) :: j, nc

    if (x%is_dev()) call x%sync()
    nc = psb_size(x%v,2_psb_ipk_)
    allocate(res(nc))
    do j=1,nc
      res(j) =  dnrm2(n,x%v(:,j),1)
    end do

  end function d_base_mlv_nrm2

  !
  !> Function  base_mlv_amax
  !! \memberof  psb_d_base_multivect_type
  !! \brief infinity-norm |x(1:n)|_\infty
  !! \param n  how many entries to consider
  function d_base_mlv_amax(n,x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_), allocatable    :: res(:)
    integer(psb_ipk_) :: j, nc

    if (x%is_dev()) call x%sync()
    nc = psb_size(x%v,2_psb_ipk_)
    allocate(res(nc))
    do j=1,nc
      res(j) =  maxval(abs(x%v(1:n,j)))
    end do

  end function d_base_mlv_amax

  !
  !> Function  base_mlv_asum
  !! \memberof  psb_d_base_multivect_type
  !! \brief 1-norm |x(1:n)|_1
  !! \param n  how many entries to consider
  function d_base_mlv_asum(n,x) result(res)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)           :: n
    real(psb_dpk_), allocatable    :: res(:)
    integer(psb_ipk_) :: j, nc

    if (x%is_dev()) call x%sync()
    nc = psb_size(x%v,2_psb_ipk_)
    allocate(res(nc))
    do j=1,nc
      res(j) =  sum(abs(x%v(1:n,j)))
    end do

  end function d_base_mlv_asum
  !
  ! Overwrite with absolute value
  !
  !
  !> Function  base_absval1
  !! \memberof  psb_d_base_vect_type
  !! \brief  Set all entries to their respective absolute values.
  !!
  subroutine d_base_mlv_absval1(x)
    implicit none
    class(psb_d_base_multivect_type), intent(inout)  :: x

    if (allocated(x%v)) then
      if (x%is_dev()) call x%sync()
      x%v =  abs(x%v)
      call x%set_host()
    end if

  end subroutine d_base_mlv_absval1

  subroutine d_base_mlv_absval2(x,y)
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    class(psb_d_base_multivect_type), intent(inout) :: y
    integer(psb_ipk_) :: info

    if (x%is_dev()) call x%sync()
    if (allocated(x%v)) then
      call y%axpby(min(x%get_nrows(),y%get_nrows()),done,x,dzero,info)
      call y%absval()
    end if

  end subroutine d_base_mlv_absval2


  function d_base_mlv_use_buffer() result(res)
    implicit none
    logical :: res

    res = .true.
  end function d_base_mlv_use_buffer

  subroutine d_base_mlv_new_buffer(n,x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    integer(psb_ipk_)               :: nc
    nc = x%get_ncols()
    call psb_realloc(n*nc,x%combuf,info)
  end subroutine d_base_mlv_new_buffer

  subroutine d_base_mlv_new_comid(n,x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(in)              :: n
    integer(psb_ipk_), intent(out)             :: info

    call psb_realloc(n,2_psb_ipk_,x%comid,info)
  end subroutine d_base_mlv_new_comid


  subroutine d_base_mlv_maybe_free_buffer(x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info


    info = 0
    if (psb_get_maybe_free_buffer())&
         &  call x%free_buffer(info)

  end subroutine d_base_mlv_maybe_free_buffer

  subroutine d_base_mlv_free_buffer(x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%combuf)) &
         &  deallocate(x%combuf,stat=info)
  end subroutine d_base_mlv_free_buffer

  subroutine d_base_mlv_free_comid(x,info)
    use psb_realloc_mod
    implicit none
    class(psb_d_base_multivect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out)             :: info

    if (allocated(x%comid)) &
         &  deallocate(x%comid,stat=info)
  end subroutine d_base_mlv_free_comid


  !
  ! Gather: Y = beta * Y + alpha * X(IDX(:))
  !
  !
  !> Function  base_mlv_gthab
  !! \memberof  psb_d_base_multivect_type
  !! \brief gather into an array
  !!    Y = beta * Y + alpha * X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param alpha
  !! \param beta
  subroutine d_base_mlv_gthab(n,idx,alpha,x,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: n, idx(:)
    real(psb_dpk_) :: alpha, beta, y(:)
    class(psb_d_base_multivect_type) :: x
    integer(psb_ipk_) :: nc

    if (x%is_dev()) call x%sync()
    if (.not.allocated(x%v)) then
      return
    end if
    nc = psb_size(x%v,2_psb_ipk_)
    call psi_gth(n,nc,idx,alpha,x%v,beta,y)

  end subroutine d_base_mlv_gthab
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_d_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  subroutine d_base_mlv_gthzv_x(i,n,idx,x,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: i,n
    class(psb_i_base_vect_type) :: idx
    real(psb_dpk_) ::  y(:)
    class(psb_d_base_multivect_type) :: x

    if (x%is_dev()) call x%sync()
    call x%gth(n,idx%v(i:),y)

  end subroutine d_base_mlv_gthzv_x

  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_d_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  subroutine d_base_mlv_gthzv(n,idx,x,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: n, idx(:)
    real(psb_dpk_) ::  y(:)
    class(psb_d_base_multivect_type) :: x
    integer(psb_ipk_) :: nc

    if (x%is_dev()) call x%sync()
    if (.not.allocated(x%v)) then
      return
    end if
    nc = psb_size(x%v,2_psb_ipk_)

    call psi_gth(n,nc,idx,x%v,y)

  end subroutine d_base_mlv_gthzv
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_d_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  subroutine d_base_mlv_gthzm(n,idx,x,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: n, idx(:)
    real(psb_dpk_) ::  y(:,:)
    class(psb_d_base_multivect_type) :: x
    integer(psb_ipk_) :: nc

    if (x%is_dev()) call x%sync()
    if (.not.allocated(x%v)) then
      return
    end if
    nc = psb_size(x%v,2_psb_ipk_)

    call psi_gth(n,nc,idx,x%v,y)

  end subroutine d_base_mlv_gthzm

  !
  ! New comm internals impl.
  !
  subroutine d_base_mlv_gthzbuf(i,ixb,n,idx,x)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: i, ixb, n
    class(psb_i_base_vect_type) :: idx
    class(psb_d_base_multivect_type) :: x
    integer(psb_ipk_) :: nc

    if (.not.allocated(x%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'gthzbuf')
      return
    end if
    if (idx%is_dev()) call idx%sync()
    if (x%is_dev()) call x%sync()
    nc = x%get_ncols()
    call x%gth(n,idx%v(i:),x%combuf(ixb:))

  end subroutine d_base_mlv_gthzbuf

  !
  ! Scatter:
  ! Y(IDX(:),:) = beta*Y(IDX(:),:) + X(:)
  !
  !
  !> Function  base_mlv_sctb
  !! \memberof  psb_d_base_multivect_type
  !! \brief scatter into a class(base_mlv_vect)
  !!    Y(IDX(:)) = beta * Y(IDX(:)) +  X(:)
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param beta
  !! \param x(:)
  subroutine d_base_mlv_sctb(n,idx,x,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: n, idx(:)
    real(psb_dpk_) :: beta, x(:)
    class(psb_d_base_multivect_type) :: y
    integer(psb_ipk_) :: nc

    if (y%is_dev()) call y%sync()
    nc = psb_size(y%v,2_psb_ipk_)
    call psi_sct(n,nc,idx,x,beta,y%v)
    call y%set_host()

  end subroutine d_base_mlv_sctb

  subroutine d_base_mlv_sctbr2(n,idx,x,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: n, idx(:)
    real(psb_dpk_) :: beta, x(:,:)
    class(psb_d_base_multivect_type) :: y
    integer(psb_ipk_) :: nc

    if (y%is_dev()) call y%sync()
    nc = y%get_ncols()
    call psi_sct(n,nc,idx,x,beta,y%v)
    call y%set_host()

  end subroutine d_base_mlv_sctbr2

  subroutine d_base_mlv_sctb_x(i,n,idx,x,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: i, n
    class(psb_i_base_vect_type) :: idx
    real( psb_dpk_) :: beta, x(:)
    class(psb_d_base_multivect_type) :: y

    call y%sct(n,idx%v(i:),x,beta)

  end subroutine d_base_mlv_sctb_x

  subroutine d_base_mlv_sctb_buf(i,iyb,n,idx,beta,y)
    use psi_serial_mod
    implicit none
    integer(psb_ipk_) :: i, iyb, n
    class(psb_i_base_vect_type) :: idx
    real(psb_dpk_) :: beta
    class(psb_d_base_multivect_type) :: y
    integer(psb_ipk_) :: nc

    if (.not.allocated(y%combuf)) then
      call psb_errpush(psb_err_alloc_dealloc_,'sctb_buf')
      return
    end if
    if (y%is_dev()) call y%sync()
    if (idx%is_dev()) call idx%sync()
    nc = y%get_ncols()
    call y%sct(n,idx%v(i:),y%combuf(iyb:),beta)
    call y%set_host()

  end subroutine d_base_mlv_sctb_buf

  !
  !> Function  base_device_wait:
  !! \memberof  psb_d_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!
  !
  subroutine d_base_mlv_device_wait()
    implicit none

  end subroutine d_base_mlv_device_wait

end module psb_d_base_multivect_mod
