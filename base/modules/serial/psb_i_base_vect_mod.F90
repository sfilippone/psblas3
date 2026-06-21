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
!
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
    integer(psb_mpk_), allocatable :: comid(:,:)
    !> vector bldstate:
    !!    null:   pristine;
    !!    build:  it's being filled with entries;
    !!    assembled: ready to use in computations;
    !!    update: accepts coefficients but only
    !!            in already existing entries.
    !!    The transitions among the states are detailed in
    !!            psb_T_vect_mod.
    integer(psb_ipk_), private :: bldstate = psb_vect_null_
    integer(psb_ipk_), private :: dupl     = psb_dupl_null_
    integer(psb_ipk_), private :: ncfs     = 0
    integer(psb_ipk_), allocatable :: iv(:)
  contains
    !
    !  Constructors/allocators
    !
    procedure, pass(x) :: bld_x    => i_base_bld_x
    procedure, pass(x) :: bld_mn   => i_base_bld_mn
    procedure, pass(x) :: bld_en   => i_base_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
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
    procedure, pass(x) :: asb_m    => i_base_asb_m
    procedure, pass(x) :: asb_e    => i_base_asb_e
    generic, public    :: asb      => asb_m, asb_e
    procedure, pass(x) :: free     => i_base_free
    procedure, pass(x) :: reinit      => i_base_reinit
    procedure, pass(x) :: set_ncfs => i_base_set_ncfs 
    procedure, pass(x) :: get_ncfs => i_base_get_ncfs
    procedure, pass(x) :: set_dupl => i_base_set_dupl 
    procedure, pass(x) :: get_dupl => i_base_get_dupl
    procedure, pass(x) :: set_state => i_base_set_state
    procedure, pass(x) :: set_null  => i_base_set_null
    procedure, pass(x) :: set_bld   => i_base_set_bld
    procedure, pass(x) :: set_upd   => i_base_set_upd
    procedure, pass(x) :: set_asb   => i_base_set_asb
    procedure, pass(x) :: get_state => i_base_get_state
    procedure, pass(x) :: is_null   => i_base_is_null
    procedure, pass(x) :: is_bld    => i_base_is_bld
    procedure, pass(x) :: is_upd    => i_base_is_upd
    procedure, pass(x) :: is_asb    => i_base_is_asb
    procedure, pass(x) :: base_cpy  => i_base_cpy
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
    procedure, pass(x) :: maybe_free_buffer  => i_base_maybe_free_buffer
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

    procedure, pass(x) :: check_addr  => i_base_check_addr






  end type psb_i_base_vect_type

  public  :: psb_i_base_vect, psb_i_base_vect_type
  private :: constructor, size_const
  interface psb_i_base_vect
    module procedure constructor, size_const
  end interface psb_i_base_vect

  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  interface 
    module subroutine i_base_bld_x(x,this,scratch)
      integer(psb_ipk_), intent(in) :: this(:)
      class(psb_i_base_vect_type), intent(inout) :: x
      logical, intent(in), optional        :: scratch
    end subroutine i_base_bld_x
  end interface

  !
  ! Create with size, but no initialization
  !

  !> Function  bld_mn:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  interface 
    module subroutine i_base_bld_mn(x,n,scratch)
      integer(psb_mpk_), intent(in) :: n
      class(psb_i_base_vect_type), intent(inout) :: x
      logical, intent(in), optional        :: scratch
    end subroutine i_base_bld_mn
  end interface
  
  !> Function  bld_en:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  interface 
    module subroutine i_base_bld_en(x,n,scratch)
      integer(psb_epk_), intent(in) :: n
      class(psb_i_base_vect_type), intent(inout) :: x
      logical, intent(in), optional        :: scratch
    end subroutine i_base_bld_en
  end interface

  !> Function  base_all:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated.
  !!  \param info  return code
  !!
  interface 
    module subroutine i_base_all(n, x, info)
      integer(psb_ipk_), intent(in)               :: n
      class(psb_i_base_vect_type), intent(out)    :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_all
  end interface

  !> Function  base_mold:
  !! \memberof  psb_i_base_vect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  interface 
    module subroutine i_base_mold(x, y, info)
      class(psb_i_base_vect_type), intent(in)   :: x
      class(psb_i_base_vect_type), intent(out), allocatable :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_mold
  end interface

  interface 
    module subroutine i_base_reinit(x, info,clear)
      class(psb_i_base_vect_type), intent(inout)    :: x
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in), optional               :: clear
    end subroutine i_base_reinit
  end interface

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
  interface 
    module subroutine i_base_ins_a(n,irl,val,dupl,x,maxr,info)
      class(psb_i_base_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, dupl, maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      integer(psb_ipk_), intent(in)        :: val(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_ins_a
  end interface


  interface 
    module subroutine i_base_ins_v(n,irl,val,dupl,x,maxr,info)
      class(psb_i_base_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, dupl, maxr
      class(psb_i_base_vect_type), intent(inout)  :: irl
      class(psb_i_base_vect_type), intent(inout)  :: val
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_ins_v
  end interface



  !
  !> Function  base_zero
  !! \memberof  psb_i_base_vect_type
  !! \brief Zero out contents
  !!
  !
  interface 
    module subroutine i_base_zero(x)
      class(psb_i_base_vect_type), intent(inout)    :: x
    end subroutine i_base_zero
  end interface



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

  interface 
    module subroutine i_base_asb_m(n, x, info, scratch)
      integer(psb_mpk_), intent(in)              :: n
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: scratch
    end subroutine i_base_asb_m
  end interface


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

  interface 
    module subroutine i_base_asb_e(n, x, info, scratch)
      integer(psb_epk_), intent(in)              :: n
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: scratch
    end subroutine i_base_asb_e
  end interface


  !
  !> Function  base_free:
  !! \memberof  psb_i_base_vect_type
  !! \brief Free vector
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine i_base_free(x, info)
      class(psb_i_base_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_free
  end interface


  !
  !> Function  base_free_buffer:
  !! \memberof  psb_i_base_vect_type
  !! \brief Free aux buffer
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine i_base_free_buffer(x,info)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_free_buffer
  end interface


  !
  !> Function  base_maybe_free_buffer:
  !! \memberof  psb_i_base_vect_type
  !! \brief Conditionally Free aux buffer.
  !!        In some derived classes, e.g. GPU,
  !!        does not really frees to avoid  runtime
  !!        costs
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine i_base_maybe_free_buffer(x,info)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_maybe_free_buffer
  end interface


  !
  !> Function  base_free_comid:
  !! \memberof  psb_i_base_vect_type
  !! \brief Free aux MPI communication id buffer
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine i_base_free_comid(x,info)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_free_comid
  end interface


  interface 
    module function i_base_get_ncfs(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_get_ncfs
  end interface

  interface 
    module function i_base_get_dupl(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_get_dupl
  end interface

  interface 
    module function i_base_get_state(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_get_state
  end interface

  interface 
    module function i_base_is_null(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      logical :: res
    end function i_base_is_null
  end interface

  interface 
    module function i_base_is_bld(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      logical :: res
    end function i_base_is_bld
  end interface

  interface 
    module function i_base_is_upd(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      logical :: res
    end function i_base_is_upd
  end interface

  interface 
    module function i_base_is_asb(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      logical :: res
    end function i_base_is_asb
  end interface

  interface 
    module subroutine  i_base_set_ncfs(n,x)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine i_base_set_ncfs
  end interface

  interface 
    module subroutine  i_base_set_dupl(n,x)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine i_base_set_dupl
  end interface

  interface 
    module subroutine  i_base_set_state(n,x)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine i_base_set_state
  end interface

  interface 
    module subroutine  i_base_set_null(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_set_null
  end interface

  interface 
    module subroutine  i_base_set_bld(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_set_bld
  end interface

  interface 
    module subroutine  i_base_set_upd(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_set_upd
  end interface

  interface 
    module subroutine  i_base_set_asb(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_set_asb
  end interface

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
  interface 
    module subroutine i_base_sync(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_sync
  end interface

  !
  !> Function  base_set_host:
  !! \memberof  psb_i_base_vect_type
  !! \brief Set_host: base version is a no-op.
  !!
  !
  interface 
    module subroutine i_base_set_host(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_set_host
  end interface

  !
  !> Function  base_set_dev:
  !! \memberof  psb_i_base_vect_type
  !! \brief Set_dev: base version is a no-op.
  !!
  !
  interface 
    module subroutine i_base_set_dev(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_set_dev
  end interface

  !
  !> Function  base_set_sync:
  !! \memberof  psb_i_base_vect_type
  !! \brief Set_sync: base version is a no-op.
  !!
  !
  interface 
    module subroutine i_base_set_sync(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_set_sync
  end interface

  !
  !> Function  base_is_dev:
  !! \memberof  psb_i_base_vect_type
  !! \brief Is  vector on external device    .
  !!
  !
  interface 
    module function i_base_is_dev(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      logical  :: res
    end function i_base_is_dev
  end interface

  !
  !> Function  base_is_host
  !! \memberof  psb_i_base_vect_type
  !! \brief Is  vector on standard memory    .
  !!
  !
  interface 
    module function i_base_is_host(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      logical  :: res
    end function i_base_is_host
  end interface

  !
  !> Function  base_is_sync
  !! \memberof  psb_i_base_vect_type
  !! \brief Is  vector on sync               .
  !!
  !
  interface 
    module function i_base_is_sync(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      logical  :: res
    end function i_base_is_sync
  end interface

  !> Function  base_cpy:
  !! \memberof  psb_d_base_vect_type
  !! \brief     base_cpy: copy base contents
  !!  \param    y returned variable
  !!
  interface 
    module subroutine i_base_cpy(x, y)
      class(psb_i_base_vect_type), intent(in)   :: x
      class(psb_i_base_vect_type), intent(out)  :: y
    end subroutine i_base_cpy
  end interface

  !
  ! Size info.
  !
  !
  !> Function  base_get_nrows
  !! \memberof  psb_i_base_vect_type
  !! \brief  Number of entries
  !!
  !
  interface 
    module function i_base_get_nrows(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_get_nrows
  end interface

  !
  !> Function  base_get_sizeof
  !! \memberof  psb_i_base_vect_type
  !! \brief  Size in bytes
  !!
  !
  interface 
    module function i_base_sizeof(x) result(res)
      class(psb_i_base_vect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function i_base_sizeof
  end interface

  !
  !> Function  base_get_fmt
  !! \memberof  psb_i_base_vect_type
  !! \brief  Format
  !!
  !
  interface 
    module function i_base_get_fmt() result(res)
      character(len=5) :: res
    end function i_base_get_fmt
  end interface

  !
  !
  !
  !> Function  base_get_vect
  !! \memberof  psb_i_base_vect_type
  !! \brief  Extract a copy of the contents
  !!
  !
  interface 
    module function  i_base_get_vect(x,n) result(res)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), allocatable                 :: res(:)
      integer(psb_ipk_), optional :: n
    end function i_base_get_vect
  end interface
  
  !
  ! Reset all values
  !
  !
  !> Function  base_set_scal
  !! \memberof  psb_i_base_vect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  interface 
    module subroutine i_base_set_scal(x,val,first,last)
      class(psb_i_base_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in) :: val
      integer(psb_ipk_), optional :: first, last
    end subroutine i_base_set_scal
  end interface

  !
  !> Function  base_set_vect
  !! \memberof  psb_i_base_vect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in
  !!
  interface 
    module subroutine i_base_set_vect(x,val,first,last)
      class(psb_i_base_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in) :: val(:)
      integer(psb_ipk_), optional :: first, last
    end subroutine i_base_set_vect
  end interface
  
  interface 
    module subroutine i_base_check_addr(x)
      class(psb_i_base_vect_type), intent(inout) :: x
    end subroutine i_base_check_addr
  end interface
  

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
  interface 
    module subroutine i_base_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: alpha, beta, y(:)
      class(psb_i_base_vect_type) :: x
    end subroutine i_base_gthab
  end interface
  
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_gthzv
  !! \memberof  psb_i_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine i_base_gthzv_x(i,n,idx,x,y)
      integer(psb_ipk_) :: i
      integer(psb_mpk_) :: n
      class(psb_i_base_vect_type) :: idx
      integer(psb_ipk_) ::  y(:)
      class(psb_i_base_vect_type) :: x
    end subroutine i_base_gthzv_x
  end interface
  
  !
  ! New comm internals impl.
  !
  interface 
    module subroutine i_base_gthzbuf(i,n,idx,x)
      integer(psb_ipk_) :: i
      integer(psb_mpk_) :: n
      class(psb_i_base_vect_type) :: idx
      class(psb_i_base_vect_type) :: x
    end subroutine i_base_gthzbuf
  end interface
  
  !
  !> Function  base_device_wait:
  !! \memberof  psb_i_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!
  !
  interface 
    module subroutine i_base_device_wait()
    end subroutine i_base_device_wait
  end interface

  interface 
    module function i_base_use_buffer() result(res)
      logical :: res
    end function i_base_use_buffer
  end interface
  
  interface 
    module subroutine i_base_new_buffer(n,x,info)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: n
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_new_buffer
  end interface

  interface 
    module subroutine i_base_new_comid(n,x,info)
      class(psb_i_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: n
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_new_comid
  end interface

  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_gthzv
  !! \memberof  psb_i_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine i_base_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) ::  y(:)
      class(psb_i_base_vect_type) :: x
    end subroutine i_base_gthzv
  end interface

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
  interface 
    module subroutine i_base_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: beta, x(:)
      class(psb_i_base_vect_type) :: y
    end subroutine i_base_sctb
  end interface
  
  interface 
    module subroutine i_base_sctb_x(i,n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      integer(psb_ipk_) :: beta, x(:)
      class(psb_i_base_vect_type) :: y
    end subroutine i_base_sctb_x
  end interface
  
  interface 
    module subroutine i_base_sctb_buf(i,n,idx,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      integer(psb_ipk_) :: beta
      class(psb_i_base_vect_type) :: y
    end subroutine i_base_sctb_buf
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

end module psb_i_base_vect_mod


module psb_i_base_multivect_mod

  use psb_error_mod
  use psb_realloc_mod
  use psb_i_base_vect_mod

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

  type psb_i_base_multivect_type
    !> Values.
    integer(psb_ipk_), allocatable :: v(:,:)
    integer(psb_ipk_), allocatable :: combuf(:)
    integer(psb_mpk_), allocatable :: comid(:,:)
    !> vector bldstate:
    !!    null:   pristine;
    !!    build:  it's being filled with entries;
    !!    assembled: ready to use in computations;
    !!    update: accepts coefficients but only
    !!            in already existing entries.
    !!    The transitions among the states are detailed in
    !!            psb_T_vect_mod.
    integer(psb_ipk_), private :: bldstate = psb_vect_null_
    integer(psb_ipk_), private :: dupl     = psb_dupl_null_
    integer(psb_ipk_), private :: ncfs     = 0
    integer(psb_ipk_), allocatable :: iv(:)
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
    procedure, pass(x) :: reinit      => i_base_mlv_reinit
    procedure, pass(x) :: set_ncfs => i_base_mlv_set_ncfs 
    procedure, pass(x) :: get_ncfs => i_base_mlv_get_ncfs
    procedure, pass(x) :: set_dupl => i_base_mlv_set_dupl 
    procedure, pass(x) :: get_dupl => i_base_mlv_get_dupl
    procedure, pass(x) :: set_state => i_base_mlv_set_state
    procedure, pass(x) :: set_null  => i_base_mlv_set_null
    procedure, pass(x) :: set_bld   => i_base_mlv_set_bld
    procedure, pass(x) :: set_upd   => i_base_mlv_set_upd
    procedure, pass(x) :: set_asb   => i_base_mlv_set_asb
    procedure, pass(x) :: get_state => i_base_mlv_get_state
    procedure, pass(x) :: is_null   => i_base_mlv_is_null
    procedure, pass(x) :: is_bld    => i_base_mlv_is_bld
    procedure, pass(x) :: is_upd    => i_base_mlv_is_upd
    procedure, pass(x) :: is_asb    => i_base_mlv_is_asb
    procedure, pass(x) :: base_cpy  => i_base_mlv_cpy
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


    !
    ! These are for handling gather/scatter in new
    ! comm internals implementation.
    !
    procedure, nopass  :: use_buffer   => i_base_mlv_use_buffer
    procedure, pass(x) :: new_buffer   => i_base_mlv_new_buffer
    procedure, nopass  :: device_wait  => i_base_mlv_device_wait
    procedure, pass(x) :: maybe_free_buffer  => i_base_mlv_maybe_free_buffer
    procedure, pass(x) :: free_buffer  => i_base_mlv_free_buffer
    procedure, pass(x) :: new_comid    => i_base_mlv_new_comid
    procedure, pass(x) :: free_comid   => i_base_mlv_free_comid

    !
    ! Gather/scatter. These are needed for MPI interfacing.
    ! May have to be reworked.
    !
    procedure, pass(x) :: gthab    => i_base_mlv_gthab
    procedure, pass(x) :: gthzv    => i_base_mlv_gthzv
    procedure, pass(x) :: gthzm    => i_base_mlv_gthzm
    procedure, pass(x) :: gthzv_x  => i_base_mlv_gthzv_x
    procedure, pass(x) :: gthzbuf  => i_base_mlv_gthzbuf
    generic, public    :: gth      => gthab, gthzv, gthzm, gthzv_x, gthzbuf
    procedure, pass(y) :: sctb     => i_base_mlv_sctb
    procedure, pass(y) :: sctbr2   => i_base_mlv_sctbr2
    procedure, pass(y) :: sctb_x   => i_base_mlv_sctb_x
    procedure, pass(y) :: sctb_buf => i_base_mlv_sctb_buf
    generic, public    :: sct      => sctb, sctbr2, sctb_x, sctb_buf
  end type psb_i_base_multivect_type

  public  :: psb_i_base_multivect, psb_i_base_multivect_type

  interface psb_i_base_multivect
    module procedure constructor, size_const
  end interface psb_i_base_multivect

  private

  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_i_base_multivect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  interface 
    module subroutine i_base_mlv_bld_x(x,this)
      integer(psb_ipk_), intent(in) :: this(:,:)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_) :: info
    end subroutine i_base_mlv_bld_x
  end interface
  
  !
  ! Create with size, but no initialization
  !

  !> Function  bld_n:
  !! \memberof  psb_i_base_multivect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  interface 
    module subroutine i_base_mlv_bld_n(x,m,n,scratch)
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_) :: info
      logical, intent(in), optional        :: scratch
    end subroutine i_base_mlv_bld_n
  end interface
  
  !> Function  base_mlv_all:
  !! \memberof  psb_i_base_multivect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated.
  !!  \param info  return code
  !!
  interface 
    module subroutine i_base_mlv_all(m,n, x, info)
      integer(psb_ipk_), intent(in)               :: m,n
      class(psb_i_base_multivect_type), intent(out) :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_mlv_all
  end interface

  !> Function  base_mlv_mold:
  !! \memberof  psb_i_base_multivect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  interface 
    module subroutine i_base_mlv_mold(x, y, info)
      class(psb_i_base_multivect_type), intent(in)   :: x
      class(psb_i_base_multivect_type), intent(out), allocatable :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_mlv_mold
  end interface

  interface 
    module subroutine i_base_mlv_reinit(x, info)
      class(psb_i_base_multivect_type), intent(out)    :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_mlv_reinit
  end interface

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
  interface 
    module subroutine i_base_mlv_ins(n,irl,val,dupl,x,maxr,info)
      class(psb_i_base_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, dupl,maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      integer(psb_ipk_), intent(in)        :: val(:,:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_mlv_ins
  end interface

  !
  !> Function  base_mlv_zero
  !! \memberof  psb_i_base_multivect_type
  !! \brief Zero out contents
  !!
  !
  interface 
    module subroutine i_base_mlv_zero(x)
      class(psb_i_base_multivect_type), intent(inout)    :: x
    end subroutine i_base_mlv_zero
  end interface

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
  
  interface 
    module subroutine i_base_mlv_asb(m,n, x, info, scratch)
      integer(psb_ipk_), intent(in)              :: m,n
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: scratch
    end subroutine i_base_mlv_asb
  end interface

  !
  !> Function  base_mlv_free:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Free vector
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine i_base_mlv_free(x, info)
      class(psb_i_base_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine i_base_mlv_free
  end interface

  interface 
    module function i_base_mlv_get_ncfs(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_mlv_get_ncfs
  end interface

  interface 
    module function i_base_mlv_get_dupl(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_mlv_get_dupl
  end interface

  interface 
    module function i_base_mlv_get_state(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_mlv_get_state
  end interface

  interface 
    module function i_base_mlv_is_null(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      logical :: res
    end function i_base_mlv_is_null
  end interface

  interface 
    module function i_base_mlv_is_bld(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      logical :: res
    end function i_base_mlv_is_bld
  end interface

  interface 
    module function i_base_mlv_is_upd(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      logical :: res
    end function i_base_mlv_is_upd
  end interface

  interface 
    module function i_base_mlv_is_asb(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      logical :: res
    end function i_base_mlv_is_asb
  end interface

  interface 
    module subroutine  i_base_mlv_set_ncfs(n,x)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine i_base_mlv_set_ncfs
  end interface

  interface 
    module subroutine  i_base_mlv_set_dupl(n,x)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine i_base_mlv_set_dupl
  end interface

  interface 
    module subroutine  i_base_mlv_set_state(n,x)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine i_base_mlv_set_state
  end interface

  interface 
    module subroutine  i_base_mlv_set_null(x)
      class(psb_i_base_multivect_type), intent(inout) :: x
    end subroutine i_base_mlv_set_null
  end interface

  interface 
    module subroutine  i_base_mlv_set_bld(x)
      class(psb_i_base_multivect_type), intent(inout) :: x
    end subroutine i_base_mlv_set_bld
  end interface

  interface 
    module subroutine  i_base_mlv_set_upd(x)
      class(psb_i_base_multivect_type), intent(inout) :: x
    end subroutine i_base_mlv_set_upd
  end interface

  interface 
    module subroutine  i_base_mlv_set_asb(x)
      class(psb_i_base_multivect_type), intent(inout) :: x
    end subroutine i_base_mlv_set_asb
  end interface

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
  interface 
    module subroutine i_base_mlv_sync(x)
      class(psb_i_base_multivect_type), intent(inout) :: x
    end subroutine i_base_mlv_sync
  end interface

  !
  !> Function  base_mlv_set_host:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Set_host: base version is a no-op.
  !!
  !
  interface 
    module subroutine i_base_mlv_set_host(x)
      class(psb_i_base_multivect_type), intent(inout) :: x
    end subroutine i_base_mlv_set_host
  end interface

  !
  !> Function  base_mlv_set_dev:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Set_dev: base version is a no-op.
  !!
  !
  interface 
    module subroutine i_base_mlv_set_dev(x)
      class(psb_i_base_multivect_type), intent(inout) :: x
    end subroutine i_base_mlv_set_dev
  end interface

  !
  !> Function  base_mlv_set_sync:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Set_sync: base version is a no-op.
  !!
  !
  interface 
    module subroutine i_base_mlv_set_sync(x)
      class(psb_i_base_multivect_type), intent(inout) :: x
    end subroutine i_base_mlv_set_sync
  end interface

  !
  !> Function  base_mlv_is_dev:
  !! \memberof  psb_i_base_multivect_type
  !! \brief Is  vector on external device    .
  !!
  !
  interface 
    module function i_base_mlv_is_dev(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      logical  :: res
    end function i_base_mlv_is_dev
  end interface

  !
  !> Function  base_mlv_is_host
  !! \memberof  psb_i_base_multivect_type
  !! \brief Is  vector on standard memory    .
  !!
  !
  interface 
    module function i_base_mlv_is_host(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      logical  :: res
    end function i_base_mlv_is_host
  end interface

  !
  !> Function  base_mlv_is_sync
  !! \memberof  psb_i_base_multivect_type
  !! \brief Is  vector on sync               .
  !!
  !
  interface 
    module function i_base_mlv_is_sync(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      logical  :: res
    end function i_base_mlv_is_sync
  end interface

  !> Function  base_cpy:
  !! \memberof  psb_d_base_vect_type
  !! \brief     base_cpy: copy base contents
  !!  \param    y returned variable
  !!
  interface 
    module subroutine i_base_mlv_cpy(x, y)
      class(psb_i_base_multivect_type), intent(in)   :: x
      class(psb_i_base_multivect_type), intent(out)  :: y
    end subroutine i_base_mlv_cpy
  end interface

  ! Size info.
  !
  !
  !> Function  base_mlv_get_nrows
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Number of entries
  !!
  !
  interface 
    module function i_base_mlv_get_nrows(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_mlv_get_nrows
  end interface

  interface 
    module function i_base_mlv_get_ncols(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function i_base_mlv_get_ncols
  end interface

  !
  !> Function  base_mlv_get_sizeof
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Size in bytesa
  !!
  !
  interface 
    module function i_base_mlv_sizeof(x) result(res)
      class(psb_i_base_multivect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function i_base_mlv_sizeof
  end interface

  !
  !> Function  base_mlv_get_fmt
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Format
  !!
  !
  interface 
    module function i_base_mlv_get_fmt() result(res)
      character(len=5) :: res
    end function i_base_mlv_get_fmt
  end interface

  !
  !
  !
  !> Function  base_mlv_get_vect
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Extract a copy of the contents
  !!
  !
  interface 
    module function  i_base_mlv_get_vect(x) result(res)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), allocatable                 :: res(:,:)
    end function i_base_mlv_get_vect
  end interface

  !
  ! Reset all values
  !
  !
  !> Function  base_mlv_set_scal
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  interface 
    module subroutine i_base_mlv_set_scal(x,val)
      class(psb_i_base_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in) :: val
    end subroutine i_base_mlv_set_scal
  end interface
  
  !
  !> Function  base_mlv_set_vect
  !! \memberof  psb_i_base_multivect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in
  !!
  interface 
    module subroutine i_base_mlv_set_vect(x,val)
      class(psb_i_base_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in) :: val(:,:)
    end subroutine i_base_mlv_set_vect
  end interface


  interface 
    module function i_base_mlv_use_buffer() result(res)
      logical :: res
    end function i_base_mlv_use_buffer
  end interface
  
  interface 
    module subroutine i_base_mlv_new_buffer(n,x,info)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: n
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_mlv_new_buffer
  end interface
  
  interface 
    module subroutine i_base_mlv_new_comid(n,x,info)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: n
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_mlv_new_comid
  end interface
  
  interface 
    module subroutine i_base_mlv_maybe_free_buffer(x,info)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_mlv_maybe_free_buffer
  end interface
  
  interface 
    module subroutine i_base_mlv_free_buffer(x,info)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_mlv_free_buffer
  end interface

  interface 
    module subroutine i_base_mlv_free_comid(x,info)
      class(psb_i_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine i_base_mlv_free_comid
  end interface

  !
  ! Gather: Y = beta * Y + alpha * X(IDX(:))
  !
  !
  !> Function  base_mlv_gthab
  !! \memberof  psb_i_base_multivect_type
  !! \brief gather into an array
  !!    Y = beta * Y + alpha * X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param alpha
  !! \param beta
  interface 
    module subroutine i_base_mlv_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: alpha, beta, y(:)
      class(psb_i_base_multivect_type) :: x
    end subroutine i_base_mlv_gthab
  end interface
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_i_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine i_base_mlv_gthzv_x(i,n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      integer(psb_ipk_) ::  y(:)
      class(psb_i_base_multivect_type) :: x
    end subroutine i_base_mlv_gthzv_x
  end interface

  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_i_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine i_base_mlv_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) ::  y(:)
      class(psb_i_base_multivect_type) :: x
    end subroutine i_base_mlv_gthzv
  end interface
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_i_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine i_base_mlv_gthzm(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) ::  y(:,:)
      class(psb_i_base_multivect_type) :: x
    end subroutine i_base_mlv_gthzm
  end interface

  !
  ! New comm internals impl.
  !
  interface 
    module subroutine i_base_mlv_gthzbuf(i,ixb,n,idx,x)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i, ixb
      class(psb_i_base_vect_type) :: idx
      class(psb_i_base_multivect_type) :: x
    end subroutine i_base_mlv_gthzbuf
  end interface
  
  !
  ! Scatter:
  ! Y(IDX(:),:) = beta*Y(IDX(:),:) + X(:)
  !
  !
  !> Function  base_mlv_sctb
  !! \memberof  psb_i_base_multivect_type
  !! \brief scatter into a class(base_mlv_vect)
  !!    Y(IDX(:)) = beta * Y(IDX(:)) +  X(:)
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param beta
  !! \param x(:)
  interface 
    module subroutine i_base_mlv_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: beta, x(:)
      class(psb_i_base_multivect_type) :: y
    end subroutine i_base_mlv_sctb
  end interface

  interface 
    module subroutine i_base_mlv_sctbr2(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      integer(psb_ipk_) :: beta, x(:,:)
      class(psb_i_base_multivect_type) :: y
    end subroutine i_base_mlv_sctbr2
  end interface
  
  interface 
    module subroutine i_base_mlv_sctb_x(i,n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      integer( psb_ipk_) :: beta, x(:)
      class(psb_i_base_multivect_type) :: y
    end subroutine i_base_mlv_sctb_x
  end interface

  interface 
    module subroutine i_base_mlv_sctb_buf(i,iyb,n,idx,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i, iyb
      class(psb_i_base_vect_type) :: idx
      integer(psb_ipk_) :: beta
      class(psb_i_base_multivect_type) :: y
    end subroutine i_base_mlv_sctb_buf
  end interface

  !
  !> Function  base_device_wait:
  !! \memberof  psb_i_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!
  !
  interface 
    module subroutine i_base_mlv_device_wait()
    end subroutine i_base_mlv_device_wait
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
    call this%asb(size(x,dim=1,kind=psb_ipk_),&
         & size(x,dim=2,kind=psb_ipk_),info)
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

end module psb_i_base_multivect_mod
