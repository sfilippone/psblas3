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
! package: psb_s_base_vect_mod
!
! This module contains the definition of the psb_s_base_vect type which
! is a container for dense vectors.
!  This is encapsulated instead of being just a simple array to allow for
!  more complicated situations, such as GPU programming, where the memory
!  area we are interested in is not easily accessible from the host/Fortran
!  side. It is also meant to be encapsulated in an outer type, to allow
!  runtime switching as per the STATE design pattern, similar to the
!  sparse matrix types.
!
!
module psb_s_base_vect_mod

  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_i_base_vect_mod
  use psb_l_base_vect_mod

  !> \namespace  psb_base_mod  \class psb_s_base_vect_type
  !! The psb_s_base_vect_type
  !! defines a middle level  real(psb_spk_) encapsulated dense vector.
  !! The encapsulation is needed, in place of a simple array, to allow
  !! for complicated situations, such as GPU programming, where the memory
  !!  area we are interested in is not easily accessible from the host/Fortran
  !!  side. It is also meant to be encapsulated in an outer type, to allow
  !!  runtime switching as per the STATE design pattern, similar to the
  !!  sparse matrix types.
  !!
  type psb_s_base_vect_type
    !> Values.
    real(psb_spk_), allocatable :: v(:)
    real(psb_spk_), allocatable :: combuf(:)
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
    procedure, pass(x) :: bld_x    => s_base_bld_x
    procedure, pass(x) :: bld_mn   => s_base_bld_mn
    procedure, pass(x) :: bld_en   => s_base_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
    procedure, pass(x) :: all      => s_base_all
    procedure, pass(x) :: mold     => s_base_mold
    !
    ! Insert/set. Assembly and free.
    ! Assembly does almost nothing here, but is important
    ! in derived classes.
    !
    procedure, pass(x) :: ins_a    => s_base_ins_a
    procedure, pass(x) :: ins_v    => s_base_ins_v
    generic, public    :: ins      => ins_a, ins_v
    procedure, pass(x) :: zero     => s_base_zero
    procedure, pass(x) :: asb_m    => s_base_asb_m
    procedure, pass(x) :: asb_e    => s_base_asb_e
    generic, public    :: asb      => asb_m, asb_e
    procedure, pass(x) :: free     => s_base_free
    procedure, pass(x) :: reinit      => s_base_reinit
    procedure, pass(x) :: set_ncfs => s_base_set_ncfs 
    procedure, pass(x) :: get_ncfs => s_base_get_ncfs
    procedure, pass(x) :: set_dupl => s_base_set_dupl 
    procedure, pass(x) :: get_dupl => s_base_get_dupl
    procedure, pass(x) :: set_state => s_base_set_state
    procedure, pass(x) :: set_null  => s_base_set_null
    procedure, pass(x) :: set_bld   => s_base_set_bld
    procedure, pass(x) :: set_upd   => s_base_set_upd
    procedure, pass(x) :: set_asb   => s_base_set_asb
    procedure, pass(x) :: get_state => s_base_get_state
    procedure, pass(x) :: is_null   => s_base_is_null
    procedure, pass(x) :: is_bld    => s_base_is_bld
    procedure, pass(x) :: is_upd    => s_base_is_upd
    procedure, pass(x) :: is_asb    => s_base_is_asb
    procedure, pass(x) :: base_cpy  => s_base_cpy
    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder.
    !
    procedure, pass(x) :: sync     => s_base_sync
    procedure, pass(x) :: is_host  => s_base_is_host
    procedure, pass(x) :: is_dev   => s_base_is_dev
    procedure, pass(x) :: is_sync  => s_base_is_sync
    procedure, pass(x) :: set_host => s_base_set_host
    procedure, pass(x) :: set_dev  => s_base_set_dev
    procedure, pass(x) :: set_sync => s_base_set_sync

    !
    ! These are for handling gather/scatter in new
    ! comm internals implementation.
    !
    procedure, nopass  :: use_buffer   => s_base_use_buffer
    procedure, pass(x) :: new_buffer   => s_base_new_buffer
    procedure, nopass  :: device_wait  => s_base_device_wait
    procedure, pass(x) :: maybe_free_buffer  => s_base_maybe_free_buffer
    procedure, pass(x) :: free_buffer  => s_base_free_buffer
    procedure, pass(x) :: new_comid    => s_base_new_comid
    procedure, pass(x) :: free_comid   => s_base_free_comid

    !
    ! Basic info
    procedure, pass(x) :: get_nrows => s_base_get_nrows
    procedure, pass(x) :: sizeof    => s_base_sizeof
    procedure, nopass  :: get_fmt   => s_base_get_fmt
    !
    ! Set/get data from/to an external array; also
    ! overload assignment.
    !
    procedure, pass(x) :: get_vect => s_base_get_vect
    procedure, pass(x) :: set_scal => s_base_set_scal
    procedure, pass(x) :: set_vect => s_base_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: get_entry=> s_base_get_entry
    procedure, pass(x) :: set_entry=> s_base_set_entry
    !
    ! Gather/scatter. These are needed for MPI interfacing.
    ! May have to be reworked.
    !
    procedure, pass(x) :: gthab    => s_base_gthab
    procedure, pass(x) :: gthzv    => s_base_gthzv
    procedure, pass(x) :: gthzv_x  => s_base_gthzv_x
    procedure, pass(x) :: gthzbuf  => s_base_gthzbuf
    generic, public    :: gth      => gthab, gthzv, gthzv_x, gthzbuf
    procedure, pass(y) :: sctb     => s_base_sctb
    procedure, pass(y) :: sctb_x   => s_base_sctb_x
    procedure, pass(y) :: sctb_buf => s_base_sctb_buf
    generic, public    :: sct      => sctb, sctb_x, sctb_buf

    procedure, pass(x) :: check_addr  => s_base_check_addr



    !
    ! Dot product and AXPBY
    !
    procedure, pass(x) :: dot_v    => s_base_dot_v
    procedure, pass(x) :: dot_a    => s_base_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => s_base_axpby_v
    procedure, pass(y) :: axpby_a  => s_base_axpby_a
    procedure, pass(z) :: axpby_v2  => s_base_axpby_v2
    procedure, pass(z) :: axpby_a2  => s_base_axpby_a2
    generic, public    :: axpby    => axpby_v, axpby_a, axpby_v2, axpby_a2
    procedure, pass(z) :: upd_xyz  => s_base_upd_xyz
    procedure, pass(w) :: xyzw     => s_base_xyzw
    
    !
    ! Vector by vector multiplication. Need all variants
    ! to handle multiple requirements from preconditioners
    !
    procedure, pass(y) :: mlt_v    => s_base_mlt_v
    procedure, pass(y) :: mlt_a    => s_base_mlt_a
    procedure, pass(z) :: mlt_a_2  => s_base_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => s_base_mlt_v_2
    procedure, pass(z) :: mlt_va   => s_base_mlt_va
    procedure, pass(z) :: mlt_av   => s_base_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2, mlt_v_2, mlt_av, &
                                        mlt_va
    !
    ! Vector-Vector operations
    !
    procedure, pass(y) :: div_v         => s_base_div_v
    procedure, pass(y) :: div_a         => s_base_div_a
    procedure, pass(y) :: div_v_check   => s_base_div_v_check
    procedure, pass(z) :: div_v2         => s_base_div_v2
    procedure, pass(z) :: div_v2_check   => s_base_div_v2_check
    procedure, pass(z) :: div_a2        => s_base_div_a2
    procedure, pass(z) :: div_a2_check  => s_base_div_a2_check
    generic, public    :: div           => div_v, div_v2, div_v_check, div_a, &
                                            div_v2_check, div_a2, div_a2_check
    procedure, pass(y) :: inv_v    => s_base_inv_v
    procedure, pass(y) :: inv_v_check => s_base_inv_v_check
    procedure, pass(y) :: inv_a2   => s_base_inv_a2
    procedure, pass(y) :: inv_a2_check => s_base_inv_a2_check
    generic, public    :: inv      => inv_v, inv_v_check, inv_a2, inv_a2_check
    !
    ! Scaling and norms
    !
    procedure, pass(x) :: scal     => s_base_scal
    procedure, pass(x) :: absval1  => s_base_absval1
    procedure, pass(x) :: absval2  => s_base_absval2
    generic, public    :: absval   => absval1, absval2
    procedure, pass(x) :: nrm2     => s_base_nrm2
    procedure, pass(x) :: amax     => s_base_amax
    procedure, pass(x) :: asum     => s_base_asum

    !
    ! Comparison and mask operation
    !
    procedure, pass(z) :: acmp_a2   => s_base_acmp_a2
    procedure, pass(z) :: acmp_v2   => s_base_acmp_v2
    generic, public    :: acmp      => acmp_a2,acmp_v2
    !
    ! Add constant value to all entry of a vector
    !
    procedure, pass(z) :: addconst_a2   => s_base_addconst_a2
    procedure, pass(z) :: addconst_v2   => s_base_addconst_v2
    generic, public    :: addconst      => addconst_a2,addconst_v2

    procedure, pass(x) :: minreal    => s_base_min
    procedure, pass(m) :: mask_v => s_base_mask_v
    procedure, pass(m) :: mask_a => s_base_mask_a
    generic, public    :: mask => mask_a, mask_v
    procedure, pass(x) :: minquotient_v  => s_base_minquotient_v
    procedure, pass(x) :: minquotient_a2 => s_base_minquotient_a2
    generic, public    :: minquotient    => minquotient_v, minquotient_a2

  end type psb_s_base_vect_type

  public  :: psb_s_base_vect, psb_s_base_vect_type
  private :: constructor, size_const
  interface psb_s_base_vect
    module procedure constructor, size_const
  end interface psb_s_base_vect

  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_s_base_vect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  interface 
    module subroutine s_base_bld_x(x,this,scratch)
      real(psb_spk_), intent(in) :: this(:)
      class(psb_s_base_vect_type), intent(inout) :: x
      logical, intent(in), optional        :: scratch
    end subroutine s_base_bld_x
  end interface

  !
  ! Create with size, but no initialization
  !

  !> Function  bld_mn:
  !! \memberof  psb_s_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  interface 
    module subroutine s_base_bld_mn(x,n,scratch)
      integer(psb_mpk_), intent(in) :: n
      class(psb_s_base_vect_type), intent(inout) :: x
      logical, intent(in), optional        :: scratch
    end subroutine s_base_bld_mn
  end interface
  
  !> Function  bld_en:
  !! \memberof  psb_s_base_vect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  interface 
    module subroutine s_base_bld_en(x,n,scratch)
      integer(psb_epk_), intent(in) :: n
      class(psb_s_base_vect_type), intent(inout) :: x
      logical, intent(in), optional        :: scratch
    end subroutine s_base_bld_en
  end interface

  !> Function  base_all:
  !! \memberof  psb_s_base_vect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated.
  !!  \param info  return code
  !!
  interface 
    module subroutine s_base_all(n, x, info)
      integer(psb_ipk_), intent(in)               :: n
      class(psb_s_base_vect_type), intent(out)    :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_all
  end interface

  !> Function  base_mold:
  !! \memberof  psb_s_base_vect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  interface 
    module subroutine s_base_mold(x, y, info)
      class(psb_s_base_vect_type), intent(in)   :: x
      class(psb_s_base_vect_type), intent(out), allocatable :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mold
  end interface

  interface 
    module subroutine s_base_reinit(x, info,clear)
      class(psb_s_base_vect_type), intent(inout)    :: x
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in), optional               :: clear
    end subroutine s_base_reinit
  end interface

  !
  ! Insert a bunch of values at specified positions.
  !
  !> Function  base_ins:
  !! \memberof  psb_s_base_vect_type
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
    module subroutine s_base_ins_a(n,irl,val,dupl,x,maxr,info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, dupl, maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      real(psb_spk_), intent(in)        :: val(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_ins_a
  end interface


  interface 
    module subroutine s_base_ins_v(n,irl,val,dupl,x,maxr,info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, dupl, maxr
      class(psb_i_base_vect_type), intent(inout)  :: irl
      class(psb_s_base_vect_type), intent(inout)  :: val
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_ins_v
  end interface



  !
  !> Function  base_zero
  !! \memberof  psb_s_base_vect_type
  !! \brief Zero out contents
  !!
  !
  interface 
    module subroutine s_base_zero(x)
      class(psb_s_base_vect_type), intent(inout)    :: x
    end subroutine s_base_zero
  end interface



  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_asb:
  !! \memberof  psb_s_base_vect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !

  interface 
    module subroutine s_base_asb_m(n, x, info, scratch)
      integer(psb_mpk_), intent(in)              :: n
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: scratch
    end subroutine s_base_asb_m
  end interface


  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_asb:
  !! \memberof  psb_s_base_vect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !

  interface 
    module subroutine s_base_asb_e(n, x, info, scratch)
      integer(psb_epk_), intent(in)              :: n
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: scratch
    end subroutine s_base_asb_e
  end interface


  !
  !> Function  base_free:
  !! \memberof  psb_s_base_vect_type
  !! \brief Free vector
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine s_base_free(x, info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_free
  end interface


  !
  !> Function  base_free_buffer:
  !! \memberof  psb_s_base_vect_type
  !! \brief Free aux buffer
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine s_base_free_buffer(x,info)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_free_buffer
  end interface


  !
  !> Function  base_maybe_free_buffer:
  !! \memberof  psb_s_base_vect_type
  !! \brief Conditionally Free aux buffer.
  !!        In some derived classes, e.g. GPU,
  !!        does not really frees to avoid  runtime
  !!        costs
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine s_base_maybe_free_buffer(x,info)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_maybe_free_buffer
  end interface


  !
  !> Function  base_free_comid:
  !! \memberof  psb_s_base_vect_type
  !! \brief Free aux MPI communication id buffer
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine s_base_free_comid(x,info)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_free_comid
  end interface


  interface 
    module function s_base_get_ncfs(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_get_ncfs
  end interface

  interface 
    module function s_base_get_dupl(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_get_dupl
  end interface

  interface 
    module function s_base_get_state(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_get_state
  end interface

  interface 
    module function s_base_is_null(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      logical :: res
    end function s_base_is_null
  end interface

  interface 
    module function s_base_is_bld(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      logical :: res
    end function s_base_is_bld
  end interface

  interface 
    module function s_base_is_upd(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      logical :: res
    end function s_base_is_upd
  end interface

  interface 
    module function s_base_is_asb(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      logical :: res
    end function s_base_is_asb
  end interface

  interface 
    module subroutine  s_base_set_ncfs(n,x)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine s_base_set_ncfs
  end interface

  interface 
    module subroutine  s_base_set_dupl(n,x)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine s_base_set_dupl
  end interface

  interface 
    module subroutine  s_base_set_state(n,x)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine s_base_set_state
  end interface

  interface 
    module subroutine  s_base_set_null(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_set_null
  end interface

  interface 
    module subroutine  s_base_set_bld(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_set_bld
  end interface

  interface 
    module subroutine  s_base_set_upd(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_set_upd
  end interface

  interface 
    module subroutine  s_base_set_asb(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_set_asb
  end interface

  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  !
  !
  !> Function  base_sync:
  !! \memberof  psb_s_base_vect_type
  !! \brief Sync: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_sync(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_sync
  end interface

  !
  !> Function  base_set_host:
  !! \memberof  psb_s_base_vect_type
  !! \brief Set_host: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_set_host(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_set_host
  end interface

  !
  !> Function  base_set_dev:
  !! \memberof  psb_s_base_vect_type
  !! \brief Set_dev: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_set_dev(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_set_dev
  end interface

  !
  !> Function  base_set_sync:
  !! \memberof  psb_s_base_vect_type
  !! \brief Set_sync: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_set_sync(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_set_sync
  end interface

  !
  !> Function  base_is_dev:
  !! \memberof  psb_s_base_vect_type
  !! \brief Is  vector on external device    .
  !!
  !
  interface 
    module function s_base_is_dev(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      logical  :: res
    end function s_base_is_dev
  end interface

  !
  !> Function  base_is_host
  !! \memberof  psb_s_base_vect_type
  !! \brief Is  vector on standard memory    .
  !!
  !
  interface 
    module function s_base_is_host(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      logical  :: res
    end function s_base_is_host
  end interface

  !
  !> Function  base_is_sync
  !! \memberof  psb_s_base_vect_type
  !! \brief Is  vector on sync               .
  !!
  !
  interface 
    module function s_base_is_sync(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      logical  :: res
    end function s_base_is_sync
  end interface

  !> Function  base_cpy:
  !! \memberof  psb_d_base_vect_type
  !! \brief     base_cpy: copy base contents
  !!  \param    y returned variable
  !!
  interface 
    module subroutine s_base_cpy(x, y)
      class(psb_s_base_vect_type), intent(in)   :: x
      class(psb_s_base_vect_type), intent(out)  :: y
    end subroutine s_base_cpy
  end interface

  !
  ! Size info.
  !
  !
  !> Function  base_get_nrows
  !! \memberof  psb_s_base_vect_type
  !! \brief  Number of entries
  !!
  !
  interface 
    module function s_base_get_nrows(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_get_nrows
  end interface

  !
  !> Function  base_get_sizeof
  !! \memberof  psb_s_base_vect_type
  !! \brief  Size in bytes
  !!
  !
  interface 
    module function s_base_sizeof(x) result(res)
      class(psb_s_base_vect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function s_base_sizeof
  end interface

  !
  !> Function  base_get_fmt
  !! \memberof  psb_s_base_vect_type
  !! \brief  Format
  !!
  !
  interface 
    module function s_base_get_fmt() result(res)
      character(len=5) :: res
    end function s_base_get_fmt
  end interface

  !
  !
  !
  !> Function  base_get_vect
  !! \memberof  psb_s_base_vect_type
  !! \brief  Extract a copy of the contents
  !!
  !
  interface 
    module function  s_base_get_vect(x,n) result(res)
      class(psb_s_base_vect_type), intent(inout) :: x
      real(psb_spk_), allocatable                 :: res(:)
      integer(psb_ipk_), optional :: n
    end function s_base_get_vect
  end interface
  
  !
  ! Reset all values
  !
  !
  !> Function  base_set_scal
  !! \memberof  psb_s_base_vect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  interface 
    module subroutine s_base_set_scal(x,val,first,last)
      class(psb_s_base_vect_type), intent(inout)  :: x
      real(psb_spk_), intent(in) :: val
      integer(psb_ipk_), optional :: first, last
    end subroutine s_base_set_scal
  end interface

  !
  !> Function  base_set_vect
  !! \memberof  psb_s_base_vect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in
  !!
  interface 
    module subroutine s_base_set_vect(x,val,first,last)
      class(psb_s_base_vect_type), intent(inout)  :: x
      real(psb_spk_), intent(in) :: val(:)
      integer(psb_ipk_), optional :: first, last
    end subroutine s_base_set_vect
  end interface
  
  interface 
    module subroutine s_base_check_addr(x)
      class(psb_s_base_vect_type), intent(inout) :: x
    end subroutine s_base_check_addr
  end interface
  
  !
  ! Get entry.
  !
  !
  !> Function  base_get_entry
  !! \memberof  psb_s_base_vect_type
  !! \brief  Get one entry from the vector
  !!
  !
  interface 
    module function s_base_get_entry(x, index) result(res)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: index
      real(psb_spk_)                           :: res
    end function s_base_get_entry
  end interface
  
  interface 
    module subroutine s_base_set_entry(x, index, val)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)             :: index
      real(psb_spk_)                           :: val
    end subroutine s_base_set_entry
  end interface
  
  !
  ! Overwrite with absolute value
  !
  !
  !> Function  base_absval1
  !! \memberof  psb_s_base_vect_type
  !! \brief  Set all entries to their respective absolute values.
  !!
  interface 
    module subroutine s_base_absval1(x)
      class(psb_s_base_vect_type), intent(inout)  :: x
    end subroutine s_base_absval1
  end interface
  
  interface 
    module subroutine s_base_absval2(x,y)
      class(psb_s_base_vect_type), intent(inout) :: x
      class(psb_s_base_vect_type), intent(inout) :: y
    end subroutine s_base_absval2
  end interface
  
  !
  ! Dot products
  !
  !
  !> Function  base_dot_v
  !! \memberof  psb_s_base_vect_type
  !! \brief  Dot product by another base_vector
  !! \param n    Number of entries to be considered
  !! \param y    The other (base_vect) to be multiplied by
  !!
  interface 
    module function s_base_dot_v(n,x,y) result(res)
      class(psb_s_base_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function
  end interface

  !
  ! Base workhorse is good old BLAS1
  !
  !
  !> Function  base_dot_a
  !! \memberof  psb_s_base_vect_type
  !! \brief  Dot product by a normal array
  !! \param n    Number of entries to be considered
  !! \param y(:) The array to be multiplied by
  !!
  interface 
    module function s_base_dot_a(n,x,y) result(res)
      class(psb_s_base_vect_type), intent(inout) :: x
      real(psb_spk_), intent(in)    :: y(:)
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function s_base_dot_a
  end interface

  !
  ! AXPBY is invoked via Y, hence the structure below.
  !
  !
  !
  !> Function  base_axpby_v
  !! \memberof  psb_s_base_vect_type
  !! \brief AXPBY  by a (base_vect) y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x     The class(base_vect) to be added
  !! \param beta scalar beta
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_axpby_v(m,alpha, x, beta, y, info)
      integer(psb_ipk_), intent(in)               :: m
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      real(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_axpby_v
  end interface

  !
  ! AXPBY is invoked via Z, hence the structure below.
  !
  !
  !
  !> Function  base_axpby_v2
  !! \memberof  psb_s_base_vect_type
  !! \brief AXPBY  by a (base_vect) z=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x     The class(base_vect) to be added
  !! \param beta scalar beta
  !! \param y     The class(base_vect) to be added
  !! \param z     The class(base_vect) to be returned
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_axpby_v2(m,alpha, x, beta, y, z, info)
      integer(psb_ipk_), intent(in)               :: m
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      class(psb_s_base_vect_type), intent(inout)  :: z
      real(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_axpby_v2
  end interface

  !
  ! AXPBY is invoked via Y, hence the structure below.
  !
  !
  !> Function  base_axpby_a
  !! \memberof  psb_s_base_vect_type
  !! \brief AXPBY  by a normal array y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x(:) The array to be added
  !! \param beta scalar beta
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_axpby_a(m,alpha, x, beta, y, info)
      integer(psb_ipk_), intent(in)               :: m
      real(psb_spk_), intent(in)        :: x(:)
      class(psb_s_base_vect_type), intent(inout)  :: y
      real(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_axpby_a
  end interface

  !
  ! AXPBY is invoked via Z, hence the structure below.
  !
  !
  !> Function  base_axpby_a2
  !! \memberof  psb_s_base_vect_type
  !! \brief AXPBY  by a normal array y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x(:) The array to be added
  !! \param beta scalar beta
  !! \param y(:) The array to be added
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_axpby_a2(m,alpha, x, beta, y, z, info)
      integer(psb_ipk_), intent(in)                :: m
      real(psb_spk_), intent(in)                  :: x(:)
      real(psb_spk_), intent(in)                  :: y(:)
      class(psb_s_base_vect_type), intent(inout) :: z
      real(psb_spk_), intent (in)                 :: alpha, beta
      integer(psb_ipk_), intent(out)               :: info
    end subroutine s_base_axpby_a2
  end interface

  !
  ! UPD_XYZ is invoked via Z, hence the structure below.
  !
  !
  !> Function  base_upd_xyz
  !! \memberof  psb_s_base_vect_type
  !! \brief UPD_XYZ combines two AXPBYS y=alpha*x+beta*y, z=gamma*y+delta*zeta
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param beta scalar beta 
  !! \param gamma scalar gamma
  !! \param delta scalar delta
  !! \param x     The class(base_vect) to be added
  !! \param y     The class(base_vect) to be added
  !! \param z     The class(base_vect) to be added
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_upd_xyz(m,alpha, beta, gamma,delta,x, y, z, info)
      integer(psb_ipk_), intent(in)               :: m
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      class(psb_s_base_vect_type), intent(inout)  :: z
      real(psb_spk_), intent (in)       :: alpha, beta, gamma, delta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_upd_xyz
  end interface

  interface 
    module subroutine s_base_xyzw(m,a,b,c,d,e,f,x, y, z, w,info)
      integer(psb_ipk_), intent(in)               :: m
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      class(psb_s_base_vect_type), intent(inout)  :: z
      class(psb_s_base_vect_type), intent(inout)  :: w
      real(psb_spk_), intent (in)                :: a,b,c,d,e,f
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_xyzw
  end interface
  
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
  !! \memberof  psb_s_base_vect_type
  !! \brief Vector entry-by-entry multiply  by a base_vect array y=x*y
  !! \param x   The class(base_vect) to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlt_v(x, y, info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlt_v
  end interface
  
  !
  !> Function  base_mlt_a
  !! \memberof  psb_s_base_vect_type
  !! \brief Vector entry-by-entry multiply  by a normal array y=x*y
  !! \param x(:) The array to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlt_a(x, y, info)
      real(psb_spk_), intent(in)        :: x(:)
      class(psb_s_base_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlt_a
  end interface

  !
  !> Function  base_mlt_a_2
  !! \memberof  psb_s_base_vect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by normal arrays
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x(:) The array to be multiplied b
  !! \param y(:) The array to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlt_a_2(alpha,x,y,beta,z,info)
      real(psb_spk_), intent(in)        :: alpha,beta
      real(psb_spk_), intent(in)        :: y(:)
      real(psb_spk_), intent(in)        :: x(:)
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlt_a_2
  end interface
  
  !
  !> Function  base_mlt_v_2
  !! \memberof  psb_s_base_vect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by class(base_vect)
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x  The class(base_vect) to be multiplied b
  !! \param y  The class(base_vect) to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
      real(psb_spk_), intent(in)        :: alpha,beta
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional     :: conjgx, conjgy
    end subroutine s_base_mlt_v_2
  end interface
  
  interface 
    module subroutine s_base_mlt_av(alpha,x,y,beta,z,info)
      real(psb_spk_), intent(in)        :: alpha,beta
      real(psb_spk_), intent(in)        :: x(:)
      class(psb_s_base_vect_type), intent(inout)  :: y
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlt_av
  end interface

  interface 
    module subroutine s_base_mlt_va(alpha,x,y,beta,z,info)
      real(psb_spk_), intent(in)        :: alpha,beta
      real(psb_spk_), intent(in)        :: y(:)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlt_va
  end interface
  
  !
  !> Function  base_div_v
  !! \memberof  psb_s_base_vect_type
  !! \brief Vector entry-by-entry divide  by a vector x=x/y
  !! \param y The array to be divided by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_div_v(x, y, info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_div_v
  end interface
  
  interface 
    module subroutine s_base_div_a(x, y, info)
      real(psb_spk_), intent(in)        :: x(:)
      class(psb_s_base_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_div_a
  end interface
  
  !
  !> Function  base_div_v2
  !! \memberof  psb_s_base_vect_type
  !! \brief Vector entry-by-entry divide  by a vector z=x/y
  !! \param y The array to be divided by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_div_v2(x, y, z, info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_div_v2
  end interface
  
  !
  !> Function  base_div_v_check
  !! \memberof  psb_s_base_vect_type
  !! \brief Vector entry-by-entry divide  by a vector x=x/y
  !! \param y The array to be divided by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_div_v_check(x, y, info, flag)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine s_base_div_v_check
  end interface
  
  !
  !> Function  base_div_v2_check
  !! \memberof  psb_s_base_vect_type
  !! \brief Vector entry-by-entry divide  by a vector z=x/y
  !! \param y The array to be divided by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_div_v2_check(x, y, z, info, flag)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine s_base_div_v2_check
  end interface
  
  !
  !> Function  base_div_a2
  !! \memberof  psb_s_base_vect_type
  !! \brief Entry-by-entry divide  between normal array z=x/y
  !! \param y(:) The array to be divided by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_div_a2(x, y, z, info)
      class(psb_s_base_vect_type), intent(inout)  :: z
      real(psb_spk_), intent(in) :: x(:)
      real(psb_spk_), intent(in)    :: y(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_div_a2
  end interface
  
  !
  !> Function  base_div_a2_check
  !! \memberof  psb_s_base_vect_type
  !! \brief Entry-by-entry divide  between normal array x=x/y and check if y(i)
  !! is different from zero
  !! \param y(:) The array to be dived by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_div_a2_check(x, y, z, info, flag)
      class(psb_s_base_vect_type), intent(inout)  :: z
      real(psb_spk_), intent(in) :: x(:)
      real(psb_spk_), intent(in)    :: y(:)
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in)              :: flag
    end subroutine s_base_div_a2_check
  end interface
  
  !
  !> Function  base_inv_v
  !! \memberof  psb_s_base_vect_type
  !! \brief Compute the entry-by-entry inverse of x and put it in y
  !! \param x The vector to be inverted
  !! \param y The vector containing the inverted vector
  !! \param info return code
  interface 
    module subroutine s_base_inv_v(x, y, info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_inv_v
  end interface
  
  !
  !> Function  base_inv_v_check
  !! \memberof  psb_s_base_vect_type
  !! \brief Compute the entry-by-entry inverse of x and put it in y, with 0 check
  !! \param x The vector to be inverted
  !! \param y The vector containing the inverted vector
  !! \param info return code
  !! \param flag if true does the check, otherwise call base_inv_v
  interface 
    module subroutine s_base_inv_v_check(x, y, info, flag)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_) :: i, n
      logical, intent(in) :: flag
    end subroutine s_base_inv_v_check
  end interface
  
  !
  !> Function  base_inv_a2
  !! \memberof  psb_s_base_vect_type
  !! \brief Compute the entry-by-entry inverse of x and put it in y,
  !! \param x(:) The array to be inverted
  !! \param y The vector containing the inverted vector
  !! \param info return code
  !
  interface 
    module subroutine s_base_inv_a2(x, y, info)
      class(psb_s_base_vect_type), intent(inout)  :: y
      real(psb_spk_), intent(in) :: x(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_inv_a2
  end interface
  
  !
  !> Function  base_inv_a2_check
  !! \memberof  psb_s_base_vect_type
  !! \brief Compute the entry-by-entry inverse of x and put it in y, with 0 check
  !! \param x(:) The array to be inverted
  !! \param y The vector containing the inverted vector
  !! \param info return code
  !! \param flag if true does the check, otherwise call base_inv_v
  !
  interface 
    module subroutine s_base_inv_a2_check(x, y, info, flag)
      class(psb_s_base_vect_type), intent(inout)  :: y
      real(psb_spk_), intent(inout) :: x(:)
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in)              :: flag
    end subroutine s_base_inv_a2_check
  end interface

  !
  !> Function  base_inv_a2_check
  !! \memberof  psb_s_base_vect_type
  !! \brief Compare entry-by-entry the vector x with the scalar c
  !! \param x The array to be compared
  !! \param z The vector containing in position i 1 if |x(i)| > c, 0 otherwise
  !! \param c The comparison term
  !! \param info return code
  !
  interface 
    module subroutine s_base_acmp_a2(x,c,z,info)
      real(psb_spk_), intent(in)             :: c
      real(psb_spk_), intent(inout)           :: x(:)
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine s_base_acmp_a2
  end interface
  
  !
  !> Function  base_cmp_v2
  !! \memberof  psb_s_base_vect_type
  !! \brief Compare entry-by-entry the vector x with the scalar c
  !! \param x The vector to be compared
  !! \param z The vector containing in position i 1 if |x(i)| > c, 0 otherwise
  !! \param c The comparison term
  !! \param info return code
  !
  interface 
    module subroutine s_base_acmp_v2(x,c,z,info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      real(psb_spk_), intent(in)                  :: c
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine s_base_acmp_v2
  end interface

  !
  ! Simple scaling
  !
  !> Function  base_scal
  !! \memberof  psb_s_base_vect_type
  !! \brief Scale all entries  x = alpha*x
  !! \param alpha   The multiplier
  !!
  interface 
    module subroutine s_base_scal(alpha, x)
      class(psb_s_base_vect_type), intent(inout)  :: x
      real(psb_spk_), intent (in)       :: alpha
    end subroutine s_base_scal
  end interface
  
  !
  ! Norms 1, 2 and infinity
  !
  !> Function  base_nrm2
  !! \memberof  psb_s_base_vect_type
  !! \brief 2-norm |x(1:n)|_2
  !! \param n  how many entries to consider
  interface 
    module function s_base_nrm2(n,x) result(res)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function s_base_nrm2
  end interface

  !
  !> Function  base_amax
  !! \memberof  psb_s_base_vect_type
  !! \brief infinity-norm |x(1:n)|_\infty
  !! \param n  how many entries to consider
  interface 
    module function s_base_amax(n,x) result(res)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function s_base_amax
  end interface

  !
  !> Function  base_min
  !! \memberof  psb_s_base_vect_type
  !! \brief min x(1:n)
  !! \param n  how many entries to consider
  interface 
    module function s_base_min(n,x) result(res)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function s_base_min
  end interface

  !
  !> Function  base_minquotient_v
  !! \memberof  psb_s_base_vect_type
  !! \brief Minimum entry of the vector entry-by-entry divide  x/y
  !! \param x The numerator vector
  !! \param y The denumerator vector
  !! \param info   return code
  !!
  interface 
    module function s_base_minquotient_v(x, y, info) result(z)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: y
      real(psb_spk_)                              :: z
      integer(psb_ipk_), intent(out)                :: info
    end function s_base_minquotient_v
  end interface

  !
  !> Function  base_minquotient_a2
  !! \memberof  psb_s_base_vect_type
  !! \brief Minimum entry of the array entry-by-entry divide  x/y
  !! \param x The numerator array
  !! \param y The denumerator array
  !! \param info   return code
  !!
  interface 
    module function s_base_minquotient_a2(x, y, info) result(z)
      class(psb_s_base_vect_type), intent(inout) :: x
      real(psb_spk_), intent(in)                 :: y(:)
      real(psb_spk_)                             :: z
      integer(psb_ipk_), intent(out)               :: info
    end function s_base_minquotient_a2
  end interface
  
  !
  !> Function  base_asum
  !! \memberof  psb_s_base_vect_type
  !! \brief 1-norm |x(1:n)|_1
  !! \param n  how many entries to consider
  interface 
    module function s_base_asum(n,x) result(res)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function s_base_asum
  end interface

  !
  ! Gather: Y = beta * Y + alpha * X(IDX(:))
  !
  !
  !> Function  base_gthab
  !! \memberof  psb_s_base_vect_type
  !! \brief gather into an array
  !!    Y = beta * Y + alpha * X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param alpha
  !! \param beta
  interface 
    module subroutine s_base_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_) :: alpha, beta, y(:)
      class(psb_s_base_vect_type) :: x
    end subroutine s_base_gthab
  end interface
  
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_gthzv
  !! \memberof  psb_s_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine s_base_gthzv_x(i,n,idx,x,y)
      integer(psb_ipk_) :: i
      integer(psb_mpk_) :: n
      class(psb_i_base_vect_type) :: idx
      real(psb_spk_) ::  y(:)
      class(psb_s_base_vect_type) :: x
    end subroutine s_base_gthzv_x
  end interface
  
  !
  ! New comm internals impl.
  !
  interface 
    module subroutine s_base_gthzbuf(i,n,idx,x)
      integer(psb_ipk_) :: i
      integer(psb_mpk_) :: n
      class(psb_i_base_vect_type) :: idx
      class(psb_s_base_vect_type) :: x
    end subroutine s_base_gthzbuf
  end interface
  
  !
  !> Function  base_device_wait:
  !! \memberof  psb_s_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_device_wait()
    end subroutine s_base_device_wait
  end interface

  interface 
    module function s_base_use_buffer() result(res)
      logical :: res
    end function s_base_use_buffer
  end interface
  
  interface 
    module subroutine s_base_new_buffer(n,x,info)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: n
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_new_buffer
  end interface

  interface 
    module subroutine s_base_new_comid(n,x,info)
      class(psb_s_base_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: n
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_new_comid
  end interface

  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_gthzv
  !! \memberof  psb_s_base_vect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine s_base_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_) ::  y(:)
      class(psb_s_base_vect_type) :: x
    end subroutine s_base_gthzv
  end interface

  !
  ! Scatter:
  ! Y(IDX(:)) = beta*Y(IDX(:)) + X(:)
  !
  !
  !> Function  base_sctb
  !! \memberof  psb_s_base_vect_type
  !! \brief scatter into a class(base_vect)
  !!    Y(IDX(:)) = beta * Y(IDX(:)) +  X(:)
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param beta
  !! \param x(:)
  interface 
    module subroutine s_base_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_) :: beta, x(:)
      class(psb_s_base_vect_type) :: y
    end subroutine s_base_sctb
  end interface
  
  interface 
    module subroutine s_base_sctb_x(i,n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      real(psb_spk_) :: beta, x(:)
      class(psb_s_base_vect_type) :: y
    end subroutine s_base_sctb_x
  end interface
  
  interface 
    module subroutine s_base_sctb_buf(i,n,idx,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      real(psb_spk_) :: beta
      class(psb_s_base_vect_type) :: y
    end subroutine s_base_sctb_buf
  end interface
  
  !
  !> Function  base_mask_a
  !! \memberof  psb_s_base_vect_type
  !! \brief Peform constraint tests looking at the value of c
  !! \param x The array to be compared
  !! \param c The array containing the information on the type of test to be
  !! performed, if c(i) = 2 ">0", if c(i) = 1 ">=0", if c(i) = 0 no test, if
  !! c(i) =-1 "<=0", if c(i) = -2 "< 0"
  !! \param m The vector containing the result of the comparison 1.0 for a
  !! failed test, and 0.0 for a passed one.
  !! \param t logical resulting from an and operation on all the tests
  !! \param info return code
  !
  interface 
    module subroutine s_base_mask_a(c,x,m,t,info)
      real(psb_spk_), intent(inout)               :: c(:)
      real(psb_spk_), intent(inout)               :: x(:)
      class(psb_s_base_vect_type), intent(inout)  :: m
      integer(psb_ipk_), intent(out)                :: info
      logical, intent(out)                          :: t
    end subroutine s_base_mask_a
  end interface
  
  !
  !> Function  base_mask_v
  !! \memberof  psb_s_base_vect_type
  !! \brief Peform constraint tests looking at the value of c
  !! \param x The vector to be compared
  !! \param c The vector containing the information on the type of test to be
  !! performed, if c(i) = 2 ">0", if c(i) = 1 ">=0", if c(i) = 0 no test, if
  !! c(i) =-1 "<=0", if c(i) = -2 "< 0"
  !! \param m The vector containing the result of the comparison 1.0 for a
  !! failed test, and 0.0 for a passed one.
  !! \param t logical resulting from an and operation on all the tests
  !! \param info return code
  !
  interface 
    module subroutine s_base_mask_v(c,x,m,t,info)
      class(psb_s_base_vect_type), intent(inout)  :: c
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_vect_type), intent(inout)  :: m
      integer(psb_ipk_), intent(out)                :: info
      logical, intent(out)                          :: t
    end subroutine s_base_mask_v
  end interface

  !
  !> Function  _base_addconst_a2
  !! \memberof  psb_s_base_vect_type
  !! \brief Add the constant b to every entry of the array x
  !! \param x The input array
  !! \param z The vector containing the x(i) + b
  !! \param b The added term
  !! \param info return code
  !
  interface 
    module subroutine s_base_addconst_a2(x,b,z,info)
      real(psb_spk_), intent(in)                  :: b
      real(psb_spk_), intent(inout)                :: x(:)
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)                :: info
    end subroutine s_base_addconst_a2
  end interface
  
  !
  !> Function  _base_addconst_v2
  !! \memberof  psb_s_base_vect_type
  !! \briefAdd the constant b to every entry of the vector x
  !! \param x The input vector
  !! \param z The vector containing the x(i) + b
  !! \param b The added term
  !! \param info return code
  !
  interface 
    module subroutine s_base_addconst_v2(x,b,z,info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      real(psb_spk_), intent(in)                  :: b
      class(psb_s_base_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)                :: info
    end subroutine s_base_addconst_v2
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
    real(psb_spk_)   :: x(:)
    type(psb_s_base_vect_type) :: this
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
    type(psb_s_base_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%asb(n,info)

  end function size_const

end module psb_s_base_vect_mod


module psb_s_base_multivect_mod

  use psb_error_mod
  use psb_realloc_mod
  use psb_s_base_vect_mod

  !> \namespace  psb_base_mod  \class psb_s_base_vect_type
  !! The psb_s_base_vect_type
  !! defines a middle level  integer(psb_ipk_) encapsulated dense vector.
  !! The encapsulation is needed, in place of a simple array, to allow
  !! for complicated situations, such as GPU programming, where the memory
  !!  area we are interested in is not easily accessible from the host/Fortran
  !!  side. It is also meant to be encapsulated in an outer type, to allow
  !!  runtime switching as per the STATE design pattern, similar to the
  !!  sparse matrix types.
  !!

  type psb_s_base_multivect_type
    !> Values.
    real(psb_spk_), allocatable :: v(:,:)
    real(psb_spk_), allocatable :: combuf(:)
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
    procedure, pass(x) :: bld_x    => s_base_mlv_bld_x
    procedure, pass(x) :: bld_n    => s_base_mlv_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: all      => s_base_mlv_all
    procedure, pass(x) :: mold     => s_base_mlv_mold
    !
    ! Insert/set. Assembly and free.
    ! Assembly does almost nothing here, but is important
    ! in derived classes.
    !
    procedure, pass(x) :: ins      => s_base_mlv_ins
    procedure, pass(x) :: zero     => s_base_mlv_zero
    procedure, pass(x) :: asb      => s_base_mlv_asb
    procedure, pass(x) :: free     => s_base_mlv_free
    procedure, pass(x) :: reinit      => s_base_mlv_reinit
    procedure, pass(x) :: set_ncfs => s_base_mlv_set_ncfs 
    procedure, pass(x) :: get_ncfs => s_base_mlv_get_ncfs
    procedure, pass(x) :: set_dupl => s_base_mlv_set_dupl 
    procedure, pass(x) :: get_dupl => s_base_mlv_get_dupl
    procedure, pass(x) :: set_state => s_base_mlv_set_state
    procedure, pass(x) :: set_null  => s_base_mlv_set_null
    procedure, pass(x) :: set_bld   => s_base_mlv_set_bld
    procedure, pass(x) :: set_upd   => s_base_mlv_set_upd
    procedure, pass(x) :: set_asb   => s_base_mlv_set_asb
    procedure, pass(x) :: get_state => s_base_mlv_get_state
    procedure, pass(x) :: is_null   => s_base_mlv_is_null
    procedure, pass(x) :: is_bld    => s_base_mlv_is_bld
    procedure, pass(x) :: is_upd    => s_base_mlv_is_upd
    procedure, pass(x) :: is_asb    => s_base_mlv_is_asb
    procedure, pass(x) :: base_cpy  => s_base_mlv_cpy
    !
    ! Sync: centerpiece of handling of external storage.
    ! Any derived class having extra storage upon sync
    ! will guarantee that both fortran/host side and
    ! external side contain the same data. The base
    ! version is only a placeholder.
    !
    procedure, pass(x) :: sync     => s_base_mlv_sync
    procedure, pass(x) :: is_host  => s_base_mlv_is_host
    procedure, pass(x) :: is_dev   => s_base_mlv_is_dev
    procedure, pass(x) :: is_sync  => s_base_mlv_is_sync
    procedure, pass(x) :: set_host => s_base_mlv_set_host
    procedure, pass(x) :: set_dev  => s_base_mlv_set_dev
    procedure, pass(x) :: set_sync => s_base_mlv_set_sync

    !
    ! Basic info
    procedure, pass(x) :: get_nrows => s_base_mlv_get_nrows
    procedure, pass(x) :: get_ncols => s_base_mlv_get_ncols
    procedure, pass(x) :: sizeof    => s_base_mlv_sizeof
    procedure, nopass  :: get_fmt   => s_base_mlv_get_fmt
    !
    ! Set/get data from/to an external array; also
    ! overload assignment.
    !
    procedure, pass(x) :: get_vect => s_base_mlv_get_vect
    procedure, pass(x) :: set_scal => s_base_mlv_set_scal
    procedure, pass(x) :: set_vect => s_base_mlv_set_vect
    generic, public    :: set      => set_vect, set_scal

    !
    ! Dot product and AXPBY
    !
    procedure, pass(x) :: dot_v    => s_base_mlv_dot_v
    procedure, pass(x) :: dot_a    => s_base_mlv_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => s_base_mlv_axpby_v
    procedure, pass(y) :: axpby_a  => s_base_mlv_axpby_a
    generic, public    :: axpby    => axpby_v, axpby_a
    !
    ! MultiVector by vector/multivector multiplication. Need all variants
    ! to handle multiple requirements from preconditioners
    !
    procedure, pass(y) :: mlt_mv   => s_base_mlv_mlt_mv
    procedure, pass(y) :: mlt_mv_v => s_base_mlv_mlt_mv_v
    procedure, pass(y) :: mlt_ar1  => s_base_mlv_mlt_ar1
    procedure, pass(y) :: mlt_ar2  => s_base_mlv_mlt_ar2
    procedure, pass(z) :: mlt_a_2  => s_base_mlv_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => s_base_mlv_mlt_v_2
!!$    procedure, pass(z) :: mlt_va   => s_base_mlv_mlt_va
!!$    procedure, pass(z) :: mlt_av   => s_base_mlv_mlt_av
    generic, public    :: mlt      => mlt_mv, mlt_mv_v, mlt_ar1, mlt_ar2, &
         & mlt_a_2, mlt_v_2 !, mlt_av, mlt_va
    !
    ! Scaling and norms
    !
    procedure, pass(x) :: scal     => s_base_mlv_scal
    procedure, pass(x) :: nrm2     => s_base_mlv_nrm2
    procedure, pass(x) :: amax     => s_base_mlv_amax
    procedure, pass(x) :: asum     => s_base_mlv_asum
    procedure, pass(x) :: absval1  => s_base_mlv_absval1
    procedure, pass(x) :: absval2  => s_base_mlv_absval2
    generic, public    :: absval   => absval1, absval2

    !
    ! These are for handling gather/scatter in new
    ! comm internals implementation.
    !
    procedure, nopass  :: use_buffer   => s_base_mlv_use_buffer
    procedure, pass(x) :: new_buffer   => s_base_mlv_new_buffer
    procedure, nopass  :: device_wait  => s_base_mlv_device_wait
    procedure, pass(x) :: maybe_free_buffer  => s_base_mlv_maybe_free_buffer
    procedure, pass(x) :: free_buffer  => s_base_mlv_free_buffer
    procedure, pass(x) :: new_comid    => s_base_mlv_new_comid
    procedure, pass(x) :: free_comid   => s_base_mlv_free_comid

    !
    ! Gather/scatter. These are needed for MPI interfacing.
    ! May have to be reworked.
    !
    procedure, pass(x) :: gthab    => s_base_mlv_gthab
    procedure, pass(x) :: gthzv    => s_base_mlv_gthzv
    procedure, pass(x) :: gthzm    => s_base_mlv_gthzm
    procedure, pass(x) :: gthzv_x  => s_base_mlv_gthzv_x
    procedure, pass(x) :: gthzbuf  => s_base_mlv_gthzbuf
    generic, public    :: gth      => gthab, gthzv, gthzm, gthzv_x, gthzbuf
    procedure, pass(y) :: sctb     => s_base_mlv_sctb
    procedure, pass(y) :: sctbr2   => s_base_mlv_sctbr2
    procedure, pass(y) :: sctb_x   => s_base_mlv_sctb_x
    procedure, pass(y) :: sctb_buf => s_base_mlv_sctb_buf
    generic, public    :: sct      => sctb, sctbr2, sctb_x, sctb_buf
  end type psb_s_base_multivect_type

  public  :: psb_s_base_multivect, psb_s_base_multivect_type

  interface psb_s_base_multivect
    module procedure constructor, size_const
  end interface psb_s_base_multivect

  private

  !
  ! Build from a sample
  !

  !> Function  bld_x:
  !! \memberof  psb_s_base_multivect_type
  !! \brief     Build method from an array
  !!  \param   x(:)  input array to be copied
  !!
  interface 
    module subroutine s_base_mlv_bld_x(x,this)
      real(psb_spk_), intent(in) :: this(:,:)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_) :: info
    end subroutine s_base_mlv_bld_x
  end interface
  
  !
  ! Create with size, but no initialization
  !

  !> Function  bld_n:
  !! \memberof  psb_s_base_multivect_type
  !! \brief     Build method with size (uninitialized data)
  !!  \param    n    size to be allocated.
  !!
  interface 
    module subroutine s_base_mlv_bld_n(x,m,n,scratch)
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_) :: info
      logical, intent(in), optional        :: scratch
    end subroutine s_base_mlv_bld_n
  end interface
  
  !> Function  base_mlv_all:
  !! \memberof  psb_s_base_multivect_type
  !! \brief     Build method with size (uninitialized data) and
  !!            allocation return code.
  !!  \param    n    size to be allocated.
  !!  \param info  return code
  !!
  interface 
    module subroutine s_base_mlv_all(m,n, x, info)
      integer(psb_ipk_), intent(in)               :: m,n
      class(psb_s_base_multivect_type), intent(out) :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_all
  end interface

  !> Function  base_mlv_mold:
  !! \memberof  psb_s_base_multivect_type
  !! \brief     Mold method: return a variable with the same dynamic type
  !!  \param    y returned variable
  !!  \param info  return code
  !!
  interface 
    module subroutine s_base_mlv_mold(x, y, info)
      class(psb_s_base_multivect_type), intent(in)   :: x
      class(psb_s_base_multivect_type), intent(out), allocatable :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_mold
  end interface

  interface 
    module subroutine s_base_mlv_reinit(x, info)
      class(psb_s_base_multivect_type), intent(out)    :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_reinit
  end interface

  !
  ! Insert a bunch of values at specified positions.
  !
  !> Function  base_mlv_ins:
  !! \memberof  psb_s_base_multivect_type
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
    module subroutine s_base_mlv_ins(n,irl,val,dupl,x,maxr,info)
      class(psb_s_base_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, dupl,maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      real(psb_spk_), intent(in)        :: val(:,:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_ins
  end interface

  !
  !> Function  base_mlv_zero
  !! \memberof  psb_s_base_multivect_type
  !! \brief Zero out contents
  !!
  !
  interface 
    module subroutine s_base_mlv_zero(x)
      class(psb_s_base_multivect_type), intent(inout)    :: x
    end subroutine s_base_mlv_zero
  end interface

  !
  ! Assembly.
  ! For derived classes: after this the vector
  ! storage is supposed to be in sync.
  !
  !> Function  base_mlv_asb:
  !! \memberof  psb_s_base_multivect_type
  !! \brief Assemble vector: reallocate as necessary.
  !!
  !!  \param n     final size
  !!  \param info  return code
  !!
  !
  
  interface 
    module subroutine s_base_mlv_asb(m,n, x, info, scratch)
      integer(psb_ipk_), intent(in)              :: m,n
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
      logical, intent(in), optional        :: scratch
    end subroutine s_base_mlv_asb
  end interface

  !
  !> Function  base_mlv_free:
  !! \memberof  psb_s_base_multivect_type
  !! \brief Free vector
  !!
  !!  \param info  return code
  !!
  !
  interface 
    module subroutine s_base_mlv_free(x, info)
      class(psb_s_base_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_free
  end interface

  interface 
    module function s_base_mlv_get_ncfs(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_mlv_get_ncfs
  end interface

  interface 
    module function s_base_mlv_get_dupl(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_mlv_get_dupl
  end interface

  interface 
    module function s_base_mlv_get_state(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_mlv_get_state
  end interface

  interface 
    module function s_base_mlv_is_null(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      logical :: res
    end function s_base_mlv_is_null
  end interface

  interface 
    module function s_base_mlv_is_bld(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      logical :: res
    end function s_base_mlv_is_bld
  end interface

  interface 
    module function s_base_mlv_is_upd(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      logical :: res
    end function s_base_mlv_is_upd
  end interface

  interface 
    module function s_base_mlv_is_asb(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      logical :: res
    end function s_base_mlv_is_asb
  end interface

  interface 
    module subroutine  s_base_mlv_set_ncfs(n,x)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine s_base_mlv_set_ncfs
  end interface

  interface 
    module subroutine  s_base_mlv_set_dupl(n,x)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine s_base_mlv_set_dupl
  end interface

  interface 
    module subroutine  s_base_mlv_set_state(n,x)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine s_base_mlv_set_state
  end interface

  interface 
    module subroutine  s_base_mlv_set_null(x)
      class(psb_s_base_multivect_type), intent(inout) :: x
    end subroutine s_base_mlv_set_null
  end interface

  interface 
    module subroutine  s_base_mlv_set_bld(x)
      class(psb_s_base_multivect_type), intent(inout) :: x
    end subroutine s_base_mlv_set_bld
  end interface

  interface 
    module subroutine  s_base_mlv_set_upd(x)
      class(psb_s_base_multivect_type), intent(inout) :: x
    end subroutine s_base_mlv_set_upd
  end interface

  interface 
    module subroutine  s_base_mlv_set_asb(x)
      class(psb_s_base_multivect_type), intent(inout) :: x
    end subroutine s_base_mlv_set_asb
  end interface

  !
  ! The base version of SYNC & friends does nothing, it's just
  ! a placeholder.
  !
  !
  !> Function  base_mlv_sync:
  !! \memberof  psb_s_base_multivect_type
  !! \brief Sync: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_mlv_sync(x)
      class(psb_s_base_multivect_type), intent(inout) :: x
    end subroutine s_base_mlv_sync
  end interface

  !
  !> Function  base_mlv_set_host:
  !! \memberof  psb_s_base_multivect_type
  !! \brief Set_host: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_mlv_set_host(x)
      class(psb_s_base_multivect_type), intent(inout) :: x
    end subroutine s_base_mlv_set_host
  end interface

  !
  !> Function  base_mlv_set_dev:
  !! \memberof  psb_s_base_multivect_type
  !! \brief Set_dev: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_mlv_set_dev(x)
      class(psb_s_base_multivect_type), intent(inout) :: x
    end subroutine s_base_mlv_set_dev
  end interface

  !
  !> Function  base_mlv_set_sync:
  !! \memberof  psb_s_base_multivect_type
  !! \brief Set_sync: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_mlv_set_sync(x)
      class(psb_s_base_multivect_type), intent(inout) :: x
    end subroutine s_base_mlv_set_sync
  end interface

  !
  !> Function  base_mlv_is_dev:
  !! \memberof  psb_s_base_multivect_type
  !! \brief Is  vector on external device    .
  !!
  !
  interface 
    module function s_base_mlv_is_dev(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      logical  :: res
    end function s_base_mlv_is_dev
  end interface

  !
  !> Function  base_mlv_is_host
  !! \memberof  psb_s_base_multivect_type
  !! \brief Is  vector on standard memory    .
  !!
  !
  interface 
    module function s_base_mlv_is_host(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      logical  :: res
    end function s_base_mlv_is_host
  end interface

  !
  !> Function  base_mlv_is_sync
  !! \memberof  psb_s_base_multivect_type
  !! \brief Is  vector on sync               .
  !!
  !
  interface 
    module function s_base_mlv_is_sync(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      logical  :: res
    end function s_base_mlv_is_sync
  end interface

  !> Function  base_cpy:
  !! \memberof  psb_d_base_vect_type
  !! \brief     base_cpy: copy base contents
  !!  \param    y returned variable
  !!
  interface 
    module subroutine s_base_mlv_cpy(x, y)
      class(psb_s_base_multivect_type), intent(in)   :: x
      class(psb_s_base_multivect_type), intent(out)  :: y
    end subroutine s_base_mlv_cpy
  end interface

  ! Size info.
  !
  !
  !> Function  base_mlv_get_nrows
  !! \memberof  psb_s_base_multivect_type
  !! \brief  Number of entries
  !!
  !
  interface 
    module function s_base_mlv_get_nrows(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_mlv_get_nrows
  end interface

  interface 
    module function s_base_mlv_get_ncols(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function s_base_mlv_get_ncols
  end interface

  !
  !> Function  base_mlv_get_sizeof
  !! \memberof  psb_s_base_multivect_type
  !! \brief  Size in bytesa
  !!
  !
  interface 
    module function s_base_mlv_sizeof(x) result(res)
      class(psb_s_base_multivect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function s_base_mlv_sizeof
  end interface

  !
  !> Function  base_mlv_get_fmt
  !! \memberof  psb_s_base_multivect_type
  !! \brief  Format
  !!
  !
  interface 
    module function s_base_mlv_get_fmt() result(res)
      character(len=5) :: res
    end function s_base_mlv_get_fmt
  end interface

  !
  !
  !
  !> Function  base_mlv_get_vect
  !! \memberof  psb_s_base_multivect_type
  !! \brief  Extract a copy of the contents
  !!
  !
  interface 
    module function  s_base_mlv_get_vect(x) result(res)
      class(psb_s_base_multivect_type), intent(inout) :: x
      real(psb_spk_), allocatable                 :: res(:,:)
    end function s_base_mlv_get_vect
  end interface

  !
  ! Reset all values
  !
  !
  !> Function  base_mlv_set_scal
  !! \memberof  psb_s_base_multivect_type
  !! \brief  Set all entries
  !! \param val   The value to set
  !!
  interface 
    module subroutine s_base_mlv_set_scal(x,val)
      class(psb_s_base_multivect_type), intent(inout)  :: x
      real(psb_spk_), intent(in) :: val
    end subroutine s_base_mlv_set_scal
  end interface
  
  !
  !> Function  base_mlv_set_vect
  !! \memberof  psb_s_base_multivect_type
  !! \brief  Set all entries
  !! \param val(:)  The vector to be copied in
  !!
  interface 
    module subroutine s_base_mlv_set_vect(x,val)
      class(psb_s_base_multivect_type), intent(inout)  :: x
      real(psb_spk_), intent(in) :: val(:,:)
    end subroutine s_base_mlv_set_vect
  end interface

  !
  ! Dot products
  !
  !
  !> Function  base_mlv_dot_v
  !! \memberof  psb_s_base_multivect_type
  !! \brief  Dot product by another base_mlv_vector
  !! \param n    Number of entries to be considered
  !! \param y    The other (base_mlv_vect) to be multiplied by
  !!
  interface 
    module function s_base_mlv_dot_v(n,x,y) result(res)
      class(psb_s_base_multivect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_), allocatable   :: res(:)
    end function s_base_mlv_dot_v
  end interface

  !
  ! Base workhorse is good old BLAS1
  !
  !
  !> Function  base_mlv_dot_a
  !! \memberof  psb_s_base_multivect_type
  !! \brief  Dot product by a normal array
  !! \param n    Number of entries to be considered
  !! \param y(:) The array to be multiplied by
  !!
  interface 
    module function s_base_mlv_dot_a(n,x,y) result(res)
      class(psb_s_base_multivect_type), intent(inout) :: x
      real(psb_spk_), intent(in)    :: y(:,:)
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_), allocatable     :: res(:)
    end function s_base_mlv_dot_a
  end interface

  !
  ! AXPBY is invoked via Y, hence the structure below.
  !
  !
  !
  !> Function  base_mlv_axpby_v
  !! \memberof  psb_s_base_multivect_type
  !! \brief AXPBY  by a (base_mlv_vect) y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x     The class(base_mlv_vect) to be added
  !! \param beta scalar alpha
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlv_axpby_v(m,alpha, x, beta, y, info, n)
      integer(psb_ipk_), intent(in)               :: m
      class(psb_s_base_multivect_type), intent(inout)  :: x
      class(psb_s_base_multivect_type), intent(inout)  :: y
      real(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), intent(in), optional     :: n
    end subroutine s_base_mlv_axpby_v
  end interface

  !
  ! AXPBY is invoked via Y, hence the structure below.
  !
  !
  !> Function  base_mlv_axpby_a
  !! \memberof  psb_s_base_multivect_type
  !! \brief AXPBY  by a normal array y=alpha*x+beta*y
  !! \param m    Number of entries to be considered
  !! \param alpha scalar alpha
  !! \param x(:) The array to be added
  !! \param beta scalar alpha
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlv_axpby_a(m,alpha, x, beta, y, info,n)
      integer(psb_ipk_), intent(in)               :: m
      real(psb_spk_), intent(in)        :: x(:,:)
      class(psb_s_base_multivect_type), intent(inout)  :: y
      real(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
      integer(psb_ipk_), intent(in), optional     :: n
    end subroutine s_base_mlv_axpby_a
  end interface

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
  !! \memberof  psb_s_base_multivect_type
  !! \brief Multivector entry-by-entry multiply  by a base_mlv_multivect  y=x*y
  !! \param x   The class(base_mlv_vect) to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlv_mlt_mv(x, y, info)
      class(psb_s_base_multivect_type), intent(inout)  :: x
      class(psb_s_base_multivect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_mlt_mv
  end interface

  interface 
    module subroutine s_base_mlv_mlt_mv_v(x, y, info)
      class(psb_s_base_vect_type), intent(inout)  :: x
      class(psb_s_base_multivect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_mlt_mv_v
  end interface

  !
  !> Function  base_mlv_mlt_ar1
  !! \memberof  psb_s_base_multivect_type
  !! \brief MultiVector entry-by-entry multiply  by a rank 1 array y=x*y
  !! \param x(:) The array to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlv_mlt_ar1(x, y, info)
      real(psb_spk_), intent(in)        :: x(:)
      class(psb_s_base_multivect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_mlt_ar1
  end interface
  
  !
  !> Function  base_mlv_mlt_ar2
  !! \memberof  psb_s_base_multivect_type
  !! \brief MultiVector entry-by-entry multiply  by a rank 2 array y=x*y
  !! \param x(:,:) The array to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlv_mlt_ar2(x, y, info)
      real(psb_spk_), intent(in)        :: x(:,:)
      class(psb_s_base_multivect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_mlt_ar2
  end interface

  !
  !> Function  base_mlv_mlt_a_2
  !! \memberof  psb_s_base_multivect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by normal arrays
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x(:) The array to be multiplied b
  !! \param y(:) The array to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlv_mlt_a_2(alpha,x,y,beta,z,info)
      real(psb_spk_), intent(in)        :: alpha,beta
      real(psb_spk_), intent(in)        :: y(:,:)
      real(psb_spk_), intent(in)        :: x(:,:)
      class(psb_s_base_multivect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine s_base_mlv_mlt_a_2
  end interface
  
  !
  !> Function  base_mlv_mlt_v_2
  !! \memberof  psb_s_base_multivect_type
  !! \brief AXPBY-like Vector entry-by-entry multiply  by class(base_mlv_vect)
  !!        z=beta*z+alpha*x*y
  !! \param alpha
  !! \param beta
  !! \param x  The class(base_mlv_vect) to be multiplied b
  !! \param y  The class(base_mlv_vect) to be multiplied by
  !! \param info   return code
  !!
  interface 
    module subroutine s_base_mlv_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
      real(psb_spk_), intent(in)        :: alpha,beta
      class(psb_s_base_multivect_type), intent(inout)  :: x
      class(psb_s_base_multivect_type), intent(inout)  :: y
      class(psb_s_base_multivect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), intent(in), optional     :: conjgx, conjgy
    end subroutine s_base_mlv_mlt_v_2
  end interface
!!$
!!$  subroutine s_base_mlv_mlt_av(alpha,x,y,beta,z,info)
!!$    real(psb_spk_), intent(in)        :: alpha,beta
!!$    real(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_s_base_multivect_type), intent(inout)  :: y
!!$    class(psb_s_base_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    call z%mlt(alpha,x,y%v,beta,info)
!!$
!!$  end subroutine s_base_mlv_mlt_av
!!$
!!$  subroutine s_base_mlv_mlt_va(alpha,x,y,beta,z,info)
!!$    real(psb_spk_), intent(in)        :: alpha,beta
!!$    real(psb_spk_), intent(in)        :: y(:)
!!$    class(psb_s_base_multivect_type), intent(inout)  :: x
!!$    class(psb_s_base_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    call z%mlt(alpha,y,x,beta,info)
!!$
!!$  end subroutine s_base_mlv_mlt_va
!!$
!!$
  !
  ! Simple scaling
  !
  !> Function  base_mlv_scal
  !! \memberof  psb_s_base_multivect_type
  !! \brief Scale all entries  x = alpha*x
  !! \param alpha   The multiplier
  !!
  interface 
    module subroutine s_base_mlv_scal(alpha, x)
      class(psb_s_base_multivect_type), intent(inout)  :: x
      real(psb_spk_), intent (in)       :: alpha
    end subroutine s_base_mlv_scal
  end interface
  
  !
  ! Norms 1, 2 and infinity
  !
  !> Function  base_mlv_nrm2
  !! \memberof  psb_s_base_multivect_type
  !! \brief 2-norm |x(1:n)|_2
  !! \param n  how many entries to consider
  interface 
    module function s_base_mlv_nrm2(n,x) result(res)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_), allocatable    :: res(:)
    end function s_base_mlv_nrm2
  end interface

  !
  !> Function  base_mlv_amax
  !! \memberof  psb_s_base_multivect_type
  !! \brief infinity-norm |x(1:n)|_\infty
  !! \param n  how many entries to consider
  interface 
    module function s_base_mlv_amax(n,x) result(res)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_), allocatable    :: res(:)
    end function s_base_mlv_amax
  end interface
  
  !
  !> Function  base_mlv_asum
  !! \memberof  psb_s_base_multivect_type
  !! \brief 1-norm |x(1:n)|_1
  !! \param n  how many entries to consider
  interface 
    module function s_base_mlv_asum(n,x) result(res)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_), allocatable    :: res(:)
    end function s_base_mlv_asum
  end interface
  
  !
  ! Overwrite with absolute value
  !
  !
  !> Function  base_absval1
  !! \memberof  psb_s_base_vect_type
  !! \brief  Set all entries to their respective absolute values.
  !!
  interface 
    module subroutine s_base_mlv_absval1(x)
      class(psb_s_base_multivect_type), intent(inout)  :: x
    end subroutine s_base_mlv_absval1
  end interface
  
  interface 
    module subroutine s_base_mlv_absval2(x,y)
      class(psb_s_base_multivect_type), intent(inout) :: x
      class(psb_s_base_multivect_type), intent(inout) :: y
      integer(psb_ipk_) :: info
    end subroutine s_base_mlv_absval2
  end interface


  interface 
    module function s_base_mlv_use_buffer() result(res)
      logical :: res
    end function s_base_mlv_use_buffer
  end interface
  
  interface 
    module subroutine s_base_mlv_new_buffer(n,x,info)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: n
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_mlv_new_buffer
  end interface
  
  interface 
    module subroutine s_base_mlv_new_comid(n,x,info)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)              :: n
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_mlv_new_comid
  end interface
  
  interface 
    module subroutine s_base_mlv_maybe_free_buffer(x,info)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_mlv_maybe_free_buffer
  end interface
  
  interface 
    module subroutine s_base_mlv_free_buffer(x,info)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_mlv_free_buffer
  end interface

  interface 
    module subroutine s_base_mlv_free_comid(x,info)
      class(psb_s_base_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine s_base_mlv_free_comid
  end interface

  !
  ! Gather: Y = beta * Y + alpha * X(IDX(:))
  !
  !
  !> Function  base_mlv_gthab
  !! \memberof  psb_s_base_multivect_type
  !! \brief gather into an array
  !!    Y = beta * Y + alpha * X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param alpha
  !! \param beta
  interface 
    module subroutine s_base_mlv_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_) :: alpha, beta, y(:)
      class(psb_s_base_multivect_type) :: x
    end subroutine s_base_mlv_gthab
  end interface
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_s_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine s_base_mlv_gthzv_x(i,n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      real(psb_spk_) ::  y(:)
      class(psb_s_base_multivect_type) :: x
    end subroutine s_base_mlv_gthzv_x
  end interface

  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_s_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine s_base_mlv_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_) ::  y(:)
      class(psb_s_base_multivect_type) :: x
    end subroutine s_base_mlv_gthzv
  end interface
  !
  ! shortcut alpha=1 beta=0
  !
  !> Function  base_mlv_gthzv
  !! \memberof  psb_s_base_multivect_type
  !! \brief gather into an array special alpha=1 beta=0
  !!    Y = X(IDX(:))
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  interface 
    module subroutine s_base_mlv_gthzm(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_) ::  y(:,:)
      class(psb_s_base_multivect_type) :: x
    end subroutine s_base_mlv_gthzm
  end interface

  !
  ! New comm internals impl.
  !
  interface 
    module subroutine s_base_mlv_gthzbuf(i,ixb,n,idx,x)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i, ixb
      class(psb_i_base_vect_type) :: idx
      class(psb_s_base_multivect_type) :: x
    end subroutine s_base_mlv_gthzbuf
  end interface
  
  !
  ! Scatter:
  ! Y(IDX(:),:) = beta*Y(IDX(:),:) + X(:)
  !
  !
  !> Function  base_mlv_sctb
  !! \memberof  psb_s_base_multivect_type
  !! \brief scatter into a class(base_mlv_vect)
  !!    Y(IDX(:)) = beta * Y(IDX(:)) +  X(:)
  !! \param n  how many entries to consider
  !! \param idx(:) indices
  !! \param beta
  !! \param x(:)
  interface 
    module subroutine s_base_mlv_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_) :: beta, x(:)
      class(psb_s_base_multivect_type) :: y
    end subroutine s_base_mlv_sctb
  end interface

  interface 
    module subroutine s_base_mlv_sctbr2(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_spk_) :: beta, x(:,:)
      class(psb_s_base_multivect_type) :: y
    end subroutine s_base_mlv_sctbr2
  end interface
  
  interface 
    module subroutine s_base_mlv_sctb_x(i,n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      real( psb_spk_) :: beta, x(:)
      class(psb_s_base_multivect_type) :: y
    end subroutine s_base_mlv_sctb_x
  end interface

  interface 
    module subroutine s_base_mlv_sctb_buf(i,iyb,n,idx,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i, iyb
      class(psb_i_base_vect_type) :: idx
      real(psb_spk_) :: beta
      class(psb_s_base_multivect_type) :: y
    end subroutine s_base_mlv_sctb_buf
  end interface

  !
  !> Function  base_device_wait:
  !! \memberof  psb_s_base_vect_type
  !! \brief device_wait: base version is a no-op.
  !!
  !
  interface 
    module subroutine s_base_mlv_device_wait()
    end subroutine s_base_mlv_device_wait
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
    real(psb_spk_)   :: x(:,:)
    type(psb_s_base_multivect_type) :: this
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
    type(psb_s_base_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%asb(m,n,info)

  end function size_const

end module psb_s_base_multivect_mod
