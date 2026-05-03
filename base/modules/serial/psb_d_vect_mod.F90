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
! package: psb_d_vect_mod
!
! This module contains the definition of the psb_d_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component.
!
module psb_d_vect_mod

  use psb_realloc_mod
  use psb_d_base_vect_mod
  use psb_i_vect_mod

  type psb_d_vect_type
    class(psb_d_base_vect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build=psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    real(psb_dpk_), allocatable :: rmtv(:)
    integer(psb_lpk_), allocatable :: rmidx(:)
  contains
    procedure, pass(x) :: get_nrows => d_vect_get_nrows
    procedure, pass(x) :: sizeof   => d_vect_sizeof
    procedure, pass(x) :: get_fmt  => d_vect_get_fmt
    procedure, pass(x) :: is_remote_build => d_vect_is_remote_build
    procedure, pass(x) :: set_remote_build => d_vect_set_remote_build
    procedure, pass(x) :: get_nrmv => d_vect_get_nrmv
    procedure, pass(x) :: set_nrmv => d_vect_set_nrmv
    procedure, pass(x) :: all      => d_vect_all
    procedure, pass(x) :: reall    => d_vect_reall
    procedure, pass(x) :: zero     => d_vect_zero
    procedure, pass(x) :: asb      => d_vect_asb
    procedure, pass(x) :: set_dupl => d_vect_set_dupl 
    procedure, pass(x) :: get_dupl => d_vect_get_dupl
    procedure, pass(x) :: set_ncfs => d_vect_set_ncfs 
    procedure, pass(x) :: get_ncfs => d_vect_get_ncfs
    procedure, pass(x) :: set_state => d_vect_set_state
    procedure, pass(x) :: set_null  => d_vect_set_null
    procedure, pass(x) :: set_bld   => d_vect_set_bld
    procedure, pass(x) :: set_upd   => d_vect_set_upd
    procedure, pass(x) :: set_asb   => d_vect_set_asb
    procedure, pass(x) :: get_state => d_vect_get_state
    procedure, pass(x) :: is_null   => d_vect_is_null
    procedure, pass(x) :: is_bld    => d_vect_is_bld
    procedure, pass(x) :: is_upd    => d_vect_is_upd
    procedure, pass(x) :: is_asb    => d_vect_is_asb
    procedure, pass(x) :: reinit    => d_vect_reinit

    procedure, pass(x) :: gthab    => d_vect_gthab
    procedure, pass(x) :: gthzv    => d_vect_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => d_vect_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => d_vect_free
    procedure, pass(x) :: ins_a    => d_vect_ins_a
    procedure, pass(x) :: ins_v    => d_vect_ins_v
    generic, public    :: ins      => ins_v, ins_a
    procedure, pass(x) :: bld_x    => d_vect_bld_x
    procedure, pass(x) :: bld_mn   => d_vect_bld_mn
    procedure, pass(x) :: bld_en   => d_vect_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
    procedure, pass(x) :: get_vect => d_vect_get_vect
    procedure, pass(x) :: cnv      => d_vect_cnv
    procedure, pass(x) :: set_scal => d_vect_set_scal
    procedure, pass(x) :: set_vect => d_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => d_vect_clone

    procedure, pass(x) :: sync     => d_vect_sync
    procedure, pass(x) :: is_host  => d_vect_is_host
    procedure, pass(x) :: is_dev   => d_vect_is_dev
    procedure, pass(x) :: is_sync  => d_vect_is_sync
    procedure, pass(x) :: set_host => d_vect_set_host
    procedure, pass(x) :: set_dev  => d_vect_set_dev
    procedure, pass(x) :: set_sync => d_vect_set_sync
    procedure, pass(x) :: check_addr => d_vect_check_addr

    procedure, pass(x) :: get_entry => d_vect_get_entry
    procedure, pass(x) :: set_entry => d_vect_set_entry

    procedure, pass(x) :: dot_v    => d_vect_dot_v
    procedure, pass(x) :: dot_a    => d_vect_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => d_vect_axpby_v
    procedure, pass(y) :: axpby_a  => d_vect_axpby_a
    procedure, pass(z) :: axpby_v2  => d_vect_axpby_v2
    procedure, pass(z) :: axpby_a2  => d_vect_axpby_a2
    generic, public    :: axpby    => axpby_v, axpby_a, axpby_v2, axpby_a2
    procedure, pass(z) :: upd_xyz  => d_vect_upd_xyz
    procedure, pass(z) :: xyzw     => d_vect_xyzw
    procedure, pass(y) :: mlt_v    => d_vect_mlt_v
    procedure, pass(y) :: mlt_a    => d_vect_mlt_a
    procedure, pass(z) :: mlt_a_2  => d_vect_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => d_vect_mlt_v_2
    procedure, pass(z) :: mlt_va   => d_vect_mlt_va
    procedure, pass(z) :: mlt_av   => d_vect_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
         & mlt_v_2, mlt_av, mlt_va
    procedure, pass(x) :: div_v    => d_vect_div_v
    procedure, pass(z) :: div_v2    => d_vect_div_v2
    procedure, pass(x) :: div_v_check => d_vect_div_v_check
    procedure, pass(x) :: div_v2_check => d_vect_div_v2_check
    procedure, pass(z) :: div_a2   => d_vect_div_a2
    procedure, pass(z) :: div_a2_check => d_vect_div_a2_check
    generic, public    :: div      => div_v, div_v2, div_v_check, &
                                        div_v2_check, div_a2, div_a2_check
    procedure, pass(y) :: inv_v    => d_vect_inv_v
    procedure, pass(y) :: inv_v_check => d_vect_inv_v_check
    procedure, pass(y) :: inv_a2   => d_vect_inv_a2
    procedure, pass(y) :: inv_a2_check => d_vect_inv_a2_check
    generic, public    :: inv      => inv_v, inv_v_check, inv_a2, inv_a2_check
    procedure, pass(x) :: scal     => d_vect_scal
    procedure, pass(x) :: absval1  => d_vect_absval1
    procedure, pass(x) :: absval2  => d_vect_absval2
    generic, public    :: absval   => absval1, absval2
    procedure, pass(x) :: nrm2std  => d_vect_nrm2
    procedure, pass(x) :: nrm2weight => d_vect_nrm2_weight
    procedure, pass(x) :: nrm2weightmask => d_vect_nrm2_weight_mask
    generic, public    :: nrm2     => nrm2std, nrm2weight, nrm2weightmask
    procedure, pass(x) :: amax     => d_vect_amax
    procedure, pass(x) :: asum     => d_vect_asum
    procedure, pass(z) :: acmp_a2   => d_vect_acmp_a2
    procedure, pass(z) :: acmp_v2   => d_vect_acmp_v2
    generic, public    :: acmp      => acmp_a2, acmp_v2
    procedure, pass(z) :: addconst_a2   => d_vect_addconst_a2
    procedure, pass(z) :: addconst_v2   => d_vect_addconst_v2
    generic, public    :: addconst      => addconst_a2, addconst_v2

    procedure, pass(x)   :: minreal   => d_vect_min
    procedure, pass(m) :: mask_v => d_vect_mask_v
    procedure, pass(m) :: mask_a => d_vect_mask_a
    generic, public    :: mask => mask_a, mask_v
    procedure, pass(x) :: minquotient_v  => d_vect_minquotient_v
    procedure, pass(x) :: minquotient_a2 => d_vect_minquotient_a2
    generic, public    :: minquotient    => minquotient_v, minquotient_a2

  end type psb_d_vect_type

  public  :: psb_d_vect
  private :: constructor, size_const
  interface psb_d_vect
    module procedure constructor, size_const
  end interface psb_d_vect

  private :: d_vect_get_nrows, d_vect_sizeof, d_vect_get_fmt, &
       & d_vect_all, d_vect_reall, d_vect_zero,  d_vect_asb, &
       & d_vect_gthab, d_vect_gthzv, d_vect_sctb, &
       & d_vect_free, d_vect_ins_a, d_vect_ins_v, d_vect_bld_x, &
       & d_vect_bld_mn, d_vect_bld_en, d_vect_get_vect, &
       & d_vect_cnv, d_vect_set_scal, &
       & d_vect_set_vect, d_vect_clone, d_vect_sync, d_vect_is_host, &
       & d_vect_is_dev, d_vect_is_sync, d_vect_set_host, &
       & d_vect_set_dev, d_vect_set_sync, &
       & d_vect_set_remote_build, d_is_remote_build, &
       & d_vect_set_dupl, d_get_dupl, d_vect_set_nrmv, d_get_nrmv

  private ::  d_vect_dot_v, d_vect_dot_a, d_vect_axpby_v, d_vect_axpby_a, &
       & d_vect_mlt_v, d_vect_mlt_a, d_vect_mlt_a_2, d_vect_mlt_v_2, &
       & d_vect_mlt_va, d_vect_mlt_av, d_vect_scal, d_vect_absval1, &
       & d_vect_absval2, d_vect_nrm2, d_vect_amax, d_vect_asum


  class(psb_d_base_vect_type), allocatable, target,&
       & save, private :: psb_d_base_vect_default

  
  interface 
    module function d_vect_get_dupl(x) result(res)
      implicit none
      class(psb_d_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function d_vect_get_dupl
  end interface

  interface 
    module subroutine d_vect_set_dupl(x,val)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine d_vect_set_dupl
  end interface

  interface 
    module function d_vect_get_ncfs(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function d_vect_get_ncfs
  end interface

  interface 
    module subroutine d_vect_set_ncfs(x,val)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine d_vect_set_ncfs
  end interface

  interface 
    module function d_vect_get_state(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function d_vect_get_state
  end interface

  interface 
    module function d_vect_is_null(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      logical :: res
    end function d_vect_is_null
  end interface

  interface 
    module function d_vect_is_bld(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      logical :: res
    end function d_vect_is_bld
  end interface

  interface 
    module function d_vect_is_upd(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      logical :: res
    end function d_vect_is_upd
  end interface

  interface 
    module function d_vect_is_asb(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      logical :: res
    end function d_vect_is_asb
  end interface

  interface 
    module subroutine  d_vect_set_state(n,x)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine d_vect_set_state
  end interface

  interface 
    module subroutine  d_vect_set_null(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_set_null
  end interface

  interface 
    module subroutine  d_vect_set_bld(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_set_bld
  end interface

  interface 
    module subroutine  d_vect_set_upd(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_set_upd
  end interface

  interface 
    module subroutine  d_vect_set_asb(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_set_asb
  end interface

  interface 
    module function d_vect_get_nrmv(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function d_vect_get_nrmv
  end interface

  interface 
    module subroutine d_vect_set_nrmv(x,val)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: val
    end subroutine d_vect_set_nrmv
  end interface

  interface 
    module function d_vect_is_remote_build(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      logical :: res
    end function d_vect_is_remote_build
  end interface
  
  interface 
    module subroutine d_vect_set_remote_build(x,val)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine d_vect_set_remote_build
  end interface
        
  interface psb_set_vect_default
    module subroutine  psb_d_set_vect_default(v)
      class(psb_d_base_vect_type), intent(in) :: v
    end subroutine psb_d_set_vect_default
  end interface

  interface psb_get_vect_default
    module function psb_d_get_vect_default(v) result(res)
      class(psb_d_vect_type), intent(in) :: v
      class(psb_d_base_vect_type), pointer :: res
    end function psb_d_get_vect_default
  end interface

  interface 
    module subroutine  psb_d_clear_vect_default()
    end subroutine psb_d_clear_vect_default
  end interface

  interface 
    module function psb_d_get_base_vect_default() result(res)
      class(psb_d_base_vect_type), pointer :: res
    end function psb_d_get_base_vect_default
  end interface

  interface 
    module subroutine d_vect_clone(x,y,info)
      class(psb_d_vect_type), intent(inout) :: x
      class(psb_d_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)        :: info
    end subroutine d_vect_clone
  end interface

  interface 
    module subroutine d_vect_bld_x(x,invect,mold,scratch)
      real(psb_dpk_), intent(in)          :: invect(:)
      class(psb_d_vect_type), intent(inout) :: x
      class(psb_d_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine d_vect_bld_x
  end interface

  interface 
    module subroutine d_vect_bld_mn(x,n,mold,scratch)
      integer(psb_mpk_), intent(in) :: n
      class(psb_d_vect_type), intent(inout) :: x
      class(psb_d_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine d_vect_bld_mn
  end interface

  interface 
    module subroutine d_vect_bld_en(x,n,mold,scratch)
      integer(psb_epk_), intent(in) :: n
      class(psb_d_vect_type), intent(inout) :: x
      class(psb_d_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine d_vect_bld_en
  end interface

  interface 
    module function  d_vect_get_vect(x,n) result(res)
      class(psb_d_vect_type), intent(inout)  :: x
      real(psb_dpk_), allocatable                 :: res(:)
      integer(psb_ipk_), optional :: n
    end function d_vect_get_vect
  end interface

  interface 
    module subroutine d_vect_set_scal(x,val,first,last)
      class(psb_d_vect_type), intent(inout)  :: x
      real(psb_dpk_), intent(in) :: val
      integer(psb_ipk_), optional :: first, last
    end subroutine d_vect_set_scal
  end interface

  interface 
    module subroutine d_vect_set_vect(x,val,first,last)
      class(psb_d_vect_type), intent(inout) :: x
      real(psb_dpk_), intent(in)         :: val(:)
      integer(psb_ipk_), optional :: first, last
    end subroutine d_vect_set_vect
  end interface

  interface 
    module subroutine d_vect_check_addr(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_check_addr
  end interface

  interface 
    module function d_vect_get_nrows(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function d_vect_get_nrows
  end interface

  interface 
    module function d_vect_sizeof(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function d_vect_sizeof
  end interface

  interface 
    module function d_vect_get_fmt(x) result(res)
      class(psb_d_vect_type), intent(in) :: x
      character(len=5) :: res
    end function d_vect_get_fmt
  end interface

  interface 
    module subroutine d_vect_all(n, x, info, mold)
      integer(psb_ipk_), intent(in)           :: n
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)      :: info
      class(psb_d_base_vect_type), intent(in), optional :: mold
    end subroutine d_vect_all
  end interface

  interface 
    module subroutine d_vect_reinit(x, info, clear)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: clear
    end subroutine d_vect_reinit
  end interface
 
  interface 
    module subroutine d_vect_reall(n, x, info)
      integer(psb_ipk_), intent(in)         :: n
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
    end subroutine d_vect_reall
  end interface

  interface 
    module subroutine d_vect_zero(x)
      class(psb_d_vect_type), intent(inout)    :: x
    end subroutine d_vect_zero
  end interface

  interface 
    module subroutine d_vect_asb(n, x, info, scratch)
      integer(psb_ipk_), intent(in)         :: n
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional        :: scratch
    end subroutine d_vect_asb
  end interface

  interface 
    module subroutine d_vect_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_dpk_) :: alpha, beta, y(:)
      class(psb_d_vect_type) :: x
    end subroutine d_vect_gthab
  end interface

  interface 
    module subroutine d_vect_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_dpk_) ::  y(:)
      class(psb_d_vect_type) :: x
    end subroutine d_vect_gthzv
  end interface

  interface 
    module subroutine d_vect_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_dpk_) :: beta, x(:)
      class(psb_d_vect_type) :: y
    end subroutine d_vect_sctb
  end interface

  interface 
    module subroutine d_vect_free(x, info)
      class(psb_d_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_free
  end interface
  
  interface 
    module subroutine d_vect_ins_a(n,irl,val,x,maxr,info)
      class(psb_d_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      real(psb_dpk_), intent(in)        :: val(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_ins_a
  end interface

  interface 
    module subroutine d_vect_ins_v(n,irl,val,x,maxr,info)
      class(psb_d_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, maxr
      class(psb_i_vect_type), intent(inout)       :: irl
      class(psb_d_vect_type), intent(inout)       :: val
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_ins_v
  end interface

  interface 
    module subroutine d_vect_cnv(x,mold)
      class(psb_d_vect_type), intent(inout) :: x
      class(psb_d_base_vect_type), intent(in), optional :: mold
      class(psb_d_base_vect_type), allocatable :: tmp
    end subroutine d_vect_cnv
  end interface
  
  interface 
    module subroutine d_vect_sync(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_sync
  end interface

  interface 
    module subroutine d_vect_set_sync(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_set_sync
  end interface

  interface 
    module subroutine d_vect_set_host(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_set_host
  end interface

  interface 
    module subroutine d_vect_set_dev(x)
      class(psb_d_vect_type), intent(inout) :: x
    end subroutine d_vect_set_dev
  end interface

  interface 
    module function d_vect_is_sync(x) result(res)
      logical :: res
      class(psb_d_vect_type), intent(inout) :: x
    end function d_vect_is_sync
  end interface

  interface 
    module function d_vect_is_host(x) result(res)
      logical :: res
      class(psb_d_vect_type), intent(inout) :: x
    end function d_vect_is_host
  end interface

  interface 
    module function d_vect_is_dev(x) result(res)
      logical :: res
      class(psb_d_vect_type), intent(inout) :: x
    end function d_vect_is_dev
  end interface


  interface 
    module function d_vect_get_entry(x,index) result(res)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)        :: index
      real(psb_dpk_) :: res
    end function d_vect_get_entry
  end interface

  interface 
    module subroutine d_vect_set_entry(x,index,val) 
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)        :: index
      real(psb_dpk_) :: val
    end subroutine d_vect_set_entry
  end interface

  interface 
    module function d_vect_dot_v(n,x,y) result(res)
      class(psb_d_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)           :: n
      real(psb_dpk_)                :: res
    end function d_vect_dot_v
  end interface

  interface 
    module function d_vect_dot_a(n,x,y) result(res)
      class(psb_d_vect_type), intent(inout) :: x
      real(psb_dpk_), intent(in)    :: y(:)
      integer(psb_ipk_), intent(in)           :: n
      real(psb_dpk_)                :: res
    end function d_vect_dot_a
  end interface

  interface 
    module subroutine d_vect_axpby_v(m,alpha, x, beta, y, info)
      integer(psb_ipk_), intent(in)               :: m
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_axpby_v
  end interface

  interface 
    module subroutine d_vect_axpby_v2(m,alpha, x, beta, y, z, info)
      integer(psb_ipk_), intent(in)            :: m
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      class(psb_d_vect_type), intent(inout)  :: z
      real(psb_dpk_), intent (in)             :: alpha, beta
      integer(psb_ipk_), intent(out)           :: info
    end subroutine d_vect_axpby_v2
  end interface

  interface 
    module subroutine d_vect_axpby_a(m,alpha, x, beta, y, info)
      integer(psb_ipk_), intent(in)               :: m
      real(psb_dpk_), intent(in)        :: x(:)
      class(psb_d_vect_type), intent(inout)  :: y
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_axpby_a
  end interface

  interface 
    module subroutine d_vect_axpby_a2(m,alpha, x, beta, y, z, info)
      integer(psb_ipk_), intent(in)               :: m
      real(psb_dpk_), intent(in)                 :: x(:)
      real(psb_dpk_), intent(in)                 :: y(:)
      class(psb_d_vect_type), intent(inout)     :: z
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_axpby_a2
  end interface

  interface 
    module subroutine d_vect_upd_xyz(m,alpha,beta,gamma,delta,x, y, z, info)
      integer(psb_ipk_), intent(in)            :: m
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      class(psb_d_vect_type), intent(inout)  :: z
      real(psb_dpk_), intent (in)     :: alpha, beta, gamma, delta
      integer(psb_ipk_), intent(out)   :: info
    end subroutine d_vect_upd_xyz
  end interface

  interface 
    module subroutine d_vect_xyzw(m,a,b,c,d,e,f,x, y, z, w, info)
      integer(psb_ipk_), intent(in)            :: m
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      class(psb_d_vect_type), intent(inout)  :: z
      class(psb_d_vect_type), intent(inout)  :: w
      real(psb_dpk_), intent (in)     :: a, b, c, d, e, f
      integer(psb_ipk_), intent(out)   :: info
    end subroutine d_vect_xyzw
  end interface

  interface 
    module subroutine d_vect_mlt_v(x, y, info)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_mlt_v
  end interface

  interface 
    module subroutine d_vect_mlt_a(x, y, info)
      real(psb_dpk_), intent(in)        :: x(:)
      class(psb_d_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_mlt_a
  end interface

  interface 
    module subroutine d_vect_mlt_a_2(alpha,x,y,beta,z,info)
      real(psb_dpk_), intent(in)         :: alpha,beta
      real(psb_dpk_), intent(in)         :: y(:)
      real(psb_dpk_), intent(in)         :: x(:)
      class(psb_d_vect_type), intent(inout) :: z
      integer(psb_ipk_), intent(out)                  :: info
    end subroutine d_vect_mlt_a_2
  end interface

  interface 
    module subroutine d_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
      real(psb_dpk_), intent(in)          :: alpha,beta
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)                   :: info
      character(len=1), intent(in), optional :: conjgx, conjgy
    end subroutine d_vect_mlt_v_2
  end interface

  interface 
    module subroutine d_vect_mlt_av(alpha,x,y,beta,z,info)
      real(psb_dpk_), intent(in)        :: alpha,beta
      real(psb_dpk_), intent(in)        :: x(:)
      class(psb_d_vect_type), intent(inout)  :: y
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_mlt_av
  end interface

  interface 
    module subroutine d_vect_mlt_va(alpha,x,y,beta,z,info)
      real(psb_dpk_), intent(in)        :: alpha,beta
      real(psb_dpk_), intent(in)        :: y(:)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_mlt_va
  end interface

  interface 
    module subroutine d_vect_div_v(x, y, info)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_div_v
  end interface

  interface 
    module subroutine d_vect_div_v2( x, y, z, info)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_div_v2
  end interface

  interface 
    module subroutine d_vect_div_v_check(x, y, info, flag)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine d_vect_div_v_check
  end interface

  interface 
    module subroutine d_vect_div_v2_check(x, y, z, info, flag)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine d_vect_div_v2_check
  end interface

  interface 
    module subroutine d_vect_div_a2(x, y, z, info)
      real(psb_dpk_), intent(in) :: x(:)
      real(psb_dpk_), intent(in)    :: y(:)
      class(psb_d_vect_type), intent(inout) :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_div_a2
  end interface

  interface 
    module subroutine d_vect_div_a2_check(x, y, z, info,flag)
      real(psb_dpk_), intent(in) :: x(:)
      real(psb_dpk_), intent(in)    :: y(:)
      class(psb_d_vect_type), intent(inout) :: z
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine d_vect_div_a2_check
  end interface

  interface 
    module subroutine d_vect_inv_v(x, y, info)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_inv_v
  end interface

  interface 
    module subroutine d_vect_inv_v_check(x, y, info, flag)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine d_vect_inv_v_check
  end interface

  interface 
    module subroutine d_vect_inv_a2(x, y, info)
      real(psb_dpk_), intent(inout) :: x(:)
      class(psb_d_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_vect_inv_a2
  end interface

  interface 
    module subroutine d_vect_inv_a2_check(x, y, info,flag)
      real(psb_dpk_), intent(inout) :: x(:)
      class(psb_d_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine d_vect_inv_a2_check
  end interface

  interface 
    module subroutine d_vect_acmp_a2(x,c,z,info)
      real(psb_dpk_), intent(in)              :: c
      real(psb_dpk_), intent(inout)           :: x(:)
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine d_vect_acmp_a2
  end interface

  interface 
    module subroutine d_vect_acmp_v2(x,c,z,info)
      real(psb_dpk_), intent(in)              :: c
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine d_vect_acmp_v2
  end interface

  interface 
    module subroutine d_vect_scal(alpha, x)
      class(psb_d_vect_type), intent(inout)  :: x
      real(psb_dpk_), intent (in)       :: alpha
    end subroutine d_vect_scal
  end interface

  interface 
    module subroutine d_vect_absval1(x)
      class(psb_d_vect_type), intent(inout)  :: x
    end subroutine d_vect_absval1
  end interface

  interface 
    module subroutine d_vect_absval2(x,y)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
    end subroutine d_vect_absval2
  end interface

  interface 
    module function d_vect_nrm2(n,x) result(res)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_dpk_)                :: res
    end function d_vect_nrm2
  end interface

  interface 
    module function d_vect_nrm2_weight(n,x,w,aux) result(res)
      class(psb_d_vect_type), intent(inout) :: x
      class(psb_d_vect_type), intent(inout) :: w
      class(psb_d_vect_type), intent(inout), optional :: aux
      integer(psb_ipk_), intent(in)           :: n
      real(psb_dpk_)                        :: res
    end function d_vect_nrm2_weight
  end interface
 
  interface 
    module function d_vect_nrm2_weight_mask(n,x,w,id,info,aux) result(res)
      class(psb_d_vect_type), intent(inout) :: x
      class(psb_d_vect_type), intent(inout) :: w
      class(psb_d_vect_type), intent(inout) :: id
      integer(psb_ipk_), intent(in)           :: n
      real(psb_dpk_)                        :: res
      integer(psb_ipk_), intent(out)          :: info
      class(psb_d_vect_type), intent(inout), optional :: aux
    end function d_vect_nrm2_weight_mask
  end interface

  interface 
    module function d_vect_amax(n,x) result(res)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_dpk_)                :: res
    end function d_vect_amax
  end interface

  interface 
    module function d_vect_min(n,x) result(res)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_dpk_)                :: res
    end function d_vect_min
  end interface

  interface 
    module function d_vect_asum(n,x) result(res)
      class(psb_d_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_dpk_)                :: res
    end function d_vect_asum
  end interface


  interface 
    module subroutine d_vect_mask_a(c,x,m,t,info)
      real(psb_dpk_), intent(inout)          :: c(:)
      real(psb_dpk_), intent(inout)          :: x(:)
      logical, intent(out)                     :: t;
      class(psb_d_vect_type), intent(inout)  :: m
      integer(psb_ipk_), intent(out)           :: info
    end subroutine d_vect_mask_a
  end interface

  interface 
    module subroutine d_vect_mask_v(c,x,m,t,info)
      class(psb_d_vect_type), intent(inout)  :: c
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: m
      logical, intent(out)                     :: t;
      integer(psb_ipk_), intent(out)           :: info
    end subroutine d_vect_mask_v
  end interface

  interface 
    module function d_vect_minquotient_v(x, y, info) result(z)
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: y
      real(psb_dpk_)                         :: z
      integer(psb_ipk_), intent(out)           :: info
    end function d_vect_minquotient_v
  end interface

  interface 
    module function d_vect_minquotient_a2(x, y, info) result(z)
      class(psb_d_vect_type), intent(inout)  :: x
      real(psb_dpk_), intent(inout)           :: y(:)
      integer(psb_ipk_), intent(out)           :: info
      real(psb_dpk_)                         :: z
    end function d_vect_minquotient_a2
  end interface



  interface 
    module subroutine d_vect_addconst_a2(x,b,z,info)
      real(psb_dpk_), intent(in)             :: b
      real(psb_dpk_), intent(inout)           :: x(:)
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine d_vect_addconst_a2
  end interface

  interface 
    module subroutine d_vect_addconst_v2(x,b,z,info)
      real(psb_dpk_), intent(in)             :: b
      class(psb_d_vect_type), intent(inout)  :: x
      class(psb_d_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine d_vect_addconst_v2
  end interface

contains

  function constructor(x) result(this)
    real(psb_dpk_)   :: x(:)
    type(psb_d_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,kind=psb_ipk_),info)

  end function constructor


  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_d_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(n)
    call this%asb(n,info)

  end function size_const

end module psb_d_vect_mod

module psb_d_multivect_mod

  use psb_d_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod


  !private

  type psb_d_multivect_type
    class(psb_d_base_multivect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build=psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    real(psb_dpk_), allocatable :: rmtv(:,:)
  contains
    procedure, pass(x) :: get_nrows => d_mvect_get_nrows
    procedure, pass(x) :: get_ncols => d_mvect_get_ncols
    procedure, pass(x) :: sizeof   => d_mvect_sizeof
    procedure, pass(x) :: get_fmt  => d_mvect_get_fmt
    procedure, pass(x) :: is_remote_build => d_mvect_is_remote_build
    procedure, pass(x) :: set_remote_build => d_mvect_set_remote_build
    procedure, pass(x) :: get_dupl => d_mvect_get_dupl
    procedure, pass(x) :: set_dupl => d_mvect_set_dupl

    procedure, pass(x) :: all      => d_mvect_all
    procedure, pass(x) :: reall    => d_mvect_reall
    procedure, pass(x) :: zero     => d_mvect_zero
    procedure, pass(x) :: asb      => d_mvect_asb
    procedure, pass(x) :: sync     => d_mvect_sync
    procedure, pass(x) :: free     => d_mvect_free
    procedure, pass(x) :: ins      => d_mvect_ins
    procedure, pass(x) :: bld_x    => d_mvect_bld_x
    procedure, pass(x) :: bld_n    => d_mvect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: get_vect => d_mvect_get_vect
    procedure, pass(x) :: cnv      => d_mvect_cnv
    procedure, pass(x) :: set_scal => d_mvect_set_scal
    procedure, pass(x) :: set_vect => d_mvect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => d_mvect_clone
    procedure, pass(x) :: gthab    => d_mvect_gthab
    procedure, pass(x) :: gthzv    => d_mvect_gthzv
    procedure, pass(x) :: gthzv_x  => d_mvect_gthzv_x
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => d_mvect_sctb
    procedure, pass(y) :: sctb_x   => d_mvect_sctb_x
    generic, public    :: sct      => sctb, sctb_x
!!$    procedure, pass(x) :: dot_v    => d_mvect_dot_v
!!$    procedure, pass(x) :: dot_a    => d_mvect_dot_a
!!$    generic, public    :: dot      => dot_v, dot_a
!!$    procedure, pass(y) :: axpby_v  => d_mvect_axpby_v
!!$    procedure, pass(y) :: axpby_a  => d_mvect_axpby_a
!!$    generic, public    :: axpby    => axpby_v, axpby_a
!!$    procedure, pass(y) :: mlt_v    => d_mvect_mlt_v
!!$    procedure, pass(y) :: mlt_a    => d_mvect_mlt_a
!!$    procedure, pass(z) :: mlt_a_2  => d_mvect_mlt_a_2
!!$    procedure, pass(z) :: mlt_v_2  => d_mvect_mlt_v_2
!!$    procedure, pass(z) :: mlt_va   => d_mvect_mlt_va
!!$    procedure, pass(z) :: mlt_av   => d_mvect_mlt_av
!!$    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
!!$         & mlt_v_2, mlt_av, mlt_va
!!$    procedure, pass(x) :: scal     => d_mvect_scal
!!$    procedure, pass(x) :: nrm2     => d_mvect_nrm2
!!$    procedure, pass(x) :: amax     => d_mvect_amax
!!$    procedure, pass(x) :: asum     => d_mvect_asum
  end type psb_d_multivect_type

  public  :: psb_d_multivect, psb_d_multivect_type,&
       & psb_set_multivect_default, psb_get_multivect_default, &
       & psb_d_base_multivect_type

  private
  interface psb_d_multivect
    module procedure constructor, size_const
  end interface psb_d_multivect

  class(psb_d_base_multivect_type), allocatable, target,&
       & save, private :: psb_d_base_multivect_default

  
  interface 
    module function d_mvect_get_dupl(x) result(res)
      class(psb_d_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function d_mvect_get_dupl
  end interface

  interface 
    module subroutine d_mvect_set_dupl(x,val)
      class(psb_d_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine d_mvect_set_dupl
  end interface

  interface 
    module function d_mvect_is_remote_build(x) result(res)
      class(psb_d_multivect_type), intent(in) :: x
      logical :: res
    end function d_mvect_is_remote_build
  end interface

  interface 
    module subroutine d_mvect_set_remote_build(x,val)
      class(psb_d_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine d_mvect_set_remote_build
  end interface

  interface 
    module subroutine  psb_d_set_multivect_default(v)
      class(psb_d_base_multivect_type), intent(in) :: v
    end subroutine psb_d_set_multivect_default
  end interface

  interface 
    module function psb_d_get_multivect_default(v) result(res)
      class(psb_d_multivect_type), intent(in) :: v
      class(psb_d_base_multivect_type), pointer :: res
    end function psb_d_get_multivect_default
  end interface

  interface 
    module function psb_d_get_base_multivect_default() result(res)
      class(psb_d_base_multivect_type), pointer :: res
    end function psb_d_get_base_multivect_default
  end interface

  interface 
    module subroutine d_mvect_clone(x,y,info)
      class(psb_d_multivect_type), intent(inout) :: x
      class(psb_d_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)        :: info
    end subroutine d_mvect_clone
  end interface

  interface 
    module subroutine d_mvect_bld_x(x,invect,mold)
      real(psb_dpk_), intent(in)          :: invect(:,:)
      class(psb_d_multivect_type), intent(out) :: x
      class(psb_d_base_multivect_type), intent(in), optional :: mold
    end subroutine d_mvect_bld_x
  end interface


  interface 
    module subroutine d_mvect_bld_n(x,m,n,mold,scratch)
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_d_multivect_type), intent(out) :: x
      class(psb_d_base_multivect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine d_mvect_bld_n
  end interface

  interface 
    module function  d_mvect_get_vect(x) result(res)
      class(psb_d_multivect_type), intent(inout)  :: x
      real(psb_dpk_), allocatable                 :: res(:,:)
    end function d_mvect_get_vect
  end interface

  interface 
    module subroutine d_mvect_set_scal(x,val)
      class(psb_d_multivect_type), intent(inout)  :: x
      real(psb_dpk_), intent(in) :: val
    end subroutine d_mvect_set_scal
  end interface

  interface 
    module subroutine d_mvect_set_vect(x,val)
      class(psb_d_multivect_type), intent(inout) :: x
      real(psb_dpk_), intent(in)         :: val(:,:)
    end subroutine d_mvect_set_vect
  end interface

  interface 
    module function d_mvect_get_nrows(x) result(res)
      class(psb_d_multivect_type), intent(in) :: x
      integer(psb_ipk_)  :: res
    end function d_mvect_get_nrows
  end interface

  interface 
    module function d_mvect_get_ncols(x) result(res)
      class(psb_d_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function d_mvect_get_ncols
  end interface

  interface 
    module function d_mvect_sizeof(x) result(res)
      class(psb_d_multivect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function d_mvect_sizeof
  end interface

  interface 
    module function d_mvect_get_fmt(x) result(res)
      class(psb_d_multivect_type), intent(in) :: x
      character(len=5) :: res
    end function d_mvect_get_fmt
  end interface

  interface 
    module subroutine d_mvect_all(m,n, x, info, mold)
      integer(psb_ipk_), intent(in)       :: m,n
      class(psb_d_multivect_type), intent(out) :: x
      class(psb_d_base_multivect_type), intent(in), optional :: mold
      integer(psb_ipk_), intent(out)      :: info
    end subroutine d_mvect_all
  end interface

  interface 
    module subroutine d_mvect_reall(m,n, x, info)
      integer(psb_ipk_), intent(in)         :: m,n
      class(psb_d_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
    end subroutine d_mvect_reall
  end interface

  interface 
    module subroutine d_mvect_zero(x)
      class(psb_d_multivect_type), intent(inout)    :: x
    end subroutine d_mvect_zero
  end interface

  interface 
    module subroutine d_mvect_asb(m,n, x, info)
      integer(psb_ipk_), intent(in)              :: m,n
      class(psb_d_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine d_mvect_asb
  end interface

  interface 
    module subroutine d_mvect_sync(x)
      class(psb_d_multivect_type), intent(inout) :: x
    end subroutine d_mvect_sync
  end interface

  interface 
    module subroutine d_mvect_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_dpk_) :: alpha, beta, y(:)
      class(psb_d_multivect_type) :: x
    end subroutine d_mvect_gthab
  end interface

  interface 
    module subroutine d_mvect_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_dpk_) ::  y(:)
      class(psb_d_multivect_type) :: x
    end subroutine d_mvect_gthzv
  end interface

  interface 
    module subroutine d_mvect_gthzv_x(i,n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      real(psb_dpk_) ::  y(:)
      class(psb_d_multivect_type) :: x
    end subroutine d_mvect_gthzv_x
  end interface

  interface 
    module subroutine d_mvect_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      real(psb_dpk_) :: beta, x(:)
      class(psb_d_multivect_type) :: y
    end subroutine d_mvect_sctb
  end interface

  interface 
    module subroutine d_mvect_sctb_x(i,n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      real(psb_dpk_) :: beta, x(:)
      class(psb_d_multivect_type) :: y
    end subroutine d_mvect_sctb_x
  end interface

  interface 
    module subroutine d_mvect_free(x, info)
      class(psb_d_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_mvect_free
  end interface

  interface 
    module subroutine d_mvect_ins(n,irl,val,x,maxr,info)
      class(psb_d_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n,maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      real(psb_dpk_), intent(in)        :: val(:,:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine d_mvect_ins
  end interface

  interface 
    module subroutine d_mvect_cnv(x,mold)
      class(psb_d_multivect_type), intent(inout) :: x
      class(psb_d_base_multivect_type), intent(in), optional :: mold
    end subroutine d_mvect_cnv
  end interface


!!$  module function d_mvect_dot_v(n,x,y) result(res)
!!$    implicit none
!!$    class(psb_d_multivect_type), intent(inout) :: x, y
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_dpk_)                :: res
!!$
!!$    res = dzero
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & res = x%v%dot(n,y%v)
!!$
!!$  end function d_mvect_dot_v
!!$
!!$  module function d_mvect_dot_a(n,x,y) result(res)
!!$    implicit none
!!$    class(psb_d_multivect_type), intent(inout) :: x
!!$    real(psb_dpk_), intent(in)    :: y(:)
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_dpk_)                :: res
!!$
!!$    res = dzero
!!$    if (allocated(x%v)) &
!!$         & res = x%v%dot(n,y)
!!$
!!$  end function d_mvect_dot_a
!!$
!!$  module subroutine d_mvect_axpby_v(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    class(psb_d_multivect_type), intent(inout)  :: x
!!$    class(psb_d_multivect_type), intent(inout)  :: y
!!$    real(psb_dpk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(x%v).and.allocated(y%v)) then
!!$      call y%v%axpby(m,alpha,x%v,beta,info)
!!$    else
!!$      info = psb_err_invalid_mvect_state_
!!$    end if
!!$
!!$  end subroutine d_mvect_axpby_v
!!$
!!$  module subroutine d_mvect_axpby_a(m,alpha, x, beta, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    real(psb_dpk_), intent(in)        :: x(:)
!!$    class(psb_d_multivect_type), intent(inout)  :: y
!!$    real(psb_dpk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(y%v)) &
!!$         & call y%v%axpby(m,alpha,x,beta,info)
!!$
!!$  end subroutine d_mvect_axpby_a
!!$
!!$
!!$  module subroutine d_mvect_mlt_v(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_d_multivect_type), intent(inout)  :: x
!!$    class(psb_d_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%mlt(x%v,info)
!!$
!!$  end subroutine d_mvect_mlt_v
!!$
!!$  module subroutine d_mvect_mlt_a(x, y, info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_dpk_), intent(in)        :: x(:)
!!$    class(psb_d_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%mlt(x,info)
!!$
!!$  end subroutine d_mvect_mlt_a
!!$
!!$
!!$  module subroutine d_mvect_mlt_a_2(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_dpk_), intent(in)         :: alpha,beta
!!$    real(psb_dpk_), intent(in)         :: y(:)
!!$    real(psb_dpk_), intent(in)         :: x(:)
!!$    class(psb_d_multivect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)                  :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%mlt(alpha,x,y,beta,info)
!!$
!!$  end subroutine d_mvect_mlt_a_2
!!$
!!$  module subroutine d_mvect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_dpk_), intent(in)          :: alpha,beta
!!$    class(psb_d_multivect_type), intent(inout)  :: x
!!$    class(psb_d_multivect_type), intent(inout)  :: y
!!$    class(psb_d_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)                   :: info
!!$    character(len=1), intent(in), optional :: conjgx, conjgy
!!$
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v).and.&
!!$         & allocated(z%v)) &
!!$         & call z%v%mlt(alpha,x%v,y%v,beta,info,conjgx,conjgy)
!!$
!!$  end subroutine d_mvect_mlt_v_2

!!$  module subroutine d_mvect_mlt_av(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_dpk_), intent(in)        :: alpha,beta
!!$    real(psb_dpk_), intent(in)        :: x(:)
!!$    class(psb_d_multivect_type), intent(inout)  :: y
!!$    class(psb_d_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v).and.allocated(y%v)) &
!!$         & call z%v%mlt(alpha,x,y%v,beta,info)
!!$
!!$  end subroutine d_mvect_mlt_av
!!$
!!$  module subroutine d_mvect_mlt_va(alpha,x,y,beta,z,info)
!!$    use psi_serial_mod
!!$    implicit none
!!$    real(psb_dpk_), intent(in)        :: alpha,beta
!!$    real(psb_dpk_), intent(in)        :: y(:)
!!$    class(psb_d_multivect_type), intent(inout)  :: x
!!$    class(psb_d_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    if (allocated(z%v).and.allocated(x%v)) &
!!$         & call z%v%mlt(alpha,x%v,y,beta,info)
!!$
!!$  end subroutine d_mvect_mlt_va
!!$
!!$  module subroutine d_mvect_scal(alpha, x)
!!$    use psi_serial_mod
!!$    implicit none
!!$    class(psb_d_multivect_type), intent(inout)  :: x
!!$    real(psb_dpk_), intent (in)       :: alpha
!!$
!!$    if (allocated(x%v)) call x%v%scal(alpha)
!!$
!!$  end subroutine d_mvect_scal
!!$
!!$
!!$  module function d_mvect_nrm2(n,x) result(res)
!!$    implicit none
!!$    class(psb_d_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_dpk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%nrm2(n)
!!$    else
!!$      res = dzero
!!$    end if
!!$
!!$  end function d_mvect_nrm2
!!$
!!$  module function d_mvect_amax(n,x) result(res)
!!$    implicit none
!!$    class(psb_d_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_dpk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%amax(n)
!!$    else
!!$      res = dzero
!!$    end if
!!$
!!$  end function d_mvect_amax
!!$
!!$  module function d_mvect_asum(n,x) result(res)
!!$    implicit none
!!$    class(psb_d_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_dpk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%asum(n)
!!$    else
!!$      res = dzero
!!$    end if
!!$
!!$  end function d_mvect_asum

contains
  
 function constructor(x) result(this)
    real(psb_dpk_)   :: x(:,:)
    type(psb_d_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld_x(x)
    call this%asb(size(x,dim=1,kind=psb_ipk_),size(x,dim=2,kind=psb_ipk_),info)

  end function constructor

  function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_d_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld_n(m,n)
    call this%asb(m,n,info)

  end function size_const

end module psb_d_multivect_mod
