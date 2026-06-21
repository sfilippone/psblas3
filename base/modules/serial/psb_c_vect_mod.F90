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
! package: psb_c_vect_mod
!
! This module contains the definition of the psb_c_vect type which
! is the outer container for dense vectors.
! Therefore all methods simply invoke the corresponding methods of the
! inner component.
!
module psb_c_vect_mod

  use psb_realloc_mod
  use psb_c_base_vect_mod
  use psb_i_vect_mod

  type psb_c_vect_type
    class(psb_c_base_vect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build=psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    complex(psb_spk_), allocatable :: rmtv(:)
    integer(psb_lpk_), allocatable :: rmidx(:)
  contains
    procedure, pass(x) :: get_nrows => c_vect_get_nrows
    procedure, pass(x) :: sizeof   => c_vect_sizeof
    procedure, pass(x) :: get_fmt  => c_vect_get_fmt
    procedure, pass(x) :: is_remote_build => c_vect_is_remote_build
    procedure, pass(x) :: set_remote_build => c_vect_set_remote_build
    procedure, pass(x) :: get_nrmv => c_vect_get_nrmv
    procedure, pass(x) :: set_nrmv => c_vect_set_nrmv
    procedure, pass(x) :: all      => c_vect_all
    procedure, pass(x) :: reall    => c_vect_reall
    procedure, pass(x) :: zero     => c_vect_zero
    procedure, pass(x) :: asb      => c_vect_asb
    procedure, pass(x) :: set_dupl => c_vect_set_dupl 
    procedure, pass(x) :: get_dupl => c_vect_get_dupl
    procedure, pass(x) :: set_ncfs => c_vect_set_ncfs 
    procedure, pass(x) :: get_ncfs => c_vect_get_ncfs
    procedure, pass(x) :: set_state => c_vect_set_state
    procedure, pass(x) :: set_null  => c_vect_set_null
    procedure, pass(x) :: set_bld   => c_vect_set_bld
    procedure, pass(x) :: set_upd   => c_vect_set_upd
    procedure, pass(x) :: set_asb   => c_vect_set_asb
    procedure, pass(x) :: get_state => c_vect_get_state
    procedure, pass(x) :: is_null   => c_vect_is_null
    procedure, pass(x) :: is_bld    => c_vect_is_bld
    procedure, pass(x) :: is_upd    => c_vect_is_upd
    procedure, pass(x) :: is_asb    => c_vect_is_asb
    procedure, pass(x) :: reinit    => c_vect_reinit

    procedure, pass(x) :: gthab    => c_vect_gthab
    procedure, pass(x) :: gthzv    => c_vect_gthzv
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => c_vect_sctb
    generic, public    :: sct      => sctb
    procedure, pass(x) :: free     => c_vect_free
    procedure, pass(x) :: ins_a    => c_vect_ins_a
    procedure, pass(x) :: ins_v    => c_vect_ins_v
    generic, public    :: ins      => ins_v, ins_a
    procedure, pass(x) :: bld_x    => c_vect_bld_x
    procedure, pass(x) :: bld_mn   => c_vect_bld_mn
    procedure, pass(x) :: bld_en   => c_vect_bld_en
    generic, public    :: bld      => bld_x, bld_mn, bld_en
    procedure, pass(x) :: get_vect => c_vect_get_vect
    procedure, pass(x) :: cnv      => c_vect_cnv
    procedure, pass(x) :: set_scal => c_vect_set_scal
    procedure, pass(x) :: set_vect => c_vect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => c_vect_clone

    procedure, pass(x) :: sync     => c_vect_sync
    procedure, pass(x) :: is_host  => c_vect_is_host
    procedure, pass(x) :: is_dev   => c_vect_is_dev
    procedure, pass(x) :: is_sync  => c_vect_is_sync
    procedure, pass(x) :: set_host => c_vect_set_host
    procedure, pass(x) :: set_dev  => c_vect_set_dev
    procedure, pass(x) :: set_sync => c_vect_set_sync
    procedure, pass(x) :: check_addr => c_vect_check_addr

    procedure, pass(x) :: get_entry => c_vect_get_entry
    procedure, pass(x) :: set_entry => c_vect_set_entry

    procedure, pass(x) :: dot_v    => c_vect_dot_v
    procedure, pass(x) :: dot_a    => c_vect_dot_a
    generic, public    :: dot      => dot_v, dot_a
    procedure, pass(y) :: axpby_v  => c_vect_axpby_v
    procedure, pass(y) :: axpby_a  => c_vect_axpby_a
    procedure, pass(z) :: axpby_v2  => c_vect_axpby_v2
    procedure, pass(z) :: axpby_a2  => c_vect_axpby_a2
    generic, public    :: axpby    => axpby_v, axpby_a, axpby_v2, axpby_a2
    procedure, pass(z) :: upd_xyz  => c_vect_upd_xyz
    procedure, pass(z) :: xyzw     => c_vect_xyzw
    procedure, pass(y) :: mlt_v    => c_vect_mlt_v
    procedure, pass(y) :: mlt_a    => c_vect_mlt_a
    procedure, pass(z) :: mlt_a_2  => c_vect_mlt_a_2
    procedure, pass(z) :: mlt_v_2  => c_vect_mlt_v_2
    procedure, pass(z) :: mlt_va   => c_vect_mlt_va
    procedure, pass(z) :: mlt_av   => c_vect_mlt_av
    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
         & mlt_v_2, mlt_av, mlt_va
    procedure, pass(y) :: div_v    => c_vect_div_v
    procedure, pass(z) :: div_v2    => c_vect_div_v2
    procedure, pass(y) :: div_v_check => c_vect_div_v_check
    procedure, pass(y) :: div_v2_check => c_vect_div_v2_check
    procedure, pass(z) :: div_a2   => c_vect_div_a2
    procedure, pass(z) :: div_a2_check => c_vect_div_a2_check
    generic, public    :: div      => div_v, div_v2, div_v_check, &
                                        div_v2_check, div_a2, div_a2_check
    procedure, pass(y) :: inv_v    => c_vect_inv_v
    procedure, pass(y) :: inv_v_check => c_vect_inv_v_check
    procedure, pass(y) :: inv_a2   => c_vect_inv_a2
    procedure, pass(y) :: inv_a2_check => c_vect_inv_a2_check
    generic, public    :: inv      => inv_v, inv_v_check, inv_a2, inv_a2_check
    procedure, pass(x) :: scal     => c_vect_scal
    procedure, pass(x) :: absval1  => c_vect_absval1
    procedure, pass(x) :: absval2  => c_vect_absval2
    generic, public    :: absval   => absval1, absval2
    procedure, pass(x) :: nrm2std  => c_vect_nrm2
    procedure, pass(x) :: nrm2weight => c_vect_nrm2_weight
    procedure, pass(x) :: nrm2weightmask => c_vect_nrm2_weight_mask
    generic, public    :: nrm2     => nrm2std, nrm2weight, nrm2weightmask
    procedure, pass(x) :: amax     => c_vect_amax
    procedure, pass(x) :: asum     => c_vect_asum
    procedure, pass(z) :: acmp_a2   => c_vect_acmp_a2
    procedure, pass(z) :: acmp_v2   => c_vect_acmp_v2
    generic, public    :: acmp      => acmp_a2, acmp_v2
    procedure, pass(z) :: addconst_a2   => c_vect_addconst_a2
    procedure, pass(z) :: addconst_v2   => c_vect_addconst_v2
    generic, public    :: addconst      => addconst_a2, addconst_v2


  end type psb_c_vect_type

  public  :: psb_c_vect, psb_c_vect_type,&
       & psb_c_set_vect_default, psb_c_get_vect_default, &
       & psb_c_clear_vect_default, psb_c_base_vect_type
  
  private :: constructor, size_const
  interface psb_c_vect
    module procedure constructor, size_const
  end interface psb_c_vect

  private :: c_vect_get_nrows, c_vect_sizeof, c_vect_get_fmt, &
       & c_vect_all, c_vect_reall, c_vect_zero,  c_vect_asb, &
       & c_vect_gthab, c_vect_gthzv, c_vect_sctb, &
       & c_vect_free, c_vect_ins_a, c_vect_ins_v, c_vect_bld_x, &
       & c_vect_bld_mn, c_vect_bld_en, c_vect_get_vect, &
       & c_vect_cnv, c_vect_set_scal, &
       & c_vect_set_vect, c_vect_clone, c_vect_sync, c_vect_is_host, &
       & c_vect_is_dev, c_vect_is_sync, c_vect_set_host, &
       & c_vect_set_dev, c_vect_set_sync, &
       & c_vect_set_remote_build, c_is_remote_build, &
       & c_vect_set_dupl, c_get_dupl, c_vect_set_nrmv, c_get_nrmv

  private ::  c_vect_dot_v, c_vect_dot_a, c_vect_axpby_v, c_vect_axpby_a, &
       & c_vect_mlt_v, c_vect_mlt_a, c_vect_mlt_a_2, c_vect_mlt_v_2, &
       & c_vect_mlt_va, c_vect_mlt_av, c_vect_scal, c_vect_absval1, &
       & c_vect_absval2, c_vect_nrm2, c_vect_amax, c_vect_asum


  class(psb_c_base_vect_type), allocatable, target,&
       & save, private :: psb_c_base_vect_default


  interface 
    module function c_vect_get_dupl(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function c_vect_get_dupl
  end interface
  
  interface 
    module subroutine c_vect_set_dupl(x,val)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine c_vect_set_dupl
  end interface
  
  interface 
    module function c_vect_get_ncfs(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function c_vect_get_ncfs
  end interface

  interface 
    module subroutine c_vect_set_ncfs(x,val)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine c_vect_set_ncfs
  end interface

  interface 
    module function c_vect_get_state(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function c_vect_get_state
  end interface

  interface 
    module function c_vect_is_null(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      logical :: res
    end function c_vect_is_null
  end interface

  interface 
    module function c_vect_is_bld(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      logical :: res
    end function c_vect_is_bld
  end interface

  interface 
    module function c_vect_is_upd(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      logical :: res
    end function c_vect_is_upd
  end interface

  interface 
    module function c_vect_is_asb(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      logical :: res
    end function c_vect_is_asb
  end interface

  interface 
    module subroutine  c_vect_set_state(n,x)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: n
    end subroutine c_vect_set_state
  end interface

  interface 
    module subroutine  c_vect_set_null(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_set_null
  end interface

  interface 
    module subroutine  c_vect_set_bld(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_set_bld
  end interface

  interface 
    module subroutine  c_vect_set_upd(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_set_upd
  end interface

  interface 
    module subroutine  c_vect_set_asb(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_set_asb
  end interface

  interface 
    module function c_vect_get_nrmv(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function c_vect_get_nrmv
  end interface

  interface 
    module subroutine c_vect_set_nrmv(x,val)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in) :: val
    end subroutine c_vect_set_nrmv
  end interface

  interface 
    module function c_vect_is_remote_build(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      logical :: res
    end function c_vect_is_remote_build
  end interface

  interface 
    module subroutine c_vect_set_remote_build(x,val)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine c_vect_set_remote_build
  end interface

  interface 
    module subroutine c_vect_clone(x,y,info)
      class(psb_c_vect_type), intent(inout) :: x
      class(psb_c_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)        :: info
    end subroutine c_vect_clone
  end interface

  interface 
    module subroutine c_vect_bld_x(x,invect,mold,scratch)
      complex(psb_spk_), intent(in)          :: invect(:)
      class(psb_c_vect_type), intent(inout) :: x
      class(psb_c_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine c_vect_bld_x
  end interface
  
  interface 
    module subroutine c_vect_bld_mn(x,n,mold,scratch)
      integer(psb_mpk_), intent(in) :: n
      class(psb_c_vect_type), intent(inout) :: x
      class(psb_c_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine c_vect_bld_mn
  end interface
  
  interface 
    module subroutine c_vect_bld_en(x,n,mold,scratch)
      integer(psb_epk_), intent(in) :: n
      class(psb_c_vect_type), intent(inout) :: x
      class(psb_c_base_vect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine c_vect_bld_en
  end interface

  interface 
    module function  c_vect_get_vect(x,n) result(res)
      class(psb_c_vect_type), intent(inout)  :: x
      complex(psb_spk_), allocatable                 :: res(:)
      integer(psb_ipk_), optional :: n
    end function c_vect_get_vect
  end interface
  
  interface 
    module subroutine c_vect_set_scal(x,val,first,last)
      class(psb_c_vect_type), intent(inout)  :: x
      complex(psb_spk_), intent(in) :: val
      integer(psb_ipk_), optional :: first, last
    end subroutine c_vect_set_scal
  end interface

  interface 
    module subroutine c_vect_set_vect(x,val,first,last)
      class(psb_c_vect_type), intent(inout) :: x
      complex(psb_spk_), intent(in)         :: val(:)
      integer(psb_ipk_), optional :: first, last
    end subroutine c_vect_set_vect
  end interface

  interface 
    module subroutine c_vect_check_addr(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_check_addr
  end interface

  interface 
    module function c_vect_get_nrows(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function c_vect_get_nrows
  end interface

  interface 
    module function c_vect_sizeof(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function c_vect_sizeof
  end interface

  interface 
    module function c_vect_get_fmt(x) result(res)
      class(psb_c_vect_type), intent(in) :: x
      character(len=5) :: res
    end function c_vect_get_fmt
  end interface

  interface 
    module subroutine c_vect_all(n, x, info, mold)
      integer(psb_ipk_), intent(in)           :: n
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)      :: info
      class(psb_c_base_vect_type), intent(in), optional :: mold
    end subroutine c_vect_all
  end interface

  interface 
    module subroutine c_vect_reinit(x, info, clear)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)      :: info
      logical, intent(in), optional       :: clear
    end subroutine c_vect_reinit
  end interface

  interface 
    module subroutine c_vect_reall(n, x, info)
      integer(psb_ipk_), intent(in)         :: n
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
    end subroutine c_vect_reall
  end interface

  interface 
    module subroutine c_vect_zero(x)
      class(psb_c_vect_type), intent(inout)    :: x
    end subroutine c_vect_zero
  end interface

  interface 
    module subroutine c_vect_asb(n, x, info, scratch)
      integer(psb_ipk_), intent(in)         :: n
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
      logical, intent(in), optional        :: scratch
    end subroutine c_vect_asb
  end interface

  interface 
    module subroutine c_vect_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_spk_) :: alpha, beta, y(:)
      class(psb_c_vect_type) :: x
    end subroutine c_vect_gthab
  end interface
  
  interface 
    module subroutine c_vect_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_spk_) ::  y(:)
      class(psb_c_vect_type) :: x
    end subroutine c_vect_gthzv
  end interface

  interface 
    module subroutine c_vect_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_spk_) :: beta, x(:)
      class(psb_c_vect_type) :: y
    end subroutine c_vect_sctb
  end interface

  interface 
    module subroutine c_vect_free(x, info)
      class(psb_c_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_free
  end interface

  interface 
    module subroutine c_vect_ins_a(n,irl,val,x,maxr,info)
      class(psb_c_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      complex(psb_spk_), intent(in)        :: val(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_ins_a
  end interface
  
  interface 
    module subroutine c_vect_ins_v(n,irl,val,x,maxr,info)
      class(psb_c_vect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n, maxr
      class(psb_i_vect_type), intent(inout)       :: irl
      class(psb_c_vect_type), intent(inout)       :: val
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_ins_v
  end interface
  
  interface 
    module subroutine c_vect_cnv(x,mold)
      class(psb_c_vect_type), intent(inout) :: x
      class(psb_c_base_vect_type), intent(in), optional :: mold
    end subroutine c_vect_cnv
  end interface
  
  interface 
    module subroutine c_vect_sync(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_sync
  end interface
  
  interface 
    module subroutine c_vect_set_sync(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_set_sync
  end interface
  
  interface 
    module subroutine c_vect_set_host(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_set_host
  end interface
  
  interface 
    module subroutine c_vect_set_dev(x)
      class(psb_c_vect_type), intent(inout) :: x
    end subroutine c_vect_set_dev
  end interface
  
  interface 
    module function c_vect_is_sync(x) result(res)
      logical :: res
      class(psb_c_vect_type), intent(inout) :: x
    end function c_vect_is_sync
  end interface
  
  interface 
    module function c_vect_is_host(x) result(res)
      logical :: res
      class(psb_c_vect_type), intent(inout) :: x
    end function c_vect_is_host
  end interface
  
  interface 
    module function c_vect_is_dev(x) result(res)
      logical :: res
      class(psb_c_vect_type), intent(inout) :: x
    end function c_vect_is_dev
  end interface


  interface 
    module function c_vect_get_entry(x,index) result(res)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)        :: index
      complex(psb_spk_) :: res
    end function c_vect_get_entry
  end interface
  
  interface 
    module subroutine c_vect_set_entry(x,index,val)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)        :: index
      complex(psb_spk_) :: val
    end subroutine c_vect_set_entry
  end interface
  
  interface 
    module function c_vect_dot_v(n,x,y) result(res)
      class(psb_c_vect_type), intent(inout) :: x, y
      integer(psb_ipk_), intent(in)           :: n
      complex(psb_spk_)                :: res
    end function c_vect_dot_v
  end interface
  
  interface 
    module function c_vect_dot_a(n,x,y) result(res)
      class(psb_c_vect_type), intent(inout) :: x
      complex(psb_spk_), intent(in)    :: y(:)
      integer(psb_ipk_), intent(in)           :: n
      complex(psb_spk_)                :: res
    end function c_vect_dot_a
  end interface
  
  interface 
    module subroutine c_vect_axpby_v(m,alpha, x, beta, y, info)
      integer(psb_ipk_), intent(in)               :: m
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      complex(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_axpby_v
  end interface
  
  interface 
    module subroutine c_vect_axpby_v2(m,alpha, x, beta, y, z, info)
      integer(psb_ipk_), intent(in)            :: m
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      class(psb_c_vect_type), intent(inout)  :: z
      complex(psb_spk_), intent (in)             :: alpha, beta
      integer(psb_ipk_), intent(out)           :: info
    end subroutine c_vect_axpby_v2
  end interface

  interface 
    module subroutine c_vect_axpby_a(m,alpha, x, beta, y, info)
      integer(psb_ipk_), intent(in)               :: m
      complex(psb_spk_), intent(in)        :: x(:)
      class(psb_c_vect_type), intent(inout)  :: y
      complex(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_axpby_a
  end interface
 
  interface 
    module subroutine c_vect_axpby_a2(m,alpha, x, beta, y, z, info)
      integer(psb_ipk_), intent(in)               :: m
      complex(psb_spk_), intent(in)                 :: x(:)
      complex(psb_spk_), intent(in)                 :: y(:)
      class(psb_c_vect_type), intent(inout)     :: z
      complex(psb_spk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_axpby_a2
  end interface
  
  interface 
    module subroutine c_vect_upd_xyz(m,alpha,beta,gamma,delta,x, y, z, info)
      integer(psb_ipk_), intent(in)            :: m
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      class(psb_c_vect_type), intent(inout)  :: z
      complex(psb_spk_), intent (in)     :: alpha, beta, gamma, delta
      integer(psb_ipk_), intent(out)   :: info
    end subroutine c_vect_upd_xyz
  end interface
  
  interface 
    module subroutine c_vect_xyzw(m,a,b,c,d,e,f,x, y, z, w, info)
      integer(psb_ipk_), intent(in)            :: m
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      class(psb_c_vect_type), intent(inout)  :: z
      class(psb_c_vect_type), intent(inout)  :: w
      complex(psb_spk_), intent (in)     :: a, b, c, d, e, f
      integer(psb_ipk_), intent(out)   :: info
    end subroutine c_vect_xyzw
  end interface

  interface 
    module subroutine c_vect_mlt_v(x, y, info)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_mlt_v
  end interface
  
  interface 
    module subroutine c_vect_mlt_a(x, y, info)
      complex(psb_spk_), intent(in)        :: x(:)
      class(psb_c_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_mlt_a
  end interface
  
  interface 
    module subroutine c_vect_mlt_a_2(alpha,x,y,beta,z,info)
      complex(psb_spk_), intent(in)         :: alpha,beta
      complex(psb_spk_), intent(in)         :: y(:)
      complex(psb_spk_), intent(in)         :: x(:)
      class(psb_c_vect_type), intent(inout) :: z
      integer(psb_ipk_), intent(out)                  :: info
    end subroutine c_vect_mlt_a_2
  end interface
  
  interface 
    module subroutine c_vect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
      complex(psb_spk_), intent(in)          :: alpha,beta
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)                   :: info
      character(len=1), intent(in), optional :: conjgx, conjgy
    end subroutine c_vect_mlt_v_2
  end interface

  interface 
    module subroutine c_vect_mlt_av(alpha,x,y,beta,z,info)
      complex(psb_spk_), intent(in)        :: alpha,beta
      complex(psb_spk_), intent(in)        :: x(:)
      class(psb_c_vect_type), intent(inout)  :: y
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_mlt_av
  end interface
  
  interface 
    module subroutine c_vect_mlt_va(alpha,x,y,beta,z,info)
      complex(psb_spk_), intent(in)        :: alpha,beta
      complex(psb_spk_), intent(in)        :: y(:)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_mlt_va
  end interface
  
  interface 
    module subroutine c_vect_div_v(x, y, info)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_div_v
  end interface
  
  interface 
    module subroutine c_vect_div_v2( x, y, z, info)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_div_v2
  end interface
  
  interface 
    module subroutine c_vect_div_v_check(x, y, info, flag)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine c_vect_div_v_check
  end interface
  
  interface 
    module subroutine c_vect_div_v2_check(x, y, z, info, flag)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine c_vect_div_v2_check
  end interface
  
  interface 
    module subroutine c_vect_div_a2(x, y, z, info)
      complex(psb_spk_), intent(in) :: x(:)
      complex(psb_spk_), intent(in)    :: y(:)
      class(psb_c_vect_type), intent(inout) :: z
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_div_a2
  end interface
  
  interface 
    module subroutine c_vect_div_a2_check(x, y, z, info,flag)
      complex(psb_spk_), intent(in) :: x(:)
      complex(psb_spk_), intent(in)    :: y(:)
      class(psb_c_vect_type), intent(inout) :: z
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine c_vect_div_a2_check
  end interface

  interface 
    module subroutine c_vect_inv_v(x, y, info)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_inv_v
  end interface
  
  interface 
    module subroutine c_vect_inv_v_check(x, y, info, flag)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine c_vect_inv_v_check
  end interface
  
  interface 
    module subroutine c_vect_inv_a2(x, y, info)
      complex(psb_spk_), intent(inout) :: x(:)
      class(psb_c_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_vect_inv_a2
  end interface
  
  interface 
    module subroutine c_vect_inv_a2_check(x, y, info,flag)
      complex(psb_spk_), intent(inout) :: x(:)
      class(psb_c_vect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)              :: info
      logical, intent(in) :: flag
    end subroutine c_vect_inv_a2_check
  end interface
  
  interface 
    module subroutine c_vect_acmp_a2(x,c,z,info)
      real(psb_spk_), intent(in)              :: c
      complex(psb_spk_), intent(inout)           :: x(:)
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine c_vect_acmp_a2
  end interface
  
  interface 
    module subroutine c_vect_acmp_v2(x,c,z,info)
      real(psb_spk_), intent(in)              :: c
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine c_vect_acmp_v2
  end interface
  
  interface 
    module subroutine c_vect_scal(alpha, x)
      class(psb_c_vect_type), intent(inout)  :: x
      complex(psb_spk_), intent (in)       :: alpha
    end subroutine c_vect_scal
  end interface
  
  interface 
    module subroutine c_vect_absval1(x)
      class(psb_c_vect_type), intent(inout)  :: x
    end subroutine c_vect_absval1
  end interface
  
  interface 
    module subroutine c_vect_absval2(x,y)
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: y
    end subroutine c_vect_absval2
  end interface
  
  interface 
    module function c_vect_nrm2(n,x) result(res)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function c_vect_nrm2
  end interface
  
  interface 
    module function c_vect_nrm2_weight(n,x,w,aux) result(res)
      class(psb_c_vect_type), intent(inout) :: x
      class(psb_c_vect_type), intent(inout) :: w
      class(psb_c_vect_type), intent(inout), optional :: aux
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                        :: res
    end function c_vect_nrm2_weight
  end interface
  
  interface 
    module function c_vect_nrm2_weight_mask(n,x,w,id,info,aux) result(res)
      class(psb_c_vect_type), intent(inout) :: x
      class(psb_c_vect_type), intent(inout) :: w
      class(psb_c_vect_type), intent(inout) :: id
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                        :: res
      integer(psb_ipk_), intent(out)          :: info
      class(psb_c_vect_type), intent(inout), optional :: aux
    end function c_vect_nrm2_weight_mask
  end interface
  
  interface 
    module function c_vect_amax(n,x) result(res)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function c_vect_amax
  end interface
  

  interface 
    module function c_vect_asum(n,x) result(res)
      class(psb_c_vect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in)           :: n
      real(psb_spk_)                :: res
    end function c_vect_asum
  end interface
  

  interface 
    module subroutine c_vect_addconst_a2(x,b,z,info)
      real(psb_spk_), intent(in)             :: b
      complex(psb_spk_), intent(inout)           :: x(:)
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine c_vect_addconst_a2
  end interface
  
  interface 
    module subroutine c_vect_addconst_v2(x,b,z,info)
      real(psb_spk_), intent(in)             :: b
      class(psb_c_vect_type), intent(inout)  :: x
      class(psb_c_vect_type), intent(inout)  :: z
      integer(psb_ipk_), intent(out)           :: info
    end subroutine c_vect_addconst_v2
  end interface

contains
  
  subroutine  psb_c_set_vect_default(v)
    class(psb_c_base_vect_type), intent(in) :: v

    if (allocated(psb_c_base_vect_default)) then
      deallocate(psb_c_base_vect_default)
    end if
    allocate(psb_c_base_vect_default, mold=v)

  end subroutine psb_c_set_vect_default

  function psb_c_get_vect_default(v) result(res)
    class(psb_c_vect_type), intent(in) :: v
    class(psb_c_base_vect_type), pointer :: res

    res => psb_c_get_base_vect_default()

  end function psb_c_get_vect_default

  subroutine  psb_c_clear_vect_default()

    if (allocated(psb_c_base_vect_default)) then
      deallocate(psb_c_base_vect_default)
    end if

  end subroutine psb_c_clear_vect_default

  function psb_c_get_base_vect_default() result(res)
    class(psb_c_base_vect_type), pointer :: res

    if (.not.allocated(psb_c_base_vect_default)) then
      allocate(psb_c_base_vect_type :: psb_c_base_vect_default)
    end if

    res => psb_c_base_vect_default

  end function psb_c_get_base_vect_default

  function constructor(x) result(this)
    complex(psb_spk_)   :: x(:)
    type(psb_c_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(x)
    call this%asb(size(x,kind=psb_ipk_),info)

  end function constructor


  function size_const(n) result(this)
    integer(psb_ipk_), intent(in) :: n
    type(psb_c_vect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld(n)
    call this%asb(n,info)

  end function size_const

end module psb_c_vect_mod


module psb_c_multivect_mod

  use psb_c_base_multivect_mod
  use psb_const_mod
  use psb_i_vect_mod

  !private

  type psb_c_multivect_type
    class(psb_c_base_multivect_type), allocatable :: v
    integer(psb_ipk_) :: nrmv = 0
    integer(psb_ipk_) :: remote_build=psb_matbld_noremote_
    integer(psb_ipk_) :: dupl = psb_dupl_add_
    complex(psb_spk_), allocatable :: rmtv(:,:)
  contains
    procedure, pass(x) :: get_nrows => c_mvect_get_nrows
    procedure, pass(x) :: get_ncols => c_mvect_get_ncols
    procedure, pass(x) :: sizeof   => c_mvect_sizeof
    procedure, pass(x) :: get_fmt  => c_mvect_get_fmt
    procedure, pass(x) :: is_remote_build => c_mvect_is_remote_build
    procedure, pass(x) :: set_remote_build => c_mvect_set_remote_build
    procedure, pass(x) :: get_dupl => c_mvect_get_dupl
    procedure, pass(x) :: set_dupl => c_mvect_set_dupl

    procedure, pass(x) :: all      => c_mvect_all
    procedure, pass(x) :: reall    => c_mvect_reall
    procedure, pass(x) :: zero     => c_mvect_zero
    procedure, pass(x) :: asb      => c_mvect_asb
    procedure, pass(x) :: sync     => c_mvect_sync
    procedure, pass(x) :: free     => c_mvect_free
    procedure, pass(x) :: ins      => c_mvect_ins
    procedure, pass(x) :: bld_x    => c_mvect_bld_x
    procedure, pass(x) :: bld_n    => c_mvect_bld_n
    generic, public    :: bld      => bld_x, bld_n
    procedure, pass(x) :: get_vect => c_mvect_get_vect
    procedure, pass(x) :: cnv      => c_mvect_cnv
    procedure, pass(x) :: set_scal => c_mvect_set_scal
    procedure, pass(x) :: set_vect => c_mvect_set_vect
    generic, public    :: set      => set_vect, set_scal
    procedure, pass(x) :: clone    => c_mvect_clone
    procedure, pass(x) :: gthab    => c_mvect_gthab
    procedure, pass(x) :: gthzv    => c_mvect_gthzv
    procedure, pass(x) :: gthzv_x  => c_mvect_gthzv_x
    generic, public    :: gth      => gthab, gthzv
    procedure, pass(y) :: sctb     => c_mvect_sctb
    procedure, pass(y) :: sctb_x   => c_mvect_sctb_x
    generic, public    :: sct      => sctb, sctb_x
!!$    procedure, pass(x) :: dot_v    => c_mvect_dot_v
!!$    procedure, pass(x) :: dot_a    => c_mvect_dot_a
!!$    generic, public    :: dot      => dot_v, dot_a
!!$    procedure, pass(y) :: axpby_v  => c_mvect_axpby_v
!!$    procedure, pass(y) :: axpby_a  => c_mvect_axpby_a
!!$    generic, public    :: axpby    => axpby_v, axpby_a
!!$    procedure, pass(y) :: mlt_v    => c_mvect_mlt_v
!!$    procedure, pass(y) :: mlt_a    => c_mvect_mlt_a
!!$    procedure, pass(z) :: mlt_a_2  => c_mvect_mlt_a_2
!!$    procedure, pass(z) :: mlt_v_2  => c_mvect_mlt_v_2
!!$    procedure, pass(z) :: mlt_va   => c_mvect_mlt_va
!!$    procedure, pass(z) :: mlt_av   => c_mvect_mlt_av
!!$    generic, public    :: mlt      => mlt_v, mlt_a, mlt_a_2,&
!!$         & mlt_v_2, mlt_av, mlt_va
!!$    procedure, pass(x) :: scal     => c_mvect_scal
!!$    procedure, pass(x) :: nrm2     => c_mvect_nrm2
!!$    procedure, pass(x) :: amax     => c_mvect_amax
!!$    procedure, pass(x) :: asum     => c_mvect_asum
  end type psb_c_multivect_type

  public  :: psb_c_multivect, psb_c_multivect_type,&
       & psb_c_set_multivect_default, psb_c_get_base_multivect_default, &
       & psb_c_clear_multivect_default, psb_c_base_multivect_type

  interface psb_c_multivect
    module procedure constructor, size_const
  end interface psb_c_multivect

  private

  class(psb_c_base_multivect_type), allocatable, target,&
       & save, private :: psb_c_base_multivect_default

  interface 
    module function c_mvect_get_dupl(x) result(res)
      class(psb_c_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function c_mvect_get_dupl
  end interface

  interface 
    module subroutine c_mvect_set_dupl(x,val)
      class(psb_c_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine c_mvect_set_dupl
  end interface

  interface 
    module function c_mvect_is_remote_build(x) result(res)
      class(psb_c_multivect_type), intent(in) :: x
      logical :: res
    end function c_mvect_is_remote_build
  end interface

  interface 
    module subroutine c_mvect_set_remote_build(x,val)
      class(psb_c_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(in), optional :: val
    end subroutine c_mvect_set_remote_build
  end interface

  interface 
    module subroutine c_mvect_clone(x,y,info)
      class(psb_c_multivect_type), intent(inout) :: x
      class(psb_c_multivect_type), intent(inout) :: y
      integer(psb_ipk_), intent(out)        :: info
    end subroutine c_mvect_clone
  end interface

  interface 
    module subroutine c_mvect_bld_x(x,invect,mold)
      complex(psb_spk_), intent(in)          :: invect(:,:)
      class(psb_c_multivect_type), intent(out) :: x
      class(psb_c_base_multivect_type), intent(in), optional :: mold
    end subroutine c_mvect_bld_x
  end interface

  interface 
    module subroutine c_mvect_bld_n(x,m,n,mold,scratch)
      integer(psb_ipk_), intent(in) :: m,n
      class(psb_c_multivect_type), intent(out) :: x
      class(psb_c_base_multivect_type), intent(in), optional :: mold
      logical, intent(in), optional        :: scratch
    end subroutine c_mvect_bld_n
  end interface

  interface
    module function  c_mvect_get_vect(x) result(res)
      class(psb_c_multivect_type), intent(inout)  :: x
      complex(psb_spk_), allocatable                 :: res(:,:)
    end function c_mvect_get_vect
  end interface
  
  interface 
    module subroutine c_mvect_set_scal(x,val)
      class(psb_c_multivect_type), intent(inout)  :: x
      complex(psb_spk_), intent(in) :: val
    end subroutine c_mvect_set_scal
  end interface
  
  interface 
    module subroutine c_mvect_set_vect(x,val)
      class(psb_c_multivect_type), intent(inout) :: x
      complex(psb_spk_), intent(in)         :: val(:,:)
    end subroutine c_mvect_set_vect
  end interface
  
  interface 
    module function c_mvect_get_nrows(x) result(res)
      class(psb_c_multivect_type), intent(in) :: x
      integer(psb_ipk_)  :: res
    end function c_mvect_get_nrows
  end interface
  
  interface 
    module function c_mvect_get_ncols(x) result(res)
      class(psb_c_multivect_type), intent(in) :: x
      integer(psb_ipk_) :: res
    end function c_mvect_get_ncols
  end interface
  
  interface 
    module function c_mvect_sizeof(x) result(res)
      class(psb_c_multivect_type), intent(in) :: x
      integer(psb_epk_) :: res
    end function c_mvect_sizeof
  end interface
  
  interface 
    module function c_mvect_get_fmt(x) result(res)
      class(psb_c_multivect_type), intent(in) :: x
      character(len=5) :: res
    end function c_mvect_get_fmt
  end interface
  
  interface 
    module subroutine c_mvect_all(m,n, x, info, mold)
      integer(psb_ipk_), intent(in)       :: m,n
      class(psb_c_multivect_type), intent(out) :: x
      class(psb_c_base_multivect_type), intent(in), optional :: mold
      integer(psb_ipk_), intent(out)      :: info
    end subroutine c_mvect_all
  end interface
  
  interface 
    module subroutine c_mvect_reall(m,n, x, info)
      integer(psb_ipk_), intent(in)         :: m,n
      class(psb_c_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)        :: info
    end subroutine c_mvect_reall
  end interface
  
  interface 
    module subroutine c_mvect_zero(x)
      class(psb_c_multivect_type), intent(inout)    :: x
    end subroutine c_mvect_zero
  end interface
  
  interface 
    module subroutine c_mvect_asb(m,n, x, info)
      integer(psb_ipk_), intent(in)              :: m,n
      class(psb_c_multivect_type), intent(inout) :: x
      integer(psb_ipk_), intent(out)             :: info
    end subroutine c_mvect_asb
  end interface
  
  interface 
    module subroutine c_mvect_sync(x)
      class(psb_c_multivect_type), intent(inout) :: x
    end subroutine c_mvect_sync
  end interface
  
  interface 
    module subroutine c_mvect_gthab(n,idx,alpha,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_spk_) :: alpha, beta, y(:)
      class(psb_c_multivect_type) :: x
    end subroutine c_mvect_gthab
  end interface
  
  interface 
    module subroutine c_mvect_gthzv(n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_spk_) ::  y(:)
      class(psb_c_multivect_type) :: x
    end subroutine c_mvect_gthzv
  end interface
  
  interface 
    module subroutine c_mvect_gthzv_x(i,n,idx,x,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      complex(psb_spk_) ::  y(:)
      class(psb_c_multivect_type) :: x
    end subroutine c_mvect_gthzv_x
  end interface
  
  interface 
    module subroutine c_mvect_sctb(n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: idx(:)
      complex(psb_spk_) :: beta, x(:)
      class(psb_c_multivect_type) :: y
    end subroutine c_mvect_sctb
  end interface

  interface 
    module subroutine c_mvect_sctb_x(i,n,idx,x,beta,y)
      integer(psb_mpk_) :: n
      integer(psb_ipk_) :: i
      class(psb_i_base_vect_type) :: idx
      complex(psb_spk_) :: beta, x(:)
      class(psb_c_multivect_type) :: y
    end subroutine c_mvect_sctb_x
  end interface
  
  interface 
    module subroutine c_mvect_free(x, info)
      class(psb_c_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_mvect_free
  end interface
  
  interface 
    module subroutine c_mvect_ins(n,irl,val,x,maxr,info)
      class(psb_c_multivect_type), intent(inout)  :: x
      integer(psb_ipk_), intent(in)               :: n,maxr
      integer(psb_ipk_), intent(in)               :: irl(:)
      complex(psb_spk_), intent(in)        :: val(:,:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine c_mvect_ins
  end interface
  
  interface 
    module subroutine c_mvect_cnv(x,mold)
      class(psb_c_multivect_type), intent(inout) :: x
      class(psb_c_base_multivect_type), intent(in), optional :: mold
      integer(psb_ipk_) :: info
    end subroutine c_mvect_cnv
  end interface


!!$  module function c_mvect_dot_v(n,x,y) result(res)
!!$    class(psb_c_multivect_type), intent(inout) :: x, y
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    complex(psb_spk_)                :: res
!!$
!!$    res = czero
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & res = x%v%dot(n,y%v)
!!$
!!$  end function c_mvect_dot_v
!!$
!!$  function c_mvect_dot_a(n,x,y) result(res)
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    complex(psb_spk_), intent(in)    :: y(:)
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    complex(psb_spk_)                :: res
!!$
!!$    res = czero
!!$    if (allocated(x%v)) &
!!$         & res = x%v%dot(n,y)
!!$
!!$  end function c_mvect_dot_a
!!$
!!$  module subroutine c_mvect_axpby_v(m,alpha, x, beta, y, info)
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    complex(psb_spk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(x%v).and.allocated(y%v)) then
!!$      call y%v%axpby(m,alpha,x%v,beta,info)
!!$    else
!!$      info = psb_err_invalid_mvect_state_
!!$    end if
!!$
!!$  end subroutine c_mvect_axpby_v
!!$
!!$  subroutine c_mvect_axpby_a(m,alpha, x, beta, y, info)
!!$    integer(psb_ipk_), intent(in)               :: m
!!$    complex(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    complex(psb_spk_), intent (in)       :: alpha, beta
!!$    integer(psb_ipk_), intent(out)              :: info
!!$
!!$    if (allocated(y%v)) &
!!$         & call y%v%axpby(m,alpha,x,beta,info)
!!$
!!$  end subroutine c_mvect_axpby_a
!!$
!!$
!!$  subroutine c_mvect_mlt_v(x, y, info)
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(x%v).and.allocated(y%v)) &
!!$         & call y%v%mlt(x%v,info)
!!$
!!$  end subroutine c_mvect_mlt_v
!!$
!!$  subroutine c_mvect_mlt_a(x, y, info)
!!$    complex(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$
!!$    info = 0
!!$    if (allocated(y%v)) &
!!$         & call y%v%mlt(x,info)
!!$
!!$  end subroutine c_mvect_mlt_a
!!$
!!$
!!$  subroutine c_mvect_mlt_a_2(alpha,x,y,beta,z,info)
!!$    complex(psb_spk_), intent(in)         :: alpha,beta
!!$    complex(psb_spk_), intent(in)         :: y(:)
!!$    complex(psb_spk_), intent(in)         :: x(:)
!!$    class(psb_c_multivect_type), intent(inout) :: z
!!$    integer(psb_ipk_), intent(out)                  :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v)) &
!!$         & call z%v%mlt(alpha,x,y,beta,info)
!!$
!!$  end subroutine c_mvect_mlt_a_2
!!$
!!$  subroutine c_mvect_mlt_v_2(alpha,x,y,beta,z,info,conjgx,conjgy)
!!$    complex(psb_spk_), intent(in)          :: alpha,beta
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    class(psb_c_multivect_type), intent(inout)  :: z
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
!!$  end subroutine c_mvect_mlt_v_2
!!$
!!$  subroutine c_mvect_mlt_av(alpha,x,y,beta,z,info)
!!$    complex(psb_spk_), intent(in)        :: alpha,beta
!!$    complex(psb_spk_), intent(in)        :: x(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: y
!!$    class(psb_c_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$    if (allocated(z%v).and.allocated(y%v)) &
!!$         & call z%v%mlt(alpha,x,y%v,beta,info)
!!$
!!$  end subroutine c_mvect_mlt_av
!!$
!!$  subroutine c_mvect_mlt_va(alpha,x,y,beta,z,info)
!!$    complex(psb_spk_), intent(in)        :: alpha,beta
!!$    complex(psb_spk_), intent(in)        :: y(:)
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    class(psb_c_multivect_type), intent(inout)  :: z
!!$    integer(psb_ipk_), intent(out)              :: info
!!$    integer(psb_ipk_) :: i, n
!!$
!!$    info = 0
!!$
!!$    if (allocated(z%v).and.allocated(x%v)) &
!!$         & call z%v%mlt(alpha,x%v,y,beta,info)
!!$
!!$  end subroutine c_mvect_mlt_va
!!$
!!$  subroutine c_mvect_scal(alpha, x)
!!$    class(psb_c_multivect_type), intent(inout)  :: x
!!$    complex(psb_spk_), intent (in)       :: alpha
!!$
!!$    if (allocated(x%v)) call x%v%scal(alpha)
!!$
!!$  end subroutine c_mvect_scal
!!$
!!$
!!$  function c_mvect_nrm2(n,x) result(res)
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%nrm2(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function c_mvect_nrm2
!!$
!!$  function c_mvect_amax(n,x) result(res)
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%amax(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function c_mvect_amax
!!$
!!$  function c_mvect_asum(n,x) result(res)
!!$    class(psb_c_multivect_type), intent(inout) :: x
!!$    integer(psb_ipk_), intent(in)           :: n
!!$    real(psb_spk_)                :: res
!!$
!!$    if (allocated(x%v)) then
!!$      res = x%v%asum(n)
!!$    else
!!$      res = szero
!!$    end if
!!$
!!$  end function c_mvect_asum

contains

  subroutine  psb_c_set_multivect_default(v)
    class(psb_c_base_multivect_type), intent(in) :: v

    if (allocated(psb_c_base_multivect_default)) then
      deallocate(psb_c_base_multivect_default)
    end if
    allocate(psb_c_base_multivect_default, mold=v)

  end subroutine psb_c_set_multivect_default

!!$  function psb_c_get_multivect_default(v) result(res)
!!$    class(psb_c_multivect_type), intent(in) :: v
!!$    class(psb_c_base_multivect_type), pointer :: res
!!$
!!$    res => psb_c_get_base_multivect_default()
!!$
!!$  end function psb_c_get_multivect_default
!!$

  function psb_c_get_base_multivect_default() result(res)
    class(psb_c_base_multivect_type), pointer :: res

    if (.not.allocated(psb_c_base_multivect_default)) then
      allocate(psb_c_base_multivect_type :: psb_c_base_multivect_default)
    end if

    res => psb_c_base_multivect_default

  end function psb_c_get_base_multivect_default

  function constructor(x) result(this)
    complex(psb_spk_)   :: x(:,:)
    type(psb_c_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld_x(x)
    call this%asb(size(x,dim=1,kind=psb_ipk_),size(x,dim=2,kind=psb_ipk_),info)

  end function constructor

  function size_const(m,n) result(this)
    integer(psb_ipk_), intent(in) :: m,n
    type(psb_c_multivect_type) :: this
    integer(psb_ipk_) :: info

    call this%bld_n(m,n)
    call this%asb(m,n,info)

  end function size_const

  
  subroutine  psb_c_clear_multivect_default()

    if (allocated(psb_c_base_multivect_default)) then
      deallocate(psb_c_base_multivect_default)
    end if

  end subroutine psb_c_clear_multivect_default

end module psb_c_multivect_mod
