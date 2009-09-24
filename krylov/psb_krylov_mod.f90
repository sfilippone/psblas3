!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
! File: psb_krylov_mod.f90
!  Interfaces for Krylov subspace iterative methods.
!
Module psb_krylov_mod
  
  use psb_const_mod
  public 

  interface psb_krylov
    module procedure psb_skrylov, psb_dkrylov, psb_ckrylov, psb_zkrylov
  end interface

  interface psb_cg
    subroutine psb_scg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop,cond)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_mod, only : psb_sprec_type
      type(psb_s_sparse_mat), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_spk_), intent(in)       :: b(:)
      real(psb_spk_), intent(inout)    :: x(:)
      real(psb_spk_), intent(in)       :: eps
      type(psb_sprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_spk_), optional, intent(out) :: err,cond
    end subroutine psb_scg
    subroutine psb_dcg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop,cond)
      use psb_base_mod, only  : psb_desc_type, psb_d_sparse_mat, psb_dpk_
      use psb_prec_mod, only : psb_dprec_type
      type(psb_d_sparse_mat), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_dpk_), intent(in)       :: b(:)
      real(psb_dpk_), intent(inout)    :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      type(psb_dprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err,cond
    end subroutine psb_dcg
    subroutine psb_ccg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_cspmat_type, psb_spk_
      use psb_prec_mod, only : psb_cprec_type
      type(psb_cspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_spk_), intent(in)    :: b(:)
      complex(psb_spk_), intent(inout) :: x(:)
      real(psb_spk_), intent(in)       :: eps
      type(psb_cprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_spk_), optional, intent(out) :: err
    end subroutine psb_ccg
    subroutine psb_zcg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_dpk_), intent(in)    :: b(:)
      complex(psb_dpk_), intent(inout) :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      type(psb_zprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_zcg
  end interface

  interface psb_bicg
    subroutine psb_sbicg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_mod, only : psb_sprec_type
      type(psb_s_sparse_mat), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_spk_), intent(in)       :: b(:)
      real(psb_spk_), intent(inout)    :: x(:)
      real(psb_spk_), intent(in)       :: eps
      type(psb_sprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_spk_), optional, intent(out) :: err
    end subroutine psb_sbicg
    subroutine psb_dbicg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_d_sparse_mat, psb_dpk_
      use psb_prec_mod, only : psb_dprec_type
      type(psb_d_sparse_mat), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_dpk_), intent(in)         :: b(:)
      real(psb_dpk_), intent(inout)      :: x(:)
      real(psb_dpk_), intent(in)         :: eps
      type(psb_dprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_dbicg
    subroutine psb_cbicg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_cspmat_type, psb_spk_
      use psb_prec_mod, only : psb_cprec_type
      type(psb_cspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_spk_), intent(in)      :: b(:)
      complex(psb_spk_), intent(inout)   :: x(:)
      real(psb_spk_), intent(in)         :: eps
      type(psb_cprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_spk_), optional, intent(out) :: err
    end subroutine psb_cbicg
    subroutine psb_zbicg(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_dpk_), intent(in)      :: b(:)
      complex(psb_dpk_), intent(inout)   :: x(:)
      real(psb_dpk_), intent(in)         :: eps
      type(psb_zprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_zbicg
  end interface

  interface psb_bicgstab
    subroutine psb_scgstab(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_mod, only : psb_sprec_type
      type(psb_s_sparse_mat), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_spk_), intent(in)       :: b(:)
      real(psb_spk_), intent(inout)    :: x(:)
      real(psb_spk_), intent(in)       :: eps
      type(psb_sprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_spk_), optional, intent(out) :: err
    end subroutine psb_scgstab
    subroutine psb_dcgstab(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_d_sparse_mat, psb_dpk_
      use psb_prec_mod, only : psb_dprec_type
      type(psb_d_sparse_mat), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      real(psb_dpk_), intent(in)       :: b(:)
      real(psb_dpk_), intent(inout)    :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      type(psb_dprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_dcgstab
    subroutine psb_ccgstab(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_cspmat_type, psb_spk_
      use psb_prec_mod, only : psb_cprec_type
      type(psb_cspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_spk_), intent(in)       :: b(:)
      complex(psb_spk_), intent(inout)    :: x(:)
      real(psb_spk_), intent(in)       :: eps
      type(psb_cprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_spk_), optional, intent(out) :: err
    end subroutine psb_ccgstab
    subroutine psb_zcgstab(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_dpk_), intent(in)       :: b(:)
      complex(psb_dpk_), intent(inout)    :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      type(psb_zprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_zcgstab
  end interface

  interface psb_bicgstabl
    Subroutine psb_scgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err, itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_mod, only : psb_sprec_type
      Type(psb_s_sparse_mat), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_sprec_type), intent(in)   :: prec
      Real(psb_spk_), Intent(in)       :: b(:)
      Real(psb_spk_), Intent(inout)    :: x(:)
      Real(psb_spk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_spk_), Optional, Intent(out) :: err
    end subroutine psb_scgstabl
    Subroutine psb_dcgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err, itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_d_sparse_mat, psb_dpk_
      use psb_prec_mod, only : psb_dprec_type
      type(psb_d_sparse_mat), intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_dprec_type), intent(in)   :: prec
      Real(psb_dpk_), Intent(in)       :: b(:)
      Real(psb_dpk_), Intent(inout)    :: x(:)
      Real(psb_dpk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err
    end subroutine psb_dcgstabl
    Subroutine psb_ccgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_cspmat_type, psb_spk_
      use psb_prec_mod, only : psb_cprec_type
      Type(psb_cspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_cprec_type), intent(in)   :: prec 
      complex(psb_spk_), Intent(in)    :: b(:)
      complex(psb_spk_), Intent(inout) :: x(:)
      Real(psb_spk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_spk_), Optional, Intent(out) :: err
    end subroutine psb_ccgstabl
    Subroutine psb_zcgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      Type(psb_zspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_zprec_type), intent(in)   :: prec 
      complex(psb_dpk_), Intent(in)    :: b(:)
      complex(psb_dpk_), Intent(inout) :: x(:)
      Real(psb_dpk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err
    end subroutine psb_zcgstabl
  end interface

  interface psb_rgmres
    Subroutine psb_srgmres(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_mod, only : psb_sprec_type
      Type(psb_s_sparse_mat), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_sprec_type), intent(in)   :: prec 
      Real(psb_spk_), Intent(in)       :: b(:)
      Real(psb_spk_), Intent(inout)    :: x(:)
      Real(psb_spk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_spk_), Optional, Intent(out) :: err
    end subroutine psb_srgmres
    Subroutine psb_drgmres(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_d_sparse_mat, psb_dpk_
      use psb_prec_mod, only : psb_dprec_type
      type(psb_d_sparse_mat), intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_dprec_type), intent(in)   :: prec 
      Real(psb_dpk_), Intent(in)       :: b(:)
      Real(psb_dpk_), Intent(inout)    :: x(:)
      Real(psb_dpk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err
    end subroutine psb_drgmres
    Subroutine psb_crgmres(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_cspmat_type, psb_spk_
      use psb_prec_mod, only : psb_cprec_type
      Type(psb_cspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_cprec_type), intent(in)   :: prec 
      complex(psb_spk_), Intent(in)    :: b(:)
      complex(psb_spk_), Intent(inout) :: x(:)
      Real(psb_spk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_spk_), Optional, Intent(out) :: err
    end subroutine psb_crgmres
    Subroutine psb_zrgmres(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      Type(psb_zspmat_type), Intent(in)  :: a
      Type(psb_desc_type), Intent(in)    :: desc_a
      type(psb_zprec_type), intent(in)   :: prec 
      complex(psb_dpk_), Intent(in)    :: b(:)
      complex(psb_dpk_), Intent(inout) :: x(:)
      Real(psb_dpk_), Intent(in)       :: eps
      integer, intent(out)               :: info
      Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
      Integer, Optional, Intent(out)     :: iter
      Real(psb_dpk_), Optional, Intent(out) :: err
    end subroutine psb_zrgmres
  end interface

  interface psb_cgs
    subroutine psb_scgs(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_mod, only : psb_sprec_type
      type(psb_s_sparse_mat), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a 
      type(psb_sprec_type), intent(in)   :: prec 
      real(psb_spk_), intent(in)       :: b(:)
      real(psb_spk_), intent(inout)    :: x(:)
      real(psb_spk_), intent(in)       :: eps
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_spk_), optional, intent(out) :: err
    end subroutine psb_scgs
    subroutine psb_dcgs(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_d_sparse_mat, psb_dpk_
      use psb_prec_mod, only : psb_dprec_type
      type(psb_d_sparse_mat), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a 
      type(psb_dprec_type), intent(in)   :: prec 
      real(psb_dpk_), intent(in)       :: b(:)
      real(psb_dpk_), intent(inout)    :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_dcgs
    subroutine psb_ccgs(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_cspmat_type, psb_spk_
      use psb_prec_mod, only : psb_cprec_type
      type(psb_cspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_spk_), intent(in)       :: b(:)
      complex(psb_spk_), intent(inout)    :: x(:)
      real(psb_spk_), intent(in)       :: eps
      type(psb_cprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_spk_), optional, intent(out) :: err
    end subroutine psb_ccgs
    subroutine psb_zcgs(a,prec,b,x,eps,&
         & desc_a,info,itmax,iter,err,itrace,istop)
      use psb_base_mod, only  : psb_desc_type, psb_zspmat_type, psb_dpk_
      use psb_prec_mod, only : psb_zprec_type
      type(psb_zspmat_type), intent(in)  :: a
      type(psb_desc_type), intent(in)    :: desc_a
      complex(psb_dpk_), intent(in)       :: b(:)
      complex(psb_dpk_), intent(inout)    :: x(:)
      real(psb_dpk_), intent(in)       :: eps
      type(psb_zprec_type), intent(in)   :: prec
      integer, intent(out)               :: info
      integer, optional, intent(in)      :: itmax, itrace,istop
      integer, optional, intent(out)     :: iter
      real(psb_dpk_), optional, intent(out) :: err
    end subroutine psb_zcgs
  end interface



  interface psb_init_conv
    module procedure psb_s_init_conv, psb_d_init_conv,&
         & psb_c_init_conv, psb_z_init_conv
  end interface

  interface psb_check_conv
    module procedure psb_s_check_conv, psb_d_check_conv,&
         & psb_c_check_conv, psb_z_check_conv
  end interface

  interface psb_end_conv
    module procedure psb_d_end_conv
  end interface

  integer, parameter, private :: bni_=1, rni_=2, ani_=3, xni_=4, bn2_=5, xn2_=6
  integer, parameter, private :: errnum_=7, errden_=8, eps_=9, rn2_=10
  integer, parameter, private :: stopc_=1, trace_=2, itmax_=3
  integer, parameter, private :: ivsz_=16
  type psb_itconv_type
    private
    integer        :: controls(ivsz_)
    real(psb_dpk_) :: values(ivsz_)
  end type psb_itconv_type

contains

  ! Subroutine: psb_skrylov
  ! 
  !    Front-end for the Krylov subspace iterations, real version
  !    
  ! Arguments:
  !
  !    methd  -  character                    The specific method; can take the values:
  !                                           CG
  !                                           CGS
  !                                           BICG
  !                                           BICGSTAB
  !                                           BICGSTABL
  !                                           RGMRES
  !                                           
  !    a      -  type(psb_s_sparse_mat)      Input: sparse matrix containing A.
  !    prec   -  type(psb_sprec_type)       Input: preconditioner
  !    b      -  real,dimension(:)            Input: vector containing the
  !                                           right hand side B
  !    x      -  real,dimension(:)            Input/Output: vector containing the
  !                                           initial guess and final solution X.
  !    eps    -  real                         Input: Stopping tolerance; the iteration is
  !                                           stopped when the error
  !                                           estimate  |err| <= eps
  !    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
  !    info   -  integer.                     Output: Return code
  !
  !    itmax  -  integer(optional)            Input: maximum number of iterations to be
  !                                           performed.
  !    iter   -  integer(optional)            Output: how many iterations have been
  !                                           performed.
  !    err    -  real   (optional)            Output: error estimate on exit
  !    itrace -  integer(optional)            Input: print an informational message
  !                                           with the error estimate every itrace
  !                                           iterations
  !    irst   -  integer(optional)            Input: restart parameter for RGMRES and 
  !                                           BICGSTAB(L) methods
  !    istop  -  integer(optional)            Input: stopping criterion, or how
  !                                           to estimate the error. 
  !                                           1: err =  |r|/|b|
  !                                           2: err =  |r|/(|a||x|+|b|)
  !                                           where r is the (preconditioned, recursive
  !                                           estimate of) residual 
  ! 

  Subroutine psb_skrylov(method,a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop,cond)

    use psb_base_mod
    use psb_prec_mod,only : psb_sprec_type, psb_dprec_type, psb_cprec_type, psb_zprec_type

    character(len=*)                   :: method
    Type(psb_s_sparse_mat), Intent(in)  :: a
    Type(psb_desc_type), Intent(in)    :: desc_a
    type(psb_sprec_type), intent(in)   :: prec 
    Real(psb_spk_), Intent(in)       :: b(:)
    Real(psb_spk_), Intent(inout)    :: x(:)
    Real(psb_spk_), Intent(in)       :: eps
    integer, intent(out)               :: info
    Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
    Integer, Optional, Intent(out)     :: iter
    Real(psb_spk_), Optional, Intent(out) :: err,cond

    integer                            :: ictxt,me,np,err_act
    character(len=20)             :: name

    info = 0
    name = 'psb_krylov'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)

    select case(psb_toupper(method))
    case('CG') 
      call  psb_cg(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop,cond)
    case('CGS') 
      call  psb_cgs(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    case('BICG') 
      call  psb_bicg(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    case('BICGSTAB') 
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    case('RGMRES')
      call  psb_rgmres(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,irst,istop)
    case('BICGSTABL')
      call  psb_bicgstabl(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,irst,istop)
    case default
      if (me==0) write(0,*) trim(name),': Warning: Unknown method  ',method,&
           & ', defaulting to BiCGSTAB'
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    end select

    if(info/=0) then
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if

  end subroutine psb_skrylov


  ! Subroutine: psb_dkrylov
  ! 
  !    Front-end for the Krylov subspace iterations, real version
  !    
  ! Arguments:
  !
  !    methd  -  character                    The specific method; can take the values:
  !                                           CG
  !                                           CGS
  !                                           BICG
  !                                           BICGSTAB
  !                                           BICGSTABL
  !                                           RGMRES
  !                                           
  !    a      -  type(psb_d_sparse_mat)      Input: sparse matrix containing A.
  !    prec   -  type(psb_dprec_type)       Input: preconditioner
  !    b      -  real,dimension(:)            Input: vector containing the
  !                                           right hand side B
  !    x      -  real,dimension(:)            Input/Output: vector containing the
  !                                           initial guess and final solution X.
  !    eps    -  real                         Input: Stopping tolerance; the iteration is
  !                                           stopped when the error
  !                                           estimate  |err| <= eps
  !    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
  !    info   -  integer.                     Output: Return code
  !
  !    itmax  -  integer(optional)            Input: maximum number of iterations to be
  !                                           performed.
  !    iter   -  integer(optional)            Output: how many iterations have been
  !                                           performed.
  !    err    -  real   (optional)            Output: error estimate on exit
  !    itrace -  integer(optional)            Input: print an informational message
  !                                           with the error estimate every itrace
  !                                           iterations
  !    irst   -  integer(optional)            Input: restart parameter for RGMRES and 
  !                                           BICGSTAB(L) methods
  !    istop  -  integer(optional)            Input: stopping criterion, or how
  !                                           to estimate the error. 
  !                                           1: err =  |r|/|b|
  !                                           2: err =  |r|/(|a||x|+|b|)
  !                                           where r is the (preconditioned, recursive
  !                                           estimate of) residual 
  ! 

  Subroutine psb_dkrylov(method,a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop,cond)

    use psb_base_mod
    use psb_prec_mod,only : psb_sprec_type, psb_dprec_type, psb_cprec_type, psb_zprec_type

    type(psb_d_sparse_mat), intent(in)  :: a

    character(len=*)                   :: method
    Type(psb_desc_type), Intent(in)    :: desc_a
    type(psb_dprec_type), intent(in)   :: prec 
    Real(psb_dpk_), Intent(in)       :: b(:)
    Real(psb_dpk_), Intent(inout)    :: x(:)
    Real(psb_dpk_), Intent(in)       :: eps
    integer, intent(out)               :: info
    Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
    Integer, Optional, Intent(out)     :: iter
    Real(psb_dpk_), Optional, Intent(out) :: err,cond

    integer                            :: ictxt,me,np,err_act
    character(len=20)             :: name

    info = 0
    name = 'psb_krylov'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)

    select case(psb_toupper(method))
    case('CG') 
      call  psb_cg(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop,cond)
    case('CGS') 
      call  psb_cgs(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    case('BICG') 
      call  psb_bicg(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    case('BICGSTAB') 
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    case('RGMRES')
      call  psb_rgmres(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,irst,istop)
    case('BICGSTABL')
      call  psb_bicgstabl(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,irst,istop)
    case default
      if (me==0) write(0,*) trim(name),': Warning: Unknown method  ',method,&
           & ', defaulting to BiCGSTAB'
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    end select

    if(info/=0) then
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if

  end subroutine psb_dkrylov


  !
  ! Subroutine: psb_ckrylov
  ! 
  !    Front-end for the Krylov subspace iterations, complexversion
  !    
  ! Arguments:
  !
  !    methd  -  character                    The specific method; can take the values:
  !                                           CG
  !                                           CGS
  !                                           BICG
  !                                           BICGSTAB
  !                                           BICGSTABL
  !                                           RGMRES
  !                                           
  !    a      -  type(psb_cspmat_type)      Input: sparse matrix containing A.
  !    prec   -  type(psb_cprec_type)       Input: preconditioner
  !    b      -  complex,dimension(:)         Input: vector containing the
  !                                           right hand side B
  !    x      -  complex,dimension(:)         Input/Output: vector containing the
  !                                           initial guess and final solution X.
  !    eps    -  real                         Input: Stopping tolerance; the iteration is
  !                                           stopped when the error
  !                                           estimate |err| <= eps
  !                                           
  !    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
  !    info   -  integer.                     Output: Return code
  !
  !    itmax  -  integer(optional)            Input: maximum number of iterations to be
  !                                           performed.
  !    iter   -  integer(optional)            Output: how many iterations have been
  !                                           performed.
  !    err    -  real   (optional)            Output: error estimate on exit
  !    itrace -  integer(optional)            Input: print an informational message
  !                                           with the error estimate every itrace
  !                                           iterations
  !    irst   -  integer(optional)            Input: restart parameter for RGMRES and 
  !                                           BICGSTAB(L) methods
  !    istop  -  integer(optional)            Input: stopping criterion, or how
  !                                           to estimate the error. 
  !                                           1: err =  |r|/|b|
  !                                           2: err =  |r|/(|a||x|+|b|)
  !                                           where r is the (preconditioned, recursive
  !                                           estimate of) residual 
  ! 
  Subroutine psb_ckrylov(method,a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop)
    use psb_base_mod
    use psb_prec_mod,only : psb_sprec_type, psb_dprec_type, psb_cprec_type, psb_zprec_type
    character(len=*)                   :: method
    Type(psb_cspmat_type), Intent(in)  :: a
    Type(psb_desc_type), Intent(in)    :: desc_a
    type(psb_cprec_type), intent(in)   :: prec 
    complex(psb_spk_), Intent(in)      :: b(:)
    complex(psb_spk_), Intent(inout)   :: x(:)
    Real(psb_spk_), Intent(in)         :: eps
    integer, intent(out)               :: info
    Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
    Integer, Optional, Intent(out)     :: iter
    Real(psb_spk_), Optional, Intent(out) :: err

    integer                            :: ictxt,me,np,err_act
    character(len=20)             :: name

    info = 0
    name = 'psb_krylov'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)


    select case(psb_toupper(method))
    case('CG') 
      call  psb_cg(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    case('CGS') 
      call  psb_cgs(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    case('BICG') 
      call  psb_bicg(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    case('BICGSTAB') 
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           & itmax,iter,err,itrace,istop)
    case('RGMRES')
      call  psb_rgmres(a,prec,b,x,eps,desc_a,info,&
           & itmax,iter,err,itrace,irst,istop)
    case('BICGSTABL')
      call  psb_bicgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
    case default
      if (me==0) write(0,*) trim(name),': Warning: Unknown method  ',method,&
           & ', defaulting to BiCGSTAB'
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    end select

    if(info/=0) then
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if

  end subroutine psb_ckrylov

  !
  ! Subroutine: psb_zkrylov
  ! 
  !    Front-end for the Krylov subspace iterations, complexversion
  !    
  ! Arguments:
  !
  !    methd  -  character                    The specific method; can take the values:
  !                                           CG
  !                                           CGS
  !                                           BICG
  !                                           BICGSTAB
  !                                           BICGSTABL
  !                                           RGMRES
  !                                           
  !    a      -  type(psb_zspmat_type)      Input: sparse matrix containing A.
  !    prec   -  type(psb_zprec_type)       Input: preconditioner
  !    b      -  complex,dimension(:)         Input: vector containing the
  !                                           right hand side B
  !    x      -  complex,dimension(:)         Input/Output: vector containing the
  !                                           initial guess and final solution X.
  !    eps    -  real                         Input: Stopping tolerance; the iteration is
  !                                           stopped when the error
  !                                           estimate |err| <= eps
  !                                           
  !    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
  !    info   -  integer.                     Output: Return code
  !
  !    itmax  -  integer(optional)            Input: maximum number of iterations to be
  !                                           performed.
  !    iter   -  integer(optional)            Output: how many iterations have been
  !                                           performed.
  !    err    -  real   (optional)            Output: error estimate on exit
  !    itrace -  integer(optional)            Input: print an informational message
  !                                           with the error estimate every itrace
  !                                           iterations
  !    irst   -  integer(optional)            Input: restart parameter for RGMRES and 
  !                                           BICGSTAB(L) methods
  !    istop  -  integer(optional)            Input: stopping criterion, or how
  !                                           to estimate the error. 
  !                                           1: err =  |r|/|b|
  !                                           2: err =  |r|/(|a||x|+|b|)
  !                                           where r is the (preconditioned, recursive
  !                                           estimate of) residual 
  ! 
  Subroutine psb_zkrylov(method,a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop)
    use psb_base_mod
    use psb_prec_mod,only : psb_sprec_type, psb_dprec_type, psb_cprec_type, psb_zprec_type
    character(len=*)                   :: method
    Type(psb_zspmat_type), Intent(in)  :: a
    Type(psb_desc_type), Intent(in)    :: desc_a
    type(psb_zprec_type), intent(in)   :: prec 
    complex(psb_dpk_), Intent(in)    :: b(:)
    complex(psb_dpk_), Intent(inout) :: x(:)
    Real(psb_dpk_), Intent(in)       :: eps
    integer, intent(out)               :: info
    Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
    Integer, Optional, Intent(out)     :: iter
    Real(psb_dpk_), Optional, Intent(out) :: err

    integer                            :: ictxt,me,np,err_act
    character(len=20)             :: name

    info = 0
    name = 'psb_krylov'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)


    select case(psb_toupper(method))
    case('CG') 
      call  psb_cg(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    case('CGS') 
      call  psb_cgs(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    case('BICG') 
      call  psb_bicg(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,istop)
    case('BICGSTAB') 
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           & itmax,iter,err,itrace,istop)
    case('RGMRES')
      call  psb_rgmres(a,prec,b,x,eps,desc_a,info,&
           & itmax,iter,err,itrace,irst,istop)
    case('BICGSTABL')
      call  psb_bicgstabl(a,prec,b,x,eps,desc_a,info,&
         &itmax,iter,err,itrace,irst,istop)
    case default
      if (me==0) write(0,*) trim(name),': Warning: Unknown method  ',method,&
           & ', defaulting to BiCGSTAB'
      call  psb_bicgstab(a,prec,b,x,eps,desc_a,info,&
           &itmax,iter,err,itrace,istop)
    end select

    if(info/=0) then
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if

  end subroutine psb_zkrylov

  subroutine log_header(methdname)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)   :: methdname
    character(len=*), parameter    :: fmt='(a18,1x,a4,3(2x,a15))'
    integer, parameter             :: outlen=18 
    character(len=len(methdname))  :: mname
    character(len=outlen)          :: outname
    
    mname = adjustl(trim(methdname))
    write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
    write(*,fmt) adjustl(outname),'Iteration','Error Estimate','Tolerance'
    
  end subroutine log_header


  subroutine log_conv(methdname,me,itx,itrace,errnum,errden,eps)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)  :: methdname
    integer, intent(in)           :: me, itx, itrace
    real(psb_dpk_), intent(in)    :: errnum, errden, eps
    character(len=*), parameter   :: fmt='(a18,1x,i4,3(2x,es15.9))'
    integer, parameter            :: outlen=18 
    character(len=len(methdname)) :: mname
    character(len=outlen)         :: outname

    if ((mod(itx,itrace) == 0).and.(me==0)) then 
      mname = adjustl(trim(methdname))
      write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
      if (errden > dzero ) then 
        write(*,fmt) adjustl(outname),itx,errnum/errden,eps
      else
        write(*,fmt) adjustl(outname),itx,errnum,eps
      end if
    endif
    
  end subroutine log_conv

  subroutine log_end(methdname,me,it,errnum,errden,eps,err,iter)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in) :: methdname
    integer, intent(in)          :: me, it
    real(psb_dpk_), intent(in) :: errnum, errden, eps
    real(psb_dpk_), optional, intent(out) :: err
    integer, optional, intent(out)  :: iter

    character(len=*), parameter  :: fmt='(a,2x,es15.9,1x,a,1x,i4,1x,a)'
    character(len=*), parameter  :: fmt1='(a,3(2x,es15.9))'
    
    if (errden == dzero) then 
      if (errnum > eps) then         
        if (me==0) then 
          write(*,fmt) trim(methdname)//' failed to converge to ',eps,&
               & ' in ',it,' iterations. '
          write(*,fmt1) 'Last iteration error estimate: ',&
               & errnum
        end if
      end if
      if (present(err)) err=errnum
    else
      if (errnum/errden > eps) then         
        if (me==0) then 
          write(*,fmt) trim(methdname)//' failed to converge to ',eps,&
               & ' in ',it,' iterations. '
          write(*,fmt1) 'Last iteration error estimate: ',&
               & errnum/errden
        end if
      endif
      if (present(err)) err=errnum/errden
    end if
    if (present(iter)) iter = it

  end subroutine log_end
  
  subroutine psb_s_init_conv(methdname,stopc,trace,itmax,a,b,eps,desc_a,stopdat,info)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)      :: methdname
    integer, intent(in)               :: stopc, trace,itmax
    type(psb_s_sparse_mat), intent(in) :: a
    real(psb_spk_), intent(in)        :: b(:), eps
    type(psb_desc_type), intent(in)   :: desc_a
    type(psb_itconv_type)             :: stopdat
    integer, intent(out)              :: info
    
    integer                           :: ictxt, me, np, err_act
    character(len=20)                 :: name

    info = 0
    name = 'psb_init_conv'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)
    
    stopdat%controls(:) = 0
    stopdat%values(:)   = 0.0d0

    stopdat%controls(stopc_) = stopc
    stopdat%controls(trace_) = trace
    stopdat%controls(itmax_) = itmax
    
    select case(stopdat%controls(stopc_))
    case (1) 
      stopdat%values(ani_) = psb_spnrmi(a,desc_a,info)
      if (info == 0) stopdat%values(bni_) = psb_geamax(b,desc_a,info)

    case (2) 
      stopdat%values(bn2_) = psb_genrm2(b,desc_a,info)

    case default
      info=5001
      call psb_errpush(info,name,i_err=(/stopc,0,0,0,0/))
      goto 9999      
    end select
    if (info /= 0) then
      call psb_errpush(4001,name,a_err="Init conv check data")
      goto 9999
    end if
    
    stopdat%values(eps_)    = eps
    stopdat%values(errnum_) = dzero
    stopdat%values(errden_) = done
    
    if ((stopdat%controls(trace_) > 0).and. (me == 0))&
         &  call log_header(methdname) 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
   
  end subroutine psb_s_init_conv
  
  subroutine psb_d_init_conv(methdname,stopc,trace,itmax,a,b,eps,desc_a,stopdat,info)
    use psb_base_mod
    implicit none 
    type(psb_d_sparse_mat), intent(in)  :: a
    character(len=*), intent(in)      :: methdname
    integer, intent(in)               :: stopc, trace,itmax
    real(psb_dpk_), intent(in)      :: b(:), eps
    type(psb_desc_type), intent(in)   :: desc_a
    type(psb_itconv_type)             :: stopdat
    integer, intent(out)              :: info
    
    integer                           :: ictxt, me, np, err_act
    character(len=20)                 :: name

    info = 0
    name = 'psb_init_conv'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)
    
    stopdat%controls(:) = 0
    stopdat%values(:)   = 0.0d0

    stopdat%controls(stopc_) = stopc
    stopdat%controls(trace_) = trace
    stopdat%controls(itmax_) = itmax
    
    select case(stopdat%controls(stopc_))
    case (1) 
      stopdat%values(ani_) = psb_spnrmi(a,desc_a,info)
      if (info == 0) stopdat%values(bni_) = psb_geamax(b,desc_a,info)

    case (2) 
      stopdat%values(bn2_) = psb_genrm2(b,desc_a,info)

    case default
      info=5001
      call psb_errpush(info,name,i_err=(/stopc,0,0,0,0/))
      goto 9999      
    end select
    if (info /= 0) then
      call psb_errpush(4001,name,a_err="Init conv check data")
      goto 9999
    end if
    
    stopdat%values(eps_)    = eps
    stopdat%values(errnum_) = dzero
    stopdat%values(errden_) = done
    
    if ((stopdat%controls(trace_) > 0).and. (me == 0))&
         &  call log_header(methdname) 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
   
  end subroutine psb_d_init_conv
  
  subroutine psb_c_init_conv(methdname,stopc,trace,itmax,a,b,eps,desc_a,stopdat,info)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)      :: methdname
    integer, intent(in)               :: stopc, trace, itmax
    type(psb_cspmat_type), intent(in) :: a
    complex(psb_spk_), intent(in)   :: b(:)
    real(psb_spk_), intent(in)      :: eps
    type(psb_desc_type), intent(in)   :: desc_a
    type(psb_itconv_type)             :: stopdat
    integer, intent(out)              :: info
    
    integer                           :: ictxt, me, np, err_act
    character(len=20)                 :: name

    info = 0
    name = 'psb_init_conv'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)
    
    stopdat%controls(:) = 0
    stopdat%values(:)   = 0.0d0

    stopdat%controls(stopc_) = stopc
    stopdat%controls(trace_) = trace
    stopdat%controls(itmax_) = itmax
    
    select case(stopdat%controls(stopc_))
    case (1) 
      stopdat%values(ani_) = psb_spnrmi(a,desc_a,info)
      if (info == 0) stopdat%values(bni_) = psb_geamax(b,desc_a,info)

    case (2) 
      stopdat%values(bn2_) = psb_genrm2(b,desc_a,info)

    case default
      info=5001
      call psb_errpush(info,name,i_err=(/stopc,0,0,0,0/))
      goto 9999      
    end select
    if (info /= 0) then
      call psb_errpush(4001,name,a_err="Init conv check data")
      goto 9999
    end if
    
    stopdat%values(eps_)    = eps
    stopdat%values(errnum_) = dzero
    stopdat%values(errden_) = done

    if ((stopdat%controls(trace_) > 0).and. (me == 0))&
         &  call log_header(methdname) 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
   
  end subroutine psb_c_init_conv
  
  subroutine psb_z_init_conv(methdname,stopc,trace,itmax,a,b,eps,desc_a,stopdat,info)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)      :: methdname
    integer, intent(in)               :: stopc, trace, itmax
    type(psb_zspmat_type), intent(in) :: a
    complex(psb_dpk_), intent(in)   :: b(:)
    real(psb_dpk_), intent(in)      :: eps
    type(psb_desc_type), intent(in)   :: desc_a
    type(psb_itconv_type)             :: stopdat
    integer, intent(out)              :: info
    
    integer                           :: ictxt, me, np, err_act
    character(len=20)                 :: name

    info = 0
    name = 'psb_init_conv'
    call psb_erractionsave(err_act)


    ictxt=psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)
    
    stopdat%controls(:) = 0
    stopdat%values(:)   = 0.0d0

    stopdat%controls(stopc_) = stopc
    stopdat%controls(trace_) = trace
    stopdat%controls(itmax_) = itmax
    
    select case(stopdat%controls(stopc_))
    case (1) 
      stopdat%values(ani_) = psb_spnrmi(a,desc_a,info)
      if (info == 0) stopdat%values(bni_) = psb_geamax(b,desc_a,info)

    case (2) 
      stopdat%values(bn2_) = psb_genrm2(b,desc_a,info)

    case default
      info=5001
      call psb_errpush(info,name,i_err=(/stopc,0,0,0,0/))
      goto 9999      
    end select
    if (info /= 0) then
      call psb_errpush(4001,name,a_err="Init conv check data")
      goto 9999
    end if
    
    stopdat%values(eps_)    = eps
    stopdat%values(errnum_) = dzero
    stopdat%values(errden_) = done

    if ((stopdat%controls(trace_) > 0).and. (me == 0))&
         &  call log_header(methdname) 

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
   
  end subroutine psb_z_init_conv
    

  function psb_s_check_conv(methdname,it,x,r,desc_a,stopdat,info)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)    :: methdname
    integer, intent(in)             :: it
    real(psb_spk_), intent(in)      :: x(:), r(:)
    type(psb_desc_type), intent(in) :: desc_a
    type(psb_itconv_type)           :: stopdat
    logical                         :: psb_s_check_conv
    integer, intent(out)            :: info

    integer                         :: ictxt, me, np, err_act
    character(len=20)               :: name

    info = 0
    name = 'psb_check_conv'
    call psb_erractionsave(err_act)

    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt,me,np)

    psb_s_check_conv = .false. 
    
    select case(stopdat%controls(stopc_)) 
    case(1)
      stopdat%values(rni_) = psb_geamax(r,desc_a,info)
      if (info == 0) stopdat%values(xni_) = psb_geamax(x,desc_a,info)
      stopdat%values(errnum_) = stopdat%values(rni_)
      stopdat%values(errden_) = &
           & (stopdat%values(ani_)*stopdat%values(xni_)+stopdat%values(bni_))
    case(2)
      stopdat%values(rn2_)    = psb_genrm2(r,desc_a,info)
      stopdat%values(errnum_) = stopdat%values(rn2_)
      stopdat%values(errden_) = stopdat%values(bn2_)

    case default
      info=4001
      call psb_errpush(info,name,a_err="Control data in stopdat messed up!")
      goto 9999      
    end select
    if (info /= 0) then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if
    
    if (stopdat%values(errden_) == dzero) then 
      psb_s_check_conv = (stopdat%values(errnum_) <= stopdat%values(eps_))
    else
      psb_s_check_conv = &
           & (stopdat%values(errnum_) <= stopdat%values(eps_)*stopdat%values(errden_))
    end if

    psb_s_check_conv = (psb_s_check_conv.or.(stopdat%controls(itmax_) <= it))
    
    if ( (stopdat%controls(trace_) > 0).and.&
         & ((mod(it,stopdat%controls(trace_))==0).or.psb_s_check_conv)) then 
      call log_conv(methdname,me,it,1,stopdat%values(errnum_),&
           & stopdat%values(errden_),stopdat%values(eps_))
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    
  end function psb_s_check_conv

  function psb_d_check_conv(methdname,it,x,r,desc_a,stopdat,info)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)    :: methdname
    integer, intent(in)             :: it
    real(psb_dpk_), intent(in)    :: x(:), r(:)
    type(psb_desc_type), intent(in) :: desc_a
    type(psb_itconv_type)           :: stopdat
    logical                         :: psb_d_check_conv
    integer, intent(out)            :: info

    integer                         :: ictxt, me, np, err_act
    character(len=20)               :: name

    info = 0
    name = 'psb_check_conv'
    call psb_erractionsave(err_act)

    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt,me,np)

    psb_d_check_conv = .false. 
    
    select case(stopdat%controls(stopc_)) 
    case(1)
      stopdat%values(rni_) = psb_geamax(r,desc_a,info)
      if (info == 0) stopdat%values(xni_) = psb_geamax(x,desc_a,info)
      stopdat%values(errnum_) = stopdat%values(rni_)
      stopdat%values(errden_) = &
           & (stopdat%values(ani_)*stopdat%values(xni_)+stopdat%values(bni_))
    case(2)
      stopdat%values(rn2_)    = psb_genrm2(r,desc_a,info)
      stopdat%values(errnum_) = stopdat%values(rn2_)
      stopdat%values(errden_) = stopdat%values(bn2_)

    case default
      info=4001
      call psb_errpush(info,name,a_err="Control data in stopdat messed up!")
      goto 9999      
    end select
    if (info /= 0) then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if
    
    if (stopdat%values(errden_) == dzero) then 
      psb_d_check_conv = (stopdat%values(errnum_) <= stopdat%values(eps_))
    else
      psb_d_check_conv = &
           & (stopdat%values(errnum_) <= stopdat%values(eps_)*stopdat%values(errden_))
    end if

    psb_d_check_conv = (psb_d_check_conv.or.(stopdat%controls(itmax_) <= it))
    
    if ( (stopdat%controls(trace_) > 0).and.&
         & ((mod(it,stopdat%controls(trace_))==0).or.psb_d_check_conv)) then 
      call log_conv(methdname,me,it,1,stopdat%values(errnum_),&
           & stopdat%values(errden_),stopdat%values(eps_))
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    
  end function psb_d_check_conv


  function psb_c_check_conv(methdname,it,x,r,desc_a,stopdat,info)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)    :: methdname
    integer, intent(in)             :: it
    complex(psb_spk_), intent(in)   :: x(:), r(:)
    type(psb_desc_type), intent(in) :: desc_a
    type(psb_itconv_type)           :: stopdat
    logical                         :: psb_c_check_conv
    integer, intent(out)            :: info

    integer                         :: ictxt, me, np, err_act
    character(len=20)               :: name

    info = 0
    name = 'psb_check_conv'
    call psb_erractionsave(err_act)

    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt,me,np)
    psb_c_check_conv = .false. 
    
    select case(stopdat%controls(stopc_)) 
    case(1)
      stopdat%values(rni_) = psb_geamax(r,desc_a,info)
      if (info == 0) stopdat%values(xni_) = psb_geamax(x,desc_a,info)
      stopdat%values(errnum_) = stopdat%values(rni_)
      stopdat%values(errden_) = &
           & (stopdat%values(ani_)*stopdat%values(xni_)+stopdat%values(bni_))
    case(2)
      stopdat%values(rn2_)  = psb_genrm2(r,desc_a,info)
      stopdat%values(errnum_) = stopdat%values(rn2_)
      stopdat%values(errden_) = stopdat%values(bn2_)

    case default
      info=4001
      call psb_errpush(info,name,a_err="Control data in stopdat messed up!")
      goto 9999      
    end select
    if (info /= 0) then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if
    
    if (stopdat%values(errden_) == dzero) then 
      psb_c_check_conv = (stopdat%values(errnum_) <= stopdat%values(eps_))
    else
      psb_c_check_conv = &
           & (stopdat%values(errnum_) <= stopdat%values(eps_)*stopdat%values(errden_))
    end if

    psb_c_check_conv = (psb_c_check_conv.or.(stopdat%controls(itmax_) <= it))
    
    if ( (stopdat%controls(trace_) > 0).and.&
         & ((mod(it,stopdat%controls(trace_))==0).or.psb_c_check_conv)) then 
      call log_conv(methdname,me,it,1,stopdat%values(errnum_),&
           & stopdat%values(errden_),stopdat%values(eps_))
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    
  end function psb_c_check_conv

  function psb_z_check_conv(methdname,it,x,r,desc_a,stopdat,info)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)    :: methdname
    integer, intent(in)             :: it
    complex(psb_dpk_), intent(in) :: x(:), r(:)
    type(psb_desc_type), intent(in) :: desc_a
    type(psb_itconv_type)           :: stopdat
    logical                         :: psb_z_check_conv
    integer, intent(out)            :: info

    integer                         :: ictxt, me, np, err_act
    character(len=20)               :: name

    info = 0
    name = 'psb_check_conv'
    call psb_erractionsave(err_act)

    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt,me,np)
    psb_z_check_conv = .false. 
    
    select case(stopdat%controls(stopc_)) 
    case(1)
      stopdat%values(rni_) = psb_geamax(r,desc_a,info)
      if (info == 0) stopdat%values(xni_) = psb_geamax(x,desc_a,info)
      stopdat%values(errnum_) = stopdat%values(rni_)
      stopdat%values(errden_) = &
           & (stopdat%values(ani_)*stopdat%values(xni_)+stopdat%values(bni_))
    case(2)
      stopdat%values(rn2_)  = psb_genrm2(r,desc_a,info)
      stopdat%values(errnum_) = stopdat%values(rn2_)
      stopdat%values(errden_) = stopdat%values(bn2_)

    case default
      info=4001
      call psb_errpush(info,name,a_err="Control data in stopdat messed up!")
      goto 9999      
    end select
    if (info /= 0) then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if
    
    if (stopdat%values(errden_) == dzero) then 
      psb_z_check_conv = (stopdat%values(errnum_) <= stopdat%values(eps_))
    else
      psb_z_check_conv = &
           & (stopdat%values(errnum_) <= stopdat%values(eps_)*stopdat%values(errden_))
    end if

    psb_z_check_conv = (psb_z_check_conv.or.(stopdat%controls(itmax_) <= it))
    
    if ( (stopdat%controls(trace_) > 0).and.&
         & ((mod(it,stopdat%controls(trace_))==0).or.psb_z_check_conv)) then 
      call log_conv(methdname,me,it,1,stopdat%values(errnum_),&
           & stopdat%values(errden_),stopdat%values(eps_))
    end if

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    
  end function psb_z_check_conv

  subroutine psb_d_end_conv(methdname,it,desc_a,stopdat,info,err,iter)
    use psb_base_mod
    implicit none 
    character(len=*), intent(in)    :: methdname
    integer, intent(in)             :: it
    type(psb_desc_type), intent(in) :: desc_a
    type(psb_itconv_type)           :: stopdat
    integer, intent(out)            :: info
    real(psb_dpk_), optional, intent(out) :: err
    integer, optional, intent(out)  :: iter

    integer                         :: ictxt, me, np, err_act
    real(psb_dpk_)                  :: errnum, errden, eps
    character(len=20)               :: name

    info = 0
    name = 'psb_end_conv'

    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt,me,np)
    
    errnum = stopdat%values(errnum_) 
    errden = stopdat%values(errden_) 
    eps    = stopdat%values(eps_) 
    call log_end(methdname,me,it,errnum,errden,eps,err,iter)

  end subroutine psb_d_end_conv

end module psb_krylov_mod


  
