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

module psb_prec_mod
  use psb_prec_type

  interface psb_precbld
    subroutine psb_sprecbld(a,desc_a,prec,info,upd)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_type, only : psb_sprec_type
      implicit none
      type(psb_s_sparse_mat), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_sprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
    end subroutine psb_sprecbld
    subroutine psb_dprecbld(a,desc_a,prec,info,upd)
      use psb_base_mod, only  : psb_desc_type, psb_d_sparse_mat, psb_dpk_
      use psb_prec_type, only : psb_dprec_type
      implicit none
      type(psb_d_sparse_mat), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_dprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
    end subroutine psb_dprecbld
    subroutine psb_cprecbld(a,desc_a,prec,info,upd)
      use psb_base_mod, only  : psb_desc_type, psb_c_sparse_mat, psb_spk_
      use psb_prec_type, only : psb_cprec_type
      implicit none
      type(psb_c_sparse_mat), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_cprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
    end subroutine psb_cprecbld
    subroutine psb_zprecbld(a,desc_a,prec,info,upd)
      use psb_base_mod, only  : psb_desc_type, psb_z_sparse_mat, psb_dpk_
      use psb_prec_type, only : psb_zprec_type
      implicit none
      type(psb_z_sparse_mat), intent(in), target  :: a
      type(psb_desc_type), intent(in), target    :: desc_a
      type(psb_zprec_type), intent(inout)        :: prec
      integer, intent(out)                       :: info
      character, intent(in),optional             :: upd
    end subroutine psb_zprecbld
  end interface

  interface psb_precinit
    subroutine psb_sprecinit(prec,ptype,info)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_type, only : psb_sprec_type
      implicit none
      type(psb_sprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
    end subroutine psb_sprecinit
    subroutine psb_dprecinit(prec,ptype,info)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      use psb_prec_type, only : psb_dprec_type
      implicit none
      type(psb_dprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
    end subroutine psb_dprecinit
    subroutine psb_cprecinit(prec,ptype,info)
      use psb_base_mod
      use psb_prec_type
      implicit none
      type(psb_cprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
    end subroutine psb_cprecinit
    subroutine psb_zprecinit(prec,ptype,info)
      use psb_base_mod
      use psb_prec_type
      implicit none
      type(psb_zprec_type), intent(inout)    :: prec
      character(len=*), intent(in)           :: ptype
      integer, intent(out)                   :: info
    end subroutine psb_zprecinit
  end interface

  interface psb_precset
    subroutine psb_sprecseti(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_type, only : psb_sprec_type
      implicit none
      type(psb_sprec_type), intent(inout)    :: prec
      integer                                :: what, val 
      integer, intent(out)                   :: info
    end subroutine psb_sprecseti
    subroutine psb_sprecsets(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat, psb_spk_
      use psb_prec_type, only : psb_sprec_type
      implicit none
      type(psb_sprec_type), intent(inout)    :: prec
      integer                                :: what
      real(psb_spk_)                       :: val 
      integer, intent(out)                   :: info
    end subroutine psb_sprecsets
    subroutine psb_dprecseti(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      use psb_prec_type, only : psb_dprec_type
      implicit none
      type(psb_dprec_type), intent(inout)    :: prec
      integer                                :: what, val 
      integer, intent(out)                   :: info
    end subroutine psb_dprecseti
    subroutine psb_dprecsetd(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_dpk_
      use psb_prec_type, only : psb_dprec_type
      implicit none
      type(psb_dprec_type), intent(inout)    :: prec
      integer                                :: what
      real(psb_dpk_)                       :: val 
      integer, intent(out)                   :: info
    end subroutine psb_dprecsetd
    subroutine psb_cprecseti(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_c_sparse_mat, psb_spk_
      use psb_prec_type, only : psb_cprec_type
      implicit none
      type(psb_cprec_type), intent(inout)    :: prec
      integer                                :: what, val 
      integer, intent(out)                   :: info
    end subroutine psb_cprecseti
    subroutine psb_cprecsets(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_c_sparse_mat, psb_spk_
      use psb_prec_type, only : psb_cprec_type
      implicit none
      type(psb_cprec_type), intent(inout)    :: prec
      integer                                :: what
      real(psb_spk_)                       :: val 
      integer, intent(out)                   :: info
    end subroutine psb_cprecsets
    subroutine psb_zprecseti(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_z_sparse_mat, psb_dpk_
      use psb_prec_type, only : psb_zprec_type
      implicit none
      type(psb_zprec_type), intent(inout)    :: prec
      integer                                :: what, val 
      integer, intent(out)                   :: info
    end subroutine psb_zprecseti
    subroutine psb_zprecsetd(prec,what,val,info)
      use psb_base_mod, only  : psb_desc_type, psb_z_sparse_mat, psb_dpk_
      use psb_prec_type, only : psb_zprec_type
      implicit none
      type(psb_zprec_type), intent(inout)    :: prec
      integer                                :: what
      real(psb_dpk_)                       :: val 
      integer, intent(out)                   :: info
    end subroutine psb_zprecsetd
  end interface

  interface psb_ilu_fct
    subroutine psb_silu_fct(a,l,u,d,info,blck)
      use psb_base_mod, only  : psb_desc_type, psb_s_sparse_mat,&
           & psb_s_csr_sparse_mat, psb_spk_
      integer, intent(out)                ::     info
      type(psb_s_sparse_mat),intent(in)    :: a
      type(psb_s_csr_sparse_mat),intent(inout) :: l,u
      type(psb_s_sparse_mat),intent(in), optional, target :: blck
      real(psb_spk_), intent(inout)     ::  d(:)
    end subroutine psb_silu_fct
    subroutine psb_dilu_fct(a,l,u,d,info,blck)
      use psb_base_mod, only  : psb_desc_type, psb_d_sparse_mat,&
           & psb_d_csr_sparse_mat, psb_dpk_
      integer, intent(out)                ::     info
      type(psb_d_sparse_mat),intent(in)    :: a
      type(psb_d_csr_sparse_mat),intent(inout) :: l,u
      type(psb_d_sparse_mat),intent(in), optional, target :: blck
      real(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine psb_dilu_fct
    subroutine psb_cilu_fct(a,l,u,d,info,blck)
      use psb_base_mod, only  : psb_desc_type, psb_c_sparse_mat, &
           & psb_c_csr_sparse_mat, psb_spk_
      integer, intent(out)                ::     info
      type(psb_c_sparse_mat),intent(in)    :: a
      type(psb_c_csr_sparse_mat),intent(inout) :: l,u
      type(psb_c_sparse_mat),intent(in), optional, target :: blck
      complex(psb_spk_), intent(inout)     ::  d(:)
    end subroutine psb_cilu_fct
    subroutine psb_zilu_fct(a,l,u,d,info,blck)
      use psb_base_mod, only  : psb_desc_type, psb_z_sparse_mat, &
           & psb_z_csr_sparse_mat, psb_dpk_
      integer, intent(out)                ::     info
      type(psb_z_sparse_mat),intent(in)    :: a
      type(psb_z_csr_sparse_mat),intent(inout) :: l,u
      type(psb_z_sparse_mat),intent(in), optional, target :: blck
      complex(psb_dpk_), intent(inout)     ::  d(:)
    end subroutine psb_zilu_fct
  end interface


end module psb_prec_mod
