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
module psi_i_mod
  use psb_desc_mod, only : psb_desc_type, psb_ipk_, psb_mpk_, psb_epk_, psb_lpk_
  use psi_m_comm_a_mod
  use psi_e_comm_a_mod
  use psb_i_base_vect_mod, only : psb_i_base_vect_type 
  use psb_i_base_multivect_mod, only : psb_i_base_multivect_type 
  use psi_i_comm_v_mod

  
  interface psi_compute_size
    subroutine psi_i_compute_size(desc_data,&
         & index_in, dl_lda, info)
      import 
      integer(psb_ipk_) :: info
      integer(psb_ipk_) :: dl_lda
      integer(psb_ipk_) :: desc_data(:), index_in(:)
    end subroutine psi_i_compute_size
  end interface

  interface psi_crea_bnd_elem
    subroutine psi_i_crea_bnd_elem(bndel,desc_a,info)
      import 
      integer(psb_ipk_), allocatable            :: bndel(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_i_crea_bnd_elem
  end interface

  interface psi_crea_index
    subroutine psi_i_crea_index(desc_a,index_in,index_out,glob_idx,nxch,nsnd,nrcv,info)
      import 
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), intent(out)                :: nxch,nsnd,nrcv
      integer(psb_ipk_), intent(in)                 :: index_in(:)
      integer(psb_ipk_), allocatable, intent(inout) :: index_out(:)
      logical                             :: glob_idx
      integer(psb_ipk_), intent(out)      :: info
    end subroutine psi_i_crea_index
  end interface

  interface psi_crea_ovr_elem
    subroutine psi_i_crea_ovr_elem(me,desc_overlap,ovr_elem,info)
      import 
      integer(psb_ipk_), intent(in)               :: me, desc_overlap(:)
      integer(psb_ipk_), allocatable, intent(out) :: ovr_elem(:,:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psi_i_crea_ovr_elem
  end interface

  interface psi_desc_index
    subroutine psi_i_desc_index(desc,index_in,dep_list,&
         & length_dl,nsnd,nrcv,desc_index,isglob_in,info)
      import 
      type(psb_desc_type) :: desc
      integer(psb_lpk_) :: index_in(:)
      integer(psb_ipk_) :: dep_list(:)
      integer(psb_ipk_),allocatable  :: desc_index(:)
      integer(psb_ipk_) :: length_dl,nsnd,nrcv
      logical         :: isglob_in
      integer(psb_ipk_) :: info
    end subroutine psi_i_desc_index
  end interface

  interface psi_dl_check
    subroutine psi_i_dl_check(dep_list,dl_lda,np,length_dl)
      import 
      integer(psb_ipk_) :: np,dl_lda,length_dl(0:np)
      integer(psb_ipk_) :: dep_list(dl_lda,0:np)
    end subroutine psi_i_dl_check
  end interface

  interface psi_sort_dl
    subroutine psi_i_sort_dl(dep_list,l_dep_list,np,info)
      import 
      integer(psb_ipk_) :: dep_list(:,:), l_dep_list(:)
      integer(psb_mpk_) :: np
      integer(psb_ipk_) :: info
    end subroutine psi_i_sort_dl
  end interface

  interface psi_extract_dep_list
    subroutine psi_i_extract_dep_list(ictxt,is_bld,is_upd,desc_str,dep_list,&
         & length_dl,np,dl_lda,mode,info)
      import 
      logical :: is_bld, is_upd
      integer(psb_ipk_) :: ictxt
      integer(psb_ipk_) :: dl_lda,mode
      integer(psb_lpk_) :: desc_str(*)
      integer(psb_ipk_) :: dep_list(dl_lda,0:np),length_dl(0:np)
      integer(psb_mpk_) :: np
      integer(psb_ipk_) :: info
    end subroutine psi_i_extract_dep_list
  end interface

  interface psi_fnd_owner
    subroutine psi_i_fnd_owner(nv,idx,iprc,desc,info)
      import 
      integer(psb_ipk_), intent(in) :: nv
      integer(psb_ipk_), intent(in) ::  idx(:)
      integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
      type(psb_desc_type), intent(in) :: desc
      integer(psb_ipk_), intent(out) :: info
    end subroutine psi_i_fnd_owner
  end interface psi_fnd_owner

  interface psi_bld_tmphalo
    subroutine psi_bld_tmphalo(desc,info)
      import 
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out) :: info
    end subroutine psi_bld_tmphalo
  end interface psi_bld_tmphalo


  interface psi_bld_tmpovrl
    subroutine psi_i_bld_tmpovrl(iv,desc,info)
      import 
      integer(psb_lpk_), intent(in)  :: iv(:)
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out) :: info
    end subroutine psi_i_bld_tmpovrl
  end interface psi_bld_tmpovrl

  interface psi_cnv_dsc
    subroutine psi_i_cnv_dsc(halo_in,ovrlap_in,ext_in,cdesc, info, mold)
      import 
      integer(psb_ipk_), intent(in)        :: halo_in(:), ovrlap_in(:),ext_in(:)
      type(psb_desc_type), intent(inout) :: cdesc
      integer(psb_ipk_), intent(out)               :: info
      class(psb_i_base_vect_type), optional, intent(in) :: mold
    end subroutine psi_i_cnv_dsc
  end interface psi_cnv_dsc

  interface psi_renum_index
    subroutine psi_i_renum_index(iperm,idx,info)
      import 
      integer(psb_ipk_), intent(out)   :: info
      integer(psb_ipk_), intent(in)    :: iperm(:)
      integer(psb_ipk_), intent(inout) :: idx(:)
    end subroutine psi_i_renum_index
  end interface psi_renum_index

  interface psi_inner_cnv
    subroutine psi_i_inner_cnvs(x,hashmask,hashv,glb_lc)
      import 
      integer(psb_ipk_), intent(in)    :: hashmask,hashv(0:),glb_lc(:,:)
      integer(psb_ipk_), intent(inout) :: x
    end subroutine psi_i_inner_cnvs
    subroutine psi_i_inner_cnvs2(x,y,hashmask,hashv,glb_lc)
      import 
      integer(psb_ipk_), intent(in)  :: hashmask,hashv(0:),glb_lc(:,:)
      integer(psb_ipk_), intent(in)  :: x
      integer(psb_ipk_), intent(out) :: y
    end subroutine psi_i_inner_cnvs2
    subroutine psi_i_inner_cnv1(n,x,hashmask,hashv,glb_lc,mask)
      import 
      integer(psb_ipk_), intent(in)    :: n,hashmask,hashv(0:),glb_lc(:,:)
      logical, intent(in), optional    :: mask(:)
      integer(psb_ipk_), intent(inout) :: x(:)
    end subroutine psi_i_inner_cnv1
    subroutine psi_i_inner_cnv2(n,x,y,hashmask,hashv,glb_lc,mask)
      import 
      integer(psb_ipk_), intent(in)  :: n, hashmask,hashv(0:),glb_lc(:,:)
      logical, intent(in),optional  :: mask(:)
      integer(psb_ipk_), intent(in)  :: x(:)
      integer(psb_ipk_), intent(out) :: y(:)
    end subroutine psi_i_inner_cnv2
  end interface psi_inner_cnv

  interface psi_bld_ovr_mst
    subroutine psi_i_bld_ovr_mst(me,ovrlap_elem,mst_idx,info)
      import 
      integer(psb_ipk_), intent(in)               :: me, ovrlap_elem(:,:)
      integer(psb_ipk_), allocatable, intent(out) :: mst_idx(:) 
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psi_i_bld_ovr_mst
  end interface
end module psi_i_mod

