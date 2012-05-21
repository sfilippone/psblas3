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
module psi_i_mod
  use psb_descriptor_type, only : psb_desc_type, psb_ipk_, psb_mpik_
  use psb_i_vect_mod,  only : psb_i_base_vect_type

  interface
    subroutine psi_compute_size(desc_data,&
         & index_in, dl_lda, info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_) :: info, dl_lda
      integer(psb_ipk_) :: desc_data(:), index_in(:)
    end subroutine psi_compute_size
  end interface

  interface
    subroutine psi_crea_bnd_elem(bndel,desc_a,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), allocatable            :: bndel(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_crea_bnd_elem
  end interface

  interface
    subroutine psi_crea_index(desc_a,index_in,index_out,glob_idx,nxch,nsnd,nrcv,info)
      import :: psb_desc_type, psb_ipk_
      type(psb_desc_type), intent(in)     :: desc_a
      integer(psb_ipk_), intent(out)                :: info,nxch,nsnd,nrcv
      integer(psb_ipk_), intent(in)                 :: index_in(:)
      integer(psb_ipk_), allocatable, intent(inout) :: index_out(:)
      logical                             :: glob_idx
    end subroutine psi_crea_index
  end interface

  interface
    subroutine psi_crea_ovr_elem(me,desc_overlap,ovr_elem,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)               :: me, desc_overlap(:)
      integer(psb_ipk_), allocatable, intent(out) :: ovr_elem(:,:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psi_crea_ovr_elem
  end interface

  interface
    subroutine psi_desc_index(desc,index_in,dep_list,&
         & length_dl,nsnd,nrcv,desc_index,isglob_in,info)
      import :: psb_desc_type, psb_ipk_
      type(psb_desc_type) :: desc
      integer(psb_ipk_) :: index_in(:),dep_list(:)
      integer(psb_ipk_),allocatable  :: desc_index(:)
      integer(psb_ipk_) :: length_dl,nsnd,nrcv,info
      logical         :: isglob_in
    end subroutine psi_desc_index
  end interface

  interface
    subroutine psi_dl_check(dep_list,dl_lda,np,length_dl)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_) :: np,dl_lda,length_dl(0:np)
      integer(psb_ipk_) :: dep_list(dl_lda,0:np)
    end subroutine psi_dl_check
  end interface

  interface
    subroutine psi_sort_dl(dep_list,l_dep_list,np,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_) :: np,dep_list(:,:), l_dep_list(:), info
    end subroutine psi_sort_dl
  end interface

  interface
    subroutine psi_extract_dep_list(ictxt,is_bld,is_upd,desc_str,dep_list,&
         & length_dl,np,dl_lda,mode,info)
      import :: psb_desc_type, psb_ipk_, psb_mpik_
      logical :: is_bld, is_upd
      integer(psb_mpik_) :: ictxt
      integer(psb_ipk_) :: np,dl_lda,mode, info
      integer(psb_ipk_) :: desc_str(*),dep_list(dl_lda,0:np),length_dl(0:np)
    end subroutine psi_extract_dep_list
  end interface

  interface psi_fnd_owner
    subroutine psi_fnd_owner(nv,idx,iprc,desc,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in) :: nv
      integer(psb_ipk_), intent(in) ::  idx(:)
      integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
      type(psb_desc_type), intent(in) :: desc
      integer(psb_ipk_), intent(out) :: info
    end subroutine psi_fnd_owner
  end interface

  interface psi_ldsc_pre_halo
    subroutine psi_ldsc_pre_halo(desc,ext_hv,info)
      import :: psb_desc_type, psb_ipk_
      type(psb_desc_type), intent(inout) :: desc
      logical, intent(in)  :: ext_hv
      integer(psb_ipk_), intent(out) :: info
    end subroutine psi_ldsc_pre_halo
  end interface

  interface psi_bld_tmphalo
    subroutine psi_bld_tmphalo(desc,info)
      import :: psb_desc_type, psb_ipk_
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out) :: info
    end subroutine psi_bld_tmphalo
  end interface


  interface psi_bld_tmpovrl
    subroutine psi_bld_tmpovrl(iv,desc,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: iv(:)
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out) :: info
    end subroutine psi_bld_tmpovrl
  end interface


  interface psi_idx_cnv
    subroutine psi_idx_cnv1(nv,idxin,desc,info,mask,owned)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)    :: nv
      integer(psb_ipk_), intent(inout) ::  idxin(:)
      type(psb_desc_type), intent(in) :: desc
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional :: mask(:)
      logical, intent(in), optional :: owned
    end subroutine psi_idx_cnv1
    subroutine psi_idx_cnv2(nv,idxin,idxout,desc,info,mask,owned)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: nv, idxin(:)
      integer(psb_ipk_), intent(out) :: idxout(:)
      type(psb_desc_type), intent(in) :: desc
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional :: mask(:)
      logical, intent(in), optional :: owned
    end subroutine psi_idx_cnv2
    subroutine psi_idx_cnvs(idxin,idxout,desc,info,mask,owned)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: idxin
      integer(psb_ipk_), intent(out) :: idxout
      type(psb_desc_type), intent(in) :: desc
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional :: mask
      logical, intent(in), optional :: owned
    end subroutine psi_idx_cnvs
    subroutine psi_idx_cnvs1(idxin,desc,info,mask,owned)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(inout)  :: idxin
      type(psb_desc_type), intent(in) :: desc
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional :: mask
      logical, intent(in), optional :: owned
    end subroutine psi_idx_cnvs1
  end interface

  interface psi_idx_ins_cnv
    subroutine psi_idx_ins_cnv1(nv,idxin,desc,info,mask,lidx)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)    :: nv
      integer(psb_ipk_), intent(inout) ::  idxin(:)
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional :: mask(:)
      integer, intent(in), optional :: lidx(:)
    end subroutine psi_idx_ins_cnv1
    subroutine psi_idx_ins_cnv2(nv,idxin,idxout,desc,info,mask,lidx)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: nv, idxin(:)
      integer(psb_ipk_), intent(out) :: idxout(:)
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional :: mask(:)
      integer, intent(in), optional :: lidx(:)
    end subroutine psi_idx_ins_cnv2
    subroutine psi_idx_ins_cnvs2(idxin,idxout,desc,info,mask,lidx)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: idxin
      integer(psb_ipk_), intent(out) :: idxout
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional :: mask
      integer, intent(in), optional :: lidx
    end subroutine psi_idx_ins_cnvs2
    subroutine psi_idx_ins_cnvs1(idxin,desc,info,mask,lidx)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(inout)  :: idxin
      type(psb_desc_type), intent(inout) :: desc
      integer(psb_ipk_), intent(out) :: info
      logical, intent(in), optional :: mask
      integer, intent(in), optional :: lidx
    end subroutine psi_idx_ins_cnvs1
  end interface

  interface psi_cnv_dsc
    subroutine psi_cnv_dsc(halo_in,ovrlap_in,ext_in,cdesc, info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)                :: halo_in(:), ovrlap_in(:),ext_in(:)
      type(psb_desc_type), intent(inout) :: cdesc
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psi_cnv_dsc
  end interface

  interface psi_renum_index
    subroutine psi_renum_index(iperm,idx,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(out)   :: info
      integer(psb_ipk_), intent(in)    :: iperm(:)
      integer(psb_ipk_), intent(inout) :: idx(:)
    end subroutine psi_renum_index
  end interface

  interface psi_inner_cnv
    subroutine psi_inner_cnvs(x,hashmask,hashv,glb_lc)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)    :: hashmask,hashv(0:),glb_lc(:,:)
      integer(psb_ipk_), intent(inout) :: x
    end subroutine psi_inner_cnvs
    subroutine psi_inner_cnvs2(x,y,hashmask,hashv,glb_lc)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: hashmask,hashv(0:),glb_lc(:,:)
      integer(psb_ipk_), intent(in)  :: x
      integer(psb_ipk_), intent(out) :: y
    end subroutine psi_inner_cnvs2
    subroutine psi_inner_cnv1(n,x,hashmask,hashv,glb_lc,mask)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)    :: n,hashmask,hashv(0:),glb_lc(:,:)
      logical, intent(in), optional    :: mask(:)
      integer(psb_ipk_), intent(inout) :: x(:)
    end subroutine psi_inner_cnv1
    subroutine psi_inner_cnv2(n,x,y,hashmask,hashv,glb_lc,mask)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: n, hashmask,hashv(0:),glb_lc(:,:)
      logical, intent(in),optional  :: mask(:)
      integer(psb_ipk_), intent(in)  :: x(:)
      integer(psb_ipk_), intent(out) :: y(:)
    end subroutine psi_inner_cnv2
  end interface

  interface 
    subroutine psi_bld_ovr_mst(me,ovrlap_elem,mst_idx,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)               :: me, ovrlap_elem(:,:)
      integer(psb_ipk_), allocatable, intent(out) :: mst_idx(:) 
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psi_bld_ovr_mst
  end interface


  interface psi_swapdata
    subroutine psi_iswapdatam(flag,n,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)         :: flag, n
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_) :: y(:,:), beta
      integer(psb_ipk_), target             :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_iswapdatam
    subroutine psi_iswapdatav(flag,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_) :: y(:), beta
      integer(psb_ipk_), target             :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_iswapdatav
    subroutine psi_iswapdata_vect(flag,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_,  psb_i_base_vect_type
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_vect_type) :: y
      integer(psb_ipk_)              :: beta 
      integer(psb_ipk_),target       :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_iswapdata_vect
    subroutine psi_iswapidxm(ictxt,icomm,flag,n,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: ictxt,icomm,flag, n
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_) :: y(:,:), beta
      integer(psb_ipk_),target       :: work(:)
      integer(psb_ipk_), intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_iswapidxm
    subroutine psi_iswapidxv(ictxt,icomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: ictxt,icomm,flag
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_) :: y(:), beta
      integer(psb_ipk_),target       :: work(:)
      integer(psb_ipk_), intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_iswapidxv
    subroutine psi_iswapidx_vect(ictxt,icomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_,  psb_i_base_vect_type
      integer(psb_ipk_), intent(in)   :: ictxt,icomm,flag
      integer(psb_ipk_), intent(out)  :: info
      class(psb_i_base_vect_type) :: y
      integer(psb_ipk_)        :: beta
      integer(psb_ipk_),target :: work(:)
      integer(psb_ipk_), intent(in)   :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_iswapidx_vect
  end interface


  interface psi_swaptran
    subroutine psi_iswaptranm(flag,n,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)         :: flag, n
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_) :: y(:,:), beta
      integer(psb_ipk_),target              :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_iswaptranm
    subroutine psi_iswaptranv(flag,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      integer(psb_ipk_) :: y(:), beta
      integer(psb_ipk_),target              :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_iswaptranv
    subroutine psi_iswaptran_vect(flag,beta,y,desc_a,work,info,data)
      import :: psb_desc_type, psb_ipk_,  psb_i_base_vect_type
      integer(psb_ipk_), intent(in)         :: flag
      integer(psb_ipk_), intent(out)        :: info
      class(psb_i_base_vect_type) :: y
      integer(psb_ipk_)              :: beta
      integer(psb_ipk_),target       :: work(:)
      type(psb_desc_type), target :: desc_a
      integer(psb_ipk_), optional           :: data
    end subroutine psi_iswaptran_vect
    subroutine psi_itranidxm(ictxt,icomm,flag,n,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: ictxt,icomm,flag, n
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_) :: y(:,:), beta
      integer(psb_ipk_), target      :: work(:)
      integer(psb_ipk_), intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_itranidxm
    subroutine psi_itranidxv(ictxt,icomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(in)  :: ictxt,icomm,flag
      integer(psb_ipk_), intent(out) :: info
      integer(psb_ipk_) :: y(:), beta
      integer(psb_ipk_), target      :: work(:)
      integer(psb_ipk_), intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_itranidxv
    subroutine psi_itranidx_vect(ictxt,icomm,flag,beta,y,idx,&
         & totxch,totsnd,totrcv,work,info)
      import :: psb_desc_type, psb_ipk_,  psb_i_base_vect_type
      integer(psb_ipk_), intent(in)   :: ictxt,icomm,flag
      integer(psb_ipk_), intent(out)  :: info
      class(psb_i_base_vect_type) :: y
      integer(psb_ipk_)        :: beta
      integer(psb_ipk_),target :: work(:)
      integer(psb_ipk_), intent(in)   :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_itranidx_vect
  end interface

  interface psi_ovrl_upd
    subroutine  psi_iovrl_updr1(x,desc_a,update,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(inout), target     :: x(:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(in)                :: update
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psi_iovrl_updr1
    subroutine  psi_iovrl_updr2(x,desc_a,update,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(inout), target     :: x(:,:)
      type(psb_desc_type), intent(in)    :: desc_a
      integer(psb_ipk_), intent(in)                :: update
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psi_iovrl_updr2
    subroutine  psi_iovrl_upd_vect(x,desc_a,update,info)
      import :: psb_desc_type, psb_ipk_, psb_i_base_vect_type
      class(psb_i_base_vect_type)       :: x
      type(psb_desc_type), intent(in)   :: desc_a
      integer(psb_ipk_), intent(in)               :: update
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psi_iovrl_upd_vect
  end interface

  interface psi_ovrl_save
    subroutine  psi_iovrl_saver1(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(inout)  :: x(:)
      integer(psb_ipk_), allocatable    :: xs(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_iovrl_saver1
    subroutine  psi_iovrl_saver2(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(inout)  :: x(:,:)
      integer(psb_ipk_), allocatable    :: xs(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_iovrl_saver2
    subroutine  psi_iovrl_save_vect(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_, psb_i_base_vect_type
      class(psb_i_base_vect_type)     :: x
      integer(psb_ipk_), allocatable     :: xs(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_iovrl_save_vect
  end interface

  interface psi_ovrl_restore
    subroutine  psi_iovrl_restrr1(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(inout)           :: x(:)
      integer(psb_ipk_) :: xs(:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_iovrl_restrr1
    subroutine  psi_iovrl_restrr2(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_
      integer(psb_ipk_), intent(inout)           :: x(:,:)
      integer(psb_ipk_) :: xs(:,:)
      type(psb_desc_type), intent(in)  :: desc_a
      integer(psb_ipk_), intent(out)             :: info
    end subroutine psi_iovrl_restrr2
    subroutine  psi_iovrl_restr_vect(x,xs,desc_a,info)
      import :: psb_desc_type, psb_ipk_, psb_i_base_vect_type
      class(psb_i_base_vect_type)     :: x
      integer(psb_ipk_)                  :: xs(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer(psb_ipk_), intent(out)            :: info
    end subroutine psi_iovrl_restr_vect
  end interface

end module psi_i_mod

