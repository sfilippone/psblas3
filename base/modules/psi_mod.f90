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
module psi_mod

  use psi_serial_mod


  interface
    subroutine psi_compute_size(desc_data,&
         & index_in, dl_lda, info)
      integer  :: info, dl_lda
      integer  :: desc_data(:), index_in(:)
    end subroutine psi_compute_size
  end interface

  interface
    subroutine psi_crea_bnd_elem(bndel,desc_a,info)
      use psb_descriptor_type
      integer, allocatable            :: bndel(:)
      type(psb_desc_type), intent(in) :: desc_a
      integer, intent(out)            :: info
    end subroutine psi_crea_bnd_elem
  end interface

  interface
    subroutine psi_crea_index(desc_a,index_in,index_out,glob_idx,nxch,nsnd,nrcv,info)
      use psb_descriptor_type
      type(psb_desc_type), intent(in)     :: desc_a
      integer, intent(out)                :: info,nxch,nsnd,nrcv
      integer, intent(in)                 :: index_in(:)
      integer, allocatable, intent(inout) :: index_out(:)
      logical                             :: glob_idx
    end subroutine psi_crea_index
  end interface

  interface
    subroutine psi_crea_ovr_elem(me,desc_overlap,ovr_elem,info)
      integer, intent(in)               :: me, desc_overlap(:)
      integer, allocatable, intent(out) :: ovr_elem(:,:)
      integer, intent(out)              :: info
    end subroutine psi_crea_ovr_elem
  end interface

  interface
    subroutine psi_desc_index(desc,index_in,dep_list,&
         & length_dl,nsnd,nrcv,desc_index,isglob_in,info)
      use psb_descriptor_type
      type(psb_desc_type) :: desc
      integer         :: index_in(:),dep_list(:)
      integer,allocatable, intent(inout)  :: desc_index(:)
      integer         :: length_dl,nsnd,nrcv,info
      logical         :: isglob_in
    end subroutine psi_desc_index
  end interface

  interface
    subroutine psi_dl_check(dep_list,dl_lda,np,length_dl)
      integer  :: np,dl_lda,length_dl(0:np)
      integer  :: dep_list(dl_lda,0:np)
    end subroutine psi_dl_check
  end interface

  interface
    subroutine psi_sort_dl(dep_list,l_dep_list,np,info)
      integer :: np,dep_list(:,:), l_dep_list(:), info
    end subroutine psi_sort_dl
  end interface

  interface psi_swapdata
    subroutine psi_dswapdatam(flag,n,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag, n
      integer, intent(out) :: info
      real(psb_dpk_)     :: y(:,:), beta
      real(psb_dpk_),target      :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_dswapdatam
    subroutine psi_dswapdatav(flag,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag
      integer, intent(out) :: info
      real(psb_dpk_)     :: y(:), beta 
      real(psb_dpk_),target      :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_dswapdatav
    subroutine psi_dswapidxm(ictxt,icomm,flag,n,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag, n
      integer, intent(out) :: info
      real(psb_dpk_)     :: y(:,:), beta
      real(psb_dpk_),target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dswapidxm
    subroutine psi_dswapidxv(ictxt,icomm,flag,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag
      integer, intent(out) :: info
      real(psb_dpk_)     :: y(:), beta
      real(psb_dpk_),target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dswapidxv
    subroutine psi_iswapdatam(flag,n,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag, n
      integer, intent(out) :: info
      integer              :: y(:,:), beta
      integer, target              :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_iswapdatam
    subroutine psi_iswapdatav(flag,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag
      integer, intent(out) :: info
      integer              :: y(:), beta
      integer, target              :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_iswapdatav
    subroutine psi_iswapidxm(ictxt,icomm,flag,n,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag, n
      integer, intent(out) :: info
      integer              :: y(:,:), beta
      integer,target       :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_iswapidxm
    subroutine psi_iswapidxv(ictxt,icomm,flag,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag
      integer, intent(out) :: info
      integer              :: y(:), beta
      integer,target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_iswapidxv
    subroutine psi_zswapdatam(flag,n,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag, n
      integer, intent(out) :: info
      complex(psb_dpk_)     :: y(:,:), beta
      complex(psb_dpk_),target   :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_zswapdatam
    subroutine psi_zswapdatav(flag,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag
      integer, intent(out) :: info
      complex(psb_dpk_)     :: y(:), beta
      complex(psb_dpk_),target   :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_zswapdatav
    subroutine psi_zswapidxm(ictxt,icomm,flag,n,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag, n
      integer, intent(out) :: info
      complex(psb_dpk_)     :: y(:,:), beta
      complex(psb_dpk_),target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_zswapidxm
    subroutine psi_zswapidxv(ictxt,icomm,flag,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag
      integer, intent(out) :: info
      complex(psb_dpk_)     :: y(:), beta
      complex(psb_dpk_),target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_zswapidxv
  end interface


  interface psi_swaptran
    subroutine psi_dswaptranm(flag,n,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag, n
      integer, intent(out) :: info
      real(psb_dpk_)     :: y(:,:), beta
      real(psb_dpk_),target     :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_dswaptranm
    subroutine psi_dswaptranv(flag,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag
      integer, intent(out) :: info
      real(psb_dpk_)     :: y(:), beta
      real(psb_dpk_),target     :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_dswaptranv
    subroutine psi_dtranidxm(ictxt,icomm,flag,n,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag, n
      integer, intent(out) :: info
      real(psb_dpk_)     :: y(:,:), beta
      real(psb_dpk_),target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dtranidxm
    subroutine psi_dtranidxv(ictxt,icomm,flag,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag
      integer, intent(out) :: info
      real(psb_dpk_)     :: y(:), beta
      real(psb_dpk_),target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_dtranidxv
    subroutine psi_iswaptranm(flag,n,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag, n
      integer, intent(out) :: info
      integer              :: y(:,:), beta
      integer,target               :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_iswaptranm
    subroutine psi_iswaptranv(flag,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag
      integer, intent(out) :: info
      integer              :: y(:), beta
      integer,target               :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_iswaptranv
    subroutine psi_itranidxm(ictxt,icomm,flag,n,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag, n
      integer, intent(out) :: info
      integer              :: y(:,:), beta
      integer, target      :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_itranidxm
    subroutine psi_itranidxv(ictxt,icomm,flag,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag
      integer, intent(out) :: info
      integer              :: y(:), beta
      integer, target      :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_itranidxv
    subroutine psi_zswaptranm(flag,n,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag, n
      integer, intent(out) :: info
      complex(psb_dpk_)     :: y(:,:), beta
      complex(psb_dpk_),target   :: work(:)
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_zswaptranm
    subroutine psi_zswaptranv(flag,beta,y,desc_a,work,info,data)
      use psb_descriptor_type
      integer, intent(in)  :: flag
      integer, intent(out) :: info
      complex(psb_dpk_)     :: y(:), beta
      complex(psb_dpk_),target   :: work(:)       
      type(psb_desc_type), target  :: desc_a
      integer, optional    :: data
    end subroutine psi_zswaptranv
    subroutine psi_ztranidxm(ictxt,icomm,flag,n,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag, n
      integer, intent(out) :: info
      complex(psb_dpk_)     :: y(:,:), beta
      complex(psb_dpk_), target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_ztranidxm
    subroutine psi_ztranidxv(ictxt,icomm,flag,beta,y,idx,totxch,totsnd,totrcv,work,info)
      use psb_const_mod
      integer, intent(in)  :: ictxt,icomm,flag
      integer, intent(out) :: info
      complex(psb_dpk_)     :: y(:), beta
      complex(psb_dpk_), target :: work(:)
      integer, intent(in)  :: idx(:),totxch,totsnd,totrcv
    end subroutine psi_ztranidxv
  end interface

  interface
    subroutine psi_extract_dep_list(desc_data,desc_str,dep_list,&
         & length_dl,np,dl_lda,mode,info)
      integer :: np,dl_lda,mode, info
      integer :: desc_str(*),desc_data(*),dep_list(dl_lda,0:np),length_dl(0:np)
    end subroutine psi_extract_dep_list
  end interface
  interface psi_fnd_owner
    subroutine psi_fnd_owner(nv,idx,iprc,desc,info)
      use psb_descriptor_type
      integer, intent(in) :: nv
      integer, intent(in) ::  idx(:)
      integer, allocatable, intent(out) ::  iprc(:)
      type(psb_desc_type), intent(in) :: desc
      integer, intent(out) :: info
    end subroutine psi_fnd_owner
  end interface

  interface psi_ldsc_pre_halo
    subroutine psi_ldsc_pre_halo(desc,ext_hv,info)
      use psb_descriptor_type
      type(psb_desc_type), intent(inout) :: desc
      logical, intent(in)  :: ext_hv
      integer, intent(out) :: info
    end subroutine psi_ldsc_pre_halo
  end interface

  interface psi_bld_hash
    subroutine psi_bld_hash(desc,info)
      use psb_descriptor_type
      type(psb_desc_type), intent(inout) :: desc
      integer, intent(out) :: info
    end subroutine psi_bld_hash
  end interface

  interface psi_bld_tmphalo
    subroutine psi_bld_tmphalo(desc,info)
      use psb_descriptor_type
      type(psb_desc_type), intent(inout) :: desc
      integer, intent(out) :: info
    end subroutine psi_bld_tmphalo
  end interface


  interface psi_bld_tmpovrl
    subroutine psi_bld_tmpovrl(iv,desc,info)
      use psb_descriptor_type
      integer, intent(in)  :: iv(:)
      type(psb_desc_type), intent(inout) :: desc
      integer, intent(out) :: info
    end subroutine psi_bld_tmpovrl
  end interface


  interface psi_idx_cnv
    subroutine psi_idx_cnv1(nv,idxin,desc,info,mask,owned)
      use psb_descriptor_type
      integer, intent(in)    :: nv
      integer, intent(inout) ::  idxin(:)
      type(psb_desc_type), intent(in) :: desc
      integer, intent(out) :: info
      logical, intent(in), optional, target :: mask(:)
      logical, intent(in), optional :: owned
    end subroutine psi_idx_cnv1
    subroutine psi_idx_cnv2(nv,idxin,idxout,desc,info,mask,owned)
      use psb_descriptor_type
      integer, intent(in)  :: nv, idxin(:)
      integer, intent(out) :: idxout(:)
      type(psb_desc_type), intent(in) :: desc
      integer, intent(out) :: info
      logical, intent(in), optional, target :: mask(:)
      logical, intent(in), optional :: owned
    end subroutine psi_idx_cnv2
    subroutine psi_idx_cnvs(idxin,idxout,desc,info,mask,owned)
      use psb_descriptor_type
      integer, intent(in)  :: idxin
      integer, intent(out) :: idxout
      type(psb_desc_type), intent(in) :: desc
      integer, intent(out) :: info
      logical, intent(in), optional, target :: mask
      logical, intent(in), optional :: owned
    end subroutine psi_idx_cnvs
  end interface

  interface psi_idx_ins_cnv
    subroutine psi_idx_ins_cnv1(nv,idxin,desc,info,mask)
      use psb_descriptor_type
      integer, intent(in)    :: nv
      integer, intent(inout) ::  idxin(:)
      type(psb_desc_type), intent(inout) :: desc
      integer, intent(out) :: info
      logical, intent(in), optional, target :: mask(:)
    end subroutine psi_idx_ins_cnv1
    subroutine psi_idx_ins_cnv2(nv,idxin,idxout,desc,info,mask)
      use psb_descriptor_type
      integer, intent(in)  :: nv, idxin(:)
      integer, intent(out) :: idxout(:)
      type(psb_desc_type), intent(inout) :: desc
      integer, intent(out) :: info
      logical, intent(in), optional, target :: mask(:)
    end subroutine psi_idx_ins_cnv2
    subroutine psi_idx_ins_cnvs(idxin,idxout,desc,info,mask)
      use psb_descriptor_type
      integer, intent(in)  :: idxin
      integer, intent(out) :: idxout
      type(psb_desc_type), intent(inout) :: desc
      integer, intent(out) :: info
      logical, intent(in), optional, target :: mask
    end subroutine psi_idx_ins_cnvs
  end interface

  interface psi_cnv_dsc
    module procedure psi_cnv_dsc
  end interface

  interface psi_inner_cnv
    module procedure psi_inner_cnv1, psi_inner_cnv2
  end interface

  interface psi_ovrl_upd
    module procedure psi_iovrl_updr1, psi_iovrl_updr2,&
         & psi_dovrl_updr1, psi_dovrl_updr2, &
         & psi_zovrl_updr1, psi_zovrl_updr2
  end interface

  interface psi_ovrl_save
    module procedure  psi_iovrl_saver1, psi_iovrl_saver2,&
         & psi_dovrl_saver1, psi_dovrl_saver2,&
         & psi_zovrl_saver1, psi_zovrl_saver2
  end interface

  interface psi_ovrl_restore
    module procedure  psi_iovrl_restrr1, psi_iovrl_restrr2,&
         & psi_dovrl_restrr1, psi_dovrl_restrr2,&
         & psi_zovrl_restrr1, psi_zovrl_restrr2
  end interface


contains

  subroutine psi_cnv_dsc(halo_in,ovrlap_in,ext_in,cdesc, info)

    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    use psb_descriptor_type
    use psb_realloc_mod
    implicit none

    !     ....scalars parameters....
    integer, intent(in)                :: halo_in(:), ovrlap_in(:),ext_in(:)
    type(psb_desc_type), intent(inout) :: cdesc
    integer, intent(out)               :: info

    !     ....local scalars....      
    integer  :: np,me
    integer  :: ictxt, err_act,nxch,nsnd,nrcv,j,k
    !     ...local array...
    integer, allocatable  :: idx_out(:), tmp_mst_idx(:)

    !     ...parameters
    integer :: debug_level, debug_unit
    logical, parameter :: debug=.false.
    character(len=20)  :: name

    name='psi_bld_cdesc'
    call psb_get_erraction(err_act)
    debug_level = psb_get_debug_level()
    debug_unit  = psb_get_debug_unit()

    info = 0
    ictxt = cdesc%matrix_data(psb_ctxt_)

    call psb_info(ictxt,me,np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif


    ! first the halo index
    if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on halo',&
         & size(halo_in)
    call psi_crea_index(cdesc,halo_in, idx_out,.false.,nxch,nsnd,nrcv,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_index')
      goto 9999
    end if
    call psb_transfer(idx_out,cdesc%halo_index,info)
    cdesc%matrix_data(psb_thal_xch_) = nxch
    cdesc%matrix_data(psb_thal_snd_) = nsnd
    cdesc%matrix_data(psb_thal_rcv_) = nrcv 

    if (debug_level>0) write(debug_unit,*) me,'Done crea_index on halo'
    if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ext'


    ! then ext index
    if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ext'
    call psi_crea_index(cdesc,ext_in, idx_out,.false.,nxch,nsnd,nrcv,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_index')
      goto 9999
    end if
    call psb_transfer(idx_out,cdesc%ext_index,info)
    cdesc%matrix_data(psb_text_xch_) = nxch
    cdesc%matrix_data(psb_text_snd_) = nsnd
    cdesc%matrix_data(psb_text_rcv_) = nrcv 

    if (debug_level>0) write(debug_unit,*) me,'Done crea_index on ext'
    if (debug_level>0) write(debug_unit,*) me,'Calling crea_index on ovrlap'

    ! then the overlap index
    call psi_crea_index(cdesc,ovrlap_in, idx_out,.true.,nxch,nsnd,nrcv,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_index')
      goto 9999
    end if
    call psb_transfer(idx_out,cdesc%ovrlap_index,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psb_transfer')
      goto 9999
    end if

    cdesc%matrix_data(psb_tovr_xch_) = nxch
    cdesc%matrix_data(psb_tovr_snd_) = nsnd
    cdesc%matrix_data(psb_tovr_rcv_) = nrcv 

    ! next  ovrlap_elem 
    if (debug_level>0) write(debug_unit,*) me,'Calling crea_ovr_elem'
    call psi_crea_ovr_elem(me,cdesc%ovrlap_index,cdesc%ovrlap_elem,info)
    if (debug_level>0) write(debug_unit,*) me,'Done crea_ovr_elem'
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_ovr_elem')
      goto 9999
    end if
    ! Extract ovr_mst_idx from ovrlap_elem 
    if (debug_level>0) write(debug_unit,*) me,'Calling bld_ovr_mst'
    call psi_bld_ovr_mst(me,cdesc%ovrlap_elem,tmp_mst_idx,info)
    if (info == 0) call psi_crea_index(cdesc,&
         & tmp_mst_idx,idx_out,.false.,nxch,nsnd,nrcv,info)
    if (debug_level>0) write(debug_unit,*) me,'Done crea_indx'
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_bld_ovr_mst')
      goto 9999
    end if
    call psb_transfer(idx_out,cdesc%ovr_mst_idx,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psb_transfer')
      goto 9999
    end if

    cdesc%matrix_data(psb_tmov_xch_) = nxch
    cdesc%matrix_data(psb_tmov_snd_) = nsnd
    cdesc%matrix_data(psb_tmov_rcv_) = nrcv 

    ! finally bnd_elem
    call psi_crea_bnd_elem(idx_out,cdesc,info)
    if (info == 0) call psb_transfer(idx_out,cdesc%bnd_elem,info)

    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_bnd_elem')
      goto 9999
    end if
    if (debug_level>0) write(debug_unit,*) me,'Done crea_bnd_elem'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return

  end subroutine psi_cnv_dsc



  subroutine psi_inner_cnv1(n,x,hashmask,hashv,glb_lc)
    integer, intent(in)    :: n,hashmask,hashv(0:),glb_lc(:,:)
    integer, intent(inout) :: x(:)

    integer :: i, ih, key, idx,nh,tmp,lb,ub,lm
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a binary search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !

    do i=1, n
      key = x(i) 
      ih  = iand(key,hashmask)
      idx = hashv(ih)
      nh  = hashv(ih+1) - hashv(ih) 
      if (nh > 0) then 
        tmp = -1 
        lb = idx
        ub = idx+nh-1
        do 
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key==glb_lc(lm,1)) then 
            tmp = lm
            exit
          else if (key<glb_lc(lm,1)) then 
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      else 
        tmp = -1
      end if
      if (tmp > 0) then 
        x(i) = glb_lc(tmp,2)
      else         
        x(i) = tmp 
      end if
    end do
  end subroutine psi_inner_cnv1


  subroutine psi_inner_cnv2(n,x,y,hashmask,hashv,glb_lc)
    integer, intent(in)  :: n, hashmask,hashv(0:),glb_lc(:,:)
    integer, intent(in)  :: x(:)
    integer, intent(out) :: y(:)

    integer :: i, ih, key, idx,nh,tmp,lb,ub,lm
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a binary search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !

    do i=1, n
      key = x(i) 
      ih  = iand(key,hashmask)
      if (ih > ubound(hashv,1) ) then 
        write(0,*) ' In inner cnv: ',ih,ubound(hashv)
      end if
      idx = hashv(ih)
      nh  = hashv(ih+1) - hashv(ih) 
      if (nh > 0) then 
        tmp = -1 
        lb = idx
        ub = idx+nh-1
        do 
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key==glb_lc(lm,1)) then 
            tmp = lm
            exit
          else if (key<glb_lc(lm,1)) then 
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      else 
        tmp = -1
      end if
      if (tmp > 0) then 
        y(i) = glb_lc(tmp,2)
      else         
        y(i) = tmp 
      end if
    end do
  end subroutine psi_inner_cnv2

  subroutine  psi_dovrl_updr1(x,desc_a,update,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    real(psb_dpk_), intent(inout), target :: x(:)
    type(psb_desc_type), intent(in)         :: desc_a
    integer, intent(in)                     :: update
    integer, intent(out)                    :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_dovrl_updr1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    ! switch on update type
    select case (update)
    case(psb_square_root_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx) = x(idx)/sqrt(real(ndm))
      end do
    case(psb_avg_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx) = x(idx)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        if (me /= desc_a%ovrlap_elem(i,3))&
             & x(idx) = dzero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = 70
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_dovrl_updr1


  subroutine  psi_dovrl_updr2(x,desc_a,update,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    real(psb_dpk_), intent(inout), target :: x(:,:)
    type(psb_desc_type), intent(in)         :: desc_a
    integer, intent(in)                     :: update
    integer, intent(out)                    :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_dovrl_updr2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    ! switch on update type
    select case (update)
    case(psb_square_root_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx,:) = x(idx,:)/sqrt(real(ndm))
      end do
    case(psb_avg_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx,:) = x(idx,:)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        if (me /= desc_a%ovrlap_elem(i,3))&
             & x(idx,:) = dzero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = 70
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_dovrl_updr2

  subroutine  psi_zovrl_updr1(x,desc_a,update,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    complex(psb_dpk_), intent(inout), target :: x(:)
    type(psb_desc_type), intent(in)         :: desc_a
    integer, intent(in)                     :: update
    integer, intent(out)                    :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_zovrl_updr1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    ! switch on update type
    select case (update)
    case(psb_square_root_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx) = x(idx)/sqrt(real(ndm))
      end do
    case(psb_avg_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx) = x(idx)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        if (me /= desc_a%ovrlap_elem(i,3))&
             & x(idx) = zzero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = 70
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_zovrl_updr1


  subroutine  psi_zovrl_updr2(x,desc_a,update,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    complex(psb_dpk_), intent(inout), target :: x(:,:)
    type(psb_desc_type), intent(in)         :: desc_a
    integer, intent(in)                     :: update
    integer, intent(out)                    :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_zovrl_updr2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    ! switch on update type
    select case (update)
    case(psb_square_root_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx,:) = x(idx,:)/sqrt(real(ndm))
      end do
    case(psb_avg_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx,:) = x(idx,:)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        if (me /= desc_a%ovrlap_elem(i,3))&
             & x(idx,:) = zzero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = 70
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_zovrl_updr2

  subroutine  psi_iovrl_updr1(x,desc_a,update,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    integer, intent(inout), target   :: x(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(in)              :: update
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_iovrl_updr1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    ! switch on update type
    select case (update)
      ! Square root does not make sense here
!!$    case(psb_square_root_)
!!$      do i=1,size(desc_a%ovrlap_elem,1)
!!$        idx = desc_a%ovrlap_elem(i,1)
!!$        ndm = desc_a%ovrlap_elem(i,2)
!!$        x(idx) = x(idx)/sqrt(real(ndm))
!!$      end do
    case(psb_avg_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx) = x(idx)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        if (me /= desc_a%ovrlap_elem(i,3))&
             & x(idx) = izero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = 70
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_iovrl_updr1


  subroutine  psi_iovrl_updr2(x,desc_a,update,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    integer, intent(inout), target   :: x(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(in)              :: update
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_iovrl_updr2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    ! switch on update type
    select case (update)
      ! Square root does not make sense here
!!$    case(psb_square_root_)
!!$      do i=1,size(desc_a%ovrlap_elem,1)
!!$        idx = desc_a%ovrlap_elem(i,1)
!!$        ndm = desc_a%ovrlap_elem(i,2)
!!$        x(idx,:) = x(idx,:)/sqrt(real(ndm))
!!$      end do
    case(psb_avg_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        ndm = desc_a%ovrlap_elem(i,2)
        x(idx,:) = x(idx,:)/real(ndm)
      end do
    case(psb_setzero_)
      do i=1,size(desc_a%ovrlap_elem,1)
        idx = desc_a%ovrlap_elem(i,1)
        if (me /= desc_a%ovrlap_elem(i,3))&
             & x(idx,:) = izero
      end do
    case(psb_sum_)
      ! do nothing

    case default 
      ! wrong value for choice argument
      info = 70
      call psb_errpush(info,name,i_err=(/3,update,0,0,0/))
      goto 9999
    end select

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_iovrl_updr2



  subroutine  psi_dovrl_saver1(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_penv_mod
    implicit none

    real(psb_dpk_), intent(inout)  :: x(:)
    real(psb_dpk_), allocatable    :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_dovrl_saver1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)
    call psb_realloc(isz,xs,info) 
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    endif

    do i=1, isz
      idx   = desc_a%ovrlap_elem(i,1)
      xs(i) = x(idx)
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_dovrl_saver1

  subroutine  psi_dovrl_restrr1(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    real(psb_dpk_), intent(inout)  :: x(:)
    real(psb_dpk_)                 :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_dovrl_restrr1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)

    do i=1, isz
      idx    = desc_a%ovrlap_elem(i,1)
      x(idx) = xs(i) 
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_dovrl_restrr1


  subroutine  psi_dovrl_saver2(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_penv_mod
    implicit none

    real(psb_dpk_), intent(inout)  :: x(:,:)
    real(psb_dpk_), allocatable    :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz, nc
    character(len=20) :: name, ch_err

    name='psi_dovrl_saver2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)
    nc  = size(x,2)
    call psb_realloc(isz,nc,xs,info) 
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    endif

    do i=1, isz
      idx     = desc_a%ovrlap_elem(i,1)
      xs(i,:) = x(idx,:)
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_dovrl_saver2

  subroutine  psi_dovrl_restrr2(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    real(psb_dpk_), intent(inout)  :: x(:,:)
    real(psb_dpk_)                 :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_dovrl_restrr2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif
    
    if (size(x,2) /= size(xs,2)) then 
      info = 4001
      call psb_errpush(info,name, a_err='Mismacth columns X vs XS')
      goto 9999
    endif
      
    
    isz = size(desc_a%ovrlap_elem,1)

    do i=1, isz
      idx      = desc_a%ovrlap_elem(i,1)
      x(idx,:) = xs(i,:) 
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_dovrl_restrr2



  subroutine  psi_zovrl_saver1(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_penv_mod
    implicit none

    complex(psb_dpk_), intent(inout)  :: x(:)
    complex(psb_dpk_), allocatable    :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_zovrl_saver1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)
    call psb_realloc(isz,xs,info) 
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    endif

    do i=1, isz
      idx   = desc_a%ovrlap_elem(i,1)
      xs(i) = x(idx)
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_zovrl_saver1

  subroutine  psi_zovrl_restrr1(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none
    
    complex(psb_dpk_), intent(inout)  :: x(:)
    complex(psb_dpk_)                 :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_zovrl_restrr1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)

    do i=1, isz
      idx    = desc_a%ovrlap_elem(i,1)
      x(idx) = xs(i) 
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_zovrl_restrr1


  subroutine  psi_zovrl_saver2(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_penv_mod
    implicit none

    complex(psb_dpk_), intent(inout)  :: x(:,:)
    complex(psb_dpk_), allocatable    :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz, nc
    character(len=20) :: name, ch_err

    name='psi_zovrl_saver2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)
    nc  = size(x,2)
    call psb_realloc(isz,nc,xs,info) 
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    endif

    do i=1, isz
      idx     = desc_a%ovrlap_elem(i,1)
      xs(i,:) = x(idx,:)
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_zovrl_saver2

  subroutine  psi_zovrl_restrr2(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    complex(psb_dpk_), intent(inout)  :: x(:,:)
    complex(psb_dpk_)                 :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_zovrl_restrr2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif
    
    if (size(x,2) /= size(xs,2)) then 
      info = 4001
      call psb_errpush(info,name, a_err='Mismacth columns X vs XS')
      goto 9999
    endif
      
    
    isz = size(desc_a%ovrlap_elem,1)

    do i=1, isz
      idx      = desc_a%ovrlap_elem(i,1)
      x(idx,:) = xs(i,:) 
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_zovrl_restrr2


  subroutine  psi_iovrl_saver1(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_penv_mod
    implicit none

    integer, intent(inout)  :: x(:)
    integer, allocatable    :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_iovrl_saver1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)
    call psb_realloc(isz,xs,info) 
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    endif

    do i=1, isz
      idx   = desc_a%ovrlap_elem(i,1)
      xs(i) = x(idx)
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_iovrl_saver1

  subroutine  psi_iovrl_restrr1(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    integer, intent(inout)  :: x(:)
    integer                 :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_iovrl_restrr1'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)

    do i=1, isz
      idx    = desc_a%ovrlap_elem(i,1)
      x(idx) = xs(i) 
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_iovrl_restrr1


  subroutine  psi_iovrl_saver2(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_penv_mod
    implicit none

    integer, intent(inout)  :: x(:,:)
    integer, allocatable    :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz, nc
    character(len=20) :: name, ch_err

    name='psi_iovrl_saver2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif

    isz = size(desc_a%ovrlap_elem,1)
    nc  = size(x,2)
    call psb_realloc(isz,nc,xs,info) 
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    endif

    do i=1, isz
      idx     = desc_a%ovrlap_elem(i,1)
      xs(i,:) = x(idx,:)
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_iovrl_saver2

  subroutine  psi_iovrl_restrr2(x,xs,desc_a,info)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    implicit none

    integer, intent(inout)  :: x(:,:)
    integer                 :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_iovrl_restrr2'
    if (psb_get_errstatus() /= 0) return 
    info = 0
    call psb_erractionsave(err_act)
    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif
    
    if (size(x,2) /= size(xs,2)) then 
      info = 4001
      call psb_errpush(info,name, a_err='Mismacth columns X vs XS')
      goto 9999
    endif
      
    
    isz = size(desc_a%ovrlap_elem,1)

    do i=1, isz
      idx      = desc_a%ovrlap_elem(i,1)
      x(idx,:) = xs(i,:) 
    end do

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return
  end subroutine psi_iovrl_restrr2


!!$  subroutine psi_dgthzm(n,k,idx,x,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, k, idx(:)
!!$    real(psb_dpk_) :: x(:,:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i, j, pt
!!$
!!$    pt=0
!!$    do j=1,k
!!$      do i=1,n
!!$        pt=pt+1
!!$        y(pt)=x(idx(i),j)
!!$      end do
!!$    end do
!!$
!!$  end subroutine psi_dgthzm
!!$
!!$  subroutine psi_dgthzv(n,idx,x,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, idx(:)
!!$    real(psb_dpk_) :: x(:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i
!!$
!!$    do i=1,n
!!$      y(i)=x(idx(i))
!!$    end do
!!$
!!$  end subroutine psi_dgthzv
!!$
!!$
!!$  subroutine psi_dsctm(n,k,idx,x,beta,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, k, idx(:)
!!$    real(psb_dpk_) :: beta, x(:), y(:,:)
!!$
!!$    ! Locals
!!$    integer :: i, j, pt
!!$
!!$    if (beta == dzero) then
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = x(pt)
!!$        end do
!!$      end do
!!$    else if (beta == done) then
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = y(idx(i),j)+x(pt)
!!$        end do
!!$      end do
!!$    else
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = beta*y(idx(i),j)+x(pt)
!!$        end do
!!$      end do
!!$    end if
!!$  end subroutine psi_dsctm
!!$
!!$  subroutine psi_dsctv(n,idx,x,beta,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, idx(:)
!!$    real(psb_dpk_) :: beta, x(:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i
!!$
!!$    if (beta == dzero) then
!!$      do i=1,n
!!$        y(idx(i)) = x(i)
!!$      end do
!!$    else if (beta == done) then
!!$      do i=1,n
!!$        y(idx(i)) = y(idx(i))+x(i)
!!$      end do
!!$    else
!!$      do i=1,n
!!$        y(idx(i)) = beta*y(idx(i))+x(i)
!!$      end do
!!$    end if
!!$  end subroutine psi_dsctv
!!$
!!$
!!$  subroutine psi_igthzm(n,k,idx,x,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, k, idx(:)
!!$    integer :: x(:,:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i, j, pt
!!$
!!$    pt=0
!!$    do j=1,k
!!$      do i=1,n
!!$        pt=pt+1
!!$        y(pt)=x(idx(i),j)
!!$      end do
!!$    end do
!!$
!!$  end subroutine psi_igthzm
!!$
!!$
!!$  subroutine psi_igthzv(n,idx,x,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, idx(:)
!!$    integer :: x(:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i
!!$
!!$    do i=1,n
!!$      y(i)=x(idx(i))
!!$    end do
!!$
!!$  end subroutine psi_igthzv
!!$
!!$
!!$
!!$  subroutine psi_isctm(n,k,idx,x,beta,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, k, idx(:)
!!$    integer :: beta, x(:), y(:,:)
!!$
!!$    ! Locals
!!$    integer :: i, j, pt
!!$
!!$    if (beta == izero) then
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = x(pt)
!!$        end do
!!$      end do
!!$    else if (beta == ione) then
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = y(idx(i),j)+x(pt)
!!$        end do
!!$      end do
!!$    else
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = beta*y(idx(i),j)+x(pt)
!!$        end do
!!$      end do
!!$    end if
!!$  end subroutine psi_isctm
!!$
!!$  subroutine psi_isctv(n,idx,x,beta,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, idx(:)
!!$    integer :: beta, x(:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i
!!$
!!$    if (beta == izero) then
!!$      do i=1,n
!!$        y(idx(i)) = x(i)
!!$      end do
!!$    else if (beta == ione) then
!!$      do i=1,n
!!$        y(idx(i)) = y(idx(i))+x(i)
!!$      end do
!!$    else
!!$      do i=1,n
!!$        y(idx(i)) = beta*y(idx(i))+x(i)
!!$      end do
!!$    end if
!!$  end subroutine psi_isctv
!!$
!!$
!!$  subroutine psi_zgthzm(n,k,idx,x,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, k, idx(:)
!!$    complex(psb_dpk_) :: x(:,:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i, j, pt
!!$
!!$    pt=0
!!$    do j=1,k
!!$      do i=1,n
!!$        pt=pt+1
!!$        y(pt)=x(idx(i),j)
!!$      end do
!!$    end do
!!$
!!$  end subroutine psi_zgthzm
!!$
!!$
!!$  subroutine psi_zgthzv(n,idx,x,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, idx(:)
!!$    complex(psb_dpk_) :: x(:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i
!!$
!!$    do i=1,n
!!$      y(i)=x(idx(i))
!!$    end do
!!$
!!$  end subroutine psi_zgthzv
!!$
!!$  subroutine psi_zsctm(n,k,idx,x,beta,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, k, idx(:)
!!$    complex(psb_dpk_) :: beta, x(:), y(:,:)
!!$
!!$    ! Locals
!!$    integer :: i, j, pt
!!$
!!$    if (beta == zzero) then
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = x(pt)
!!$        end do
!!$      end do
!!$    else if (beta == zone) then
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = y(idx(i),j)+x(pt)
!!$        end do
!!$      end do
!!$    else
!!$      pt=0
!!$      do j=1,k
!!$        do i=1,n
!!$          pt=pt+1
!!$          y(idx(i),j) = beta*y(idx(i),j)+x(pt)
!!$        end do
!!$      end do
!!$    end if
!!$  end subroutine psi_zsctm
!!$
!!$
!!$  subroutine psi_zsctv(n,idx,x,beta,y)
!!$
!!$    use psb_const_mod
!!$    implicit none
!!$
!!$    integer :: n, idx(:)
!!$    complex(psb_dpk_) :: beta, x(:), y(:)
!!$
!!$    ! Locals
!!$    integer :: i
!!$
!!$    if (beta == zzero) then
!!$      do i=1,n
!!$        y(idx(i)) = x(i)
!!$      end do
!!$    else if (beta == zone) then
!!$      do i=1,n
!!$        y(idx(i)) = y(idx(i))+x(i)
!!$      end do
!!$    else
!!$      do i=1,n
!!$        y(idx(i)) = beta*y(idx(i))+x(i)
!!$      end do
!!$    end if
!!$  end subroutine psi_zsctv
  
  subroutine psi_bld_ovr_mst(me,ovrlap_elem,mst_idx,info)
    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    use psb_descriptor_type
    use psb_realloc_mod
    implicit none

    !     ....scalars parameters....
    integer, intent(in)               :: me, ovrlap_elem(:,:)
    integer, allocatable, intent(out) :: mst_idx(:) 
    integer, intent(out)              :: info

    integer  :: i, j, proc, nov,isz, ip, err_act, idx
    character(len=20)  :: name

    name='psi_bld_ovr_mst'
    call psb_get_erraction(err_act)

    nov = size(ovrlap_elem,1)
    isz = 3*nov+1
    call psb_realloc(isz,mst_idx,info) 
    if (info /= 0) then
      call psb_errpush(4001,name,a_err='reallocate')
      goto 9999
    end if
    mst_idx = -1
    j = 1
    do i=1, nov
      proc = ovrlap_elem(i,3)
      if (me /= proc) then 
        idx = ovrlap_elem(i,1)
        mst_idx(j+0) = proc
        mst_idx(j+1) = 1
        mst_idx(j+2) = idx
        j = j + 3
      end if
    end do
    mst_idx(j) = -1 

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psi_bld_ovr_mst

end module psi_mod
