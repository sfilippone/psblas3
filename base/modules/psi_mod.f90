!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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
  use psi_gthsct_mod

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
     subroutine psi_crea_ovr_elem(desc_overlap,ovr_elem,info)
       integer :: desc_overlap(:)
       integer, allocatable, intent(inout) :: ovr_elem(:)
       integer, intent(out)                :: info
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
       real(kind(1.d0))     :: y(:,:), beta
       real(kind(1.d0)),target      :: work(:)
       type(psb_desc_type), target  :: desc_a
       integer, optional    :: data
     end subroutine psi_dswapdatam
     subroutine psi_dswapdatav(flag,beta,y,desc_a,work,info,data)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       real(kind(1.d0))     :: y(:), beta 
       real(kind(1.d0)),target      :: work(:)
       type(psb_desc_type), target  :: desc_a
       integer, optional    :: data
     end subroutine psi_dswapdatav
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
     subroutine psi_zswapdatam(flag,n,beta,y,desc_a,work,info,data)
       use psb_descriptor_type
       integer, intent(in)  :: flag, n
       integer, intent(out) :: info
       complex(kind(1.d0))     :: y(:,:), beta
       complex(kind(1.d0)),target   :: work(:)
       type(psb_desc_type), target  :: desc_a
       integer, optional    :: data
     end subroutine psi_zswapdatam
     subroutine psi_zswapdatav(flag,beta,y,desc_a,work,info,data)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       complex(kind(1.d0))     :: y(:), beta
       complex(kind(1.d0)),target   :: work(:)
       type(psb_desc_type), target  :: desc_a
       integer, optional    :: data
     end subroutine psi_zswapdatav
  end interface


  interface psi_swaptran
     subroutine psi_dswaptranm(flag,n,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag, n
       integer, intent(out) :: info
       real(kind(1.d0))     :: y(:,:), beta
       real(kind(1.d0)),target     :: work(:)
       type(psb_desc_type), target  :: desc_a
     end subroutine psi_dswaptranm
     subroutine psi_dswaptranv(flag,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       real(kind(1.d0))     :: y(:), beta
       real(kind(1.d0)),target     :: work(:)
       type(psb_desc_type), target  :: desc_a
     end subroutine psi_dswaptranv
     subroutine psi_iswaptranm(flag,n,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag, n
       integer, intent(out) :: info
       integer              :: y(:,:), beta
       integer,target               :: work(:)
       type(psb_desc_type), target  :: desc_a
     end subroutine psi_iswaptranm
     subroutine psi_iswaptranv(flag,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       integer              :: y(:), beta
       integer,target               :: work(:)
       type(psb_desc_type), target  :: desc_a
     end subroutine psi_iswaptranv
     subroutine psi_zswaptranm(flag,n,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag, n
       integer, intent(out) :: info
       complex(kind(1.d0))     :: y(:,:), beta
       complex(kind(1.d0)),target   :: work(:)
       type(psb_desc_type), target  :: desc_a
     end subroutine psi_zswaptranm
     subroutine psi_zswaptranv(flag,beta,y,desc_a,work,info)
       use psb_descriptor_type
       integer, intent(in)  :: flag
       integer, intent(out) :: info
       complex(kind(1.d0))     :: y(:), beta
       complex(kind(1.d0)),target   :: work(:)       
       type(psb_desc_type), target  :: desc_a
     end subroutine psi_zswaptranv
   end interface

  interface psi_cnv_dsc
    module procedure psi_cnv_dsc
  end interface

  interface psi_inner_cnv
    module procedure psi_inner_cnv1, psi_inner_cnv2
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

contains
  
  subroutine psi_cnv_dsc(halo_in,ovrlap_in,ext_in,cdesc, info)

    use psb_const_mod
    use psb_error_mod
    use psb_penv_mod
    use psb_descriptor_type
    use psb_realloc_mod
    implicit none

    !     ....scalars parameters....
    integer, intent(in)  :: halo_in(:), ovrlap_in(:),ext_in(:)
    type(psb_desc_type), intent(inout) :: cdesc
    integer, intent(out)  :: info

    !     ....local scalars....      
    integer  :: i,np,me,proc, max_index
    integer  :: ictxt, err_act,nxch,nsnd,nrcv
    !     ...local array...
    integer  :: int_err(5)
    integer, allocatable  :: idx_out(:)

    !     ...parameters
    logical, parameter :: debug=.false.
    character(len=20)  :: name

    name='psi_bld_cdesc'
    call psb_get_erraction(err_act)

    info = 0
    ictxt = cdesc%matrix_data(psb_ctxt_)

    call psb_info(ictxt,me,np)
    if (np == -1) then
      info = 2010
      call psb_errpush(info,name)
      goto 9999
    endif


    ! first the halo index
    if (debug) write(0,*) me,'Calling crea_index on halo'
    call psi_crea_index(cdesc,halo_in, idx_out,.false.,nxch,nsnd,nrcv,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_index')
      goto 9999
    end if
    call psb_transfer(idx_out,cdesc%halo_index,info)
    cdesc%matrix_data(psb_thal_xch_) = nxch
    cdesc%matrix_data(psb_thal_snd_) = nsnd
    cdesc%matrix_data(psb_thal_rcv_) = nrcv 
    
    if (debug) write(0,*) me,'Done crea_index on halo'
    if (debug) write(0,*) me,'Calling crea_index on ext'


    ! then ext index
    if (debug) write(0,*) me,'Calling crea_index on ext'
    call psi_crea_index(cdesc,ext_in, idx_out,.false.,nxch,nsnd,nrcv,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_index')
      goto 9999
    end if
    call psb_transfer(idx_out,cdesc%ext_index,info)
    cdesc%matrix_data(psb_text_xch_) = nxch
    cdesc%matrix_data(psb_text_snd_) = nsnd
    cdesc%matrix_data(psb_text_rcv_) = nrcv 
    
    if (debug) write(0,*) me,'Done crea_index on ext'
    if (debug) write(0,*) me,'Calling crea_index on ovrlap'

    ! then the overlap index

    call psi_crea_index(cdesc,ovrlap_in, idx_out,.true.,nxch,nsnd,nrcv,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_index')
      goto 9999
    end if
    call psb_transfer(idx_out,cdesc%ovrlap_index,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_transfer')
      goto 9999
    end if

    cdesc%matrix_data(psb_tovr_xch_) = nxch
    cdesc%matrix_data(psb_tovr_snd_) = nsnd
    cdesc%matrix_data(psb_tovr_rcv_) = nrcv 

    if (debug) write(0,*) me,'Calling crea_ovr_elem'
    ! next  ovrlap_elem 
    call psi_crea_ovr_elem(cdesc%ovrlap_index,cdesc%ovrlap_elem,info)
    if (debug) write(0,*) me,'Done crea_ovr_elem'
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_ovr_elem')
      goto 9999
    end if

    ! finally bnd_elem
    call psi_crea_bnd_elem(idx_out,cdesc,info)
    if (info == 0) call psb_transfer(idx_out,cdesc%bnd_elem,info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psi_crea_bnd_elem')
      goto 9999
    end if
    if (debug) write(0,*) me,'Done crea_bnd_elem'

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



  subroutine psi_inner_cnv1(n,x,hashsize,hashmask,hashv,glb_lc)
    integer, intent(in)    :: n, hashsize,hashmask,hashv(0:),glb_lc(:,:)
    integer, intent(inout) :: x(:)

    integer :: i, ih, key, idx,nh,tmp,lb,ub,lm
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


  subroutine psi_inner_cnv2(n,x,y,hashsize,hashmask,hashv,glb_lc)
    integer, intent(in)  :: n, hashsize,hashmask,hashv(0:),glb_lc(:,:)
    integer, intent(in)  :: x(:)
    integer, intent(out) :: y(:)

    integer :: i, ih, key, idx,nh,tmp,lb,ub,lm
    
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

end module psi_mod
