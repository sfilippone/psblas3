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
!
!
! package: psb_hash_map_mod
!    Defines the HASH_MAP type.
!
! This is the index map of choice for large index spaces. 
!  If the global index space is very large (larger than the threshold value
!  which may be set by the user), then it is not advisable to have a full
!  GLOB_TO_LOC array; therefore we only record the global indices that do have a 
!  local counterpart, so that the local storage will be proportional to 
!  N_COL.
!  The idea is that  glb_lc(:,1) will hold sorted global indices, and
!  glb_lc(:,2) the corresponding local indices, so that we may do a binary search.
!  To cut down  the search time we partition glb_lc into a set of lists
!  addressed by  hashv(:) based on the value of the lowest
!  PSB_HASH_BITS bits of the  global index. 
!  During the build phase glb_lc() will store the indices of the internal points,
!  i.e. local indices 1:NROW, since those are known ad CDALL time.
!  The halo indices that we encounter during the build phase are put in
!  a PSB_HASH_TYPE data structure, which implements a very simple hash; this
!  hash  will nonetheless be quite fast at low occupancy rates.
!  At assembly time, we move everything into hashv(:) and glb_lc(:,:).
!
module psb_hash_map_mod
  use psb_const_mod
  use psb_desc_const_mod
  use psb_indx_map_mod
  use psb_hash_mod 
  use psb_penv_mod
  use psb_sort_mod
  use psb_realloc_mod
  use psb_error_mod
  
  type, extends(psb_indx_map) :: psb_hash_map

    integer(psb_lpk_) :: hashvsize, hashvmask
    integer(psb_ipk_), allocatable :: hashv(:)
    integer(psb_lpk_), allocatable :: glb_lc(:,:), loc_to_glob(:)
    type(psb_hash_type)  :: hash 

  contains

    procedure, pass(idxmap)  :: init_vl    => hash_init_vl
    procedure, pass(idxmap)  :: hash_map_init => hash_init_vg

    procedure, pass(idxmap)  :: sizeof    => hash_sizeof
    procedure, pass(idxmap)  :: asb       => hash_asb
    procedure, pass(idxmap)  :: free      => hash_free
    procedure, pass(idxmap)  :: clone     => hash_clone
    procedure, pass(idxmap)  :: reinit    => hash_reinit
    procedure, nopass        :: get_fmt   => hash_get_fmt

    procedure, nopass        :: row_extendable => hash_row_extendable

    procedure, pass(idxmap)  :: ll2gs1     => hash_l2gs1
    procedure, pass(idxmap)  :: ll2gs2     => hash_l2gs2
    procedure, pass(idxmap)  :: ll2gv1     => hash_l2gv1
    procedure, pass(idxmap)  :: ll2gv2     => hash_l2gv2

    procedure, pass(idxmap)  :: lg2ls1     => hash_g2ls1
    procedure, pass(idxmap)  :: lg2ls2     => hash_g2ls2
    procedure, pass(idxmap)  :: lg2lv1     => hash_g2lv1
    procedure, pass(idxmap)  :: lg2lv2     => hash_g2lv2

    procedure, pass(idxmap)  :: lg2ls1_ins => hash_g2ls1_ins
    procedure, pass(idxmap)  :: lg2ls2_ins => hash_g2ls2_ins
    procedure, pass(idxmap)  :: lg2lv1_ins => hash_g2lv1_ins
    procedure, pass(idxmap)  :: lg2lv2_ins => hash_g2lv2_ins

    procedure, pass(idxmap)  :: bld_g2l_map => hash_bld_g2l_map

  end type psb_hash_map

  private :: hash_init_vl, hash_init_vg, hash_sizeof, hash_asb, &
       & hash_free, hash_get_fmt, hash_l2gs1, hash_l2gs2, &
       & hash_l2gv1, hash_l2gv2, hash_g2ls1, hash_g2ls2, &
       & hash_g2lv1, hash_g2lv2, hash_g2ls1_ins, hash_g2ls2_ins, &
       & hash_g2lv1_ins, hash_g2lv2_ins, hash_init_vlu, &
       & hash_bld_g2l_map, hash_inner_cnvs2, hash_inner_cnvs1, &
       & hash_inner_cnv2, hash_inner_cnv1, hash_row_extendable 

  integer(psb_ipk_), private :: psb_laddsz=500

  interface hash_inner_cnv 
    module procedure  hash_inner_cnvs2, hash_inner_cnv2,&
         &  hash_inner_cnvs1,  hash_inner_cnv1
  end interface hash_inner_cnv
  private :: hash_inner_cnv
  interface hash_srch 
#if defined(PSB_IPK4) && defined(PSB_LPK8)
    module procedure  hash_srch_ipk, hash_srch_lpk
#else
    module procedure hash_srch_ipk
#endif
  end interface hash_srch
  private :: hash_srch
  integer, parameter, private :: seqsrchmax=6

contains

  function hash_row_extendable() result(val)
    implicit none 
    logical :: val
    val = .true.
  end function hash_row_extendable

  function hash_sizeof(idxmap) result(val)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_epk_) :: val

    val = idxmap%psb_indx_map%sizeof() 
    val = val + 2 * psb_sizeof_ip
    if (allocated(idxmap%hashv)) &
         & val = val + size(idxmap%hashv)*psb_sizeof_ip
    if (allocated(idxmap%glb_lc)) &
         & val = val + size(idxmap%glb_lc)*psb_sizeof_lp
    if (allocated(idxmap%loc_to_glob)) &
         & val = val + size(idxmap%loc_to_glob)*psb_sizeof_lp
    val = val + psb_sizeof(idxmap%hash)

  end function hash_sizeof


  subroutine hash_free(idxmap)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_) :: info

    if (allocated(idxmap%hashv)) &
         & deallocate(idxmap%hashv)
    if (allocated(idxmap%glb_lc)) &
         & deallocate(idxmap%glb_lc)
    call psb_free(idxmap%hash,info) 
    call idxmap%psb_indx_map%free()

  end subroutine hash_free


  subroutine hash_l2gs1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_lpk_) :: idxv(1)
    info = 0
    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idx
    call idxmap%l2gip(idxv,info,owned=owned)
    idx = idxv(1)

  end subroutine hash_l2gs1

  subroutine hash_l2gs2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_lpk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    integer(psb_lpk_) :: idxv(1)
    info = 0
    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idxin
    call idxmap%l2gip(idxv,info,owned=owned)
    idxout = idxv(1)


  end subroutine hash_l2gs2


  subroutine hash_l2gv1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: i
    logical :: owned_
    info = 0
    if (present(mask)) then 
      if (size(mask) < size(idx)) then 
        info = -1
        return
      end if
    end if
    if (present(owned)) then 
      owned_ = owned
    else
      owned_ = .false.
    end if

    if (present(mask)) then 

      ! $ o m p parallel do default(none) schedule(dynamic) &
      ! $ o m p shared(mask,idx,idxmap,owned_) &
      ! $ o m p private(i) 
      do i=1, size(idx)
        if (mask(i)) then 
          if ((1<=idx(i)).and.(idx(i) <= idxmap%local_rows)) then
            idx(i) = idxmap%loc_to_glob(idx(i))
          else if ((idxmap%local_rows < idx(i)).and.(idx(i) <= idxmap%local_cols)&
               & .and.(.not.owned_)) then
            idx(i) = idxmap%loc_to_glob(idx(i))
          else 
            idx(i) = -1
          end if
        end if
      end do
      ! $ o m p end parallel do 
    else  if (.not.present(mask)) then 

      ! $ o m p parallel do default(none) schedule(dynamic) &
      ! $ o m p shared(idx,idxmap,owned_) &
      ! $ o m p private(i) 
      do i=1, size(idx)
        if ((1<=idx(i)).and.(idx(i) <= idxmap%local_rows)) then
          idx(i) = idxmap%loc_to_glob(idx(i))
        else if ((idxmap%local_rows < idx(i)).and.(idx(i) <= idxmap%local_cols)&
             & .and.(.not.owned_)) then
          idx(i) = idxmap%loc_to_glob(idx(i))
        else 
          idx(i) = -1
        end if
      end do
      ! $ o m p end parallel do
    end if
  end subroutine hash_l2gv1

  subroutine hash_l2gv2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_lpk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: is, im

    is = size(idxin)
    im = min(is,size(idxout))
    idxout(1:im) = idxin(1:im)
    call idxmap%l2gip(idxout(1:im),info,mask,owned)
    if (is > im) then 
      write(0,*) 'l2gv2 err -3'
      info = -3 
    end if

  end subroutine hash_l2gv2


  subroutine hash_g2ls1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_lpk_) :: idxv(1)
    info = 0

    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idx 
    call idxmap%g2lip(idxv,info,owned=owned)
    idx = idxv(1) 

  end subroutine hash_g2ls1

  subroutine hash_g2ls2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_lpk_) :: idxv(1)
    info = 0

    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idxin
    call idxmap%g2lip(idxv,info,owned=owned)
    idxout = idxv(1) 

  end subroutine hash_g2ls2


  subroutine hash_g2lv1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: i, lip, nrow, nrm, is
    integer(psb_lpk_) :: ncol, ip, tlip, mglob
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_)   :: iam, np
    logical :: owned_

    info = 0
    ctxt = idxmap%get_ctxt()
    call psb_info(ctxt,iam,np) 
    if (present(mask)) then 
      if (size(mask) < size(idx)) then 
        info = -1
        return
      end if
    end if
    if (present(owned)) then 
      owned_ = owned
    else
      owned_ = .false.
    end if

    is = size(idx)

    mglob = idxmap%get_gr()
    nrow  = idxmap%get_lr()
    ncol  = idxmap%get_lc()
    if (owned_) then 
      nrm = nrow
    else
      nrm = ncol
    end if
    if (present(mask)) then 

      if (idxmap%is_asb()) then 

        call hash_inner_cnv(is,idx,idxmap%hashvmask,&
             & idxmap%hashv,idxmap%glb_lc,mask=mask, nrm=nrm)

      else if (idxmap%is_valid()) then 

          ! $ o m p parallel do default(none) schedule(dynamic) &
          ! $ o m p shared(mask,is,idx,mglob,idxmap,nrm,ncol,nrow,owned_) &
          ! $ o m p private(i,ip,lip,tlip,info)
          do i = 1, is
            if (mask(i)) then 
              ip = idx(i) 
              if ((ip < 1 ).or.(ip>mglob)) then 
                idx(i) = -1
                cycle
              endif
              call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,&
                   & idxmap%glb_lc,nrm)
              if (lip < 0) then 
                call psb_hash_searchkey(ip,tlip,idxmap%hash,info)
                lip = tlip
                info = 0
              end if
              if (owned_) then 
                if (lip<=nrow) then 
                  idx(i) = lip
                else 
                  idx(i) = -1
                endif
              else
                idx(i) = lip
              endif
            end if
          enddo
          ! $ o m p end parallel do

!!$          call hash_inner_cnv(is,idx,idxmap%hashvmask,idxmap%hashv,&
!!$               & idxmap%glb_lc,nrm=nrm,mask=mask)
!!$
!!$          do i = 1, is
!!$            lip = idx(i)
!!$            if (lip < 0) then 
!!$              call psb_hash_searchkey(ip,tlip,idxmap%hash,info)
!!$              lip = tlip
!!$              info = 0
!!$              if (owned_) then 
!!$                if (lip<=nrow) then 
!!$                  idx(i) = lip
!!$                else 
!!$                  idx(i) = -1
!!$                endif
!!$              else
!!$                idx(i) = lip
!!$              endif
!!$            end if
!!$          enddo
!!$        
      else 
        write(0,*) 'Hash status: invalid ',idxmap%get_state()
        idx(1:is) = -1
        info = -1
      end if

    else  if (.not.present(mask)) then 

      if (idxmap%is_asb()) then 

        call hash_inner_cnv(is,idx,idxmap%hashvmask,&
             & idxmap%hashv,idxmap%glb_lc,nrm=nrm)

      else if (idxmap%is_valid()) then 

        ! $ o m p parallel do default(none) schedule(dynamic) &
        ! $ o m p shared(is,idx,mglob,idxmap,nrm,ncol,nrow,owned_) &
        ! $ o m p private(i,ip,lip,tlip,info) 
        do i = 1, is
          ip = idx(i) 
          if ((ip < 1 ).or.(ip>mglob)) then 
            idx(i) = -1
            cycle
          endif
          call hash_inner_cnv(ip,lip,idxmap%hashvmask,&
               & idxmap%hashv,idxmap%glb_lc,nrm)
          if (lip < 0) then
            call psb_hash_searchkey(ip,tlip,idxmap%hash,info)
            lip = tlip
            info = 0
          end if
          if (owned_) then 
            if (lip<=nrow) then 
              idx(i) = lip
            else 
              idx(i) = -1
            endif
          else
            idx(i) = lip
          endif
        enddo
        ! $ o m p end parallel do 
      else 
        write(0,*) 'Hash status: invalid ',idxmap%get_state()
        idx(1:is) = -1
        info = -1
      end if
    end if
  end subroutine hash_g2lv1

  subroutine hash_g2lv2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: i, lip, nrow, nrm, is, im
    integer(psb_lpk_) :: ncol, ip, tlip, mglob
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_)   :: iam, np
    logical :: owned_
    is = size(idxin)
    im = min(is,size(idxout))

    info = 0
    ctxt = idxmap%get_ctxt()
    call psb_info(ctxt,iam,np) 

    if (present(mask)) then 
      if (size(mask) < size(idxin)) then 
        info = -1
        return
      end if
    end if
    if (present(owned)) then 
      owned_ = owned
    else
      owned_ = .false.
    end if

    is = min(size(idxin), size(idxout))

    mglob = idxmap%get_gr()
    nrow  = idxmap%get_lr()
    ncol  = idxmap%get_lc()
    if (owned_) then 
      nrm = nrow
    else
      nrm = ncol
    end if
    if (present(mask)) then 

      if (idxmap%is_asb()) then 

        call hash_inner_cnv(is,idxin,idxout,idxmap%hashvmask,&
             & idxmap%hashv,idxmap%glb_lc,mask=mask, nrm=nrm)

      else if (idxmap%is_valid()) then 
        call hash_inner_cnv(is,idxin,idxout,idxmap%hashvmask,&
             & idxmap%hashv,idxmap%glb_lc,nrm=nrm,mask=mask)
        ! $ o m p parallel do default(none) schedule(dynamic) &
        ! $ o m p shared(is,idxin,idxout,mglob,idxmap,nrm,ncol,nrow,owned_) &
        ! $ o m p private(i,ip,lip,tlip,info) 
        do i = 1, is
          if (mask(i).and.(idxout(i)<0)) then 
            ip = idxin(i) 
            if ((ip < 1 ).or.(ip>mglob)) then 
              idxout(i) = -1
              cycle
            endif
            if (idxout(i) < 0) then
              call psb_hash_searchkey(ip,tlip,idxmap%hash,info)
              lip = tlip
              info = 0
              if (owned_) then 
                if (lip<=nrow) then 
                  idxout(i) = lip
                else 
                  idxout(i) = -1
                endif
              else
                idxout(i) = lip
              endif
            end if
          end if
        enddo
        ! $ o m p end parallel do 
      else 
        write(0,*) 'Hash status: invalid ',idxmap%get_state()
        idxout(1:is) = -1
        info = -1
      end if

    else  if (.not.present(mask)) then 

      if (idxmap%is_asb()) then 

        call hash_inner_cnv(is,idxin,idxout,idxmap%hashvmask,&
             & idxmap%hashv,idxmap%glb_lc,nrm=nrm)

      else if (idxmap%is_valid()) then 
        call hash_inner_cnv(is,idxin,idxout,idxmap%hashvmask,&
             & idxmap%hashv,idxmap%glb_lc,nrm=nrm)
        ! $ o m p parallel do default(none) schedule(dynamic) &
        ! $ o m p shared(is,idxin,idxout,mglob,idxmap,nrm,ncol,nrow,owned_) &
        ! $ o m p private(i,ip,lip,tlip,info) 
        do i = 1, is
          ip = idxin(i) 
          if ((ip < 1 ).or.(ip>mglob)) then 
            idxout(i) = -1
            cycle
          endif
          if (idxout(i) < 0) then
            call psb_hash_searchkey(ip,tlip,idxmap%hash,info)
            lip = tlip
            info = 0
            if (owned_) then 
              if (lip<=nrow) then 
                idxout(i) = lip
              else 
                idxout(i) = -1
              endif
            else
              idxout(i) = lip
            endif
          end if
        enddo
        ! $ o m p end parallel do 
      else 
        write(0,*) 'Hash status: invalid ',idxmap%get_state()
        idxout(1:is) = -1
        info = -1
      end if
    end if
  end subroutine hash_g2lv2



  subroutine hash_g2ls1_ins(idx,idxmap,info,mask,lidx)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx

    integer(psb_lpk_) :: idxv(1)  
    integer(psb_ipk_) :: lidxv(1)

    info = 0
    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idx
    if (present(lidx)) then 
      lidxv(1) = lidx
      call idxmap%g2lip_ins(idxv,info,lidx=lidxv)
    else
      call idxmap%g2lip_ins(idxv,info)
    end if
    idx = idxv(1) 

  end subroutine hash_g2ls1_ins

  subroutine hash_g2ls2_ins(idxin,idxout,idxmap,info,mask,lidx)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx

    integer(psb_lpk_) :: idxv(1)  
    integer(psb_ipk_) :: lidxv(1)

    info = 0
    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idxin
    if (present(lidx)) then 
      lidxv(1) = lidx
      call idxmap%g2lip_ins(idxv,info,lidx=lidxv)
    else
      call idxmap%g2lip_ins(idxv,info)
    end if
    idxout = idxv(1) 

  end subroutine hash_g2ls2_ins

  ! #################### THESIS ####################
  subroutine hash_g2lv1_ins(idx,idxmap,info,mask,lidx)
    use psb_timers_mod
#ifdef PSB_OPENMP    
    use omp_lib
#endif
    implicit none 

    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional    :: mask(:)
    integer(psb_ipk_), intent(in), optional    :: lidx(:)

    integer(psb_ipk_) :: i, is, lip, nrow, ncol,&
         & err_act
    integer(psb_lpk_) :: mglob, ip, nxt, tlip
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: me, np,ith
    character(len=20)   :: name,ch_err
    integer(psb_ipk_), allocatable   :: tidx(:)
!!$    logical :: use_openmp = .true.
#ifdef PSB_OPENMP
    integer(kind = OMP_lock_kind) :: ins_lck
#endif
    logical, volatile :: isLoopValid
    logical, parameter  :: do_timings=.true.
    integer(psb_ipk_), save  :: ins_phase1=-1, ins_phase2=-1, ins_phase3=-1, ins_phase4=-1
    integer(psb_ipk_), save  :: ins_phase11=-1, ins_phase12=-1

    info = psb_success_
    name = 'hash_g2lv1_ins'
    call psb_erractionsave(err_act)

    ctxt = idxmap%get_ctxt()
    call psb_info(ctxt, me, np)
    if ((do_timings).and.(ins_phase1==-1))       &
         & ins_phase1 = psb_get_timer_idx("HSHINS: inner_cnv ")
    if ((do_timings).and.(ins_phase2==-1))       &
         & ins_phase2 = psb_get_timer_idx("HSINS: srchins_lp")
!!$  if ((do_timings).and.(ins_phase3==-1))       &
!!$       & ins_phase3 = psb_get_timer_idx("HSHINS: csput")
!!$  if ((do_timings).and.(ins_phase4==-1))       &
!!$       & ins_phase4 = psb_get_timer_idx("HSHINS: rmt%csput")
    is = size(idx)
    call psb_realloc(is,tidx,info)
    call idxmap%lg2lv2_ins(idx,tidx,info,mask=mask,lidx=lidx)
    idx(1:is) = tidx(1:is)
  end subroutine hash_g2lv1_ins
  
  subroutine hash_g2lv2_ins(idxin,idxout,idxmap,info,mask,lidx)
    use psb_timers_mod
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)
    integer(psb_ipk_) :: is, im
    integer(psb_ipk_) :: i, lip, nrow, ncol
    integer(psb_lpk_) :: mglob, ip, nxt, tlip
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: me, np, ith, err_act
    character(len=20)   :: name,ch_err
    logical, volatile :: isLoopValid
    logical, parameter  :: do_timings=.false.

    info = psb_success_
    name = 'hash_g2lv2_ins'
    call psb_erractionsave(err_act)

    ctxt = idxmap%get_ctxt()
    call psb_info(ctxt, me, np)
    is = size(idxin)
    is = min(is,size(idxout))
    mglob = idxmap%get_gr()
    nrow  = idxmap%get_lr()
    !write(0,*)me, name, ':', present(lidx),present(mask),idxmap%is_bld()
    isLoopValid = .true.
    if (idxmap%is_bld()) then 

      if (present(lidx)) then
        if (present(mask)) then 
          ! $ o m p parallel do default(none) schedule(dynamic) &
          ! $ o m p shared(lidx,mask,name,me,is,idxin,idxout,ins_lck,mglob,idxmap,ncol,nrow,psb_laddsz) &
          ! $ o m p private(i,ip,lip,tlip,nxt,info) &
          ! $ o m p reduction(.AND.:isLoopValid)          
          do i = 1, is
            ncol  = idxmap%get_lc()
            if (mask(i)) then 
              ip = idxin(i) 
              if ((ip < 1 ).or.(ip>mglob) ) then 
                idxout(i) = -1
                cycle
              endif
              call hash_inner_cnv(ip,lip,idxmap%hashvmask,&
                   & idxmap%hashv,idxmap%glb_lc,ncol)
              if (lip < 0) then
                nxt = lidx(i)
                if (nxt <= nrow) then 
                  idxout(i) = -1
                  cycle
                endif
                call psb_hash_searchinskey(ip,tlip,nxt,idxmap%hash,info)
                lip = tlip
                if (info >=0) then 
                  if (nxt == tlip) then 
                    ncol = max(ncol,nxt)
                    call psb_ensure_size(ncol,idxmap%loc_to_glob,info,&
                         & pad=-1_psb_lpk_)
                    if (info /= psb_success_) then
                      !write(0,*) 'Error spot'
                      write(0,*)'Problem 5:',info,lip,size(idxmap%loc_to_glob)
                      info = lip
                      call psb_errpush(psb_err_from_subroutine_ai_,name,&
                           &a_err='psb_ensure_size',i_err=(/info/))
                      isLoopValid = .false.
                    end if
                    idxmap%loc_to_glob(nxt)  = ip
                    call idxmap%set_lc(ncol)
                  endif
                  info = psb_success_
                else
                  call psb_errpush(psb_err_from_subroutine_ai_,name,&
                       & a_err='SearchInsKeyVal',i_err=(/info/))
                  isLoopValid = .false.
                end if
              end if
              idxout(i) = lip
              info = psb_success_
            else
              idxout(i) = -1
            end if
          enddo
          ! $ o m p end parallel do

        else if (.not.present(mask)) then 

          ! $ o m p parallel do default(none) schedule(dynamic) &
          ! $ o m p shared(lidx,name,me,is,idxin,idxout,ins_lck,mglob,idxmap,ncol,nrow,psb_laddsz) &
          ! $ o m p private(i,ip,lip,tlip,nxt,info) &
          ! $ o m p reduction(.AND.:isLoopValid)          
          do i = 1, is
            ncol  = idxmap%get_lc()
            ip    = idxin(i) 
            if ((ip < 1 ).or.(ip>mglob)) then 
              idxout(i) = -1
              cycle
            endif
            call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,&
                 & idxmap%glb_lc,ncol)
            if (lip < 0) then 
              nxt = lidx(i)
              if (nxt <= nrow) then 
                idxout(i) = -1
                cycle
              endif
              call psb_hash_searchinskey(ip,tlip,nxt,idxmap%hash,info)
              lip = tlip

              if (info >=0) then 
                if (nxt == lip) then 
                  ncol = max(nxt,ncol)
                  call psb_ensure_size(ncol,idxmap%loc_to_glob,info,&
                       & pad=-1_psb_lpk_)
                  if (info /= psb_success_) then
                    !write(0,*) 'Error spot'                      
                    write(0,*)'Problem 6:',info,lip,size(idxmap%loc_to_glob)
                    info = lip
                    call psb_errpush(psb_err_from_subroutine_ai_,name,&
                         &a_err='psb_ensure_size',i_err=(/info/))
                    isLoopValid = .false.
                  end if
                  idxmap%loc_to_glob(nxt)  = ip
                  call idxmap%set_lc(ncol)
                endif
                info = psb_success_
              else
                call psb_errpush(psb_err_from_subroutine_ai_,name,&
                     & a_err='SearchInsKeyVal',i_err=(/info/))
                isLoopValid = .false.
              end if
            end if
            idxout(i) = lip
            info = psb_success_
          enddo
          ! $ o m p end parallel do

        end if

      else if (.not.present(lidx)) then 

        if (present(mask)) then
          ncol = idxmap%get_lc()
          call hash_inner_cnv(is,idxin,idxout,idxmap%hashvmask,idxmap%hashv,&
               & idxmap%glb_lc,nrm=ncol, mask=mask)
          !          write(0,*) me,' v2 after hash_inner_cnv ',idx(1:is)
          do i = 1, is
            if (mask(i).and.(idxout(i)<0)) then
              ncol = idxmap%get_lc()
              nxt  = ncol + 1 
              ip  = idxin(i)
              call psb_hash_searchinskey(ip,tlip,nxt,idxmap%hash,info)
              lip = tlip
              !if (i==1) write(0,*) me,' v2 isrchins:',i,lip

              if (info >=0) then 
                if (nxt == lip) then 
                  ncol = nxt
                  call psb_ensure_size(ncol,idxmap%loc_to_glob,info,&
                       & pad=-1_psb_lpk_)
                  if (info /= psb_success_) then
                    write(0,*)'Problem 7:',info,lip,size(idxmap%loc_to_glob)
                    info = lip
                    call psb_errpush(psb_err_from_subroutine_ai_,name,&
                         & a_err='psb_ensure_size',i_err=(/info/))
                    isLoopValid = .false.
                  end if
                  idxmap%loc_to_glob(nxt)  = ip
                  call idxmap%set_lc(ncol)
                endif
                info = psb_success_
              else
                call psb_errpush(psb_err_from_subroutine_ai_,name,&
                     & a_err='SearchInsKeyVal',i_err=(/info/))
                isLoopValid = .false.
              end if
              idxout(i) = lip
              info = psb_success_
            else if (.not.mask(i)) then 
              idxout(i) = -1
            end if
          enddo
          !          write(0,*) me,' v2 after cleanup ',idx(1:is)
        else if (.not.present(mask)) then 

          do i = 1, is
            ncol  = idxmap%get_lc()
            ip = idxin(i) 
            if ((ip < 1 ).or.(ip>mglob)) then 
              idxout(i) = -1
              cycle
            endif
            nxt = ncol + 1 
            call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,&
                 & idxmap%glb_lc,ncol)
            if (lip > 0) then
              idxout(i) = lip
              info = psb_success_
            else
              call psb_hash_searchinskey(ip,tlip,nxt,idxmap%hash,info)
              lip = tlip

              if (info >=0) then 
                if (nxt == lip) then 
                  ncol = nxt
                  call psb_ensure_size(ncol,idxmap%loc_to_glob,info,&
                       & pad=-1_psb_lpk_)
                  if (info /= psb_success_) then
                    write(0,*)'Problem 8:',info,lip,size(idxmap%loc_to_glob)
                    info = lip
                    ch_err='psb_ensure_size'
                    call psb_errpush(psb_err_from_subroutine_ai_,name,&
                         &a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                    isLoopValid = .false.

                  end if
                  idxmap%loc_to_glob(nxt)  = ip
                  call idxmap%set_lc(ncol)
                endif
                info = psb_success_
              else
                ch_err='SearchInsKeyVal'
                call psb_errpush(psb_err_from_subroutine_ai_,name,&
                     & a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                isLoopValid = .false.
              end if
              idxout(i) = lip
              info = psb_success_
            end if
          enddo

        end if
      end if
    else 
      ! Wrong state
      idxout(:) = -1
      info = -1
    end if
    if (.not. isLoopValid) goto 9999
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return

  end subroutine hash_g2lv2_ins
  ! ################## END THESIS #########################

  !
  ! init from VL, with checks on input.
  !
  subroutine hash_init_vl(idxmap,ctxt,vl,info)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    type(psb_ctxt_type), intent(in)    :: ctxt
    integer(psb_lpk_), intent(in)  :: vl(:)
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_ipk_) :: iam, np
    integer(psb_ipk_) ::  i,  nlu, nl, int_err(5)
    integer(psb_lpk_) ::  m, nrt
    character(len=20), parameter :: name='hash_map_init_vl'
    real(psb_dpk_) :: t0, t1, t2,t3, t4, t5
    info = 0
    call psb_info(ctxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ctxt'
      info = -1
      return
    end if
    nl = size(vl) 

    m   = maxval(vl(1:nl))
    nrt = nl
    call psb_sum(ctxt,nrt)
    call psb_max(ctxt,m)

    if ((m /= nrt).and.(iam == psb_root_))  then 
      write(psb_err_unit,*) trim(name),&
           & ' Warning: we got to hash_init_vl but there is a mismatch'
      write(psb_err_unit,*) trim(name),&
           & '        : in the global sizes!',m,nrt

    end if
    call hash_init_vlu(idxmap,ctxt,m,nl,vl,info)

  end subroutine hash_init_vl

  subroutine hash_init_vg(idxmap,ctxt,vg,info)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    type(psb_ctxt_type), intent(in)    :: ctxt
    integer(psb_ipk_), intent(in)  :: vg(:)
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_ipk_) :: iam, np
    integer(psb_ipk_) :: i, j, nl, int_err(5)
    integer(psb_lpk_) :: n
    integer(psb_lpk_), allocatable :: vlu(:)

    info = 0
    call psb_info(ctxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ctxt:'
      info = -1
      return
    end if

    n    = size(vg)
    nl   = 0
    do i=1, n
      if ((vg(i)<0).or.(vg(i)>=np)) then 
        info = psb_err_partfunc_wrong_pid_
        int_err(1) = 3
        int_err(2) = vg(i)
        int_err(3) = i
        exit
      endif
      if (vg(i) == iam) nl = nl + 1 
    end do

    allocate(vlu(nl), stat=info) 
    if (info /= 0) then 
      info = -1
      return
    end if

    j = 0
    do i=1, n
      if (vg(i) == iam) then 
        j      = j + 1 
        vlu(j) = i
      end if
    end do


    call hash_init_vlu(idxmap,ctxt,n,nl,vlu,info)    


  end subroutine hash_init_vg

  !
  ! init from VL, with no checks on input
  !
  subroutine hash_init_vlu(idxmap,ctxt,ntot,nl,vlu,info)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    type(psb_ctxt_type), intent(in)    :: ctxt
    integer(psb_lpk_), intent(in)  :: vlu(:), ntot
    integer(psb_ipk_), intent(in)  :: nl
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_ipk_) :: iam, np
    integer(psb_ipk_) :: i, j, lc2, nlu, m, nrt,int_err(5)
    character(len=20), parameter :: name='hash_map_init_vlu'

    info = 0
    call psb_info(ctxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ctxt:'
      info = -1
      return
    end if

    idxmap%global_rows  = ntot
    idxmap%global_cols  = ntot
    idxmap%local_rows   = nl
    idxmap%local_cols   = nl
    idxmap%ctxt         = ctxt
    idxmap%state        = psb_desc_bld_

    lc2 = int(1.5*nl) 
    call psb_realloc(lc2,idxmap%loc_to_glob,info) 
    if (info /= 0)  then
      info = -2
      return
    end if

    call psb_hash_init(nl,idxmap%hash,info)
    if (info /= 0) then 
      write(0,*) 'from Hash_Init inside init_vlu',info
      info = -3
      return 
    endif

    do i=1, nl
      idxmap%loc_to_glob(i) = vlu(i) 
    end do

    call hash_bld_g2l_map(idxmap,info)
    call idxmap%set_state(psb_desc_bld_)

  end subroutine hash_init_vlu


  subroutine hash_bld_g2l_map(idxmap,info)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_)   :: iam, np
    integer(psb_ipk_)   :: i, j, m, nl
    integer(psb_ipk_)   :: ih, nh, idx, nbits
    integer(psb_lpk_)   :: key, hsize, hmask
    character(len=20), parameter :: name='hash_map_init_vlu'

    info = 0
    ctxt = idxmap%get_ctxt()

    call psb_info(ctxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ctxt:'
      info = -1
      return
    end if

    nl = idxmap%get_lc()

    call psb_realloc(nl,2,idxmap%glb_lc,info) 

    nbits = psb_hash_bits
    hsize = 2**nbits
    do 
      if (hsize < 0) then 
        ! This should never happen for sane values
        ! of psb_max_hash_bits.
        write(psb_err_unit,*) &
             & 'Error: hash size overflow ',hsize,nbits
        info = -2 
        return
      end if
      if (hsize > nl) exit
      if (nbits >= psb_max_hash_bits) exit
      nbits = nbits + 1
      hsize = hsize * 2 
    end do

    hmask = hsize - 1 
    idxmap%hashvsize = hsize
    idxmap%hashvmask = hmask

    if (info == psb_success_) &
         & call psb_realloc(hsize+1,idxmap%hashv,info,lb=0_psb_lpk_)
    if (info /= psb_success_) then 
      ! !$      ch_err='psb_realloc'
      ! !$      call psb_errpush(info,name,a_err=ch_err)
      ! !$      goto 9999
      info = -4 
      return
    end if

    idxmap%hashv(:) = 0

    do i=1, nl
      key = idxmap%loc_to_glob(i) 
      ih  = iand(key,hmask) 
      idxmap%hashv(ih) = idxmap%hashv(ih) + 1
    end do

    nh = idxmap%hashv(0) 
    idx = 1

    do i=1, hsize
      idxmap%hashv(i-1) = idx
      idx = idx + nh
      nh = idxmap%hashv(i)
    end do

    do i=1, nl
      key                  = idxmap%loc_to_glob(i)
      ih                   = iand(key,hmask)
      idx                  = idxmap%hashv(ih) 
      idxmap%glb_lc(idx,1) = key
      idxmap%glb_lc(idx,2) = i
      idxmap%hashv(ih)     = idxmap%hashv(ih) + 1
    end do

    do i = hsize, 1, -1 
      idxmap%hashv(i) = idxmap%hashv(i-1)
    end do

    idxmap%hashv(0) = 1
    do i=0, hsize-1 
      idx = idxmap%hashv(i)
      nh  = idxmap%hashv(i+1) - idxmap%hashv(i) 
      if (nh > 1) then 
        call psb_msort(idxmap%glb_lc(idx:idx+nh-1,1),&
             & ix=idxmap%glb_lc(idx:idx+nh-1,2),&
             & flag=psb_sort_keep_idx_)
      end if
    end do

  end subroutine hash_bld_g2l_map


  subroutine hash_asb(idxmap,info)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(out)     :: info

    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_)   :: iam, np 
    integer(psb_ipk_)   :: nhal

    info = 0 
    ctxt = idxmap%get_ctxt()
    call psb_info(ctxt,iam,np)

    nhal = max(0,idxmap%local_cols-idxmap%local_rows)

    call hash_bld_g2l_map(idxmap,info)
    if (info /= 0) then 
      write(0,*) 'Error from bld_g2l_map', info
      return
    end if


    call psb_free(idxmap%hash,info)

    if (info /= 0) then
      write(0,*) 'Error from hash free', info
      return
    end if

    call idxmap%set_state(psb_desc_asb_)

  end subroutine hash_asb

  function hash_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'HASH'
  end function hash_get_fmt


  subroutine hash_inner_cnvs1(x,hashmask,hashv,glb_lc,nrm)
    implicit none 
    integer(psb_lpk_), intent(in)    :: hashmask,glb_lc(:,:)
    integer(psb_ipk_), intent(in)    :: hashv(0:)
    integer(psb_lpk_), intent(inout) :: x
    integer(psb_ipk_), intent(in)    :: nrm
    integer(psb_ipk_) :: idx,nh,tmp
    integer(psb_lpk_) :: key, ih
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !

    key = x
    ih  = iand(key,hashmask)
    idx = hashv(ih)
    nh  = hashv(ih+1) - hashv(ih) 
    tmp = hash_srch(key,idx,nh,glb_lc(:,1))
    if (tmp > 0) then 
      x = glb_lc(tmp,2)
      if (x > nrm) then 
        x = -1 
      end if
    else         
      x = tmp 
    end if
  end subroutine hash_inner_cnvs1

  subroutine hash_inner_cnvs2(x,y,hashmask,hashv,glb_lc,nrm)
    implicit none 
    integer(psb_ipk_), intent(in)  :: hashv(0:)
    integer(psb_lpk_), intent(in)  :: hashmask, x, glb_lc(:,:)
    integer(psb_ipk_), intent(out) :: y
    integer(psb_ipk_), intent(in)  :: nrm
    integer(psb_ipk_) :: idx,nh,tmp
    integer(psb_lpk_) :: ih, key
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !

    key = x
    ih  = iand(key,hashmask)
    idx = hashv(ih)
    nh  = hashv(ih+1) - hashv(ih) 
    tmp = hash_srch(key,idx,nh,glb_lc(:,1))
    if (tmp > 0) then 
      y = glb_lc(tmp,2)
      if (y > nrm) then 
        y = -1 
      end if
    else         
      y = tmp 
    end if
  end subroutine hash_inner_cnvs2


  subroutine hash_inner_cnv1(n,x,hashmask,hashv,glb_lc,mask,nrm)
    implicit none 
    integer(psb_ipk_), intent(in)  :: n, hashv(0:)
    integer(psb_lpk_), intent(inout) :: x(:)
    integer(psb_lpk_), intent(in)  :: glb_lc(:,:),hashmask
    logical, intent(in), optional  :: mask(:)
    integer(psb_ipk_), intent(in), optional  :: nrm

    integer(psb_ipk_) :: i, nh,tmp
    integer(psb_lpk_) :: ih, key, idx
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !
    if (present(mask)) then 
      ! $ o m p parallel do default(none) schedule(dynamic) &
      ! $ o m p shared(n,hashv,hashmask,x,glb_lc,nrm,mask) &
      ! $ o m p private(i,key,idx,ih,nh,tmp,lb,ub,lm) 
      do i=1, n
        if (mask(i)) then 
          key = x(i) 
          ih  = iand(key,hashmask)
          idx = hashv(ih)
          nh  = hashv(ih+1) - hashv(ih)
          tmp = hash_srch(key,idx,nh,glb_lc(:,1))
          if (tmp > 0) then 
            x(i) = glb_lc(tmp,2)
            if (present(nrm)) then 
              if (x(i) > nrm) then 
                x(i) = -1 
              end if
            end if
          else         
            x(i) = tmp 
          end if
        end if
      end do
      ! $ o m p end parallel do 
    else
      ! $ o m p parallel do default(none) schedule(dynamic) &
      ! $ o m p shared(n,hashv,hashmask,x,glb_lc,nrm) &
      ! $ o m p private(i,key,idx,ih,nh,tmp,lb,ub,lm) 
      do i=1, n
        key = x(i) 
        ih  = iand(key,hashmask)
        idx = hashv(ih)
        nh  = hashv(ih+1) - hashv(ih)
        tmp = hash_srch(key,idx,nh,glb_lc(:,1))
        if (tmp > 0) then 
          x(i) = glb_lc(tmp,2)
          if (present(nrm)) then 
            if (x(i) > nrm) then 
              x(i) = -1 
            end if
          end if
        else         
          x(i) = tmp 
        end if
      end do
      ! $ o m p end parallel do 
    end if
  end subroutine hash_inner_cnv1

  subroutine hash_inner_cnv2(n,x,y,hashmask,hashv,glb_lc,mask,nrm)
    implicit none 
    integer(psb_ipk_), intent(in)  :: n, hashv(0:)
    integer(psb_lpk_), intent(in)  :: hashmask,glb_lc(:,:)
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: nrm
    integer(psb_lpk_), intent(in)  :: x(:)
    integer(psb_ipk_), intent(out) :: y(:)

    integer(psb_ipk_) :: i, idx,nh,tmp
    integer(psb_lpk_) :: ih, key
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !
    if (present(mask)) then 
      ! $ o m p parallel do default(none) schedule(dynamic) &
      ! $ o m p shared(n,hashv,hashmask,x,y,glb_lc,nrm,mask,psb_err_unit) &
      ! $ o m p private(i,key,idx,ih,nh,tmp,lb,ub,lm) 
      do i=1, n
        if (mask(i)) then 
          key = x(i) 
          ih  = iand(key,hashmask)
          if (ih > ubound(hashv,1) ) then 
            write(psb_err_unit,*) ' In inner cnv: ',ih,ubound(hashv)
          end if
          idx = hashv(ih)
          nh  = hashv(ih+1) - hashv(ih)
          tmp = hash_srch(key,idx,nh,glb_lc(:,1))
          if (tmp > 0) then 
            y(i) = glb_lc(tmp,2)
            if (present(nrm)) then 
              if (y(i) > nrm) then 
                y(i) = -1 
              end if
            end if
          else         
            y(i) = tmp 
          end if
        end if
      end do
      ! $ o m p end parallel do 
    else

      ! $ o m p parallel do default(none) schedule(dynamic) &
      ! $ o m p shared(n,hashv,hashmask,x,y,glb_lc,nrm,psb_err_unit) &
      ! $ o m p private(i,key,idx,ih,nh,tmp,lb,ub,lm) 
      do i=1, n
        key = x(i) 
        ih  = iand(key,hashmask)
        if (ih > ubound(hashv,1) ) then 
          write(psb_err_unit,*) ' In inner cnv: ',ih,ubound(hashv)
        end if
        idx = hashv(ih)
        nh  = hashv(ih+1) - hashv(ih) 
        tmp = hash_srch(key,idx,nh,glb_lc(:,1))
        if (tmp > 0) then 
          y(i) = glb_lc(tmp,2)
          if (present(nrm)) then 
            if (y(i) > nrm) then 
              y(i) = -1 
            end if
          end if
        else         
          y(i) = tmp 
        end if
      end do
      ! $ o m p end parallel do 
    end if
  end subroutine hash_inner_cnv2

  function hash_srch_ipk(key,idx,nh,glb_lc) result(res)
    integer(psb_lpk_), intent(in) :: key
    integer(psb_lpk_), intent(in) :: glb_lc(:)
    integer(psb_ipk_), intent(in) :: idx
    integer(psb_ipk_), intent(in) :: nh
    integer(psb_ipk_) :: res
    !
    integer(psb_ipk_) :: lb,ub,lm
    res = -1
    if (nh > 0) then
      if (nh <= seqsrchmax) then
        !
        ! If the list is short, a sequential search is enough 
        !
        do lm=idx,idx+nh-1
          if (key == glb_lc(lm)) then 
            res = lm
            exit
          end if
        end do
      else
        !
        ! Otherwise use binary 
        !
        lb = idx
        ub = idx+nh-1
        do 
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key == glb_lc(lm)) then 
            res = lm
            exit
          else if (key<glb_lc(lm)) then 
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      end if
    end if
  end function hash_srch_ipk

#if defined(PSB_IPK4) && defined(PSB_LPK8) 
  function hash_srch_lpk(key,idx,nh,glb_lc) result(res)
    integer(psb_lpk_), intent(in) :: key
    integer(psb_lpk_), intent(in) :: glb_lc(:)
    integer(psb_lpk_), intent(in) :: idx
    integer(psb_ipk_), intent(in) :: nh
    integer(psb_ipk_) :: res
    !
    integer(psb_ipk_) :: lb,ub,lm
    res = -1
    if (nh > 0) then
      if (nh <= seqsrchmax) then
        !
        ! If the list is short, a sequential search is enough 
        !
        do lm=idx,idx+nh-1
          if (key == glb_lc(lm)) then 
            res = lm
            exit
          end if
        end do
      else
        !
        ! Otherwise use binary 
        !
        lb = idx
        ub = idx+nh-1
        do 
          if (lb>ub) exit
          lm = (lb+ub)/2
          if (key == glb_lc(lm)) then 
            res = lm
            exit
          else if (key<glb_lc(lm)) then 
            ub = lm - 1
          else
            lb = lm + 1
          end if
        end do
      end if
    end if
  end function hash_srch_lpk
#endif

  subroutine hash_clone(idxmap,outmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_hash_map), intent(inout)    :: idxmap
    class(psb_indx_map), allocatable, intent(out) :: outmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='hash_clone'
    logical, parameter :: debug=.false.

    info = psb_success_
    call psb_get_erraction(err_act)
    if (allocated(outmap)) then
      call outmap%free() 
      deallocate(outmap,stat=info)
    end if
    if (info /= 0) then 
      write(0,*) 'Error: could not cleanup output'
      info = -87
      goto 9999
    end if


    allocate(psb_hash_map :: outmap, stat=info )
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if

    select type (outmap)
    type is (psb_hash_map)
      call idxmap%psb_indx_map%cpy(outmap%psb_indx_map,info)
      if (info == psb_success_) then 
        outmap%hashvsize    = idxmap%hashvsize
        outmap%hashvmask    = idxmap%hashvmask
      end if
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%loc_to_glob,outmap%loc_to_glob,info)
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%hashv,outmap%hashv,info)
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%glb_lc,outmap%glb_lc,info)
      if (info == psb_success_)&
           &  call psb_hash_copy(idxmap%hash,outmap%hash,info)

    class default
      ! This should be impossible 
      info = -1
    end select

    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return


9999 call psb_error_handler(err_act)

    return
  end subroutine hash_clone

  subroutine hash_reinit(idxmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_hash_map), intent(inout)    :: idxmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act, nr,nc,k, nl
    integer(psb_lpk_) :: lk
    integer(psb_lpk_) :: ntot
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: me, np
    integer(psb_ipk_), allocatable :: lidx(:), tadj(:), th_own(:)
    integer(psb_lpk_), allocatable :: gidx(:)
    character(len=20)  :: name='hash_reinit'
    logical, parameter :: debug=.false.

    info = psb_success_
    call psb_get_erraction(err_act)
    ctxt = idxmap%get_ctxt()
    nr = idxmap%get_lr()
    nc = idxmap%get_lc()
    ntot = idxmap%get_gr()


    lidx = (/(k,k=1,nc)/)
    gidx = (/(lk,lk=1,nc)/)
    call idxmap%l2gip(gidx,info)
    tadj = idxmap%get_p_adjcncy()
    call idxmap%get_halo_owner(th_own,info)

    call idxmap%free()
    call hash_init_vlu(idxmap,ctxt,ntot,nr,gidx(1:nr),info) 
    if (nc>nr) then 
      call idxmap%g2lip_ins(gidx(nr+1:nc),info,lidx=lidx(nr+1:nc))
    end if
    call idxmap%set_p_adjcncy(tadj)
    call idxmap%set_halo_owner(th_own,info)

    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return


9999 call psb_error_handler(err_act)

    return
  end subroutine hash_reinit


end module psb_hash_map_mod
