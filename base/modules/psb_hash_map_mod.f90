!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
  
  type, extends(psb_indx_map) :: psb_hash_map

    integer(psb_ipk_) :: hashvsize, hashvmask
    integer(psb_ipk_), allocatable :: hashv(:), glb_lc(:,:), loc_to_glob(:)
    type(psb_hash_type)  :: hash

  contains

    procedure, pass(idxmap)  :: init_vl    => hash_init_vl
    procedure, pass(idxmap)  :: hash_map_init => hash_init_vg

    procedure, pass(idxmap)  :: sizeof    => hash_sizeof
    procedure, pass(idxmap)  :: asb       => hash_asb
    procedure, pass(idxmap)  :: free      => hash_free
    procedure, pass(idxmap)  :: clone     => hash_clone
    procedure, nopass        :: get_fmt   => hash_get_fmt

    procedure, nopass        :: row_extendable => hash_row_extendable

    procedure, pass(idxmap)  :: l2gs1     => hash_l2gs1
    procedure, pass(idxmap)  :: l2gs2     => hash_l2gs2
    procedure, pass(idxmap)  :: l2gv1     => hash_l2gv1
    procedure, pass(idxmap)  :: l2gv2     => hash_l2gv2

    procedure, pass(idxmap)  :: g2ls1     => hash_g2ls1
    procedure, pass(idxmap)  :: g2ls2     => hash_g2ls2
    procedure, pass(idxmap)  :: g2lv1     => hash_g2lv1
    procedure, pass(idxmap)  :: g2lv2     => hash_g2lv2

    procedure, pass(idxmap)  :: g2ls1_ins => hash_g2ls1_ins
    procedure, pass(idxmap)  :: g2ls2_ins => hash_g2ls2_ins
    procedure, pass(idxmap)  :: g2lv1_ins => hash_g2lv1_ins
    procedure, pass(idxmap)  :: g2lv2_ins => hash_g2lv2_ins

    procedure, pass(idxmap)  :: hash_cpy
    generic, public          :: assignment(=) => hash_cpy
    procedure, pass(idxmap)  :: bld_g2l_map => hash_bld_g2l_map

  end type psb_hash_map

  private :: hash_init_vl, hash_init_vg, hash_sizeof, hash_asb, &
       & hash_free, hash_get_fmt, hash_l2gs1, hash_l2gs2, &
       & hash_l2gv1, hash_l2gv2, hash_g2ls1, hash_g2ls2, &
       & hash_g2lv1, hash_g2lv2, hash_g2ls1_ins, hash_g2ls2_ins, &
       & hash_g2lv1_ins, hash_g2lv2_ins, hash_init_vlu, &
       & hash_bld_g2l_map,  hash_inner_cnvs1, hash_inner_cnvs2,&
       & hash_inner_cnv1, hash_inner_cnv2, hash_row_extendable 

  integer(psb_ipk_), private :: laddsz=500

  interface hash_inner_cnv 
    module procedure  hash_inner_cnvs1, hash_inner_cnvs2,&
         & hash_inner_cnv1, hash_inner_cnv2 
  end interface hash_inner_cnv
  private :: hash_inner_cnv

contains

  function hash_row_extendable() result(val)
    implicit none 
    logical :: val
    val = .true.
  end function hash_row_extendable

  function hash_sizeof(idxmap) result(val)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_long_int_k_) :: val

    val = idxmap%psb_indx_map%sizeof() 
    val = val + 2 * psb_sizeof_int
    if (allocated(idxmap%hashv)) &
         & val = val + size(idxmap%hashv)*psb_sizeof_int
    if (allocated(idxmap%glb_lc)) &
         & val = val + size(idxmap%glb_lc)*psb_sizeof_int
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
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: idxv(1)
    info = 0
    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idx
    call idxmap%l2g(idxv,info,owned=owned)
    idx = idxv(1)

  end subroutine hash_l2gs1

  subroutine hash_l2gs2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    idxout = idxin
    call idxmap%l2g(idxout,info,mask,owned)

  end subroutine hash_l2gs2


  subroutine hash_l2gv1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx(:)
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

    else  if (.not.present(mask)) then 

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

    end if

  end subroutine hash_l2gv1

  subroutine hash_l2gv2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: is, im

    is = size(idxin)
    im = min(is,size(idxout))
    idxout(1:im) = idxin(1:im)
    call idxmap%l2g(idxout(1:im),info,mask,owned)
    if (is > im) then 
      write(0,*) 'l2gv2 err -3'
      info = -3 
    end if

  end subroutine hash_l2gv2


  subroutine hash_g2ls1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: idxv(1)
    info = 0

    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idx 
    call idxmap%g2l(idxv,info,owned=owned)
    idx = idxv(1) 

  end subroutine hash_g2ls1

  subroutine hash_g2ls2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    idxout = idxin
    call idxmap%g2l(idxout,info,mask,owned)

  end subroutine hash_g2ls2


  subroutine hash_g2lv1(idx,idxmap,info,mask,owned)
    use psb_penv_mod
    use psb_sort_mod
    implicit none 
    class(psb_hash_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: i, is, mglob, ip, lip, nrow, ncol, nrm 
    integer(psb_mpik_) :: ictxt, iam, np
    logical :: owned_

    info = 0
    ictxt = idxmap%get_ctxt()
    call psb_info(ictxt,iam,np) 

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

        do i = 1, is
          if (mask(i)) then 
            ip = idx(i) 
            if ((ip < 1 ).or.(ip>mglob)) then 
              idx(i) = -1
              cycle
            endif
            call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,idxmap%glb_lc,nrm)
            if (lip < 0) &
                 &  call psb_hash_searchkey(ip,lip,idxmap%hash,info)
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

        do i = 1, is
          ip = idx(i) 
          if ((ip < 1 ).or.(ip>mglob)) then 
            idx(i) = -1
            cycle
          endif
          call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,idxmap%glb_lc,nrm)
          if (lip < 0) &
               &  call psb_hash_searchkey(ip,lip,idxmap%hash,info)
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
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: is, im

    is = size(idxin)
    im = min(is,size(idxout))
    idxout(1:im) = idxin(1:im)
    call idxmap%g2l(idxout(1:im),info,mask,owned)
    if (is > im) then 
      write(0,*) 'g2lv2 err -3'
      info = -3 
    end if

  end subroutine hash_g2lv2



  subroutine hash_g2ls1_ins(idx,idxmap,info,mask,lidx)
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer, intent(in), optional :: lidx

    integer(psb_ipk_) :: idxv(1), lidxv(1)

    info = 0
    if (present(mask)) then 
      if (.not.mask) return
    end if

    idxv(1) = idx
    if (present(lidx)) then 
      lidxv(1) = lidx
      call idxmap%g2l_ins(idxv,info,lidx=lidxv)
    else
      call idxmap%g2l_ins(idxv,info)
    end if
    idx = idxv(1) 

  end subroutine hash_g2ls1_ins

  subroutine hash_g2ls2_ins(idxin,idxout,idxmap,info,mask,lidx)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer, intent(in), optional :: lidx


    idxout = idxin
    call idxmap%g2l_ins(idxout,info,mask=mask,lidx=lidx)

  end subroutine hash_g2ls2_ins


  subroutine hash_g2lv1_ins(idx,idxmap,info,mask,lidx)
    use psb_error_mod
    use psb_realloc_mod
    use psb_sort_mod
    use psb_penv_mod
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional    :: mask(:)
    integer, intent(in), optional    :: lidx(:)

    integer(psb_ipk_) :: i, is, mglob, ip, lip, nrow, ncol, &
         & nxt, err_act
    integer(psb_mpik_) :: ictxt, me, np
    character(len=20)  :: name,ch_err

    info = psb_success_
    name = 'hash_g2l_ins'
    call psb_erractionsave(err_act)

    ictxt = idxmap%get_ctxt()
    call psb_info(ictxt, me, np)

    is = size(idx)

    if (present(mask)) then 
      if (size(mask) < size(idx)) then 
        info = -1
        return
      end if
    end if
    if (present(lidx)) then 
      if (size(lidx) < size(idx)) then 
        info = -1
        return
      end if
    end if


    mglob = idxmap%get_gr()
    nrow  = idxmap%get_lr()
    if (idxmap%is_bld()) then 

      if (present(lidx)) then
        if (present(mask)) then 
          do i = 1, is
            ncol  = idxmap%get_lc()
            if (mask(i)) then 
              ip = idx(i) 
              if ((ip < 1 ).or.(ip>mglob) ) then 
                idx(i) = -1
                cycle
              endif
              call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,idxmap%glb_lc,ncol)
              if (lip < 0) then 
                nxt = lidx(i)
                if (nxt <= nrow) then 
                  idx(i) = -1
                  cycle
                endif
                call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
                if (info >=0) then 
                  if (nxt == lip) then 
                    ncol = max(ncol,nxt)
                    call psb_ensure_size(ncol,idxmap%loc_to_glob,info,pad=-ione,addsz=laddsz)
                    if (info /= psb_success_) then
                      info=1
                      ch_err='psb_ensure_size'
                      call psb_errpush(psb_err_from_subroutine_ai_,name,&
                           &a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                      goto 9999
                    end if
                    idxmap%loc_to_glob(nxt)  = ip
                    call idxmap%set_lc(ncol)
                  endif
                  info = psb_success_
                else
                  ch_err='SearchInsKeyVal'
                  call psb_errpush(psb_err_from_subroutine_ai_,name,&
                       & a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                  goto 9999
                end if
              end if
              idx(i) = lip
              info = psb_success_
            else
              idx(i) = -1
            end if
          enddo

        else if (.not.present(mask)) then 

          do i = 1, is
            ncol  = idxmap%get_lc()
            ip    = idx(i) 
            if ((ip < 1 ).or.(ip>mglob)) then 
              idx(i) = -1
              cycle
            endif
            call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,idxmap%glb_lc,ncol)
            if (lip < 0) then 
              nxt = lidx(i)
              if (nxt <= nrow) then 
                idx(i) = -1
                cycle
              endif
              call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)

              if (info >=0) then 
                if (nxt == lip) then 
                  ncol = max(nxt,ncol)
                  call psb_ensure_size(ncol,idxmap%loc_to_glob,info,pad=-ione,addsz=laddsz)
                  if (info /= psb_success_) then
                    info=1
                    ch_err='psb_ensure_size'
                    call psb_errpush(psb_err_from_subroutine_ai_,name,&
                         &a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                    goto 9999
                  end if
                  idxmap%loc_to_glob(nxt)  = ip
                  call idxmap%set_lc(ncol)
                endif
                info = psb_success_
              else
                ch_err='SearchInsKeyVal'
                call psb_errpush(psb_err_from_subroutine_ai_,name,&
                     & a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                goto 9999
              end if
            end if
            idx(i) = lip
            info = psb_success_
          enddo

        end if

      else if (.not.present(lidx)) then 

        if (present(mask)) then 
          do i = 1, is
            ncol  = idxmap%get_lc()
            if (mask(i)) then 
              ip = idx(i) 
              if ((ip < 1 ).or.(ip>mglob)) then 
                idx(i) = -1
                cycle
              endif
              nxt = ncol + 1 
              call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,idxmap%glb_lc,ncol)
              if (lip < 0) &
                   &  call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)

              if (info >=0) then 
                if (nxt == lip) then 
                  ncol = nxt
                  call psb_ensure_size(ncol,idxmap%loc_to_glob,info,pad=-ione,addsz=laddsz)
                  if (info /= psb_success_) then
                    info=1
                    ch_err='psb_ensure_size'
                    call psb_errpush(psb_err_from_subroutine_ai_,name,&
                         &a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                    goto 9999
                  end if
                  idxmap%loc_to_glob(nxt)  = ip
                  call idxmap%set_lc(ncol)
                endif
                info = psb_success_
              else
                ch_err='SearchInsKeyVal'
                call psb_errpush(psb_err_from_subroutine_ai_,name,&
                     & a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                goto 9999
              end if
              idx(i) = lip
              info = psb_success_
            else
              idx(i) = -1
            end if
          enddo

        else if (.not.present(mask)) then 

          do i = 1, is
            ncol  = idxmap%get_lc()
            ip = idx(i) 
            if ((ip < 1 ).or.(ip>mglob)) then 
              idx(i) = -1
              cycle
            endif
            nxt = ncol + 1 
            call hash_inner_cnv(ip,lip,idxmap%hashvmask,idxmap%hashv,idxmap%glb_lc,ncol)
            if (lip < 0) &
                 &  call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)

            if (info >=0) then 
              if (nxt == lip) then 
                ncol = nxt
                call psb_ensure_size(ncol,idxmap%loc_to_glob,info,pad=-ione,addsz=laddsz)
                if (info /= psb_success_) then
                  info=1
                  ch_err='psb_ensure_size'
                  call psb_errpush(psb_err_from_subroutine_ai_,name,&
                       &a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
                  goto 9999
                end if
                idxmap%loc_to_glob(nxt)  = ip
                call idxmap%set_lc(ncol)
              endif
              info = psb_success_
            else
              ch_err='SearchInsKeyVal'
              call psb_errpush(psb_err_from_subroutine_ai_,name,&
                   & a_err=ch_err,i_err=(/info,izero,izero,izero,izero/))
              goto 9999
            end if
            idx(i) = lip
            info = psb_success_
          enddo


        end if
      end if
    else 
      ! Wrong state
      idx = -1
      info = -1
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_ret_) then
      return
    else
      call psb_error(ictxt)
    end if
    return

  end subroutine hash_g2lv1_ins

  subroutine hash_g2lv2_ins(idxin,idxout,idxmap,info,mask,lidx)
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer, intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: is, im

    is = size(idxin)
    im = min(is,size(idxout))
    idxout(1:im) = idxin(1:im)
    call idxmap%g2l_ins(idxout(1:im),info,mask=mask,lidx=lidx)
    if (is > im) then 
      write(0,*) 'g2lv2_ins err -3'
      info = -3 
    end if

  end subroutine hash_g2lv2_ins

  subroutine hash_init_vl(idxmap,ictxt,vl,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_sort_mod
    use psb_realloc_mod
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: vl(:)
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_mpik_) :: iam, np
    integer(psb_ipk_) ::  i,  nlu, nl, m, nrt,int_err(5)
    integer(psb_ipk_), allocatable :: vlu(:), ix(:)
    character(len=20), parameter :: name='hash_map_init_vl'

    info = 0
    call psb_info(ictxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ictxt:',ictxt
      info = -1
      return
    end if

    nl = size(vl) 

    m   = maxval(vl(1:nl))
    nrt = nl
    call psb_sum(ictxt,nrt)
    call psb_max(ictxt,m)

    allocate(vlu(nl), ix(nl), stat=info) 
    if (info /= 0) then 
      info = -1
      return
    end if

    do i=1,nl
      if ((vl(i)<1).or.(vl(i)>m)) then 
        info = psb_err_entry_out_of_bounds_
        int_err(1) = i
        int_err(2) = vl(i)
        int_err(3) = nl
        int_err(4) = m
        exit
      endif
      vlu(i) = vl(i) 
    end do

    if ((m /= nrt).and.(iam == psb_root_))  then 
      write(psb_err_unit,*) trim(name),&
           & ' Warning: globalcheck=.false., but there is a mismatch'
      write(psb_err_unit,*) trim(name),&
           & '        : in the global sizes!',m,nrt

    end if

    call psb_msort(vlu,ix)
    nlu = 1
    do i=2,nl
      if (vlu(i) /= vlu(nlu)) then
        nlu = nlu + 1 
        vlu(nlu) = vlu(i)
        ix(nlu) = ix(i)
      end if
    end do
    call psb_msort(ix(1:nlu),vlu(1:nlu),flag=psb_sort_keep_idx_)
    
    nlu = nl
    call hash_init_vlu(idxmap,ictxt,m,nlu,vlu,info)    

  end subroutine hash_init_vl

  subroutine hash_init_vg(idxmap,ictxt,vg,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: vg(:)
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_mpik_) :: iam, np
    integer(psb_ipk_) :: i, j, nl, n, int_err(5)
    integer(psb_ipk_), allocatable :: vlu(:)

    info = 0
    call psb_info(ictxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ictxt:',ictxt
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


    call hash_init_vlu(idxmap,ictxt,n,nl,vlu,info)    


  end subroutine hash_init_vg


  subroutine hash_init_vlu(idxmap,ictxt,ntot,nl,vlu,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_sort_mod
    use psb_realloc_mod
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_mpik_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: vlu(:), nl, ntot
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_mpik_) :: iam, np
    integer(psb_ipk_) :: i, j, lc2, nlu, m, nrt,int_err(5)
    character(len=20), parameter :: name='hash_map_init_vlu'

    info = 0
    call psb_info(ictxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ictxt:',ictxt
      info = -1
      return
    end if

    idxmap%global_rows  = ntot
    idxmap%global_cols  = ntot
    idxmap%local_rows   = nl
    idxmap%local_cols   = nl
    idxmap%ictxt        = ictxt
    idxmap%state        = psb_desc_bld_
    call psb_get_mpicomm(ictxt,idxmap%mpic)

    lc2 = int(1.5*nl) 
    allocate(idxmap%loc_to_glob(lc2),stat=info) 
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
    use psb_penv_mod
    use psb_error_mod
    use psb_sort_mod
    use psb_realloc_mod
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_mpik_) :: ictxt, iam, np
    integer(psb_ipk_) :: i, j, m, nl
    integer(psb_ipk_) :: key, ih, nh, idx, nbits, hsize, hmask
    character(len=20), parameter :: name='hash_map_init_vlu'

    info = 0
    ictxt = idxmap%get_ctxt()

    call psb_info(ictxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ictxt:',ictxt
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
         & call psb_realloc(hsize+1,idxmap%hashv,info,lb=0_psb_ipk_)
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
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_hash_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(out) :: info

    integer(psb_mpik_) :: ictxt, iam, np 
    integer(psb_ipk_) :: nhal

    info = 0 
    ictxt = idxmap%get_ctxt()
    call psb_info(ictxt,iam,np)

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

    integer(psb_ipk_), intent(in)    :: hashmask,hashv(0:),glb_lc(:,:)
    integer(psb_ipk_), intent(inout) :: x
    integer(psb_ipk_), intent(in)    :: nrm
    integer(psb_ipk_) :: ih, key, idx,nh,tmp,lb,ub,lm
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a binary search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !

    key = x
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
        if (key == glb_lc(lm,1)) then 
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
      x = glb_lc(tmp,2)
      if (x > nrm) then 
        x = -1 
      end if
    else         
      x = tmp 
    end if
  end subroutine hash_inner_cnvs1

  subroutine hash_inner_cnvs2(x,y,hashmask,hashv,glb_lc,nrm)
    integer(psb_ipk_), intent(in)  :: hashmask,hashv(0:),glb_lc(:,:)
    integer(psb_ipk_), intent(in)  :: x
    integer(psb_ipk_), intent(out) :: y
    integer(psb_ipk_), intent(in)  :: nrm
    integer(psb_ipk_) :: ih, key, idx,nh,tmp,lb,ub,lm
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a binary search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !

    key = x
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
        if (key == glb_lc(lm,1)) then 
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
      y = glb_lc(tmp,2)
      if (y > nrm) then 
        y = -1 
      end if
    else         
      y = tmp 
    end if
  end subroutine hash_inner_cnvs2


  subroutine hash_inner_cnv1(n,x,hashmask,hashv,glb_lc,mask,nrm)
    integer(psb_ipk_), intent(in)    :: n,hashmask,hashv(0:),glb_lc(:,:)
    logical, intent(in), optional  :: mask(:)
    integer(psb_ipk_), intent(in), optional  :: nrm
    integer(psb_ipk_), intent(inout) :: x(:)

    integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a binary search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !
    if (present(mask)) then 
      do i=1, n
        if (mask(i)) then 
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
              if (key == glb_lc(lm,1)) then 
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
    else
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
            if (key == glb_lc(lm,1)) then 
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
          if (present(nrm)) then 
            if (x(i) > nrm) then 
              x(i) = -1 
            end if
          end if
        else         
          x(i) = tmp 
        end if
      end do
    end if
  end subroutine hash_inner_cnv1

  subroutine hash_inner_cnv2(n,x,y,hashmask,hashv,glb_lc,mask,nrm)
    integer(psb_ipk_), intent(in)  :: n, hashmask,hashv(0:),glb_lc(:,:)
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: nrm
    integer(psb_ipk_), intent(in)  :: x(:)
    integer(psb_ipk_), intent(out) :: y(:)

    integer(psb_ipk_) :: i, ih, key, idx,nh,tmp,lb,ub,lm
    !
    ! When a large descriptor is assembled the indices 
    ! are kept in a (hashed) list of ordered lists. 
    ! Thus we first hash the index, then we do a binary search on the 
    ! ordered sublist. The hashing is based on the low-order bits 
    ! for a width of psb_hash_bits 
    !
    if (present(mask)) then 
      do i=1, n
        if (mask(i)) then 
          key = x(i) 
          ih  = iand(key,hashmask)
          if (ih > ubound(hashv,1) ) then 
            write(psb_err_unit,*) ' In inner cnv: ',ih,ubound(hashv)
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
              if (key == glb_lc(lm,1)) then 
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

    else

      do i=1, n
        key = x(i) 
        ih  = iand(key,hashmask)
        if (ih > ubound(hashv,1) ) then 
          write(psb_err_unit,*) ' In inner cnv: ',ih,ubound(hashv)
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
            if (key == glb_lc(lm,1)) then 
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
          if (present(nrm)) then 
            if (y(i) > nrm) then 
              y(i) = -1 
            end if
          end if
        else         
          y(i) = tmp 
        end if
      end do
    end if
  end subroutine hash_inner_cnv2


  subroutine hash_clone(idxmap,outmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_hash_map), intent(in)    :: idxmap
    class(psb_indx_map), allocatable, intent(out) :: outmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='hash_clone'
    logical, parameter :: debug=.false.

    info = psb_success_
    call psb_get_erraction(err_act)
    if (allocated(outmap)) then 
      write(0,*) 'Error: should not be allocated on input'
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
      if (info == psb_success_) then 
        outmap%psb_indx_map = idxmap%psb_indx_map
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
!!$      outmap = idxmap
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

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act /= psb_act_ret_) then
      call psb_error()
    end if
    return
  end subroutine hash_clone


  subroutine hash_cpy(outmap,idxmap)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    class(psb_hash_map), intent(in) :: idxmap
    type(psb_hash_map), intent(out) :: outmap
    integer(psb_ipk_) :: info

    info = psb_success_
    outmap%psb_indx_map = idxmap%psb_indx_map
    outmap%hashvsize    = idxmap%hashvsize
    outmap%hashvmask    = idxmap%hashvmask
    if (info == psb_success_)&
         &  call psb_safe_ab_cpy(idxmap%loc_to_glob,outmap%loc_to_glob,info)
    if (info == psb_success_)&
         &  call psb_safe_ab_cpy(idxmap%hashv,outmap%hashv,info)
    if (info == psb_success_)&
         &  call psb_safe_ab_cpy(idxmap%glb_lc,outmap%glb_lc,info)
    if (info == psb_success_)&
         &  call psb_hash_copy(idxmap%hash,outmap%hash,info)
  end subroutine hash_cpy

    
  
end module psb_hash_map_mod
