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
! package: psb_gen_block_map_mod
!    Defines the GEN_BLOCK_MAP type.
!
!    It is the implementation of the general BLOCK distribution,
!    i.e. process I gets the I-th block of consecutive indices.
!    It needs to store the limits of the owned block, plus the global
!    indices of the local halo.
!    The choice is to store the boundaries of ALL blocks, since in general
!    there will be few processes, compared to indices, so it is possible
!    to answer the ownership question without resorting to data exchange
!    (well, the data exchange is needed but only once at initial allocation
!    time). 
!
!
module psb_gen_block_map_mod
  use psb_const_mod
  use psb_desc_const_mod
  use psb_indx_map_mod
  use psb_hash_mod
  
  type, extends(psb_indx_map) :: psb_gen_block_map
    integer(psb_lpk_) :: min_glob_row   = -1
    integer(psb_lpk_) :: max_glob_row   = -1
    integer(psb_lpk_), allocatable :: loc_to_glob(:), srt_l2g(:,:), vnl(:)
    type(psb_hash_type)  :: hash
  contains

    procedure, pass(idxmap)  :: gen_block_map_init => block_init

    procedure, pass(idxmap)  :: sizeof    => block_sizeof
    procedure, pass(idxmap)  :: asb       => block_asb
    procedure, pass(idxmap)  :: free      => block_free
    procedure, pass(idxmap)  :: clone     => block_clone
    procedure, pass(idxmap)  :: reinit    => block_reinit
    procedure, nopass        :: get_fmt   => block_get_fmt

!!$    procedure, pass(idxmap)  :: l2gs1 => block_l2gs1
!!$    procedure, pass(idxmap)  :: l2gs2 => block_l2gs2
!!$    procedure, pass(idxmap)  :: l2gv1 => block_l2gv1
!!$    procedure, pass(idxmap)  :: l2gv2 => block_l2gv2

    procedure, pass(idxmap)  :: ll2gs1 => block_ll2gs1
    procedure, pass(idxmap)  :: ll2gs2 => block_ll2gs2
    procedure, pass(idxmap)  :: ll2gv1 => block_ll2gv1
    procedure, pass(idxmap)  :: ll2gv2 => block_ll2gv2

!!$    procedure, pass(idxmap)  :: g2ls1 => block_g2ls1
!!$    procedure, pass(idxmap)  :: g2ls2 => block_g2ls2
!!$    procedure, pass(idxmap)  :: g2lv1 => block_g2lv1
!!$    procedure, pass(idxmap)  :: g2lv2 => block_g2lv2

    procedure, pass(idxmap)  :: lg2ls1 => block_lg2ls1
    procedure, pass(idxmap)  :: lg2ls2 => block_lg2ls2
    procedure, pass(idxmap)  :: lg2lv1 => block_lg2lv1
    procedure, pass(idxmap)  :: lg2lv2 => block_lg2lv2

!!$    procedure, pass(idxmap)  :: g2ls1_ins => block_g2ls1_ins
!!$    procedure, pass(idxmap)  :: g2ls2_ins => block_g2ls2_ins
!!$    procedure, pass(idxmap)  :: g2lv1_ins => block_g2lv1_ins
!!$    procedure, pass(idxmap)  :: g2lv2_ins => block_g2lv2_ins

    procedure, pass(idxmap)  :: lg2ls1_ins => block_lg2ls1_ins
    procedure, pass(idxmap)  :: lg2ls2_ins => block_lg2ls2_ins
    procedure, pass(idxmap)  :: lg2lv1_ins => block_lg2lv1_ins
    procedure, pass(idxmap)  :: lg2lv2_ins => block_lg2lv2_ins

    procedure, pass(idxmap)  :: fnd_owner => block_fnd_owner

  end type psb_gen_block_map

  private ::  block_init, block_sizeof, block_asb, block_free,&
       & block_l2gs1, block_l2gs2, block_l2gv1, block_l2gv2, &
       & block_ll2gs1, block_ll2gs2, block_ll2gv1, block_ll2gv2, &
       & block_g2ls1, block_g2ls2, block_g2lv1, block_g2lv2, &
       & block_g2ls1_ins, block_g2ls2_ins, block_g2lv1_ins, block_g2lv2_ins, &
       & block_lg2ls1_ins, block_lg2ls2_ins, block_lg2lv1_ins, block_lg2lv2_ins, &
       & block_clone, block_reinit,&
       & block_get_fmt, gen_block_search, l_gen_block_search

  interface gen_block_search
    module procedure gen_block_search, l_gen_block_search
  end interface gen_block_search

  integer(psb_ipk_), private :: laddsz=500

contains

  
  function block_sizeof(idxmap) result(val)
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
    integer(psb_epk_) :: val
    
    val = idxmap%psb_indx_map%sizeof() 
    val = val + 2 * psb_sizeof_lp
    if (allocated(idxmap%loc_to_glob)) &
         & val = val + size(idxmap%loc_to_glob)*psb_sizeof_lp
    if (allocated(idxmap%srt_l2g)) &
         & val = val + size(idxmap%srt_l2g)*psb_sizeof_lp
    if (allocated(idxmap%vnl)) &
         & val = val + size(idxmap%vnl)*psb_sizeof_lp
    val = val + psb_sizeof(idxmap%hash)
  end function block_sizeof


  subroutine block_free(idxmap)
    implicit none 
    class(psb_gen_block_map), intent(inout) :: idxmap
    integer(psb_ipk_) :: info
    if (allocated(idxmap%loc_to_glob)) &
         & deallocate(idxmap%loc_to_glob)
    if (allocated(idxmap%srt_l2g)) &
         & deallocate(idxmap%srt_l2g)

    if (allocated(idxmap%srt_l2g)) &
         & deallocate(idxmap%vnl)
    call psb_free(idxmap%hash,info)
    call idxmap%psb_indx_map%free()

  end subroutine block_free

!!$
!!$  subroutine block_l2gs1(idx,idxmap,info,mask,owned)
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(in) :: idxmap
!!$    integer(psb_ipk_), intent(inout) :: idx
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask
!!$    logical, intent(in), optional :: owned
!!$    integer(psb_ipk_) :: idxv(1)
!!$    info = 0
!!$    if (present(mask)) then 
!!$      if (.not.mask) return
!!$    end if
!!$
!!$    idxv(1) = idx
!!$    call idxmap%l2gip(idxv,info,owned=owned)
!!$    idx = idxv(1)
!!$
!!$  end subroutine block_l2gs1
!!$
!!$  subroutine block_l2gs2(idxin,idxout,idxmap,info,mask,owned)
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(in) :: idxmap
!!$    integer(psb_ipk_), intent(in)    :: idxin
!!$    integer(psb_ipk_), intent(out)   :: idxout
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask
!!$    logical, intent(in), optional :: owned
!!$
!!$    idxout = idxin
!!$    call idxmap%l2gip(idxout,info,mask,owned)
!!$    
!!$  end subroutine block_l2gs2
!!$
!!$
!!$  subroutine block_l2gv1(idx,idxmap,info,mask,owned)
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(in) :: idxmap
!!$    integer(psb_ipk_), intent(inout) :: idx(:)
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask(:)
!!$    logical, intent(in), optional :: owned
!!$    integer(psb_ipk_) :: i
!!$    logical :: owned_
!!$    info = 0
!!$
!!$    if (present(mask)) then 
!!$      if (size(mask) < size(idx)) then 
!!$        info = -1
!!$        return
!!$      end if
!!$    end if
!!$    if (present(owned)) then 
!!$      owned_ = owned
!!$    else
!!$      owned_ = .false.
!!$    end if
!!$
!!$    if (present(mask)) then 
!!$
!!$      do i=1, size(idx)
!!$        if (mask(i)) then 
!!$          if ((1<=idx(i)).and.(idx(i) <= idxmap%local_rows)) then
!!$            idx(i) = idxmap%min_glob_row + idx(i) - 1
!!$          else if ((idxmap%local_rows < idx(i)).and.(idx(i) <= idxmap%local_cols)&
!!$               & .and.(.not.owned_)) then
!!$            idx(i) = idxmap%loc_to_glob(idx(i)-idxmap%local_rows)
!!$          else 
!!$            idx(i) = -1
!!$            info = -1
!!$          end if
!!$        end if
!!$      end do
!!$
!!$    else  if (.not.present(mask)) then 
!!$
!!$      do i=1, size(idx)
!!$        if ((1<=idx(i)).and.(idx(i) <= idxmap%local_rows)) then
!!$          idx(i) = idxmap%min_glob_row + idx(i) - 1
!!$        else if ((idxmap%local_rows < idx(i)).and.(idx(i) <= idxmap%local_cols)&
!!$             & .and.(.not.owned_)) then
!!$          idx(i) = idxmap%loc_to_glob(idx(i)-idxmap%local_rows)
!!$        else 
!!$          idx(i) = -1
!!$          info = -1
!!$        end if
!!$      end do
!!$
!!$    end if
!!$
!!$  end subroutine block_l2gv1
!!$
!!$  subroutine block_l2gv2(idxin,idxout,idxmap,info,mask,owned)
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(in) :: idxmap
!!$    integer(psb_ipk_), intent(in)    :: idxin(:)
!!$    integer(psb_ipk_), intent(out)   :: idxout(:)
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask(:)
!!$    logical, intent(in), optional :: owned
!!$    integer(psb_ipk_) :: is, im, i
!!$    logical :: owned_
!!$
!!$    info = 0
!!$
!!$    is = size(idxin)
!!$    im = min(is,size(idxout))
!!$
!!$    if (present(mask)) then 
!!$      if (size(mask) < im) then 
!!$        info = -1
!!$        return
!!$      end if
!!$    end if
!!$    if (present(owned)) then 
!!$      owned_ = owned
!!$    else
!!$      owned_ = .false.
!!$    end if
!!$
!!$    if (present(mask)) then 
!!$
!!$      do i=1, im
!!$        if (mask(i)) then 
!!$          if ((1<=idxin(i)).and.(idxin(i) <= idxmap%local_rows)) then
!!$            idxout(i) = idxmap%min_glob_row + idxin(i) - 1
!!$          else if ((idxmap%local_rows < idxin(i)).and.(idxin(i) <= idxmap%local_cols)&
!!$               & .and.(.not.owned_)) then
!!$            idxout(i) = idxmap%loc_to_glob(idxin(i)-idxmap%local_rows)
!!$          else 
!!$            idxout(i) = -1
!!$            info = -1
!!$          end if
!!$        end if
!!$      end do
!!$
!!$    else  if (.not.present(mask)) then 
!!$
!!$      do i=1, im
!!$        if ((1<=idxin(i)).and.(idxin(i) <= idxmap%local_rows)) then
!!$          idxout(i) = idxmap%min_glob_row + idxin(i) - 1
!!$        else if ((idxmap%local_rows < idxin(i)).and.(idxin(i) <= idxmap%local_cols)&
!!$             & .and.(.not.owned_)) then
!!$          idxout(i) = idxmap%loc_to_glob(idxin(i)-idxmap%local_rows)
!!$        else 
!!$          idxout(i) = -1
!!$          info = -1
!!$        end if
!!$      end do
!!$
!!$    end if
!!$
!!$    if (is > im) then 
!!$      info = -3 
!!$    end if
!!$
!!$  end subroutine block_l2gv2
!!$

  subroutine block_ll2gs1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
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

  end subroutine block_ll2gs1

  subroutine block_ll2gs2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_lpk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    idxout = idxin
    call idxmap%l2gip(idxout,info,mask,owned)
    
  end subroutine block_ll2gs2


  subroutine block_ll2gv1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
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

      do i=1, size(idx)
        if (mask(i)) then 
          if ((1<=idx(i)).and.(idx(i) <= idxmap%local_rows)) then
            idx(i) = idxmap%min_glob_row + idx(i) - 1
          else if ((idxmap%local_rows < idx(i)).and.(idx(i) <= idxmap%local_cols)&
               & .and.(.not.owned_)) then
            idx(i) = idxmap%loc_to_glob(idx(i)-idxmap%local_rows)
          else 
            idx(i) = -1
            info = -1
          end if
        end if
      end do

    else  if (.not.present(mask)) then 

      do i=1, size(idx)
        if ((1<=idx(i)).and.(idx(i) <= idxmap%local_rows)) then
          idx(i) = idxmap%min_glob_row + idx(i) - 1
        else if ((idxmap%local_rows < idx(i)).and.(idx(i) <= idxmap%local_cols)&
             & .and.(.not.owned_)) then
          idx(i) = idxmap%loc_to_glob(idx(i)-idxmap%local_rows)
        else 
          idx(i) = -1
          info = -1
        end if
      end do

    end if

  end subroutine block_ll2gv1

  subroutine block_ll2gv2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_lpk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: is, im, i
    logical :: owned_

    info = 0

    is = size(idxin)
    im = min(is,size(idxout))

    if (present(mask)) then 
      if (size(mask) < im) then 
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

      do i=1, im
        if (mask(i)) then 
          if ((1<=idxin(i)).and.(idxin(i) <= idxmap%local_rows)) then
            idxout(i) = idxmap%min_glob_row + idxin(i) - 1
          else if ((idxmap%local_rows < idxin(i)).and.(idxin(i) <= idxmap%local_cols)&
               & .and.(.not.owned_)) then
            idxout(i) = idxmap%loc_to_glob(idxin(i)-idxmap%local_rows)
          else 
            idxout(i) = -1
            info = -1
          end if
        end if
      end do

    else  if (.not.present(mask)) then 

      do i=1, im
        if ((1<=idxin(i)).and.(idxin(i) <= idxmap%local_rows)) then
          idxout(i) = idxmap%min_glob_row + idxin(i) - 1
        else if ((idxmap%local_rows < idxin(i)).and.(idxin(i) <= idxmap%local_cols)&
             & .and.(.not.owned_)) then
          idxout(i) = idxmap%loc_to_glob(idxin(i)-idxmap%local_rows)
        else 
          idxout(i) = -1
          info = -1
        end if
      end do

    end if

    if (is > im) then 
      info = -3 
    end if

  end subroutine block_ll2gv2

!!$  subroutine block_g2ls1(idx,idxmap,info,mask,owned)
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(in) :: idxmap
!!$    integer(psb_ipk_), intent(inout) :: idx
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask
!!$    logical, intent(in), optional :: owned
!!$    integer(psb_ipk_) :: idxv(1)
!!$    info = 0
!!$
!!$    if (present(mask)) then 
!!$      if (.not.mask) return
!!$    end if
!!$    
!!$    idxv(1) = idx 
!!$    call idxmap%g2lip(idxv,info,owned=owned)
!!$    idx = idxv(1) 
!!$      
!!$  end subroutine block_g2ls1
!!$
!!$  subroutine block_g2ls2(idxin,idxout,idxmap,info,mask,owned)
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(in) :: idxmap
!!$    integer(psb_ipk_), intent(in)    :: idxin
!!$    integer(psb_ipk_), intent(out)   :: idxout
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask
!!$    logical, intent(in), optional :: owned
!!$
!!$    idxout = idxin
!!$    call idxmap%g2lip(idxout,info,mask,owned)
!!$    
!!$  end subroutine block_g2ls2
!!$
!!$
!!$  subroutine block_g2lv1(idx,idxmap,info,mask,owned)
!!$    use psb_penv_mod
!!$    use psb_sort_mod
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(in) :: idxmap
!!$    integer(psb_ipk_), intent(inout) :: idx(:)
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask(:)
!!$    logical, intent(in), optional :: owned
!!$    integer(psb_ipk_) :: i, nv, is, ip, lip 
!!$    integer(psb_lpk_) :: tidx
!!$    integer(psb_mpk_) :: ictxt, iam, np
!!$    logical :: owned_
!!$
!!$    info = 0
!!$    ictxt = idxmap%get_ctxt()
!!$    call psb_info(ictxt,iam,np) 
!!$
!!$    if (present(mask)) then 
!!$      if (size(mask) < size(idx)) then 
!!$! !$        write(0,*) 'Block g2l: size of mask', size(mask),size(idx)
!!$        info = -1
!!$        return
!!$      end if
!!$    end if
!!$    if (present(owned)) then 
!!$      owned_ = owned
!!$    else
!!$      owned_ = .false.
!!$    end if
!!$
!!$    is = size(idx)
!!$    if (present(mask)) then 
!!$
!!$      if (idxmap%is_asb()) then 
!!$        do i=1, is
!!$          if (mask(i)) then 
!!$            if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
!!$              idx(i) = idx(i) - idxmap%min_glob_row + 1
!!$            else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)&
!!$                 &.and.(.not.owned_)) then
!!$              nv  = size(idxmap%srt_l2g,1)
!!$              tidx = idx(i)
!!$              idx(i) = psb_bsrch(tidx,nv,idxmap%srt_l2g(:,1))
!!$              if (idx(i) > 0) idx(i) = idxmap%srt_l2g(idx(i),2)+idxmap%local_rows
!!$            else 
!!$              idx(i) = -1
!!$            end if
!!$          end if
!!$        end do
!!$      else if (idxmap%is_valid()) then 
!!$        do i=1,is
!!$          if (mask(i)) then 
!!$            if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
!!$              idx(i) = idx(i) - idxmap%min_glob_row + 1
!!$            else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)&
!!$                 &.and.(.not.owned_)) then
!!$              ip = idx(i)
!!$              call psb_hash_searchkey(ip,lip,idxmap%hash,info)
!!$              if (lip > 0) idx(i) = lip + idxmap%local_rows
!!$            else 
!!$              idx(i) = -1
!!$            end if
!!$          end if
!!$        end do
!!$      else 
!!$! !$        write(0,*) 'Block status: invalid ',idxmap%get_state()
!!$        idx(1:is) = -1
!!$        info = -1
!!$      end if
!!$
!!$    else  if (.not.present(mask)) then 
!!$
!!$      if (idxmap%is_asb()) then 
!!$        do i=1, is
!!$          if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
!!$            idx(i) = idx(i) - idxmap%min_glob_row + 1
!!$          else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)&
!!$               &.and.(.not.owned_)) then
!!$            nv  = size(idxmap%srt_l2g,1)
!!$            tidx = idx(i)
!!$            idx(i) = psb_bsrch(tidx,nv,idxmap%srt_l2g(:,1))
!!$            if (idx(i) > 0) idx(i) = idxmap%srt_l2g(idx(i),2)+idxmap%local_rows
!!$          else 
!!$            idx(i) = -1
!!$          end if
!!$        end do
!!$
!!$      else if (idxmap%is_valid()) then 
!!$        do i=1,is
!!$          if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
!!$            idx(i) = idx(i) - idxmap%min_glob_row + 1
!!$          else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)&
!!$               &.and.(.not.owned_)) then
!!$            ip = idx(i)
!!$            call psb_hash_searchkey(ip,lip,idxmap%hash,info)
!!$            if (lip > 0) idx(i) = lip + idxmap%local_rows
!!$          else 
!!$            idx(i) = -1
!!$          end if
!!$        end do
!!$      else 
!!$! !$        write(0,*) 'Block status: invalid ',idxmap%get_state()
!!$        idx(1:is) = -1
!!$        info = -1
!!$      end if
!!$
!!$    end if
!!$
!!$  end subroutine block_g2lv1
!!$
!!$  subroutine block_g2lv2(idxin,idxout,idxmap,info,mask,owned)
!!$    use psb_penv_mod
!!$    use psb_sort_mod
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(in) :: idxmap
!!$    integer(psb_ipk_), intent(in)    :: idxin(:)
!!$    integer(psb_ipk_), intent(out)   :: idxout(:)
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask(:)
!!$    logical, intent(in), optional :: owned
!!$
!!$    integer(psb_ipk_) :: i, nv, is, ip, lip, im
!!$    integer(psb_lpk_) :: tidx
!!$    integer(psb_mpk_) :: ictxt, iam, np
!!$    logical :: owned_
!!$
!!$    info = 0
!!$    ictxt = idxmap%get_ctxt()
!!$    call psb_info(ictxt,iam,np) 
!!$    is = size(idxin)
!!$    im = min(is,size(idxout))
!!$
!!$    if (present(mask)) then 
!!$      if (size(mask) < im) then 
!!$! !$        write(0,*) 'Block g2l: size of mask', size(mask),size(idx)
!!$        info = -1
!!$        return
!!$      end if
!!$    end if
!!$    if (present(owned)) then 
!!$      owned_ = owned
!!$    else
!!$      owned_ = .false.
!!$    end if
!!$
!!$    if (present(mask)) then 
!!$
!!$      if (idxmap%is_asb()) then 
!!$        do i=1, im
!!$          if (mask(i)) then 
!!$            if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
!!$              idxout(i) = idxin(i) - idxmap%min_glob_row + 1
!!$            else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)&
!!$                 &.and.(.not.owned_)) then
!!$              nv  = size(idxmap%srt_l2g,1)
!!$              tidx = idxin(i)
!!$              idxout(i) = psb_bsrch(tidx,nv,idxmap%srt_l2g(:,1))
!!$              if (idxout(i) > 0) idxout(i) = idxmap%srt_l2g(idxout(i),2)+idxmap%local_rows
!!$            else 
!!$              idxout(i) = -1
!!$            end if
!!$          end if
!!$        end do
!!$      else if (idxmap%is_valid()) then 
!!$        do i=1,im
!!$          if (mask(i)) then 
!!$            if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
!!$              idxout(i) = idxin(i) - idxmap%min_glob_row + 1
!!$            else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)&
!!$                 &.and.(.not.owned_)) then
!!$              ip = idxin(i)
!!$              call psb_hash_searchkey(ip,lip,idxmap%hash,info)
!!$              if (lip > 0) idxout(i) = lip + idxmap%local_rows
!!$            else 
!!$              idxout(i) = -1
!!$            end if
!!$          end if
!!$        end do
!!$      else 
!!$! !$        write(0,*) 'Block status: invalid ',idxmap%get_state()
!!$        idxout(1:im) = -1
!!$        info = -1
!!$      end if
!!$
!!$    else  if (.not.present(mask)) then 
!!$
!!$      if (idxmap%is_asb()) then 
!!$        do i=1, im
!!$          if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
!!$            idxout(i) = idxin(i) - idxmap%min_glob_row + 1
!!$          else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)&
!!$               &.and.(.not.owned_)) then
!!$            nv  = size(idxmap%srt_l2g,1)
!!$            tidx = idxin(i)
!!$            idxout(i) = psb_bsrch(tidx,nv,idxmap%srt_l2g(:,1))
!!$            if (idxout(i) > 0) idxout(i) = idxmap%srt_l2g(idxout(i),2)+idxmap%local_rows
!!$          else 
!!$            idxout(i) = -1
!!$          end if
!!$        end do
!!$
!!$      else if (idxmap%is_valid()) then 
!!$        do i=1,im
!!$          if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
!!$            idxout(i) = idxin(i) - idxmap%min_glob_row + 1
!!$          else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)&
!!$               &.and.(.not.owned_)) then
!!$            ip = idxin(i)
!!$            call psb_hash_searchkey(ip,lip,idxmap%hash,info)
!!$            if (lip > 0) idxout(i) = lip + idxmap%local_rows
!!$          else 
!!$            idxout(i) = -1
!!$          end if
!!$        end do
!!$      else 
!!$! !$        write(0,*) 'Block status: invalid ',idxmap%get_state()
!!$        idxout(1:im) = -1
!!$        info = -1
!!$      end if
!!$
!!$    end if
!!$
!!$    if (is > im) info = -3 
!!$
!!$  end subroutine block_g2lv2


  subroutine block_lg2ls1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
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
      
  end subroutine block_lg2ls1

  subroutine block_lg2ls2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
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
      
  end subroutine block_lg2ls2


  subroutine block_lg2lv1(idx,idxmap,info,mask,owned)
    use psb_penv_mod
    use psb_sort_mod
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_ipk_) :: i, nv, is, ip, lip 
    integer(psb_lpk_) :: tidx
    integer(psb_mpk_) :: ictxt, iam, np
    logical :: owned_

    info = 0
    ictxt = idxmap%get_ctxt()
    call psb_info(ictxt,iam,np) 

    if (present(mask)) then 
      if (size(mask) < size(idx)) then 
!!$        write(0,*) 'Block g2l: size of mask', size(mask),size(idx)
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
    if (present(mask)) then 

      if (idxmap%is_asb()) then 
        do i=1, is
          if (mask(i)) then 
            if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
              idx(i) = idx(i) - idxmap%min_glob_row + 1
            else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)&
                 &.and.(.not.owned_)) then
              nv  = size(idxmap%srt_l2g,1)
              tidx = idx(i)
              idx(i) = psb_bsrch(tidx,nv,idxmap%srt_l2g(:,1))
              if (idx(i) > 0) idx(i) = idxmap%srt_l2g(idx(i),2)+idxmap%local_rows
            else 
              idx(i) = -1
            end if
          end if
        end do
      else if (idxmap%is_valid()) then 
        do i=1,is
          if (mask(i)) then 
            if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
              idx(i) = idx(i) - idxmap%min_glob_row + 1
            else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)&
                 &.and.(.not.owned_)) then
              ip = idx(i)
              call psb_hash_searchkey(ip,lip,idxmap%hash,info)
              if (lip > 0) idx(i) = lip + idxmap%local_rows
            else 
              idx(i) = -1
            end if
          end if
        end do
      else 
!!$        write(0,*) 'Block status: invalid ',idxmap%get_state()
        idx(1:is) = -1
        info = -1
      end if

    else  if (.not.present(mask)) then 

      if (idxmap%is_asb()) then 
        do i=1, is
          if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
            idx(i) = idx(i) - idxmap%min_glob_row + 1
          else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)&
               &.and.(.not.owned_)) then
            nv  = size(idxmap%srt_l2g,1)
            tidx = idx(i)
            idx(i) = psb_bsrch(tidx,nv,idxmap%srt_l2g(:,1))
            if (idx(i) > 0) idx(i) = idxmap%srt_l2g(idx(i),2)+idxmap%local_rows
          else 
            idx(i) = -1
          end if
        end do

      else if (idxmap%is_valid()) then 
        do i=1,is
          if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
            idx(i) = idx(i) - idxmap%min_glob_row + 1
          else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)&
               &.and.(.not.owned_)) then
            ip = idx(i)
            call psb_hash_searchkey(ip,lip,idxmap%hash,info)
            if (lip > 0) idx(i) = lip + idxmap%local_rows
          else 
            idx(i) = -1
          end if
        end do
      else 
!!$        write(0,*) 'Block status: invalid ',idxmap%get_state()
        idx(1:is) = -1
        info = -1
      end if

    end if

  end subroutine block_lg2lv1

  subroutine block_lg2lv2(idxin,idxout,idxmap,info,mask,owned)
    use psb_penv_mod
    use psb_sort_mod
    implicit none 
    class(psb_gen_block_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: i, nv, is, ip, lip, im
    integer(psb_lpk_) :: tidx
    integer(psb_mpk_) :: ictxt, iam, np
    logical :: owned_

    info = 0
    ictxt = idxmap%get_ctxt()
    call psb_info(ictxt,iam,np) 
    is = size(idxin)
    im = min(is,size(idxout))

    if (present(mask)) then 
      if (size(mask) < im) then 
!!$        write(0,*) 'Block g2l: size of mask', size(mask),size(idx)
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

      if (idxmap%is_asb()) then 
        do i=1, im
          if (mask(i)) then 
            if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
              idxout(i) = idxin(i) - idxmap%min_glob_row + 1
            else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)&
                 &.and.(.not.owned_)) then
              nv  = size(idxmap%srt_l2g,1)
              tidx = idxin(i)
              idxout(i) = psb_bsrch(tidx,nv,idxmap%srt_l2g(:,1))
              if (idxout(i) > 0) idxout(i) = idxmap%srt_l2g(idxout(i),2)+idxmap%local_rows
            else 
              idxout(i) = -1
            end if
          end if
        end do
      else if (idxmap%is_valid()) then 
        do i=1,im
          if (mask(i)) then 
            if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
              idxout(i) = idxin(i) - idxmap%min_glob_row + 1
            else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)&
                 &.and.(.not.owned_)) then
              ip = idxin(i)
              call psb_hash_searchkey(ip,lip,idxmap%hash,info)
              if (lip > 0) idxout(i) = lip + idxmap%local_rows
            else 
              idxout(i) = -1
            end if
          end if
        end do
      else 
!!$        write(0,*) 'Block status: invalid ',idxmap%get_state()
        idxout(1:im) = -1
        info = -1
      end if

    else  if (.not.present(mask)) then 

      if (idxmap%is_asb()) then 
        do i=1, im
          if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
            idxout(i) = idxin(i) - idxmap%min_glob_row + 1
          else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)&
               &.and.(.not.owned_)) then
            nv  = size(idxmap%srt_l2g,1)
            tidx = idxin(i)
            idxout(i) = psb_bsrch(tidx,nv,idxmap%srt_l2g(:,1))
            if (idxout(i) > 0) idxout(i) = idxmap%srt_l2g(idxout(i),2)+idxmap%local_rows
          else 
            idxout(i) = -1
          end if
        end do

      else if (idxmap%is_valid()) then 
        do i=1,im
          if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
            idxout(i) = idxin(i) - idxmap%min_glob_row + 1
          else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)&
               &.and.(.not.owned_)) then
            ip = idxin(i)
            call psb_hash_searchkey(ip,lip,idxmap%hash,info)
            if (lip > 0) idxout(i) = lip + idxmap%local_rows
          else 
            idxout(i) = -1
          end if
        end do
      else 
!!$        write(0,*) 'Block status: invalid ',idxmap%get_state()
        idxout(1:im) = -1
        info = -1
      end if

    end if

    if (is > im) info = -3 

  end subroutine block_lg2lv2

!!$  subroutine block_g2ls1_ins(idx,idxmap,info,mask, lidx)
!!$    use psb_realloc_mod
!!$    use psb_sort_mod
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(inout) :: idxmap
!!$    integer(psb_ipk_), intent(inout) :: idx
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask
!!$    integer(psb_ipk_), intent(in), optional :: lidx
!!$    
!!$    integer(psb_ipk_) :: idxv(1), lidxv(1)
!!$
!!$    info = 0
!!$    if (present(mask)) then 
!!$      if (.not.mask) return
!!$    end if
!!$    idxv(1) = idx
!!$    if (present(lidx)) then 
!!$      lidxv(1) = lidx
!!$      call idxmap%g2lip_ins(idxv,info,lidx=lidxv)
!!$    else
!!$      call idxmap%g2lip_ins(idxv,info)
!!$    end if
!!$    idx = idxv(1) 
!!$
!!$  end subroutine block_g2ls1_ins
!!$
!!$  subroutine block_g2ls2_ins(idxin,idxout,idxmap,info,mask,lidx)
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(inout) :: idxmap
!!$    integer(psb_ipk_), intent(in)    :: idxin
!!$    integer(psb_ipk_), intent(out)   :: idxout
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask
!!$    integer(psb_ipk_), intent(in), optional :: lidx
!!$
!!$    idxout = idxin
!!$    call idxmap%g2lip_ins(idxout,info,mask=mask,lidx=lidx)
!!$    
!!$  end subroutine block_g2ls2_ins
!!$
!!$
!!$  subroutine block_g2lv1_ins(idx,idxmap,info,mask,lidx)
!!$    use psb_realloc_mod
!!$    use psb_sort_mod
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(inout) :: idxmap
!!$    integer(psb_ipk_), intent(inout) :: idx(:)
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask(:)
!!$    integer(psb_ipk_), intent(in), optional :: lidx(:)
!!$
!!$    integer(psb_ipk_) :: i, nv, is, ix
!!$    integer(psb_ipk_) :: ip, lip, nxt
!!$
!!$
!!$    info = 0
!!$    is = size(idx)
!!$
!!$    if (present(mask)) then 
!!$      if (size(mask) < size(idx)) then 
!!$        info = -1
!!$        return
!!$      end if
!!$    end if
!!$    if (present(lidx)) then 
!!$      if (size(lidx) < size(idx)) then 
!!$        info = -1
!!$        return
!!$      end if
!!$    end if
!!$
!!$
!!$    if (idxmap%is_asb()) then 
!!$      ! State is wrong for this one ! 
!!$      idx = -1
!!$      info = -1
!!$
!!$    else if (idxmap%is_valid()) then 
!!$
!!$      if (present(lidx)) then 
!!$        if (present(mask)) then 
!!$
!!$          do i=1, is
!!$            if (mask(i)) then 
!!$              if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
!!$                idx(i) = idx(i) - idxmap%min_glob_row + 1
!!$              else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
!!$
!!$                if (lidx(i) <= idxmap%local_rows) then 
!!$                  info = -5
!!$                  return
!!$                end if
!!$                nxt = lidx(i)-idxmap%local_rows
!!$                ip  = idx(i) 
!!$                call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
!!$                if (info >= 0) then 
!!$                  if (lip == nxt) then 
!!$                    ! We have added one item
!!$                    call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
!!$                    if (info /= 0) then 
!!$                      info = -4
!!$                      return
!!$                    end if
!!$                    idxmap%local_cols       = max(lidx(i),idxmap%local_cols)
!!$                    idxmap%loc_to_glob(nxt) = idx(i)
!!$                  end if
!!$                  info = psb_success_
!!$                else
!!$                  info = -5
!!$                  return
!!$                end if
!!$                idx(i)  = lip + idxmap%local_rows
!!$              else 
!!$                idx(i) = -1
!!$                info = -1
!!$              end if
!!$            end if
!!$          end do
!!$
!!$        else if (.not.present(mask)) then 
!!$
!!$          do i=1, is
!!$
!!$            if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
!!$              idx(i) = idx(i) - idxmap%min_glob_row + 1
!!$            else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
!!$              if (lidx(i) <= idxmap%local_rows) then 
!!$                info = -5
!!$                return
!!$              end if
!!$              nxt = lidx(i)-idxmap%local_rows
!!$              ip = idx(i) 
!!$              call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
!!$
!!$              if (info >= 0) then 
!!$                if (lip == nxt) then 
!!$                  ! We have added one item
!!$                  call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
!!$                  if (info /= 0) then 
!!$                    info = -4
!!$                    return
!!$                  end if
!!$                  idxmap%local_cols       = max(lidx(i),idxmap%local_cols)
!!$                  idxmap%loc_to_glob(nxt) = idx(i)
!!$                end if
!!$                info = psb_success_
!!$              else
!!$                info = -5
!!$                return
!!$              end if
!!$              idx(i)  = lip + idxmap%local_rows
!!$            else 
!!$              idx(i) = -1
!!$              info = -1
!!$            end if
!!$          end do
!!$        end if
!!$
!!$      else if (.not.present(lidx)) then 
!!$
!!$        if (present(mask)) then 
!!$          do i=1, is
!!$            if (mask(i)) then 
!!$              if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
!!$                idx(i) = idx(i) - idxmap%min_glob_row + 1
!!$              else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
!!$                nv  = idxmap%local_cols-idxmap%local_rows
!!$                nxt = nv + 1 
!!$                ip = idx(i) 
!!$                call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
!!$                if (info >= 0) then 
!!$                  if (lip == nxt) then 
!!$                    ! We have added one item
!!$                    call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
!!$                    if (info /= 0) then 
!!$                      info = -4
!!$                      return
!!$                    end if
!!$                    idxmap%local_cols       = nxt + idxmap%local_rows
!!$                    idxmap%loc_to_glob(nxt) = idx(i)
!!$                  end if
!!$                  info = psb_success_
!!$                else
!!$                  info = -5
!!$                  return
!!$                end if
!!$                idx(i)  = lip + idxmap%local_rows
!!$              else 
!!$                idx(i) = -1
!!$                info = -1
!!$              end if
!!$            end if
!!$          end do
!!$
!!$        else if (.not.present(mask)) then 
!!$
!!$          do i=1, is
!!$
!!$            if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
!!$              idx(i) = idx(i) - idxmap%min_glob_row + 1
!!$            else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
!!$              nv  = idxmap%local_cols-idxmap%local_rows
!!$              nxt = nv + 1 
!!$              ip = idx(i) 
!!$              call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
!!$
!!$              if (info >= 0) then 
!!$                if (lip == nxt) then 
!!$                  ! We have added one item
!!$                  call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
!!$                  if (info /= 0) then 
!!$                    info = -4
!!$                    return
!!$                  end if
!!$                  idxmap%local_cols       = nxt + idxmap%local_rows
!!$                  idxmap%loc_to_glob(nxt) = idx(i)
!!$                end if
!!$                info = psb_success_
!!$              else
!!$                info = -5
!!$                return
!!$              end if
!!$              idx(i)  = lip + idxmap%local_rows
!!$            else 
!!$              idx(i) = -1
!!$              info = -1
!!$            end if
!!$          end do
!!$        end if
!!$      end if
!!$
!!$    else 
!!$      idx = -1
!!$      info = -1
!!$    end if
!!$
!!$  end subroutine block_g2lv1_ins
!!$
!!$  subroutine block_g2lv2_ins(idxin,idxout,idxmap,info,mask,lidx)
!!$    use psb_realloc_mod
!!$    use psb_sort_mod
!!$    implicit none 
!!$    class(psb_gen_block_map), intent(inout) :: idxmap
!!$    integer(psb_ipk_), intent(in)    :: idxin(:)
!!$    integer(psb_ipk_), intent(out)   :: idxout(:)
!!$    integer(psb_ipk_), intent(out)   :: info 
!!$    logical, intent(in), optional :: mask(:)
!!$    integer(psb_ipk_), intent(in), optional :: lidx(:)
!!$
!!$    integer(psb_ipk_) :: i, nv, is, ix, im
!!$    integer(psb_ipk_) :: ip, lip, nxt
!!$
!!$
!!$    info = 0
!!$    
!!$    is = size(idxin)
!!$    im = min(is,size(idxout))
!!$
!!$    if (present(mask)) then 
!!$      if (size(mask) < im) then 
!!$        info = -1
!!$        return
!!$      end if
!!$    end if
!!$    if (present(lidx)) then 
!!$      if (size(lidx) < im) then 
!!$        info = -1
!!$        return
!!$      end if
!!$    end if
!!$
!!$    if (idxmap%is_asb()) then 
!!$      ! State is wrong for this one ! 
!!$      idxout = -1
!!$      info   = -1
!!$
!!$    else if (idxmap%is_valid()) then 
!!$
!!$      if (present(lidx)) then 
!!$        if (present(mask)) then 
!!$
!!$          do i=1, im
!!$            if (mask(i)) then 
!!$              if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
!!$                idxout(i) = idxin(i) - idxmap%min_glob_row + 1
!!$              else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
!!$
!!$                if (lidx(i) <= idxmap%local_rows) then 
!!$                  info = -5
!!$                  return
!!$                end if
!!$                nxt = lidx(i)-idxmap%local_rows
!!$                ip  = idxin(i) 
!!$                call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
!!$                if (info >= 0) then 
!!$                  if (lip == nxt) then 
!!$                    ! We have added one item
!!$                    call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
!!$                    if (info /= 0) then 
!!$                      info = -4
!!$                      return
!!$                    end if
!!$                    idxmap%local_cols       = max(lidx(i),idxmap%local_cols)
!!$                    idxmap%loc_to_glob(nxt) = idxin(i)
!!$                  end if
!!$                  info = psb_success_
!!$                else
!!$                  info = -5
!!$                  return
!!$                end if
!!$                idxout(i)  = lip + idxmap%local_rows
!!$              else 
!!$                idxout(i) = -1
!!$                info      = -1
!!$              end if
!!$            end if
!!$          end do
!!$
!!$        else if (.not.present(mask)) then 
!!$
!!$          do i=1, im
!!$
!!$            if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
!!$              idxout(i) = idxin(i) - idxmap%min_glob_row + 1
!!$            else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
!!$              if (lidx(i) <= idxmap%local_rows) then 
!!$                info = -5
!!$                return
!!$              end if
!!$              nxt = lidx(i)-idxmap%local_rows
!!$              ip  = idxin(i) 
!!$              call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
!!$
!!$              if (info >= 0) then 
!!$                if (lip == nxt) then 
!!$                  ! We have added one item
!!$                  call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
!!$                  if (info /= 0) then 
!!$                    info = -4
!!$                    return
!!$                  end if
!!$                  idxmap%local_cols       = max(lidx(i),idxmap%local_cols)
!!$                  idxmap%loc_to_glob(nxt) = idxin(i)
!!$                end if
!!$                info = psb_success_
!!$              else
!!$                info = -5
!!$                return
!!$              end if
!!$              idxout(i)  = lip + idxmap%local_rows
!!$            else 
!!$              idxout(i) = -1
!!$              info      = -1
!!$            end if
!!$          end do
!!$        end if
!!$
!!$      else if (.not.present(lidx)) then 
!!$
!!$        if (present(mask)) then 
!!$          do i=1, im
!!$            if (mask(i)) then 
!!$              if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
!!$                idxout(i) = idxin(i) - idxmap%min_glob_row + 1
!!$              else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
!!$                nv  = idxmap%local_cols-idxmap%local_rows
!!$                nxt = nv + 1 
!!$                ip = idxin(i) 
!!$                call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
!!$                if (info >= 0) then 
!!$                  if (lip == nxt) then 
!!$                    ! We have added one item
!!$                    call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
!!$                    if (info /= 0) then 
!!$                      info = -4
!!$                      return
!!$                    end if
!!$                    idxmap%local_cols       = nxt + idxmap%local_rows
!!$                    idxmap%loc_to_glob(nxt) = idxin(i)
!!$                  end if
!!$                  info = psb_success_
!!$                else
!!$                  info = -5
!!$                  return
!!$                end if
!!$                idxout(i)  = lip + idxmap%local_rows
!!$              else 
!!$                idxout(i) = -1
!!$                info      = -1
!!$              end if
!!$            end if
!!$          end do
!!$
!!$        else if (.not.present(mask)) then 
!!$
!!$          do i=1, im
!!$
!!$            if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
!!$              idxout(i) = idxin(i) - idxmap%min_glob_row + 1
!!$            else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
!!$              nv  = idxmap%local_cols-idxmap%local_rows
!!$              nxt = nv + 1 
!!$              ip = idxin(i) 
!!$              call psb_hash_searchinskey(ip,lip,nxt,idxmap%hash,info)
!!$
!!$              if (info >= 0) then 
!!$                if (lip == nxt) then 
!!$                  ! We have added one item
!!$                  call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
!!$                  if (info /= 0) then 
!!$                    info = -4
!!$                    return
!!$                  end if
!!$                  idxmap%local_cols       = nxt + idxmap%local_rows
!!$                  idxmap%loc_to_glob(nxt) = idxin(i)
!!$                end if
!!$                info = psb_success_
!!$              else
!!$                info = -5
!!$                return
!!$              end if
!!$              idxout(i)  = lip + idxmap%local_rows
!!$            else 
!!$              idxout(i) = -1
!!$              info      = -1
!!$            end if
!!$          end do
!!$        end if
!!$      end if
!!$
!!$    else 
!!$      idxout = -1
!!$      info   = -1
!!$    end if
!!$
!!$    if (is > im) then 
!!$! !$      write(0,*) 'g2lv2_ins err -3'
!!$      info = -3 
!!$    end if
!!$
!!$  end subroutine block_g2lv2_ins

  subroutine block_lg2ls1_ins(idx,idxmap,info,mask, lidx)
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_gen_block_map), intent(inout) :: idxmap
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

  end subroutine block_lg2ls1_ins

  subroutine block_lg2ls2_ins(idxin,idxout,idxmap,info,mask,lidx)
    implicit none 
    class(psb_gen_block_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin
    integer(psb_ipk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    integer(psb_ipk_), intent(in), optional :: lidx
    integer(psb_lpk_) :: tidx
    tidx = idxin
    call idxmap%g2lip_ins(tidx,info,mask=mask,lidx=lidx)
    idxout = tidx
  end subroutine block_lg2ls2_ins


  subroutine block_lg2lv1_ins(idx,idxmap,info,mask,lidx)
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_gen_block_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: i, nv, is, ix
    integer(psb_lpk_) :: ip, lip, lnxt
    integer(psb_ipk_) :: nxt


    info = 0
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


    if (idxmap%is_asb()) then 
      ! State is wrong for this one ! 
      idx = -1
      info = -1

    else if (idxmap%is_valid()) then 

      if (present(lidx)) then 
        if (present(mask)) then 

          do i=1, is
            if (mask(i)) then 
              if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
                idx(i) = idx(i) - idxmap%min_glob_row + 1
              else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then

                if (lidx(i) <= idxmap%local_rows) then 
                  info = -5
                  return
                end if
                lnxt = lidx(i)-idxmap%local_rows
                ip   = idx(i) 
                call psb_hash_searchinskey(ip,lip,lnxt,idxmap%hash,info)
                nxt = lnxt
                if (info >= 0) then 
                  if (lip == nxt) then 
                    ! We have added one item
                    call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
                    if (info /= 0) then 
                      info = -4
                      return
                    end if
                    idxmap%local_cols       = max(lidx(i),idxmap%local_cols)
                    idxmap%loc_to_glob(nxt) = idx(i)
                  end if
                  info = psb_success_
                else
                  info = -5
                  return
                end if
                idx(i)  = lip + idxmap%local_rows
              else 
                idx(i) = -1
                info = -1
              end if
            end if
          end do

        else if (.not.present(mask)) then 

          do i=1, is

            if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
              idx(i) = idx(i) - idxmap%min_glob_row + 1
            else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
              if (lidx(i) <= idxmap%local_rows) then 
                info = -5
                return
              end if
              lnxt = lidx(i)-idxmap%local_rows
              ip = idx(i) 
              call psb_hash_searchinskey(ip,lip,lnxt,idxmap%hash,info)
              nxt = lnxt 
              if (info >= 0) then 
                if (lip == nxt) then 
                  ! We have added one item
                  call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
                  if (info /= 0) then 
                    info = -4
                    return
                  end if
                  idxmap%local_cols       = max(lidx(i),idxmap%local_cols)
                  idxmap%loc_to_glob(nxt) = idx(i)
                end if
                info = psb_success_
              else
                info = -5
                return
              end if
              idx(i)  = lip + idxmap%local_rows
            else 
              idx(i) = -1
              info = -1
            end if
          end do
        end if

      else if (.not.present(lidx)) then 

        if (present(mask)) then 
          do i=1, is
            if (mask(i)) then 
              if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
                idx(i) = idx(i) - idxmap%min_glob_row + 1
              else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
                nv   = idxmap%local_cols-idxmap%local_rows
                lnxt = nv + 1 
                ip  = idx(i) 
                call psb_hash_searchinskey(ip,lip,lnxt,idxmap%hash,info)
                nxt = lnxt 
                if (info >= 0) then 
                  if (lip == nxt) then 
                    ! We have added one item
                    call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
                    if (info /= 0) then 
                      info = -4
                      return
                    end if
                    idxmap%local_cols       = nxt + idxmap%local_rows
                    idxmap%loc_to_glob(nxt) = idx(i)
                  end if
                  info = psb_success_
                else
                  info = -5
                  return
                end if
                idx(i)  = lip + idxmap%local_rows
              else 
                idx(i) = -1
                info = -1
              end if
            end if
          end do

        else if (.not.present(mask)) then 

          do i=1, is

            if ((idxmap%min_glob_row <= idx(i)).and.(idx(i) <= idxmap%max_glob_row)) then
              idx(i) = idx(i) - idxmap%min_glob_row + 1
            else if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
              nv   = idxmap%local_cols-idxmap%local_rows
              lnxt = nv + 1 
              ip  = idx(i) 
              call psb_hash_searchinskey(ip,lip,lnxt,idxmap%hash,info)
              nxt = lnxt
              if (info >= 0) then 
                if (lip == nxt) then 
                  ! We have added one item
                  call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
                  if (info /= 0) then 
                    info = -4
                    return
                  end if
                  idxmap%local_cols       = nxt + idxmap%local_rows
                  idxmap%loc_to_glob(nxt) = idx(i)
                end if
                info = psb_success_
              else
                info = -5
                return
              end if
              idx(i)  = lip + idxmap%local_rows
            else 
              idx(i) = -1
              info = -1
            end if
          end do
        end if
      end if

    else 
      idx = -1
      info = -1
    end if

  end subroutine block_lg2lv1_ins

  subroutine block_lg2lv2_ins(idxin,idxout,idxmap,info,mask,lidx)
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_gen_block_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: i, nv, is, ix, im
    integer(psb_lpk_) :: ip, lip, lnxt
    integer(psb_ipk_) :: nxt


    info = 0
    
    is = size(idxin)
    im = min(is,size(idxout))

    if (present(mask)) then 
      if (size(mask) < im) then 
        info = -1
        return
      end if
    end if
    if (present(lidx)) then 
      if (size(lidx) < im) then 
        info = -1
        return
      end if
    end if

    if (idxmap%is_asb()) then 
      ! State is wrong for this one ! 
      idxout = -1
      info   = -1

    else if (idxmap%is_valid()) then 

      if (present(lidx)) then 
        if (present(mask)) then 

          do i=1, im
            if (mask(i)) then 
              if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
                idxout(i) = idxin(i) - idxmap%min_glob_row + 1
              else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then

                if (lidx(i) <= idxmap%local_rows) then 
                  info = -5
                  return
                end if
                lnxt = lidx(i)-idxmap%local_rows
                ip  = idxin(i) 
                call psb_hash_searchinskey(ip,lip,lnxt,idxmap%hash,info)
                nxt = lnxt
                if (info >= 0) then 
                  if (lip == nxt) then
                    ! We have added one item
                    call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
                    if (info /= 0) then 
                      info = -4
                      return
                    end if
                    idxmap%local_cols       = max(lidx(i),idxmap%local_cols)
                    idxmap%loc_to_glob(nxt) = idxin(i)
                  end if
                  info = psb_success_
                else
                  info = -5
                  return
                end if
                idxout(i)  = lip + idxmap%local_rows
              else 
                idxout(i) = -1
                info      = -1
              end if
            end if
          end do

        else if (.not.present(mask)) then 

          do i=1, im

            if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
              idxout(i) = idxin(i) - idxmap%min_glob_row + 1
            else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
              if (lidx(i) <= idxmap%local_rows) then 
                info = -5
                return
              end if
              lnxt = lidx(i)-idxmap%local_rows
              ip  = idxin(i) 
              call psb_hash_searchinskey(ip,lip,lnxt,idxmap%hash,info)
              nxt = lnxt
              if (info >= 0) then 
                if (lip == nxt) then 
                  ! We have added one item
                  call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
                  if (info /= 0) then 
                    info = -4
                    return
                  end if
                  idxmap%local_cols       = max(lidx(i),idxmap%local_cols)
                  idxmap%loc_to_glob(nxt) = idxin(i)
                end if
                info = psb_success_
              else
                info = -5
                return
              end if
              idxout(i)  = lip + idxmap%local_rows
            else 
              idxout(i) = -1
              info      = -1
            end if
          end do
        end if

      else if (.not.present(lidx)) then 

        if (present(mask)) then 
          do i=1, im
            if (mask(i)) then 
              if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
                idxout(i) = idxin(i) - idxmap%min_glob_row + 1
              else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
                nv  = idxmap%local_cols-idxmap%local_rows
                lnxt = nv + 1 
                ip = idxin(i) 
                call psb_hash_searchinskey(ip,lip,lnxt,idxmap%hash,info)
                nxt = lnxt
                if (info >= 0) then 
                  if (lip == nxt) then 
                    ! We have added one item
                    call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
                    if (info /= 0) then 
                      info = -4
                      return
                    end if
                    idxmap%local_cols       = nxt + idxmap%local_rows
                    idxmap%loc_to_glob(nxt) = idxin(i)
                  end if
                  info = psb_success_
                else
                  info = -5
                  return
                end if
                idxout(i)  = lip + idxmap%local_rows
              else 
                idxout(i) = -1
                info      = -1
              end if
            end if
          end do

        else if (.not.present(mask)) then 

          do i=1, im

            if ((idxmap%min_glob_row <= idxin(i)).and.(idxin(i) <= idxmap%max_glob_row)) then
              idxout(i) = idxin(i) - idxmap%min_glob_row + 1
            else if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
              nv  = idxmap%local_cols-idxmap%local_rows
              lnxt = nv + 1 
              ip = idxin(i) 
              call psb_hash_searchinskey(ip,lip,lnxt,idxmap%hash,info)
              nxt = lnxt
              if (info >= 0) then 
                if (lip == nxt) then 
                  ! We have added one item
                  call psb_ensure_size(nxt,idxmap%loc_to_glob,info,addsz=laddsz)
                  if (info /= 0) then 
                    info = -4
                    return
                  end if
                  idxmap%local_cols       = nxt + idxmap%local_rows
                  idxmap%loc_to_glob(nxt) = idxin(i)
                end if
                info = psb_success_
              else
                info = -5
                return
              end if
              idxout(i)  = lip + idxmap%local_rows
            else 
              idxout(i) = -1
              info      = -1
            end if
          end do
        end if
      end if

    else 
      idxout = -1
      info   = -1
    end if

    if (is > im) then 
!!$      write(0,*) 'g2lv2_ins err -3'
      info = -3 
    end if

  end subroutine block_lg2lv2_ins

  subroutine block_fnd_owner(idx,iprc,idxmap,info)
    use psb_penv_mod
    implicit none 
    integer(psb_lpk_), intent(in) :: idx(:)
    integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
    class(psb_gen_block_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: ictxt, iam, np, nv, ip, i
    integer(psb_lpk_) :: tidx
    
    ictxt = idxmap%get_ctxt()
    call psb_info(ictxt,iam,np)
    nv = size(idx)
    allocate(iprc(nv),stat=info) 
    if (info /= 0) then 
!!$      write(0,*) 'Memory allocation failure in repl_map_fnd-owner'
      return
    end if
    do i=1, nv
      tidx = idx(i) 
      ip = gen_block_search(tidx-1,np+1,idxmap%vnl)
      iprc(i) = ip - 1
    end do

  end subroutine block_fnd_owner



  subroutine block_init(idxmap,ictxt,nl,info)
    use psb_penv_mod
    use psb_realloc_mod
    use psb_error_mod
    implicit none 
    class(psb_gen_block_map), intent(inout) :: idxmap
    integer(psb_mpk_), intent(in)  :: ictxt
    integer(psb_ipk_), intent(in)  :: nl
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_mpk_) :: iam, np
    integer(psb_ipk_) :: i, ntot
    integer(psb_lpk_), allocatable :: vnl(:)

    info = 0
    call psb_info(ictxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ictxt:',ictxt
      info = -1
      return
    end if
    allocate(vnl(0:np),stat=info)
    if (info /= 0)  then
      info = -2
      return
    end if
    
    vnl(:)   = 0
    vnl(iam) = nl
    call psb_sum(ictxt,vnl)
    ntot = sum(vnl)
    vnl(1:np) = vnl(0:np-1)
    vnl(0) = 0
    do i=1,np
      vnl(i) = vnl(i) + vnl(i-1)
    end do
    if (ntot /= vnl(np)) then 
!!$      write(0,*) ' Mismatch in block_init ',ntot,vnl(np)
    end if
    
    idxmap%global_rows  = ntot
    idxmap%global_cols  = ntot
    idxmap%local_rows   = nl
    idxmap%local_cols   = nl
    idxmap%ictxt        = ictxt
    idxmap%state        = psb_desc_bld_
    call psb_get_mpicomm(ictxt,idxmap%mpic)
    idxmap%min_glob_row = vnl(iam)+1
    idxmap%max_glob_row = vnl(iam+1) 
    call move_alloc(vnl,idxmap%vnl)
    call psb_realloc(nl,idxmap%loc_to_glob,info) 
    if (info /= 0)  then
      info = -2
      return
    end if
    call psb_hash_init(nl,idxmap%hash,info)
    call idxmap%set_state(psb_desc_bld_)
    
  end subroutine block_init


  subroutine block_asb(idxmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_gen_block_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(out) :: info
    
    integer(psb_ipk_) :: nhal
    integer(psb_mpk_) :: ictxt, iam, np 
    
    info = 0 
    ictxt = idxmap%get_ctxt()
    call psb_info(ictxt,iam,np)

    nhal = idxmap%local_cols-idxmap%local_rows

    call psb_realloc(nhal,idxmap%loc_to_glob,info)
    call psb_realloc(nhal,2,idxmap%srt_l2g,info)
    idxmap%srt_l2g(1:nhal,1) = idxmap%loc_to_glob(1:nhal)

    call psb_msort(idxmap%srt_l2g(:,1),&
         & ix=idxmap%srt_l2g(:,2),dir=psb_sort_up_)

    call psb_free(idxmap%hash,info)
    call idxmap%set_state(psb_desc_asb_)
  end subroutine block_asb

  function block_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'BLOCK'
  end function block_get_fmt


  subroutine block_clone(idxmap,outmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_gen_block_map), intent(inout)    :: idxmap
    class(psb_indx_map), allocatable, intent(out) :: outmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='block_clone'
    logical, parameter :: debug=.false.

    info = psb_success_
    call psb_get_erraction(err_act)
    if (allocated(outmap)) then 
      write(0,*) 'Error: should not be allocated on input'
      info = -87
      goto 9999
    end if
    
    allocate(psb_gen_block_map :: outmap, stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if

    select type (outmap)
    type is (psb_gen_block_map) 
      if (info == psb_success_) then 
        outmap%psb_indx_map = idxmap%psb_indx_map
        outmap%min_glob_row = idxmap%min_glob_row
        outmap%max_glob_row = idxmap%max_glob_row
      end if
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%loc_to_glob,outmap%loc_to_glob,info)
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%vnl,outmap%vnl,info)
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%srt_l2g,outmap%srt_l2g,info)
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
  end subroutine block_clone


  subroutine block_reinit(idxmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_gen_block_map), intent(inout)    :: idxmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act, nr,nc,k, nl, ictxt
    integer(psb_ipk_), allocatable :: lidx(:)
    integer(psb_lpk_), allocatable :: idx(:)
    character(len=20)  :: name='block_reinit'
    logical, parameter :: debug=.false.

    info = psb_success_
    call psb_get_erraction(err_act)

    nr = idxmap%get_lr()
    nc = idxmap%get_lc()
    if (nc>nr) then 
      lidx = (/(k,k=nr+1,nc)/)
      idx  = (/(k,k=nr+1,nc)/)
      call idxmap%l2gip(idx,info)
    end if
    if (info /= 0) &
         & write(0,*) 'From l2gip',info

    
    call psb_hash_init(nr,idxmap%hash,info)
    if (info /= 0) &
         & write(0,*) 'From hash_init',info
    call idxmap%set_state(psb_desc_bld_)
    if (nc>nr) then 
      call idxmap%g2lip_ins(idx,info,lidx=lidx)
    end if

      
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return


9999 call psb_error_handler(err_act)

    return
  end subroutine block_reinit


  !
  ! This is a purely internal version of "binary" search
  ! specialized for gen_block usage.
  !
  function  gen_block_search(key,n,v) result(ipos)
    implicit none
    integer(psb_lpk_) :: key
    integer(psb_ipk_) :: ipos, n
    integer(psb_ipk_) :: v(:)

    integer(psb_ipk_) :: lb, ub, m

    if (n < 5) then 
      ! don't bother with binary search for very
      ! small vectors
      ipos = 0
      do
        if (ipos == n) return
        if (key < v(ipos+1)) return 
        ipos = ipos + 1 
      end do
    else
      lb = 1 
      ub = n
      ipos = -1 

      do while (lb <= ub) 
        m = (lb+ub)/2
        if (key==v(m))  then
          ipos = m 
          return
        else if (key < v(m))  then
          ub = m-1
        else 
          lb = m + 1
        end if
      enddo
      if (v(ub) > key) then
        ub = ub - 1 
      end if
      ipos = ub 
    endif
    return
  end function gen_block_search

  function  l_gen_block_search(key,n,v) result(ipos)
    implicit none
    integer(psb_ipk_) :: ipos, n
    integer(psb_lpk_) :: key
    integer(psb_lpk_) :: v(:)

    integer(psb_ipk_) :: lb, ub, m

    if (n < 5) then 
      ! don't bother with binary search for very
      ! small vectors
      ipos = 0
      do
        if (ipos == n) return
        if (key < v(ipos+1)) return 
        ipos = ipos + 1 
      end do
    else
      lb = 1 
      ub = n
      ipos = -1 

      do while (lb <= ub) 
        m = (lb+ub)/2
        if (key==v(m))  then
          ipos = m 
          return
        else if (key < v(m))  then
          ub = m-1
        else 
          lb = m + 1
        end if
      enddo
      if (v(ub) > key) then
        ub = ub - 1 
      end if
      ipos = ub 
    endif
    return
  end function l_gen_block_search


end module psb_gen_block_map_mod
