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
! package: psb_list_map_mod
!    Defines the LIST_MAP type.
!
! This is essentially the original PSBLAS index map. We assume that
! 1. We have room for GLOB_TO_LOC and LOC_TO_GLOB
! 2. There could be an overlap, so we don't store explicitly who owns an index.
!
!
module psb_list_map_mod
  use psb_const_mod
  use psb_desc_const_mod
  use psb_indx_map_mod
  
  type, extends(psb_indx_map) :: psb_list_map
    integer(psb_ipk_) :: pnt_h          = -1 
    integer(psb_lpk_), allocatable :: loc_to_glob(:)
    integer(psb_ipk_), allocatable :: glob_to_loc(:)
  contains
    procedure, pass(idxmap)  :: init_vl    => list_initlvl

    procedure, pass(idxmap)  :: sizeof    => list_sizeof
    procedure, pass(idxmap)  :: asb       => list_asb
    procedure, pass(idxmap)  :: free      => list_free
    procedure, pass(idxmap)  :: clone     => list_clone
    procedure, pass(idxmap)  :: reinit    => list_reinit
    procedure, nopass        :: get_fmt   => list_get_fmt
    procedure, nopass        :: row_extendable => list_row_extendable

    procedure, pass(idxmap)  :: ll2gs1 => list_ll2gs1
    procedure, pass(idxmap)  :: ll2gs2 => list_ll2gs2
    procedure, pass(idxmap)  :: ll2gv1 => list_ll2gv1
    procedure, pass(idxmap)  :: ll2gv2 => list_ll2gv2

    procedure, pass(idxmap)  :: lg2ls1 => list_lg2ls1
    procedure, pass(idxmap)  :: lg2ls2 => list_lg2ls2
    procedure, pass(idxmap)  :: lg2lv1 => list_lg2lv1
    procedure, pass(idxmap)  :: lg2lv2 => list_lg2lv2

    procedure, pass(idxmap)  :: lg2ls1_ins => list_lg2ls1_ins
    procedure, pass(idxmap)  :: lg2ls2_ins => list_lg2ls2_ins
    procedure, pass(idxmap)  :: lg2lv1_ins => list_lg2lv1_ins
    procedure, pass(idxmap)  :: lg2lv2_ins => list_lg2lv2_ins

  end type psb_list_map

  private :: list_initvl, list_sizeof, list_asb, list_free,&
       & list_get_fmt, list_l2gs1, list_l2gs2, list_l2gv1,&
       & list_l2gv2, list_g2ls1, list_g2ls2, list_g2lv1,&
       & list_g2lv2, list_g2ls1_ins, list_g2ls2_ins,&
       & list_g2lv1_ins, list_g2lv2_ins, list_row_extendable

  integer(psb_ipk_), private :: laddsz=500

contains
    
  function list_row_extendable() result(val)
    implicit none 
    logical :: val
    val = .true.
  end function list_row_extendable

  function list_sizeof(idxmap) result(val)
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
    integer(psb_epk_) :: val
    
    val = idxmap%psb_indx_map%sizeof()

    if (allocated(idxmap%loc_to_glob)) &
         & val = val + size(idxmap%loc_to_glob)*psb_sizeof_ip
    if (allocated(idxmap%glob_to_loc)) &
         & val = val + size(idxmap%glob_to_loc)*psb_sizeof_ip

  end function list_sizeof


  subroutine list_free(idxmap)
    implicit none 
    class(psb_list_map), intent(inout) :: idxmap
    
    if (allocated(idxmap%loc_to_glob)) &
         & deallocate(idxmap%loc_to_glob)
    if (allocated(idxmap%glob_to_loc)) &
         & deallocate(idxmap%glob_to_loc)

    call idxmap%psb_indx_map%free()

  end subroutine list_free

  subroutine list_ll2gs1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
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

  end subroutine list_ll2gs1

  subroutine list_ll2gs2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin
    integer(psb_lpk_), intent(out)   :: idxout
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask
    logical, intent(in), optional :: owned

    idxout = idxin
    call idxmap%l2gip(idxout,info,mask,owned)
    
  end subroutine list_ll2gs2


  subroutine list_ll2gv1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_lpk_) :: i
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
          if ((1<=idx(i)).and.(idx(i) <= idxmap%get_lr())) then
            idx(i) = idxmap%loc_to_glob(idx(i))
          else if ((idxmap%get_lr() < idx(i)).and.(idx(i) <= idxmap%local_cols)&
               & .and.(.not.owned_)) then
            idx(i) = idxmap%loc_to_glob(idx(i))
          else 
            idx(i) = -1
          end if
        end if
      end do

    else  if (.not.present(mask)) then 

      do i=1, size(idx)
        if ((1<=idx(i)).and.(idx(i) <= idxmap%get_lr())) then
          idx(i) = idxmap%loc_to_glob(idx(i))
        else if ((idxmap%get_lr() < idx(i)).and.(idx(i) <= idxmap%local_cols)&
             & .and.(.not.owned_)) then
          idx(i) = idxmap%loc_to_glob(idx(i))
        else 
          idx(i) = -1
        end if
      end do

    end if

  end subroutine list_ll2gv1

  subroutine list_ll2gv2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(in)    :: idxin(:)
    integer(psb_lpk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_lpk_) :: is, im
    
    is = size(idxin)
    im = min(is,size(idxout))
    idxout(1:im) = idxin(1:im)
    call idxmap%l2gip(idxout(1:im),info,mask,owned)
    if (is > im) info = -3 

  end subroutine list_ll2gv2

  subroutine list_lg2ls1(idx,idxmap,info,mask,owned)
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
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
      
  end subroutine list_lg2ls1

  subroutine list_lg2ls2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
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

    idxv = idxin
    call idxmap%g2lip(idxv,info,owned=owned)
    idxout = idxv(1)
    
  end subroutine list_lg2ls2


  subroutine list_lg2lv1(idx,idxmap,info,mask,owned)
    use psb_sort_mod
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned
    integer(psb_lpk_) :: i, is, ix
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

    is = size(idx)

    if (present(mask)) then 
      if (idxmap%is_valid()) then 
        do i=1,is
          if (mask(i)) then 
            if ((1 <= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
              ix = idxmap%glob_to_loc(idx(i))
              if ((ix > idxmap%get_lr()).and.(owned_)) ix = -1
              idx(i) = ix
            else 
              idx(i) = -1
            end if
          end if
        end do
      else 
        idx(1:is) = -1
        info = -1
      end if

    else  if (.not.present(mask)) then 

      if (idxmap%is_valid()) then 
        do i=1, is
          if ((1 <= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
            ix = idxmap%glob_to_loc(idx(i))
                if ((ix > idxmap%get_lr()).and.(owned_)) ix = -1
            idx(i) = ix
          else 
            idx(i) = -1
          end if
        end do
      else 
        idx(1:is) = -1
        info = -1
      end if
 
    end if

  end subroutine list_lg2lv1

  subroutine list_lg2lv2(idxin,idxout,idxmap,info,mask,owned)
    implicit none 
    class(psb_list_map), intent(in) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    logical, intent(in), optional :: owned

    integer(psb_ipk_) :: im
    integer(psb_lpk_) :: i, is, ix
    logical :: owned_

    info = 0

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

    if (present(mask)) then 
      if (idxmap%is_valid()) then 
        do i=1,is
          if (mask(i)) then 
            if ((1 <= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
              ix = idxmap%glob_to_loc(idxin(i))
              if ((ix > idxmap%get_lr()).and.(owned_)) ix = -1
              idxout(i) = ix
            else 
              idxout(i) = -1
            end if
          end if
        end do
      else 
        idxout(1:is) = -1
        info = -1
      end if

    else  if (.not.present(mask)) then 

      if (idxmap%is_valid()) then 
        do i=1, is
          if ((1 <= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
            ix = idxmap%glob_to_loc(idxin(i))
                if ((ix > idxmap%get_lr()).and.(owned_)) ix = -1
            idxout(i) = ix
          else 
            idxout(i) = -1
          end if
        end do
      else 
        idxout(1:is) = -1
        info = -1
      end if
 
    end if

  end subroutine list_lg2lv2

  subroutine list_lg2ls1_ins(idx,idxmap,info,mask,lidx)
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_list_map), intent(inout) :: idxmap
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

  end subroutine list_lg2ls1_ins

  subroutine list_lg2ls2_ins(idxin,idxout,idxmap,info,mask,lidx)
    implicit none 
    class(psb_list_map), intent(inout) :: idxmap
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
    
  end subroutine list_lg2ls2_ins


  subroutine list_lg2lv1_ins(idx,idxmap,info,mask,lidx)
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_list_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(inout) :: idx(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: ix, lix
    integer(psb_lpk_) :: i, is

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
              if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
                ix = idxmap%glob_to_loc(idx(i))                
                if (ix < 0) then 
                  ix = lidx(i) 
                  call psb_ensure_size(ix,idxmap%loc_to_glob,info,addsz=laddsz)
                  if ((ix <= idxmap%local_rows).or.(info /= 0)) then 
                    info = -4
                    return
                  end if
                  idxmap%local_cols          = max(ix,idxmap%local_cols)
                  idxmap%loc_to_glob(ix)     = idx(i)
                  idxmap%glob_to_loc(idx(i)) = ix
                end if
                idx(i) = ix
              else 
                idx(i) = -1
              end if
            end if
          end do

        else if (.not.present(mask)) then 

          do i=1, is
            if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
              ix = idxmap%glob_to_loc(idx(i))
              if (ix < 0) then 
                ix = lidx(i) 
                call psb_ensure_size(ix,idxmap%loc_to_glob,info,addsz=laddsz)
                if ((ix <= idxmap%local_rows).or.(info /= 0)) then 
                  info = -4
                  return
                end if
                idxmap%local_cols          = max(ix,idxmap%local_cols)
                idxmap%loc_to_glob(ix)     = idx(i)
                idxmap%glob_to_loc(idx(i)) = ix
              end if
              idx(i) = ix
            else 
              idx(i) = -1
            end if
          end do
        end if

      else if (.not.present(lidx)) then

        if (present(mask)) then 
          do i=1, is
            if (mask(i)) then 
              if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
                ix = idxmap%glob_to_loc(idx(i))
                if (ix < 0) then 
                  ix = idxmap%local_cols + 1
                  call psb_ensure_size(ix,idxmap%loc_to_glob,info,addsz=laddsz)
                  if (info /= 0) then 
                    info = -4
                    return
                  end if
                  idxmap%local_cols      = ix
                  idxmap%loc_to_glob(ix) = idx(i)
                  idxmap%glob_to_loc(idx(i)) = ix
                end if
                idx(i) = ix
              else 
                idx(i) = -1
              end if
            end if
          end do

        else if (.not.present(mask)) then 

          do i=1, is
            if ((1<= idx(i)).and.(idx(i) <= idxmap%global_rows)) then
              ix = idxmap%glob_to_loc(idx(i))
              if (ix < 0) then 
                ix = idxmap%local_cols + 1
                call psb_ensure_size(ix,idxmap%loc_to_glob,info,addsz=laddsz)
                if (info /= 0) then 
                  info = -4
                  return
                end if
                idxmap%local_cols      = ix
                idxmap%loc_to_glob(ix) = idx(i)
                idxmap%glob_to_loc(idx(i)) = ix
              end if
              idx(i) = ix
            else 
              idx(i) = -1
            end if
          end do
        end if
      end if

    else 
      idx = -1
      info = -1
    end if

  end subroutine list_lg2lv1_ins

  subroutine list_lg2lv2_ins(idxin,idxout,idxmap,info,mask,lidx)
    use psb_realloc_mod
    implicit none 
    class(psb_list_map), intent(inout) :: idxmap
    integer(psb_lpk_), intent(in)    :: idxin(:)
    integer(psb_ipk_), intent(out)   :: idxout(:)
    integer(psb_ipk_), intent(out)   :: info 
    logical, intent(in), optional :: mask(:)
    integer(psb_ipk_), intent(in), optional :: lidx(:)

    integer(psb_ipk_) :: ix, lix
    integer(psb_lpk_) :: i, is

    info = 0
    is = min(size(idxin),size(idxout))

    if (present(mask)) then 
      if (size(mask) < size(idxin)) then 
        info = -1
        return
      end if
    end if
    if (present(lidx)) then 
      if (size(lidx) < size(idxin)) then 
        info = -1
        return
      end if
    end if


    if (idxmap%is_asb()) then 
      ! State is wrong for this one ! 
      idxout = -1
      info = -1

    else if (idxmap%is_valid()) then 

      if (present(lidx)) then 
        if (present(mask)) then 
          do i=1, is
            if (mask(i)) then 
              if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
                ix = idxmap%glob_to_loc(idxin(i))                
                if (ix < 0) then 
                  ix = lidx(i) 
                  call psb_ensure_size(ix,idxmap%loc_to_glob,info,addsz=laddsz)
                  if ((ix <= idxmap%local_rows).or.(info /= 0)) then 
                    info = -4
                    return
                  end if
                  idxmap%local_cols            = max(ix,idxmap%local_cols)
                  idxmap%loc_to_glob(ix)       = idxin(i)
                  idxmap%glob_to_loc(idxin(i)) = ix
                end if
                idxout(i) = ix
              else 
                idxout(i) = -1
              end if
            end if
          end do

        else if (.not.present(mask)) then 

          do i=1, is
            if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
              ix = idxmap%glob_to_loc(idxin(i))
              if (ix < 0) then 
                ix = lidx(i) 
                call psb_ensure_size(ix,idxmap%loc_to_glob,info,addsz=laddsz)
                if ((ix <= idxmap%local_rows).or.(info /= 0)) then 
                  info = -4
                  return
                end if
                idxmap%local_cols            = max(ix,idxmap%local_cols)
                idxmap%loc_to_glob(ix)       = idxin(i)
                idxmap%glob_to_loc(idxin(i)) = ix
              end if
              idxout(i) = ix
            else 
              idxout(i) = -1
            end if
          end do
        end if

      else if (.not.present(lidx)) then

        if (present(mask)) then 
          do i=1, is
            if (mask(i)) then 
              if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
                ix = idxmap%glob_to_loc(idxin(i))
                if (ix < 0) then 
                  ix = idxmap%local_cols + 1
                  call psb_ensure_size(ix,idxmap%loc_to_glob,info,addsz=laddsz)
                  if (info /= 0) then 
                    info = -4
                    return
                  end if
                  idxmap%local_cols            = ix
                  idxmap%loc_to_glob(ix)       = idxin(i)
                  idxmap%glob_to_loc(idxin(i)) = ix
                end if
                idxout(i) = ix
              else 
                idxout(i) = -1
              end if
            end if
          end do

        else if (.not.present(mask)) then 

          do i=1, is
            if ((1<= idxin(i)).and.(idxin(i) <= idxmap%global_rows)) then
              ix = idxmap%glob_to_loc(idxin(i))
              if (ix < 0) then 
                ix = idxmap%local_cols + 1
                call psb_ensure_size(ix,idxmap%loc_to_glob,info,addsz=laddsz)
                if (info /= 0) then 
                  info = -4
                  return
                end if
                idxmap%local_cols            = ix
                idxmap%loc_to_glob(ix)       = idxin(i)
                idxmap%glob_to_loc(idxin(i)) = ix
              end if
              idxout(i) = ix
            else 
              idxout(i) = -1
            end if
          end do
        end if
      end if

    else 
      idxout = -1
      info = -1
    end if

  end subroutine list_lg2lv2_ins

  subroutine list_initvl(idxmap,ctxt,vl,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_list_map), intent(inout) :: idxmap
    type(psb_ctxt_type), intent(in) :: ctxt
    integer(psb_ipk_), intent(in)  :: vl(:)
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_lpk_) :: nl
    integer(psb_lpk_), allocatable :: lvl(:)
    integer(psb_ipk_) :: iam, np

    info = 0
    call psb_info(ctxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ctxt:'
      info = -1
      return
    end if

    nl = size(vl) 
    allocate(lvl(nl),stat=info)
    if (info /= 0) then
      info = -1
      return
    end if

    lvl(1:nl) = vl(1:nl)
    call idxmap%init_vl(ctxt,lvl,info)
   
  end subroutine list_initvl


  subroutine list_initlvl(idxmap,ctxt,vl,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_list_map), intent(inout) :: idxmap
    type(psb_ctxt_type), intent(in) :: ctxt
    integer(psb_lpk_), intent(in)  :: vl(:)
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_lpk_) ::  i, ix, nl, n, nrt
    integer(psb_ipk_) :: iam, np

    info = 0
    call psb_info(ctxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ctxt:'
      info = -1
      return
    end if

    nl = size(vl) 
    

    n   = maxval(vl(1:nl))
    nrt = nl
    call psb_sum(ctxt,nrt)
    call psb_max(ctxt,n)


    if (n /= nrt) then 
      write(psb_err_unit,*) 'Size mismatch', n, nrt
      info = -1
      return
    end if
    
    idxmap%global_rows  = n
    idxmap%global_cols  = n

    allocate(idxmap%loc_to_glob(n),idxmap%glob_to_loc(n),stat=info) 
    if (info /= 0)  then
      info = -2
      return
    end if

    idxmap%ctxt = ctxt
    idxmap%state = psb_desc_bld_
    idxmap%mpic  = psb_get_mpi_comm(ctxt)
    do i=1, n
      idxmap%glob_to_loc(i) = -1
    end do
    
    do i=1, nl 
      ix = vl(i) 
      idxmap%loc_to_glob(i)  = ix
      idxmap%glob_to_loc(ix) = i
    end do
    
    idxmap%local_rows   = nl
    idxmap%local_cols   = nl
    call idxmap%set_state(psb_desc_bld_)
   
  end subroutine list_initlvl


  subroutine list_asb(idxmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    use psb_sort_mod
    implicit none 
    class(psb_list_map), intent(inout) :: idxmap
    integer(psb_ipk_), intent(out) :: info
    
    integer(psb_ipk_) :: nhal
    type(psb_ctxt_type) :: ctxt
    integer(psb_mpk_) :: iam, np 
    
    info = 0 
    ctxt = idxmap%get_ctxt()
    call psb_info(ctxt,iam,np)

    nhal = idxmap%local_cols
    call psb_realloc(nhal,idxmap%loc_to_glob,info)

    call idxmap%set_state(psb_desc_asb_)
    
  end subroutine list_asb

  function list_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'LIST'
  end function list_get_fmt


  subroutine list_clone(idxmap,outmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_list_map), intent(inout)    :: idxmap
    class(psb_indx_map), allocatable, intent(out) :: outmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='list_clone'
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
    
    allocate(psb_list_map :: outmap, stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if

    select type (outmap)
    type is (psb_list_map) 
      call idxmap%psb_indx_map%cpy(outmap%psb_indx_map,info)
      if (info == psb_success_) then 
        outmap%pnt_h        = idxmap%pnt_h
      end if
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%loc_to_glob,outmap%loc_to_glob,info)
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%glob_to_loc,outmap%glob_to_loc,info)
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
  end subroutine list_clone


  subroutine list_reinit(idxmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_list_map), intent(inout)    :: idxmap
    integer(psb_ipk_), intent(out) :: info
    character(len=20)  :: name='list_reinit'
    logical, parameter :: debug=.false.

    info = psb_success_

    call idxmap%set_state(psb_desc_bld_)

    return

  end subroutine list_reinit


end module psb_list_map_mod
