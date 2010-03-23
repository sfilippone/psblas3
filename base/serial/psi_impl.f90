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

  subroutine psi_renum_index(iperm,idx,info)
    use psi_mod, psi_protect_name =>  psi_renum_index
    use psb_serial_mod 
    implicit none 

    integer, intent(out)   :: info
    integer, intent(in)    :: iperm(:)
    integer, intent(inout) :: idx(:)

    integer :: i,j,k,nh
    
    i=1
    k=idx(i)
    do while (k /= -1) 
      i = i+1
      nh = idx(i)
      do j = i+1, i+nh
        idx(j) = iperm(idx(j))
      enddo
      i  = i + nh + 1
      nh = idx(i)
      do j = i+1, i+nh
        idx(j) = iperm(idx(j))
      enddo
      i = i + nh + 1
      k = idx(i)
    enddo

  end subroutine psi_renum_index

  subroutine psi_renum_idxmap(nc,iperm,idxmap,info)
    use psi_mod, psi_protect_name =>  psi_renum_idxmap
    use psb_serial_mod 
    implicit none 

    integer, intent(out)   :: info
    integer, intent(in)    :: nc,iperm(:)
    type(psb_idxmap_type), intent(inout) :: idxmap
    
    integer, allocatable :: itmp(:)
    integer :: i,j,k,nh
    
    if (nc > size(iperm)) then 
      info = 2 
      return
    endif

    if (idxmap%state == psb_desc_large_) then 

      allocate(itmp(size(idxmap%loc_to_glob)), stat=i)
      if (i/=0) then
        info = 4001
        return
      end if
      do i=1,nc
        itmp(i) = idxmap%loc_to_glob(iperm(i))
      end do
      do i=1, size(idxmap%glb_lc,1)
        idxmap%glb_lc(i,2) = iperm(idxmap%glb_lc(i,2))
      end do
      do i=1, nc 
        idxmap%loc_to_glob(i) = itmp(i)
      end do
        
    else

      do i=1, nc
        idxmap%glob_to_loc(idxmap%loc_to_glob(iperm(i))) = i  
      enddo
      do i=1,size(idxmap%glob_to_loc)
        j = idxmap%glob_to_loc(i)
        if (j>0) then 
          idxmap%loc_to_glob(j) = i
        endif
      enddo
    end if
      
  end subroutine psi_renum_idxmap
  
  subroutine psi_cnv_dsc(halo_in,ovrlap_in,ext_in,cdesc, info)

    use psi_mod, psi_protect_name =>  psi_cnv_dsc
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
    call psb_move_alloc(idx_out,cdesc%halo_index,info)
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
    call psb_move_alloc(idx_out,cdesc%ext_index,info)
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
    call psb_move_alloc(idx_out,cdesc%ovrlap_index,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psb_move_alloc')
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
    call psb_move_alloc(idx_out,cdesc%ovr_mst_idx,info)
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psb_move_alloc')
      goto 9999
    end if

    cdesc%matrix_data(psb_tmov_xch_) = nxch
    cdesc%matrix_data(psb_tmov_snd_) = nsnd
    cdesc%matrix_data(psb_tmov_rcv_) = nrcv 

    ! finally bnd_elem
    call psi_crea_bnd_elem(idx_out,cdesc,info)
    if (info == 0) call psb_move_alloc(idx_out,cdesc%bnd_elem,info)

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


  subroutine psi_inner_cnvs(x,hashmask,hashv,glb_lc)
    use psi_mod, psi_protect_name => psi_inner_cnvs

    integer, intent(in)    :: hashmask,hashv(0:),glb_lc(:,:)
    integer, intent(inout) :: x

    integer :: i, ih, key, idx,nh,tmp,lb,ub,lm
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
      x = glb_lc(tmp,2)
    else         
      x = tmp 
    end if
  end subroutine psi_inner_cnvs
  
  subroutine psi_inner_cnvs2(x,y,hashmask,hashv,glb_lc)
    use psi_mod, psi_protect_name =>  psi_inner_cnvs2
    integer, intent(in)  :: hashmask,hashv(0:),glb_lc(:,:)
    integer, intent(in)  :: x
    integer, intent(out) :: y

    integer :: i, ih, key, idx,nh,tmp,lb,ub,lm
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
      y = glb_lc(tmp,2)
    else         
      y = tmp 
    end if
  end subroutine psi_inner_cnvs2


  subroutine psi_inner_cnv1(n,x,hashmask,hashv,glb_lc,mask)
    use psi_mod, psi_protect_name =>  psi_inner_cnv1
    integer, intent(in)    :: n,hashmask,hashv(0:),glb_lc(:,:)
    logical, intent(in), optional    :: mask(:)
    integer, intent(inout) :: x(:)

    integer :: i, ih, key, idx,nh,tmp,lb,ub,lm
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
    end if
  end subroutine psi_inner_cnv1

  subroutine psi_inner_cnv2(n,x,y,hashmask,hashv,glb_lc,mask)
    use psi_mod, psi_protect_name =>  psi_inner_cnv2
    integer, intent(in)  :: n, hashmask,hashv(0:),glb_lc(:,:)
    logical, intent(in),optional  :: mask(:)
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
    if (present(mask)) then 
      do i=1, n
        if (mask(i)) then 
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
        end if
      end do
    else
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
    end if
  end subroutine psi_inner_cnv2

  subroutine  psi_sovrl_updr1(x,desc_a,update,info)
    use psi_mod, psi_protect_name =>   psi_sovrl_updr1

    implicit none

    real(psb_spk_), intent(inout), target :: x(:)
    type(psb_desc_type), intent(in)         :: desc_a
    integer, intent(in)                     :: update
    integer, intent(out)                    :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_sovrl_updr1'
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
             & x(idx) = szero
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
  end subroutine psi_sovrl_updr1


  subroutine  psi_sovrl_updr2(x,desc_a,update,info)
    use psi_mod, psi_protect_name =>   psi_sovrl_updr2

    implicit none

    real(psb_spk_), intent(inout), target :: x(:,:)
    type(psb_desc_type), intent(in)         :: desc_a
    integer, intent(in)                     :: update
    integer, intent(out)                    :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_sovrl_updr2'
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
             & x(idx,:) = szero
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
  end subroutine psi_sovrl_updr2

  subroutine  psi_dovrl_updr1(x,desc_a,update,info)
    use psi_mod, psi_protect_name =>   psi_dovrl_updr1

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
    use psi_mod, psi_protect_name =>   psi_dovrl_updr2

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

  subroutine  psi_covrl_updr1(x,desc_a,update,info)
    use psi_mod, psi_protect_name =>   psi_covrl_updr1

    implicit none

    complex(psb_spk_), intent(inout), target :: x(:)
    type(psb_desc_type), intent(in)         :: desc_a
    integer, intent(in)                     :: update
    integer, intent(out)                    :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_covrl_updr1'
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
             & x(idx) = czero
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
  end subroutine psi_covrl_updr1


  subroutine  psi_covrl_updr2(x,desc_a,update,info)
    use psi_mod, psi_protect_name =>   psi_covrl_updr2

    implicit none

    complex(psb_spk_), intent(inout), target :: x(:,:)
    type(psb_desc_type), intent(in)         :: desc_a
    integer, intent(in)                     :: update
    integer, intent(out)                    :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, ndm
    character(len=20) :: name, ch_err

    name='psi_covrl_updr2'
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
             & x(idx,:) = czero
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
  end subroutine psi_covrl_updr2

  subroutine  psi_zovrl_updr1(x,desc_a,update,info)
    use psi_mod, psi_protect_name =>   psi_zovrl_updr1

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
    use psi_mod, psi_protect_name =>   psi_zovrl_updr2

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
    use psi_mod, psi_protect_name =>   psi_iovrl_updr1

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
    use psi_mod, psi_protect_name =>   psi_iovrl_updr2

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


  subroutine  psi_sovrl_saver1(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_sovrl_saver1
    use psb_realloc_mod

    implicit none

    real(psb_spk_), intent(inout)  :: x(:)
    real(psb_spk_), allocatable    :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_sovrl_saver1'
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
  end subroutine psi_sovrl_saver1

  subroutine  psi_sovrl_restrr1(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_sovrl_restrr1

    implicit none

    real(psb_spk_), intent(inout)  :: x(:)
    real(psb_spk_)                 :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_sovrl_restrr1'
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
  end subroutine psi_sovrl_restrr1


  subroutine  psi_sovrl_saver2(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_sovrl_saver2
    use psb_realloc_mod

    implicit none

    real(psb_spk_), intent(inout)  :: x(:,:)
    real(psb_spk_), allocatable    :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz, nc
    character(len=20) :: name, ch_err

    name='psi_sovrl_saver2'
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
  end subroutine psi_sovrl_saver2

  subroutine  psi_sovrl_restrr2(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_sovrl_restrr2

    implicit none

    real(psb_spk_), intent(inout)  :: x(:,:)
    real(psb_spk_)                 :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_sovrl_restrr2'
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
  end subroutine psi_sovrl_restrr2


  subroutine  psi_dovrl_saver1(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_dovrl_saver1
    use psb_realloc_mod

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
    use psi_mod, psi_protect_name =>   psi_dovrl_restrr1

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
    use psi_mod, psi_protect_name =>   psi_dovrl_saver2
    use psb_realloc_mod

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
    use psi_mod, psi_protect_name =>   psi_dovrl_restrr2

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

  subroutine  psi_covrl_saver1(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_covrl_saver1
    use psb_realloc_mod

    implicit none

    complex(psb_spk_), intent(inout)  :: x(:)
    complex(psb_spk_), allocatable    :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_covrl_saver1'
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
  end subroutine psi_covrl_saver1

  subroutine  psi_covrl_restrr1(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_covrl_restrr1

    implicit none

    complex(psb_spk_), intent(inout)  :: x(:)
    complex(psb_spk_)                 :: xs(:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_covrl_restrr1'
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
  end subroutine psi_covrl_restrr1


  subroutine  psi_covrl_saver2(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_covrl_saver2
    use psb_realloc_mod

    implicit none

    complex(psb_spk_), intent(inout)  :: x(:,:)
    complex(psb_spk_), allocatable    :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz, nc
    character(len=20) :: name, ch_err

    name='psi_covrl_saver2'
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
  end subroutine psi_covrl_saver2

  subroutine  psi_covrl_restrr2(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_covrl_restrr2

    implicit none

    complex(psb_spk_), intent(inout)  :: x(:,:)
    complex(psb_spk_)                 :: xs(:,:)
    type(psb_desc_type), intent(in)  :: desc_a
    integer, intent(out)             :: info

    ! locals
    integer           :: ictxt, np, me, err_act, i, idx, isz
    character(len=20) :: name, ch_err

    name='psi_covrl_restrr2'
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
  end subroutine psi_covrl_restrr2


  subroutine  psi_zovrl_saver1(x,xs,desc_a,info)
    use psi_mod, psi_protect_name =>   psi_zovrl_saver1

    use psb_realloc_mod

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
    use psi_mod, psi_protect_name =>   psi_zovrl_restrr1

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
    use psi_mod, psi_protect_name =>   psi_zovrl_saver2

    use psb_realloc_mod

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
    use psi_mod, psi_protect_name =>   psi_zovrl_restrr2

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
    use psi_mod, psi_protect_name =>   psi_iovrl_saver1

    use psb_realloc_mod

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
    use psi_mod, psi_protect_name =>   psi_iovrl_restrr1

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
    use psi_mod, psi_protect_name =>   psi_iovrl_saver2
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
    use psi_mod, psi_protect_name =>   psi_iovrl_restrr2

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

  subroutine psi_bld_ovr_mst(me,ovrlap_elem,mst_idx,info)
    use psi_mod, psi_protect_name =>  psi_bld_ovr_mst

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

