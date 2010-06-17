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
!
! File: psi_bld_hash.f90
!
! Subroutine: psi_bld_hash
!   Build a hashed list of ordered sublists of the indices
!   contained in loc_to_glob. 
!   
! 
! Arguments: 
!    desc     - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
!
subroutine psi_bld_g2lmap(desc,info)
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psi_mod, psb_protect_name => psi_bld_g2lmap
  implicit none
  type(psb_desc_type), intent(inout) :: desc
  integer, intent(out) :: info

  integer          ::  i,j,np,me,lhalo,nhalo,nbits,hsize,hmask,&
       & n_col, err_act,  key, ih, nh, idx, nk,icomm
  integer             :: ictxt,n_row
  character(len=20)   :: name,ch_err

  info = psb_success_
  name = 'psi_bld_g2lmap'
  call psb_erractionsave(err_act)

  ictxt = psb_cd_get_context(desc)
  icomm = psb_cd_get_mpic(desc)
  n_row = psb_cd_get_local_rows(desc)
  n_col = psb_cd_get_local_cols(desc)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif


  if (.not.(psb_is_bld_desc(desc).and.psb_is_large_desc(desc))) then 
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if


  nk = n_col
  call psb_realloc(nk,2,desc%idxmap%glb_lc,info) 
  
  nbits = psb_hash_bits
  hsize = 2**nbits
  do 
    if (hsize < 0) then 
      ! This should never happen for sane values
      ! of psb_max_hash_bits.
      write(0,*) 'Error: hash size overflow ',hsize,nbits
      info = -2 
      return
    end if
    if (hsize > nk) exit
    if (nbits >= psb_max_hash_bits) exit
    nbits = nbits + 1
    hsize = hsize * 2 
  end do
  hmask = hsize - 1 
  desc%idxmap%hashvsize = hsize
  desc%idxmap%hashvmask = hmask
  if (info == psb_success_) call psb_realloc(hsize+1,desc%idxmap%hashv,info,lb=0)
  if (info /= psb_success_) then 
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! Build a hashed table of sorted lists to search for 
  ! indices.
  desc%idxmap%hashv(0:hsize) = 0
  do i=1, nk
    key = desc%idxmap%loc_to_glob(i)
    ih  = iand(key,hmask) 
    desc%idxmap%hashv(ih) = desc%idxmap%hashv(ih) + 1
  end do
  nh = desc%idxmap%hashv(0) 
  idx = 1
  do i=1, hsize
    desc%idxmap%hashv(i-1) = idx
    idx = idx + nh
    nh = desc%idxmap%hashv(i)
  end do
  do i=1, nk
    key                       = desc%idxmap%loc_to_glob(i)
    ih                        = iand(key,hmask)
    idx                       = desc%idxmap%hashv(ih) 
    desc%idxmap%glb_lc(idx,1) = key
    desc%idxmap%glb_lc(idx,2) = i
    desc%idxmap%hashv(ih)     = desc%idxmap%hashv(ih) + 1
  end do
  do i = hsize, 1, -1 
    desc%idxmap%hashv(i) = desc%idxmap%hashv(i-1)
  end do
  desc%idxmap%hashv(0) = 1
  do i=0, hsize-1 
    idx = desc%idxmap%hashv(i)
    nh  = desc%idxmap%hashv(i+1) - desc%idxmap%hashv(i) 
    if (nh > 1) then 
      call psb_msort(desc%idxmap%glb_lc(idx:idx+nh-1,1),&
           & ix=desc%idxmap%glb_lc(idx:idx+nh-1,2),&
           & flag=psb_sort_keep_idx_)
    end if
  end do

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


end subroutine psi_bld_g2lmap
