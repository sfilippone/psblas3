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
! File: psi_bld_tmpovrl.f90
!
! Subroutine: psi_bld_tmpovrl
!   Build initial versions of overlap  exchange lists.
!   When the descriptor is for a large index space, we cannot build 
!   the data exchange lists "on-the-fly", but we also want to keep using the 
!   same format conversion routines we use in the small index space case, 
!   hence this adapter routine.
!   
! 
! Arguments:
!    iv(:)    - integer               Initial list.
!                                     index
!                                     nprocs (sharing it)
!                                     procs(1:nprocs)
!                                     End marked with -1
!                                        
!    desc     - type(psb_desc_type).  The communication descriptor.        
!    info     - integer.              return code.
!
subroutine psi_bld_tmpovrl(iv,desc,info)
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psi_mod, psb_protect_name => psi_bld_tmpovrl
  implicit none
  integer, intent(in)  :: iv(:)
  type(psb_desc_type), intent(inout) :: desc
  integer, intent(out)  :: info

  !locals
  Integer              :: counter,i,j,np,me,loc_row,err,loc_col,nprocs,&
       & l_ov_ix,l_ov_el,idx, err_act, itmpov, k, glx, icomm
  integer, allocatable  :: ov_idx(:),ov_el(:,:)

  integer             :: ictxt,n_row, debug_unit, debug_level
  character(len=20)   :: name,ch_err

  info = psb_success_
  name = 'psi_bld_tmpovrl'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt = psb_cd_get_context(desc)
  icomm = psb_cd_get_mpic(desc)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  l_ov_ix=0
  l_ov_el=0
  i = 1
  do while (iv(i) /= -1) 
    idx = iv(i)
    i       = i + 1
    nprocs  = iv(i)
    i       = i + 1
    l_ov_ix = l_ov_ix+3*(nprocs-1)
    l_ov_el = l_ov_el + 1
    i       = i + nprocs     
  enddo

  l_ov_ix = l_ov_ix+3  

  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': Ov len',l_ov_ix,l_ov_el

  allocate(ov_idx(l_ov_ix),ov_el(l_ov_el,3), stat=info)
  if (info /= psb_no_err_) then
    info=psb_err_from_subroutine_
    err=info
    call psb_errpush(err,name,a_err='psb_realloc')
    goto 9999
  end if

  l_ov_ix=0
  l_ov_el=0
  i = 1
  do while (iv(i) /= -1) 
    idx = iv(i)
    i   = i+1
    nprocs = iv(i)
    l_ov_el          = l_ov_el+1
    ov_el(l_ov_el,1) = idx                    ! Index
    ov_el(l_ov_el,2) = nprocs                 ! How many procs
    ov_el(l_ov_el,3) = minval(iv(i+1:i+nprocs))  ! master proc
    do j=1, nprocs
      if (iv(i+j) /= me) then
        ov_idx(l_ov_ix+1) = iv(i+j)
        ov_idx(l_ov_ix+2) = 1
        ov_idx(l_ov_ix+3) = idx
        l_ov_ix = l_ov_ix+3
      endif
    enddo
    i = i + nprocs + 1
  enddo
  l_ov_ix         = l_ov_ix + 1
  ov_idx(l_ov_ix) = -1
  call psb_move_alloc(ov_idx,desc%ovrlap_index,info) 
  if (info == psb_success_) call psb_move_alloc(ov_el,desc%ovrlap_elem,info)


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


end subroutine psi_bld_tmpovrl
