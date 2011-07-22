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
! File: psi_idx_ins_cnv.f90
!
! Subroutine: psi_idx_ins_cnv1
!   Converts a bunch of indices from global to local numbering. 
!   This routine is called while the descriptor is in the build state;
!   the idea is that if an index is not yet marked as local, it is a new 
!   connection to another process, i.e. a new entry into the halo. 
!   But we still need the mask, because we have to take out the column indices 
!   corresponding to row indices we do not own (see psb_cdins for how this is used). 
! 
! Arguments: 
!    nv        - integer                   Number of indices required 
!    idxin(:)  - integer                   Required indices, overwritten on output
!                                          output is negative for masked entries
!    desc      - type(psb_desc_type).    The communication descriptor.        
!    info      - integer.                  return code.
!    mask(:)   - logical, optional         Only do the conversion for specific indices.
!    
subroutine psi_idx_ins_cnv1(nv,idxin,desc,info,mask)
  use psi_mod, psb_protect_name => psi_idx_ins_cnv1
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none
  integer, intent(in)    :: nv
  integer, intent(inout) ::  idxin(:)
  type(psb_desc_type), intent(inout) :: desc
  integer, intent(out) :: info
  logical, intent(in), optional :: mask(:)
  integer :: ictxt,mglob, nglob
  integer                :: np, me
  integer                :: nrow,ncol, err_act
  integer                :: pnt_halo, nh, ip, lip,nxt,lipf,i,k,isize
  logical                :: pnt_h_ok
  integer, parameter     :: relocsz=200
  character(len=20)      :: name,ch_err

  info = psb_success_
  name = 'psb_idx_ins_cnv'
  call psb_erractionsave(err_act)

  ictxt = desc%get_context()
  mglob = desc%get_global_rows()
  nglob = desc%get_global_cols()
  nrow  = desc%get_local_rows()
  ncol  = desc%get_local_cols()

  call psb_info(ictxt, me, np)

  if ((.not.allocated(desc%indxmap)).or.&
       & (.not.psb_is_bld_desc(desc))) then 
    info =  psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (nv < 0) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (nv == 0) return


  if (size(idxin) < nv) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if


  if (present(mask)) then 
    if (size(mask) < nv) then 
      info = 1111
      call psb_errpush(info,name)
      goto 9999
    end if
  endif


  call desc%indxmap%g2l_ins(idxin(1:nv),info,mask)
  
  if (info /= 0) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='g2l_ins') 
    goto 9999      
  end if
    
!!$  desc%matrix_data(psb_n_col_) = desc%indxmap%get_lc()

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

end subroutine psi_idx_ins_cnv1
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
! Subroutine: psi_idx_ins_cnv2
!   Converts a bunch of indices from global to local numbering. 
!   This routine is called while the descriptor is in the build state;
!   the idea is that if an index is not yet marked as local, it is a new 
!   connection to another process, i.e. a new entry into the halo. 
!   But we still need the mask, because we have to take out the column indices 
!   corresponding to row indices we do not own (see psb_cdins for how this is used). 
! 
! Arguments: 
!    nv        - integer                   Number of indices required 
!    idxin(:)  - integer                   Required indices
!    idxout(:) - integer                   Output values (negative for masked entries)
!    desc      - type(psb_desc_type).    The communication descriptor.        
!    info      - integer.                  return code.
!    mask(:)   - logical, optional         Only do the conversion for specific indices.
!    
subroutine psi_idx_ins_cnv2(nv,idxin,idxout,desc,info,mask)
  use psi_mod, psb_protect_name => psi_idx_ins_cnv2
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psi_mod
  implicit none
  integer, intent(in)  :: nv, idxin(:)
  integer, intent(out) :: idxout(:)
  type(psb_desc_type), intent(inout) :: desc
  integer, intent(out) :: info
  logical, intent(in), optional :: mask(:)
  integer :: i,ictxt,k,mglob, nglob
  integer                :: np, me, isize
  integer                :: pnt_halo,nrow,ncol, nh, ip, err_act,lip,nxt,lipf
  logical                :: pnt_h_ok
  integer, parameter     :: relocsz=200
  character(len=20)      :: name,ch_err

  info = psb_success_
  name = 'psb_idx_ins_cnv'
  call psb_erractionsave(err_act)

  ictxt = desc%get_context()
  mglob = desc%get_global_rows()
  nglob = desc%get_global_cols()
  nrow  = desc%get_local_rows()
  ncol  = desc%get_local_cols()

  call psb_info(ictxt, me, np)

  if (.not.psb_is_ok_desc(desc)) then 
    info = psb_err_input_matrix_unassembled_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (nv < 0) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (nv == 0) return

  if (size(idxin) < nv) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if

  if (size(idxout) < nv) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if

  idxout(1:nv) = idxin(1:nv)
  call psi_idx_ins_cnv(nv,idxout,desc,info,mask)

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


end subroutine psi_idx_ins_cnv2
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
! Subroutine: psi_idx_ins_cnvs
!   Converts an index  from global to local numbering. 
!   This routine is called while the descriptor is in the build state;
!   the idea is that if an index is not yet marked as local, it is a new 
!   connection to another process, i.e. a new entry into the halo. 
!   But we still need the mask, because we have to take out the column indices 
!   corresponding to row indices we do not own (see psb_cdins for how this is used). 
! 
! Arguments: 
!    idxin     - integer                   Required index s
!    idxout    - integer                   Output value  (negative for masked entries)
!    desc      - type(psb_desc_type).    The communication descriptor.        
!    info      - integer.                  return code.
!    mask      - logical, optional         Only do the conversion for specific indices.
!    
subroutine psi_idx_ins_cnvs2(idxin,idxout,desc,info,mask)
  use psi_mod, psb_protect_name => psi_idx_ins_cnvs2
  use psb_descriptor_type
  integer, intent(in)  :: idxin
  integer, intent(out) :: idxout
  type(psb_desc_type), intent(inout) :: desc
  integer, intent(out) :: info
  logical, intent(in), optional :: mask
  integer  :: iout(1) 
  logical  :: mask_(1)
  
  if (present(mask)) then 
    mask_ = mask
  else
    mask_ = .true.
  end if

  iout(1) = idxin
  call psi_idx_ins_cnv(1,iout,desc,info,mask_)
  idxout  = iout(1) 
  return

end subroutine psi_idx_ins_cnvs2
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
! Subroutine: psi_idx_ins_cnvs
!   Converts an index  from global to local numbering. 
!   This routine is called while the descriptor is in the build state;
!   the idea is that if an index is not yet marked as local, it is a new 
!   connection to another process, i.e. a new entry into the halo. 
!   But we still need the mask, because we have to take out the column indices 
!   corresponding to row indices we do not own (see psb_cdins for how this is used). 
! 
! Arguments: 
!    idxin     - integer                   Required index s
!    idxout    - integer                   Output value  (negative for masked entries)
!    desc      - type(psb_desc_type).    The communication descriptor.        
!    info      - integer.                  return code.
!    mask      - logical, optional         Only do the conversion for specific indices.
!    
subroutine psi_idx_ins_cnvs1(idxin,desc,info,mask)
  use psi_mod, psb_protect_name => psi_idx_ins_cnvs1
  use psb_descriptor_type
  integer, intent(inout)  :: idxin
  type(psb_desc_type), intent(inout) :: desc
  integer, intent(out) :: info
  logical, intent(in), optional :: mask
  integer  :: iout(1) 
  logical  :: mask_(1)
  
  if (present(mask)) then 
    mask_ = mask
  else
    mask_ = .true.
  end if

  iout(1) = idxin
  call psi_idx_ins_cnv(1,iout,desc,info,mask_)
  idxin   = iout(1) 

  return

end subroutine psi_idx_ins_cnvs1
