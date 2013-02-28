!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
! File: psi_idx_cnv.f90
!
! Subroutine: psi_idx_cnv1
!   Converts a bunch of indices from global to local numbering. 
!   
! 
! Arguments: 
!    nv       - integer                   Number of indices required 
!    idxin(:) - integer                   Required indices,   overwritten on output.
!    desc     - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
!    mask(:)  - logical, optional         Only do the conversion for specific indices.
!    owned    - logical,optional          Restrict to local indices, no halo 
!                                         (default false)
subroutine psi_idx_cnv1(nv,idxin,desc,info,mask,owned)
  use psb_desc_mod
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psi_mod, psb_protect_name => psi_idx_cnv1
  implicit none
  integer(psb_ipk_), intent(in)    :: nv
  integer(psb_ipk_), intent(inout) ::  idxin(:)
  type(psb_desc_type), intent(in) :: desc
  integer(psb_ipk_), intent(out) :: info
  logical, intent(in), optional :: mask(:)
  logical, intent(in), optional :: owned
  integer(psb_ipk_) :: ictxt,mglob, nglob,ip,lip,i
  integer(psb_ipk_) :: np, me
  integer(psb_ipk_) :: nrow,ncol, err_act
  integer(psb_ipk_), allocatable   :: itmp(:)
  integer(psb_ipk_), parameter     :: relocsz=200
  character(len=20)      :: name
  logical                :: owned_

  info = psb_success_
  name = 'psb_idx_cnv'
  call psb_erractionsave(err_act)

  ictxt   = desc%get_context()
  mglob   = desc%get_global_rows()
  nglob   = desc%get_global_cols()
  nrow    = desc%get_local_rows()
  ncol    = desc%get_local_cols()

  call psb_info(ictxt, me, np)

  if (.not.allocated(desc%indxmap))then 
    info =  psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.psb_is_valid_desc(desc)) then 
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

  call desc%indxmap%g2lip(idxin(1:nv),info,mask=mask,owned=owned)

  if (info /= 0) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='g2l') 
    goto 9999      
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

end subroutine psi_idx_cnv1
!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
! Subroutine: psi_idx_cnv2
!   Converts a bunch of indices from global to local numbering. 
!   
! 
! Arguments: 
!    nv        - integer                   Number of indices required 
!    idxin(:)  - integer                   Required indices
!    idxout(:) - integer                   Output values, negative for invalid input.
!    desc      - type(psb_desc_type).    The communication descriptor.        
!    info      - integer.                  return code.
!    mask(:)   - logical, optional         Only do the conversion for specific indices.
!    owned     - logical,optional          Restrict to local indices, no halo
!                                          (default false)
subroutine psi_idx_cnv2(nv,idxin,idxout,desc,info,mask,owned)
  use psb_desc_mod
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psi_mod, psb_protect_name => psi_idx_cnv2
  implicit none
  integer(psb_ipk_), intent(in)  :: nv, idxin(:)
  integer(psb_ipk_), intent(out) :: idxout(:)
  type(psb_desc_type), intent(in) :: desc
  integer(psb_ipk_), intent(out) :: info
  logical, intent(in), optional :: mask(:)
  logical, intent(in), optional :: owned
  integer(psb_ipk_) :: i,ictxt,mglob, nglob
  integer(psb_ipk_) :: np, me
  integer(psb_ipk_) :: nrow,ncol, err_act
  integer(psb_ipk_), parameter     :: relocsz=200
  character(len=20)      :: name
  logical, pointer       :: mask_(:)
  logical                :: owned_

  info = psb_success_
  name = 'psb_idx_cnv'
  call psb_erractionsave(err_act)

  ictxt   = desc%get_context()
  mglob   = desc%get_global_rows()
  nglob   = desc%get_global_cols()
  nrow    = desc%get_local_rows()
  ncol    = desc%get_local_cols()


  call psb_info(ictxt, me, np)

  if (.not.psb_is_ok_desc(desc)) then 
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

  if (size(idxout) < nv) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if

  idxout(1:nv) = idxin(1:nv) 
  call psi_idx_cnv1(nv,idxout,desc,info,mask=mask,owned=owned)  

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


end subroutine psi_idx_cnv2
!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
! Subroutine: psi_idx_cnvs
!   Converts an index from global to local numbering. 
!   
! 
! Arguments: 
!    idxin     - integer                   Required index   
!    idxout    - integer                   Output value, negative for invalid input.
!    desc      - type(psb_desc_type).    The communication descriptor.        
!    info      - integer.                  return code.
!    mask      - logical, optional         Only do the conversion if true.
!    owned     - logical,optional          Restrict to local indices, no halo 
!                                          (default false)
subroutine psi_idx_cnvs(idxin,idxout,desc,info,mask,owned)

  use psi_mod, psb_protect_name => psi_idx_cnvs
  use psb_desc_mod
  integer(psb_ipk_), intent(in)  :: idxin
  integer(psb_ipk_), intent(out) :: idxout
  type(psb_desc_type), intent(in) :: desc
  integer(psb_ipk_), intent(out) :: info
  logical, intent(in), optional :: mask
  logical, intent(in), optional :: owned
  logical  :: mask_(1)
  integer(psb_ipk_) :: iout(1) 
  
  if (present(mask)) then
    mask_ = mask
  else
    mask_=.true.
  end if
  iout = idxin
  call psi_idx_cnv(ione,iout,desc,info,mask=mask_,owned=owned)
  idxout=iout(1)

  return

end subroutine psi_idx_cnvs
subroutine psi_idx_cnvs1(idxin,desc,info,mask,owned)

  use psi_mod, psb_protect_name => psi_idx_cnvs1
  use psb_desc_mod
  integer(psb_ipk_), intent(inout)  :: idxin
  type(psb_desc_type), intent(in) :: desc
  integer(psb_ipk_), intent(out) :: info
  logical, intent(in), optional :: mask
  logical, intent(in), optional :: owned
  logical  :: mask_(1)
  integer(psb_ipk_) :: iout(1) 
  
  if (present(mask)) then
    mask_ = mask
  else
    mask_=.true.
  end if

  iout(1) = idxin
  call psi_idx_cnv(ione,iout,desc,info,mask=mask_,owned=owned)
  idxin   = iout(1)

  return
  
end subroutine psi_idx_cnvs1
