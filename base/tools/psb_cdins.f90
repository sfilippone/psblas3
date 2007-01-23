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
! File: psb_cdins.f90
!
! Subroutine: psb_cdins
!   Takes as input a cloud of points and updates the descriptor accordingly.
! 
! Parameters: 
!    nz       - integer.                       The number of points to insert.
!    ia       - integer,dimension(:).          The row indices of the points.
!    ja       - integer,dimension(:).          The column indices of the points.
!    desc_a   - type(<psb_desc_type>).         The communication descriptor to be freed.
!    info     - integer.                       Eventually returns an error code.
subroutine psb_cdins(nz,ia,ja,desc_a,info,ila,jla)

  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psi_mod
  implicit none

  !....PARAMETERS...
  Type(psb_desc_type), intent(inout) :: desc_a
  Integer, intent(in)                :: nz,ia(:),ja(:)
  integer, intent(out)               :: info
  integer, optional, intent(out)     :: ila(:), jla(:)

  !LOCALS.....

  integer :: i,ictxt,row,k,dectype,mglob, nglob,err
  integer                :: np, me, isize
  integer                :: pnt_halo,nrow,ncol, nh, ip,jp, err_act,lip,ljp,nxt
  logical, parameter     :: debug=.false.
  integer, parameter     :: relocsz=200
  integer, allocatable   :: ila_(:), jla_(:)
  character(len=20)      :: name,ch_err

  info = 0
  name = 'psb_cdins'
  call psb_erractionsave(err_act)

  ictxt   = psb_cd_get_context(desc_a)
  dectype = psb_cd_get_dectype(desc_a)
  mglob   = psb_cd_get_global_rows(desc_a)
  nglob   = psb_cd_get_global_cols(desc_a)
  nrow    = psb_cd_get_local_rows(desc_a)
  ncol    = psb_cd_get_local_cols(desc_a)

  call psb_info(ictxt, me, np)

  if (.not.psb_is_bld_desc(desc_a)) then 
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif

  if (nz <= 0) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (size(ia) < nz) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if

  if (size(ja) < nz) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (present(ila)) then 
    if (size(ila) < nz) then 
      info = 1111
      call psb_errpush(info,name)
      goto 9999
    end if
  end if
  if (present(jla)) then 
    if (size(jla) < nz) then 
      info = 1111
      call psb_errpush(info,name)
      goto 9999
    end if
  end if

  if (present(ila).and.present(jla)) then 
    call psi_idx_cnv(nz,ia,ila,desc_a,info,owned=.true.)
    call psi_idx_ins_cnv(nz,ja,jla,desc_a,info,mask=(ila(1:nz)>0))

  else
    if (present(ila).or.present(jla)) then 
      write(0,*) 'Inconsistent call : ',present(ila),present(jla)
    endif
    allocate(ila_(nz),jla_(nz),stat=info)
    if (info /= 0) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    end if
    call psi_idx_cnv(nz,ia,ila_,desc_a,info,owned=.true.)
    call psi_idx_ins_cnv(nz,ja,jla_,desc_a,info,mask=(ila_(1:nz)>0))
    deallocate(ila_,jla_)
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_cdins

