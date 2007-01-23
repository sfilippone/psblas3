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
subroutine psi_fnd_owner(nv,idx,iprc,desc,info)
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psi_mod, only : psi_idx_cnv
  implicit none 
  integer, intent(in) :: nv
  integer, intent(in) ::  idx(:)
  integer, allocatable, intent(out) ::  iprc(:)
  type(psb_desc_type), intent(in) :: desc
  integer, intent(out) :: info


  integer,allocatable :: hsz(:),hidx(:),helem(:),hproc(:)

  integer          ::  i,j,err,n_row,n_col, err_act,ih,nh,icomm,hsize
  integer             :: ictxt,np,me
  logical, parameter  :: debug=.false., debugwrt=.false.
  character(len=20)   :: name,ch_err

  info = 0
  name = 'psi_fnd_owner'
  call psb_erractionsave(err_act)

  ictxt   = psb_cd_get_context(desc)
  n_row   = psb_cd_get_local_rows(desc)
  n_col   = psb_cd_get_local_cols(desc)
  call psb_get_mpicomm(ictxt,icomm )


  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (nv < 0 ) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.(psb_is_ok_desc(desc))) then 
    write(0,*) 'Invalid input descriptor in psi_fnd_owner'
  end if


  Allocate(hidx(np+1),hsz(np),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if
  hsz       = 0
  hsz(me+1) = nv
  call psb_amx(ictxt,hsz,info)
  hidx(1)   = 1
  do i=1, np
    hidx(i+1) = hidx(i) + hsz(i)
  end do
  hsize = hidx(np+1)-1
  Allocate(helem(hsize),hproc(hsize),stat=info)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if
  helem(:) = 0 
  ih = hidx(me+1)
  do i=1, hsz(me+1) 
    helem(ih+i-1) = idx(i) 
  end do
  call psb_amx(ictxt,helem,info)
  call psi_idx_cnv(hsize,helem,desc,info,owned=.true.)
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='psi_idx_cnv')
    goto 9999      
  end if
  
  do i=1,hsize
    if ((0< helem(i)).and. (helem(i) <= n_row))then 
      hproc(i) = me+1
    else
      hproc(i) = 0
    end if
  end do

  call psb_amx(ictxt,hproc,info)

  if (nv > 0) then 
    call psb_realloc(nv,iprc,info)
    ih = hidx(me+1)
    do i=1, hsz(me+1) 
      iprc(i) = hproc(ih+i-1) - 1 
    end do 
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

end subroutine psi_fnd_owner
