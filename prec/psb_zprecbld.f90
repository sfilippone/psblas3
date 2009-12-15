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
subroutine psb_zprecbld(a,desc_a,p,info,upd)

  use psb_sparse_mod
  use psb_prec_mod, psb_protect_name => psb_zprecbld
  Implicit None

  type(psb_z_sparse_mat), intent(in), target  :: a
  type(psb_desc_type), intent(in), target    :: desc_a
  type(psb_zprec_type),intent(inout)         :: p
  integer, intent(out)                       :: info
  character, intent(in), optional            :: upd


  ! Local scalars
  Integer      :: err, n_row, n_col,ictxt,&
       & me,np,mglob, err_act
  integer      :: int_err(5)
  character    :: upd_

  integer,parameter  :: iroot=psb_root_,iout=60,ilout=40
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  name = 'psb_precbld'

  info = 0
  int_err(1) = 0
  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)

  if (present(upd)) then 
    upd_ = psb_toupper(upd)
  else
    upd_='F'
  endif
  if ((upd_ == 'F').or.(upd_ == 'T')) then
    ! ok
  else
    upd_='F'
  endif
  n_row   = psb_cd_get_local_rows(desc_a)
  n_col   = psb_cd_get_local_cols(desc_a)
  mglob   = psb_cd_get_global_rows(desc_a)
  !
  ! Should add check to ensure all procs have the same... 
  !
  ! ALso should define symbolic names for the preconditioners. 
  !

  if (.not.allocated(p%prec)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if

  call p%prec%precbld(a,desc_a,info,upd)
  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_zprecbld

