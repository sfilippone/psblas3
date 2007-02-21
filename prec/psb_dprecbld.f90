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
subroutine psb_dprecbld(a,desc_a,p,info,upd)

  use psb_base_mod
  use psb_prec_mod, psb_protect_name => psb_dprecbld
  Implicit None

  type(psb_dspmat_type), target           :: a
  type(psb_desc_type), intent(in), target :: desc_a
  type(psb_dprec_type),intent(inout)      :: p
  integer, intent(out)                    :: info
  character, intent(in), optional         :: upd

  ! Local scalars
  Integer      :: err, nnzero, n_row, n_col,I,j,k,ictxt,&
       & me,mycol,np,npcol,mglob,lw, mtype, nrg, nzg, err_act
  real(kind(1.d0))         :: temp, real_err(5)
  real(kind(1.d0)),pointer :: gd(:), work(:)
  integer      :: int_err(5)
  character    :: iupd

  logical, parameter :: debug=.false.
  integer,parameter  :: iroot=0,iout=60,ilout=40
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  name = 'psb_precbld'

  if (debug) write(0,*) 'Entering precbld',P%prec,desc_a%matrix_data(:)
  info = 0
  int_err(1) = 0
  ictxt = psb_cd_get_context(desc_a)

  if (debug) write(0,*) 'Preconditioner psb_info'
  call psb_info(ictxt, me, np)

  if (present(upd)) then 
    if (debug) write(0,*) 'UPD ', upd
    if ((upd.eq.'F').or.(upd.eq.'T')) then
      iupd=upd
    else
      iupd='F'
    endif
  else
    iupd='F'
  endif
  n_row   = psb_cd_get_local_rows(desc_a)
  n_col   = psb_cd_get_local_cols(desc_a)
  mglob   = psb_cd_get_global_rows(desc_a)
  !
  ! Should add check to ensure all procs have the same... 
  !
  ! ALso should define symbolic names for the preconditioners. 
  !

  call psb_check_def(p%iprcparm(p_type_),'base_prec',&
       &  diag_,is_legal_prec)

  call psb_nullify_desc(p%desc_data)

  select case(p%iprcparm(p_type_)) 
  case (noprec_)
    ! Do nothing. 
    call psb_cdcpy(desc_a,p%desc_data,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_cdcpy'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case (diag_)

    call psb_diagsc_bld(a,desc_a,p,iupd,info)
    if(debug) write(0,*)me,': out of psb_diagsc_bld'
    if(info /= 0) then
      info=4010
      ch_err='psb_diagsc_bld'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case (bjac_)

    call psb_check_def(p%iprcparm(f_type_),'fact',&
         &  f_ilu_n_,is_legal_ml_fact)

    if (debug) write(0,*)me, ': Calling PSB_BJAC_BLD'
    if (debug) call psb_barrier(ictxt)
    call psb_bjac_bld(a,desc_a,p,iupd,info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_bjac_bld')
      goto 9999
    end if

  case default
    info=4010
    ch_err='Unknown p_type_'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999

  end select

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_dprecbld

