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
subroutine psb_sdiagsc_bld(a,desc_a,p,upd,info)

  use psb_base_mod
  use psb_prec_mod, psb_protect_name => psb_sdiagsc_bld
  Implicit None

  type(psb_sspmat_type), intent(in), target :: a
  type(psb_desc_type), intent(in)           :: desc_a
  type(psb_sprec_type),intent(inout)        :: p
  character, intent(in)                     :: upd
  integer, intent(out)                      :: info


  ! Local scalars
  Integer      :: err, n_row, n_col,I,ictxt,&
       & me,np,mglob, err_act
  integer      :: int_err(5)

  integer,parameter  :: iroot=psb_root_,iout=60,ilout=40
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  name = 'psb_diagsc_bld'

  info = 0
  int_err(1) = 0
  ictxt = psb_cd_get_context(desc_a)
  n_row = psb_cd_get_local_rows(desc_a)
  n_col = psb_cd_get_local_cols(desc_a)
  mglob = psb_cd_get_global_rows(desc_a)

  call psb_info(ictxt, me, np)

  ! diagonal scaling

  call psb_realloc(n_col,p%d,info)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_realloc')
    goto 9999
  end if

  !
  ! Retrieve the diagonal entries of the matrix A
  !
  call psb_sp_getdiag(a,p%d,info)
  if(info /= 0) then
    info=4010
    ch_err='psb_sp_getdiag'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  !
  ! Copy into p%desc_data the descriptor associated to A
  !
  call psb_cdcpy(desc_a,p%desc_Data,info)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_cdcpy')
    goto 9999
  end if

  !
  ! The i-th diagonal entry of the preconditioner is set to one if the
  ! corresponding entry a_ii of the sparse matrix A is zero; otherwise 
  ! it is set to one/a_ii
  !
  do i=1,n_row
    if (p%d(i) == dzero) then
      p%d(i) = done
    else
      p%d(i) = done/p%d(i)
    endif
  end do

  if (a%pl(1) /= 0) then
    !
    ! Apply the same row permutation as in the sparse matrix A
    !
    call  psb_gelp('n',a%pl,p%d,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_gelp'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_sdiagsc_bld

