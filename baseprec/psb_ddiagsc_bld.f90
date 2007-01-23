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
subroutine psb_ddiagsc_bld(a,desc_a,p,upd,info)

  use psb_base_mod
  use psb_prec_type
  Implicit None

  type(psb_dspmat_type), target           :: a
  type(psb_desc_type), intent(in)         :: desc_a
  type(psb_dbaseprc_type),intent(inout)   :: p
  character, intent(in)                   :: upd
  integer, intent(out)                    :: info


  ! Local scalars
  Integer      :: err, n_row, n_col,I,j,k,ictxt,&
       & me,np,mglob,lw, err_act
  real(kind(1.d0)),allocatable  :: gd(:), work(:)
  integer      :: int_err(5)
  character    :: iupd

  logical, parameter :: debug=.false.   
  integer,parameter  :: iroot=0,iout=60,ilout=40
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  err=0
  call psb_erractionsave(err_act)
  name = 'psb_diagsc_bld'

  if (debug) write(0,*) 'Entering diagsc_bld'
  info = 0
  int_err(1) = 0
  ictxt = psb_cd_get_context(desc_a)
  n_row = psb_cd_get_local_rows(desc_a)
  n_col = psb_cd_get_local_cols(desc_a)
  mglob = psb_cd_get_global_rows(desc_a)

  if (debug) write(0,*) 'Preconditioner Blacs_gridinfo'
  call psb_info(ictxt, me, np)

  if (debug) write(0,*) 'Precond: Diagonal scaling'
  ! diagonal scaling

  call psb_realloc(n_col,p%d,info)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_realloc')
    goto 9999
  end if

  call psb_csrws(p%d,a,info,trans='N')
  if(info /= 0) then
    info=4010
    ch_err='psb_csrws'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_cdcpy(desc_a,p%desc_Data,info)
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_cdcpy')
    goto 9999
  end if

  if (debug) write(ilout+me,*) 'VDIAG ',n_row
  do i=1,n_row
    if (p%d(i).eq.dzero) then
      p%d(i) = done
    else
      p%d(i) = done/p%d(i)
    endif

    if (debug) write(ilout+me,*) i,desc_a%loc_to_glob(i), p%d(i)
    if (p%d(i).lt.0.d0) then
      write(0,*) me,'Negative RWS? ',i,p%d(i)
    endif
  end do
  if (a%pl(1) /= 0) then
    allocate(work(n_row),stat=info)
    if (info /= 0) then
      info=4000
      call psb_errpush(info,name)
      goto 9999
    end if
    call  psb_gelp('n',a%pl,p%d,desc_a,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_dgelp'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(work)
  endif

  if (debug) then
    allocate(gd(mglob),stat=info)       
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    call  psb_gather(gd, p%d, desc_a, info, iroot=iroot)
    if(info /= 0) then
      info=4010
      ch_err='psb_dgatherm'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (me.eq.iroot) then
      write(iout+np,*) 'VDIAG CHECK ',mglob
      do i=1,mglob
        write(iout+np,*) i,gd(i)
      enddo
    endif
    deallocate(gd)
  endif
  if (debug) write(*,*) 'Preconditioner DIAG computed OK'


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_ddiagsc_bld

