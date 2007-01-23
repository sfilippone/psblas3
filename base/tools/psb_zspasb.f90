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
! File: psb_zspasb.f90
!
! Subroutine: psb_zspasb
!    Assembly sparse matrix and set psblas communications
!    structures.
! 
! Parameters: 
!    a        - type(<psb_zspmat_type>).          The sparse matrix to be allocated.      
!    desc_a   - type(<psb_desc_type>).            The communication descriptor to be updated.
!    info     - integer.                          Eventually returns an error code.
!    afmt     - character,dimension(5)(optional). The output format.
!    up       - character(optional).              ???
!    dup      - integer(optional).                ???
!
subroutine psb_zspasb(a,desc_a, info, afmt, upd, dupl)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psi_mod
  use psb_error_mod
  use psb_string_mod
  use psb_penv_mod
  implicit none


  !...Parameters....
  type(psb_zspmat_type), intent (inout)   :: a
  type(psb_desc_type), intent(in)         :: desc_a
  integer, intent(out)                    :: info
  integer,optional, intent(in)            :: dupl, upd
  character, optional, intent(in)         :: afmt*5
  !....Locals....
  integer               :: int_err(5)
  type(psb_zspmat_type) :: atemp
  integer               :: np,me,n_col,iout, err_act
  integer               :: spstate
  integer               :: upd_, dupl_
  integer               :: ictxt,n_row
  logical, parameter    :: debug=.false., debugwrt=.false.
  character(len=20)     :: name, ch_err

  info = 0
  int_err(1)=0
  name = 'psb_spasb'
  call psb_erractionsave(err_act)

  ictxt    = psb_cd_get_context(desc_a)
  n_row    = psb_cd_get_local_rows(desc_a)
  n_col    = psb_cd_get_local_cols(desc_a)

  ! check on BLACS grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_asb_desc(desc_a)) then 
    info = 600
    int_err(1) = psb_cd_get_dectype(desc_a)
    call psb_errpush(info,name)
    goto 9999
  endif

  if (debug) Write (*, *) '   Begin matrix assembly...'

  !check on errors encountered in psdspins

  spstate     = a%infoa(psb_state_) 
  if (spstate == psb_spmat_bld_) then 
    !
    ! First case: we come from a fresh build. 
    ! 

    n_row = psb_cd_get_local_rows(desc_a)
    n_col = psb_cd_get_local_cols(desc_a)

    !
    ! Second step: handle the local matrix part. 
    !
    if (present(upd)) then 
      upd_=upd
    else
      upd_ = psb_upd_dflt_
    endif

    if (present(dupl)) then
      select case(dupl) 
      case (psb_dupl_ovwrt_,psb_dupl_add_,psb_dupl_err_)
        dupl_ = dupl
      case default
        dupl_ =  psb_dupl_def_
      end select
    else
      dupl_ =  psb_dupl_def_
    endif

    a%m = n_row
    a%k = n_col

    call psb_sp_clone(a,atemp,info)
    if(info /= psb_no_err_) then
      info=4010
      ch_err='psb_sp_clone'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
      ! convert to user requested format after the temp copy
    end if

    if (present(afmt)) then
      a%fida = afmt
    else 
      a%fida = '???'
    endif

    if (debugwrt) then
      iout = 30+me
      open(iout)
      call psb_csprt(iout,atemp,head='Input mat')
      close(iout)
    endif

    ! Do the real conversion into the requested storage format
    ! result is put in A
    call psb_csdp(atemp,a,info,ifc=2,upd=upd_,dupl=dupl_)

    IF (debug) WRITE (*, *) me,'   ASB:  From ZCSDP',info,' ',A%FIDA
    if (info /= psb_no_err_) then    
      info=4010
      ch_err='psb_csdp'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif

    if (debugwrt) then
      iout = 60+me
      open(iout)
      call psb_csprt(iout,a,head='Output mat')
      close(iout)
    endif

    call psb_sp_free(atemp,info)


  else if (spstate == psb_spmat_upd_) then
    !
    ! Second  case: we come from an update loop.
    ! 


    ! Right now, almost nothing to be done, but this 
    ! may change in the future
    ! as we revise the implementation of the update routine. 
    call psb_sp_all(atemp,1,info)
    atemp%m=a%m
    atemp%k=a%k
    ! check on allocation
    if (info /= psb_no_err_) then    
      info=4010
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif

    call psb_csdp(atemp,a,info,check='R')
    ! check on error retuned by zcsdp
    if (info /= psb_no_err_) then
      info = 4010
      ch_err='psb_csdp90'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_sp_free(atemp,info)
    if (info /= psb_no_err_) then
      info = 4010
      ch_err='sp_free'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  else

    info = 600
    call psb_errpush(info,name)
    goto 9999
    if (debug) write(0,*) 'Sparse matrix state:',spstate,psb_spmat_bld_,psb_spmat_upd_

  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_zspasb
