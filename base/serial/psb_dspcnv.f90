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
! File:  psb_dspcnv.f90 
!
! Subroutine: psb_dspcnv2
!    This subroutine converts the storage format of a matrix.
!
! Arguments:
! 
!   a      -  type(psb_spmat_type), input     The input matrix to be converted.
!   b      -  type(psb_spmat_type), output    The assembled output matrix.
!   info   -  integer.                        Return code
!   afmt   -  character, optional             The desired storage format 
!
subroutine psb_dspcnv2(a, b,info,afmt,upd,dupl)
  use psb_const_mod
  use psb_error_mod
  use psb_spmat_type
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_dspcnv2
  implicit none
  !....Parameters...
  Type(psb_dspmat_type), intent(in)      :: A
  Type(psb_dspmat_type), intent(out)     :: B
  Integer, intent(out)                   :: info
  Integer, intent(in), optional          :: upd,dupl
  character(len=*), optional, intent(in) :: afmt

  !...Locals...
  real(psb_dpk_)              :: d(1)
  real(psb_dpk_), allocatable :: work(:)
  type(psb_dspmat_type)         :: temp_a
  Integer                       :: nzr, ntry, ifc_, ia1_size,&
       & ia2_size, aspk_size,size_req,n_row,n_col,upd_,dupl_
  integer                       :: err_act
  character                     :: check_,trans_,unitd_
  character(len=5)              :: afmt_
  Integer, Parameter            :: maxtry=8
  integer                       :: debug_level, debug_unit
  character(len=20)             :: name, ch_err

  name='psb_spcnv2'
  info  = 0
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ntry=0
  
  ifc_ = 2


  check_ = 'N'
  trans_ = 'N'
  unitd_ = 'U'
  allocate(work(max(size(a%ia1),size(a%ia2))+a%m+1000),stat=info)

  if (info /= 0) then
    info=2040
    call psb_errpush(info,name)
    goto 9999
  end if


  if (present(upd)) then 
    call psb_sp_setifld(upd,psb_upd_,b,info)    
  end if
  if (present(dupl)) then 
    call psb_sp_setifld(dupl,psb_dupl_,b,info)    
  end if
  if (present(afmt)) then 
    afmt_ = psb_tolower(afmt)
  else
    afmt_ = psb_fidef_
  end if

  upd_ = psb_sp_getifld(psb_upd_,b,info)
  select case(upd_)
  case(psb_upd_srch_,psb_upd_perm_)
    ! Legal value, do nothing
  case default
    ! Fix bad value
    upd_ = psb_upd_dflt_
    call psb_sp_setifld(upd_,psb_upd_,b,info)    
  end select

  dupl_ = psb_sp_getifld(psb_dupl_,b,info)
  select case(dupl_)
  case(psb_dupl_ovwrt_,psb_dupl_add_,psb_dupl_err_)
    ! Legal value, do nothing
  case default
    ! Fix bad value
    dupl_ = psb_dupl_def_
    call psb_sp_setifld(dupl_,psb_dupl_,b,info)    
  end select

  !  ...matrix conversion...
  b%m=a%m
  b%k=a%k
  b%fida=psb_tolower(afmt_)
  size_req = psb_sp_get_nnzeros(a)
  !
  n_row=b%m 
  n_col=b%k
  call psb_cest(afmt_,n_row,n_col,size_req,&
       & ia1_size, ia2_size, aspk_size, upd_,info)
  b%fida=afmt_
  
  if (info /= psb_no_err_) then    
    info=4010
    ch_err='psb_cest'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif
  if (.not.allocated(b%aspk).or.&
       &.not.allocated(b%ia1).or.&
       &.not.allocated(b%ia2).or.&
       &.not.allocated(b%pl).or.&
       &.not.allocated(b%pr)) then
    call psb_sp_reall(b,ia1_size,ia2_size,aspk_size,info)
  else if ((size(b%aspk)  < aspk_size) .or.&
       &(size(b%ia1) < ia1_size) .or.&
       &(size(b%ia2) < ia2_size) .or.&
       &(size(b%pl)  < b%m) .or.&
       &(size(b%pr) < b%k )) then 
    call psb_sp_reall(b,ia1_size,ia2_size,aspk_size,info)
  endif

  if (info /= psb_no_err_) then    
    info=4010
    ch_err='psb_sp_reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  endif

  b%pl(:)  = 0
  b%pr(:)  = 0

  b%descra = a%descra
  if (debug_level >= psb_debug_serial_) &
       & write(debug_unit,*) trim(name),': size_req 1:',&
       & size_req, trans_,upd_,dupl_,b%fida,b%descra

  select case (psb_tolower(a%fida))

  case ('csr')

    select case (psb_tolower(b%fida))

    case ('csr')

      call dcrcr(trans_, a%m, a%k, unitd_, d, a%descra, a%aspk,&
           & a%ia1, a%ia2, a%infoa, b%pl, b%descra, b%aspk, b%ia1,&
           & b%ia2, b%infoa, b%pr, size(b%aspk), size(b%ia1),&
           & size(b%ia2), work, size(work), info)


      if (info/=0) then
        info=4010
        ch_err='dcrcr'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    case ('jad')

      !...converting to JAD
      !...output matrix may not be big enough
      do

        call dcrjd(trans_, a%m, a%k, unitd_, d, a%descra, a%aspk,&
             & a%ia1, a%ia2, a%infoa, b%pl, b%descra, b%aspk, b%ia1,&
             & b%ia2, b%infoa, b%pr, size(b%aspk), size(b%ia1),&
             & size(b%ia2), work, size(work), nzr, info)
        if (info /= 0) then
          call psb_errpush(4010,name,a_err='dcrjd')
          goto 9999
        endif

        ntry = ntry + 1
        if (debug_level >= psb_debug_serial_) then 
          write(debug_unit,*) trim(name),' On out from dcrjad ',nzr,info
        end if
        if (nzr == 0) exit
        if (ntry > maxtry ) then 
          write(debug_unit,*) trim(name),&
               & ' Tried reallocating for DCRJAD for ',maxtry,': giving up now.'
          info=2040
          call psb_errpush(info,name)
          goto 9999
        endif
        call psb_cest(afmt_,n_row,n_col,nzr,&
             & ia1_size, ia2_size, aspk_size, upd_,info)
        call psb_sp_reall(b,ia1_size,ia2_size,aspk_size,info)
        if (info /= 0) then
          info=2040
          call psb_errpush(info,name)
          goto 9999
        endif

      end do

      if (info/=0) then
        call psb_errpush(info,name)
        goto 9999
      end if

    case ('coo')
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),' Calling CRCO ',a%descra
      call dcrco(trans_, a%m, a%k, unitd_, d, a%descra, a%aspk,&
           & a%ia1, a%ia2, a%infoa, b%pl, b%descra, b%aspk, b%ia1,&
           & b%ia2, b%infoa, b%pr, size(b%aspk), size(b%ia1),&
           & size(b%ia2), work, size(work), info)

      if (info/=0) then
        call psb_errpush(4010,name,a_err='dcrco')
        goto 9999
      end if
    case default
      info=4010
      call psb_errpush(info,name)
      goto 9999

    end select

  case ('coo','coi')

    select case (psb_tolower(b%fida))

    case ('csr')

      call dcocr(trans_, a%m, a%k, unitd_, d, a%descra, a%aspk,&
           & a%ia2, a%ia1, a%infoa, b%pl, b%descra, b%aspk, b%ia1,&
           & b%ia2, b%infoa, b%pr, size(b%aspk), size(b%ia1),&
           & size(b%ia2), work, 2*size(work), info)

      if (info/=0) then
        call psb_errpush(4010,name,a_err='dcocr')
        goto 9999
      end if

    case ('jad')

      call psb_sp_all(temp_a, size(b%ia1),size(b%ia2),size(b%aspk),info)
      if (info /= 0) then
        info=2040
        call psb_errpush(info,name)
        goto 9999
      endif
      temp_a%m = a%m
      temp_a%k = a%k
      temp_a%infoa=b%infoa
      !...Dirty trick: converting to CSR and then to JAD

      call dcocr(trans_, a%m, a%k, unitd_, d, a%descra, a%aspk,&
           & a%ia2, a%ia1, a%infoa, temp_a%pl, temp_a%descra, &
           & temp_a%aspk, temp_a%ia1, temp_a%ia2, temp_a%infoa, temp_a%pr, &
           & size(temp_a%aspk), size(temp_a%ia1),&
           & size(temp_a%ia2), work, 2*size(work), info)

      if (info/=0) then
        call psb_errpush(4010,name,a_err='dcocr')
        goto 9999
      end if

      do
        call dcrjd(trans_, temp_a%m, temp_a%k, unitd_, d, temp_a%descra, &
             & temp_a%aspk, temp_a%ia1, temp_a%ia2, temp_a%infoa, &
             & b%pl, b%descra, b%aspk, b%ia1, b%ia2, b%infoa, b%pr, &
             & size(b%aspk), size(b%ia1),&
             & size(b%ia2), work, size(work), nzr, info)
        if (info/=0) then
          call psb_errpush(4010,name,a_err='dcrjd')
          goto 9999
        end if

        ntry = ntry + 1
        if (debug_level >= psb_debug_serial_) then 
          write(debug_unit,*) trim(name),' On out from dcrjad ',nzr,info
        end if
        if (nzr == 0) exit
        if (ntry > maxtry ) then 
          write(debug_unit,*) trim(name),' Tried reallocating for DCRJAD for ',&
               & maxtry,': giving up now.'
          info=2040
          call psb_errpush(info,name)
          goto 9999
        endif

        call psb_sp_reall(b,nzr,info,ifc=ifc_)
        if (info /= 0) then
          info=2040
          call psb_errpush(info,name)
          goto 9999
        endif

      end do

    case ('coo')

      call dcoco(trans_, a%m, a%k, unitd_, d, a%descra, a%aspk,&
           & a%ia1, a%ia2, a%infoa, b%pl, b%descra, b%aspk, b%ia1,&
           & b%ia2, b%infoa, b%pr, size(b%aspk), size(b%ia1),&
           & size(b%ia2), work, 2*size(work), info)
      if (info/=0) then
        call psb_errpush(4010,name,a_err='dcoco')
        goto 9999
      end if

    case default
      info=4010
      call psb_errpush(info,name)
      goto 9999
    end select

  case default
    info=4010
    call psb_errpush(info,name)
    goto 9999

  end select

  if (psb_sp_getifld(psb_upd_,b,info) /= psb_upd_perm_) then 
    call psb_sp_trim(b,info)
  end if


  call psb_sp_setifld(psb_spmat_asb_,psb_state_,b,info)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dspcnv2


!
! Subroutine: psb_dspcnv1
!    This subroutine converts in place the storage format of a matrix.
!
! Arguments:
! 
!   a      -  type(psb_spmat_type), inout     The input matrix to be converted.
!   info   -  integer.                        Return code
!   afmt   -  character, optional             The desired storage format 
!
subroutine psb_dspcnv1(a, info, afmt, upd, dupl)

  use psb_spmat_type
  use psb_regen_mod
  use psb_serial_mod, psb_protect_name => psb_dspcnv1
  use psb_const_mod
  use psi_mod
  use psb_error_mod
  use psb_string_mod
  !use psb_penv_mod

  implicit none


  !...Parameters....
  type(psb_dspmat_type), intent (inout)   :: a
  integer, intent(out)                    :: info
  integer,optional, intent(in)            :: dupl, upd
  character(len=*), optional, intent(in)         :: afmt

  integer               :: int_err(5)
  type(psb_dspmat_type) :: atemp
  integer               :: err_act
  integer               :: spstate
  integer               :: upd_, dupl_
  integer               :: debug_level, debug_unit
  character(len=20)     :: name, ch_err
  logical               :: inplace

  info = 0
  int_err(1)=0
  name = 'psb_spcnv1'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if (present(upd)) then 
    call psb_sp_setifld(upd,psb_upd_,a,info)    
  end if
  if (present(dupl)) then 
    call psb_sp_setifld(dupl,psb_dupl_,a,info)    
  end if

  upd_ = psb_sp_getifld(psb_upd_,a,info)
  select case(upd_)
  case(psb_upd_srch_,psb_upd_perm_)
    ! Legal value, do nothing
  case default
    ! Fix bad value
    upd_ = psb_upd_dflt_
    call psb_sp_setifld(upd_,psb_upd_,a,info)    
  end select

  dupl_ = psb_sp_getifld(psb_dupl_,a,info)
  select case(dupl_)
  case(psb_dupl_ovwrt_,psb_dupl_add_,psb_dupl_err_)
    ! Legal value, do nothing
  case default
    ! Fix bad value
    dupl_ = psb_dupl_def_
    call psb_sp_setifld(dupl_,psb_dupl_,a,info)    
  end select


  spstate = psb_sp_getifld(psb_state_,a,info)
  if (info /= 0) then 
    goto 9999
  endif
  if (debug_level >= psb_debug_serial_) &
       & write(debug_unit,*) trim(name),': Sparse matrix state:',&
       & spstate,psb_spmat_bld_,psb_spmat_upd_
  if (spstate /= psb_spmat_upd_) then 
    ! Should we first figure out if we can do it in place? 
    if (debug_level >= psb_debug_serial_) &
         & write(debug_unit,*) trim(name),&
         & ': Update:',upd_,psb_upd_srch_,psb_upd_perm_
    inplace = .false.
    if (upd_ == psb_upd_srch_) then 
      if (present(afmt)) then 
        select case (psb_tolower(a%fida))
        case('coo')
          select case(psb_tolower(afmt))
          case('coo') 
            call psb_fixcoo(a,info)
            inplace = .true.
          case('csr') 
            call psb_ipcoo2csr(a,info)
            inplace = .true.
          case('csc')
            call psb_ipcoo2csc(a,info)
            inplace = .true.
          end select
        case('csr')
          select case(psb_tolower(afmt))
          case('coo') 
            call psb_ipcsr2coo(a,info)
            inplace = .true.
          end select
        end select
      end if
      if (inplace) then 
        if (info == 0) call psb_sp_trim(a,info)
        if (info /= 0) then
          info=4010
          ch_err='inplace cnv'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        endif
        call psb_sp_setifld(psb_spmat_asb_,psb_state_,a,info)
        call psb_erractionrestore(err_act)
        return
      end if
    end if

    call psb_sp_clone(a,atemp,info)
    if(info /= psb_no_err_) then
      info=4010
      ch_err='psb_sp_clone'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
      ! convert to user requested format after the temp copy
    end if

    ! Do the real conversion into the requested storage format
    ! result is put in A
    call psb_spcnv(atemp,a,info,afmt=afmt,upd=upd,dupl=dupl)

    if (debug_level >= psb_debug_serial_)&
         & write(debug_unit, *) trim(name),&
         & ':  From SPCNV',info,' ',a%fida

    if (info /= psb_no_err_) then    
      info=4010
      ch_err='psb_csdp'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif

    call psb_sp_free(atemp,info)

  else if (spstate == psb_spmat_upd_) then
    !
    ! Second  case: we come from an update loop.
    ! 
    select case(psb_tolower(a%fida))
    case('csr')
      call csr_regen(a,info)
    case ('coo','coi')
      call coo_regen(a,info)
    case('jad')
      call jad_regen(a,info)
    case default
      info=136
      ch_err=a%fida 
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end select
    
    if (info /= psb_no_err_) then
      info = 4010
      ch_err='xx_regen'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  else

    info = 600
    call psb_errpush(info,name)
    goto 9999
    if (debug_level >= psb_debug_serial_)&
         & write(debug_unit,*) trim(name),&
         & 'Sparse matrix state:',spstate,psb_spmat_bld_,psb_spmat_upd_

  endif

  call psb_sp_setifld(psb_spmat_asb_,psb_state_,a,info)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
  
end subroutine psb_dspcnv1
