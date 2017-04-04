!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006, 2010, 2015, 2017
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
! File: psb_zspins.f90
!
! Subroutine: psb_zspins
!    Takes a cloud of coefficients and inserts them into a sparse matrix.
!    Note: coefficients with a row index not belonging to the current process are
!    ignored. 
!    If desc_a is in the build state this routine implies a call to psb_cdins. 
! 
! Arguments: 
!    nz       - integer.                    The number of points to insert.
!    ia(:)    - integer                     The row indices of the coefficients.
!    ja(:)    - integer                     The column indices of the coefficients.
!    val(:)   - complex                     The values of the coefficients to be inserted.
!    a        - type(psb_dspmat_type).      The sparse destination matrix.      
!    desc_a   - type(psb_desc_type).        The communication descriptor.
!    info     - integer.                    Error code
!    rebuild  - logical                     Allows to reopen a matrix under
!                                           certain circumstances.
!
subroutine psb_zspins(nz,ia,ja,val,a,desc_a,info,rebuild,local)
  use psb_base_mod, psb_protect_name => psb_zspins
  use psi_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(inout)    :: desc_a
  type(psb_zspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in)                   :: nz,ia(:),ja(:)
  complex(psb_dpk_), intent(in)         :: val(:)
  integer(psb_ipk_), intent(out)                  :: info
  logical, intent(in), optional         :: rebuild, local
  !locals.....

  integer(psb_ipk_) :: nrow, err_act, ncol, spstate
  integer(psb_ipk_) :: ictxt,np,me
  logical, parameter     :: debug=.false.
  integer(psb_ipk_), parameter     :: relocsz=200
  logical                :: rebuild_, local_
  integer(psb_ipk_), allocatable   :: ila(:),jla(:)
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name

  info = psb_success_
  name = 'psb_zspins'
  call psb_erractionsave(err_act)

  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)

  if (nz < 0) then 
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
  if (size(val) < nz) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (nz == 0) return

  if (present(rebuild)) then 
    rebuild_ = rebuild
  else
    rebuild_ = .false.
  endif

  if (present(local)) then 
    local_ = local
  else
    local_ = .false.
  endif

  if (desc_a%is_bld()) then 

    if (local_) then
      info = psb_err_invalid_a_and_cd_state_
      call psb_errpush(info,name)
      goto 9999
    else      
      allocate(ila(nz),jla(nz),stat=info)
      if (info /= psb_success_) then
        ierr(1) = info
        call psb_errpush(psb_err_from_subroutine_ai_,name,&
             & a_err='allocate',i_err=ierr)
        goto 9999
      end if

      call desc_a%indxmap%g2l(ia(1:nz),ila(1:nz),info,owned=.true.)    
      if (info == 0) call desc_a%indxmap%g2l_ins(ja(1:nz),jla(1:nz),info,&
           & mask=(ila(1:nz)>0))

      if (info /= psb_success_) then
        ierr(1) = info
        call psb_errpush(psb_err_from_subroutine_ai_,name,&
             & a_err='psb_cdins',i_err=ierr)
        goto 9999
      end if
      nrow = desc_a%get_local_rows()
      ncol = desc_a%get_local_cols()

      if (a%is_bld()) then 
        call a%csput(nz,ila,jla,val,ione,nrow,ione,ncol,info)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='a%csput')
          goto 9999
        end if
      else
        info = psb_err_invalid_a_and_cd_state_
        call psb_errpush(info,name)
        goto 9999
      end if
    endif

  else if (desc_a%is_asb()) then 

    nrow = desc_a%get_local_rows()
    ncol = desc_a%get_local_cols()
    if (local_) then
      call a%csput(nz,ia,ja,val,ione,nrow,ione,ncol,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='a%csput')
        goto 9999
      end if
    else
      allocate(ila(nz),jla(nz),stat=info)
      if (info /= psb_success_) then
        ierr(1) = info
        call psb_errpush(psb_err_from_subroutine_ai_,name,&
             & a_err='allocate',i_err=ierr)
        goto 9999
      end if

      call desc_a%indxmap%g2l(ia(1:nz),ila(1:nz),info)
      if (info == 0) call desc_a%indxmap%g2l(ja(1:nz),jla(1:nz),info)

      call a%csput(nz,ila,jla,val,ione,nrow,ione,ncol,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='a%csput')
        goto 9999
      end if
    end if
  else
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_zspins


subroutine psb_zspins_2desc(nz,ia,ja,val,a,desc_ar,desc_ac,info)
  use psb_base_mod, psb_protect_name => psb_zspins_2desc
  use psi_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(in)      :: desc_ar
  type(psb_desc_type), intent(inout)   :: desc_ac
  type(psb_zspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in)                  :: nz,ia(:),ja(:)
  complex(psb_dpk_), intent(in)        :: val(:)
  integer(psb_ipk_), intent(out)                 :: info
  !locals.....

  integer(psb_ipk_) :: nrow, err_act, ncol, spstate
  integer(psb_ipk_) :: ictxt,np,me
  logical, parameter     :: debug=.false.
  integer(psb_ipk_), parameter     :: relocsz=200
  integer(psb_ipk_), allocatable   :: ila(:),jla(:)
  integer(psb_ipk_) :: ierr(5)
  character(len=20) :: name

  info = psb_success_
  if (psb_errstatus_fatal()) return 
  name = 'psb_dspins'
  call psb_erractionsave(err_act)
  if (.not.desc_ar%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if
  if (.not.desc_ac%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ictxt = desc_ar%get_context()
  call psb_info(ictxt, me, np)

  if (nz < 0) then 
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
  if (size(val) < nz) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (nz == 0) return

  if (desc_ac%is_bld()) then 

    allocate(ila(nz),jla(nz),stat=info)
    if (info /= psb_success_) then
      ierr(1) = info
      call psb_errpush(psb_err_from_subroutine_ai_,name,&
           & a_err='allocate',i_err=ierr)
      goto 9999
    end if

    call desc_ar%indxmap%g2l(ia(1:nz),ila(1:nz),info,owned=.true.)
    if (info == 0) call desc_ac%indxmap%g2l_ins(ja(1:nz),jla(1:nz),info,&
         & mask=(ila(1:nz)>0))

    if (psb_errstatus_fatal()) then
      ierr(1) = info 
      call psb_errpush(psb_err_from_subroutine_ai_,name,&
           & a_err='psb_cdins',i_err=ierr)
      goto 9999
    end if

    nrow = desc_ar%get_local_rows()
    ncol = desc_ac%get_local_cols()

    call a%csput(nz,ila,jla,val,ione,nrow,ione,ncol,info)
    if (psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='a%csput')
      goto 9999
    end if

  else if (desc_ac%is_asb()) then 

    write(psb_err_unit,*) 'Why are you calling me on an assembled desc_ac?'
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999

  else
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_zspins_2desc


subroutine psb_zspins_v(nz,ia,ja,val,a,desc_a,info,rebuild,local)
  use psb_base_mod, psb_protect_name => psb_zspins_v
  use psi_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(inout)    :: desc_a
  type(psb_zspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in)        :: nz
  type(psb_i_vect_type), intent(inout) :: ia,ja
  type(psb_z_vect_type), intent(inout) :: val
  integer(psb_ipk_), intent(out)                  :: info
  logical, intent(in), optional         :: rebuild, local
  !locals.....

  integer(psb_ipk_) :: nrow, err_act, ncol, spstate
  integer(psb_ipk_) :: ictxt,np,me
  logical, parameter     :: debug=.false.
  integer(psb_ipk_), parameter     :: relocsz=200
  logical                :: rebuild_, local_
  integer(psb_ipk_), allocatable   :: ila(:),jla(:)
  real(psb_dpk_) :: t1,t2,t3,tcnv,tcsput
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name

  info = psb_success_
  name = 'psb_zspins'
  call psb_erractionsave(err_act)

  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)

  if (nz < 0) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (ia%get_nrows() < nz) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if

  if (ja%get_nrows() < nz) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (val%get_nrows() < nz) then 
    info = 1111
    call psb_errpush(info,name)
    goto 9999
  end if
  if (nz == 0) return

  if (present(rebuild)) then 
    rebuild_ = rebuild
  else
    rebuild_ = .false.
  endif

  if (present(local)) then 
    local_ = local
  else
    local_ = .false.
  endif

  if (desc_a%is_bld()) then 

    if (local_) then
      info = psb_err_invalid_a_and_cd_state_
      call psb_errpush(info,name)
      goto 9999
    else      
      allocate(ila(nz),jla(nz),stat=info)
      if (info /= psb_success_) then
        ierr(1) = info
        call psb_errpush(psb_err_from_subroutine_ai_,name,&
             & a_err='allocate',i_err=ierr)
        goto 9999
      end if
      if (ia%is_dev()) call ia%sync()
      if (ja%is_dev()) call ja%sync()
      if (val%is_dev()) call val%sync()

      call desc_a%indxmap%g2l(ia%v%v(1:nz),ila(1:nz),info,owned=.true.)    
      call desc_a%indxmap%g2l_ins(ja%v%v(1:nz),jla(1:nz),info,mask=(ila(1:nz)>0))

      if (info /= psb_success_) then
        ierr(1) = info
        call psb_errpush(psb_err_from_subroutine_ai_,name,&
             & a_err='psb_cdins',i_err=ierr)
        goto 9999
      end if
      nrow = desc_a%get_local_rows()
      ncol = desc_a%get_local_cols()

      if (a%is_bld()) then 
        call a%csput(nz,ila,jla,val%v%v,ione,nrow,ione,ncol,info)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_
          call psb_errpush(info,name,a_err='a%csput')
          goto 9999
        end if
      else
        info = psb_err_invalid_a_and_cd_state_
        call psb_errpush(info,name)
        goto 9999
      end if
    endif

  else if (desc_a%is_asb()) then 

    nrow = desc_a%get_local_rows()
    ncol = desc_a%get_local_cols()
    if (local_) then
      call a%csput(nz,ia,ja,val,ione,nrow,ione,ncol,info)      
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='a%csput')
        goto 9999
      end if
    else
      allocate(ila(nz),jla(nz),stat=info)
      if (info /= psb_success_) then
        ierr(1) = info
        call psb_errpush(psb_err_from_subroutine_ai_,name,&
             & a_err='allocate',i_err=ierr)
        goto 9999
      end if
      if (ia%is_dev()) call ia%sync()
      if (ja%is_dev()) call ja%sync()
      if (val%is_dev()) call val%sync()

      call desc_a%indxmap%g2l(ia%v%v(1:nz),ila(1:nz),info)
      if (info == 0) call desc_a%indxmap%g2l(ja%v%v(1:nz),jla(1:nz),info)
      if (info == 0) call a%csput(nz,ila,jla,val%v%v,ione,nrow,ione,ncol,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='a%csput')
        goto 9999
      end if
    end if
  else
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_zspins_v
