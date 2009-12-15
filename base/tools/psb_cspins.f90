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
! File: psb_cspins.f90
!
! Subroutine: psb_cspins
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
!    a        - type(psb_d_sparse_mat).      The sparse destination matrix.      
!    desc_a   - type(psb_desc_type).        The communication descriptor.
!    info     - integer.                    Error code
!    rebuild  - logical                     Allows to reopen a matrix under
!                                           certain circumstances.
!
subroutine psb_cspins(nz,ia,ja,val,a,desc_a,info,rebuild)
  use psb_sparse_mod, psb_protect_name => psb_cspins
  use psi_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(inout)    :: desc_a
  type(psb_c_sparse_mat), intent(inout) :: a
  integer, intent(in)                   :: nz,ia(:),ja(:)
  complex(psb_spk_), intent(in)         :: val(:)
  integer, intent(out)                  :: info
  logical, intent(in), optional         :: rebuild
  !locals.....

  integer :: nrow, err_act, ncol, spstate
  integer                :: ictxt,np,me
  logical, parameter     :: debug=.false.
  integer, parameter     :: relocsz=200
  logical                :: rebuild_
  integer, allocatable   :: ila(:),jla(:)
  character(len=20)  :: name, ch_err

  info = 0
  name = 'psb_cspins'
  call psb_erractionsave(err_act)


  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)

  if (.not.psb_is_ok_desc(desc_a)) then 
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif

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
  if (nz==0) return

  if (present(rebuild)) then 
    rebuild_ = rebuild
  else
    rebuild_ = .false.
  endif

  if (psb_is_bld_desc(desc_a)) then 
    if (psb_is_large_desc(desc_a)) then 

      allocate(ila(nz),jla(nz),stat=info)
      if (info /= 0) then
        ch_err='allocate'
        call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
        goto 9999
      end if
      call  psb_cdins(nz,ia,ja,desc_a,info,ila=ila,jla=jla)
      if (info /= 0) then
        ch_err='psb_cdins'
        call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
        goto 9999
      end if
      nrow = psb_cd_get_local_rows(desc_a)
      ncol = psb_cd_get_local_cols(desc_a)

      if (a%is_bld()) then 
        call a%csput(nz,ila,jla,val,1,nrow,1,ncol,info)
        if (info /= 0) then
          info=4010
          ch_err='psb_coins'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
      else
        info = 1123
        call psb_errpush(info,name)
        goto 9999
      end if

    else

      call  psb_cdins(nz,ia,ja,desc_a,info)
      if (info /= 0) then
        ch_err='psb_cdins'
        call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
        goto 9999
      end if
      nrow = psb_cd_get_local_rows(desc_a)
      ncol = psb_cd_get_local_cols(desc_a)

      if (a%is_bld()) then 
        call a%csput(nz,ia,ja,val,1,nrow,1,ncol,info,gtl=desc_a%idxmap%glob_to_loc)
        if (info /= 0) then
          info=4010
          ch_err='psb_coins'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
      else
        info = 1123
        call psb_errpush(info,name)
        goto 9999
      end if

    end if
    
  else if (psb_is_asb_desc(desc_a)) then 

    if (psb_is_large_desc(desc_a)) then 

      allocate(ila(nz),jla(nz),stat=info)
      if (info /= 0) then
        ch_err='allocate'
        call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
        goto 9999
      end if

      ila(1:nz) = ia(1:nz)
      jla(1:nz) = ja(1:nz)
      call psb_glob_to_loc(ila(1:nz),desc_a,info,iact='I')
      call psb_glob_to_loc(jla(1:nz),desc_a,info,iact='I')
      nrow = psb_cd_get_local_rows(desc_a)
      ncol = psb_cd_get_local_cols(desc_a)

      call a%csput(nz,ila,jla,val,1,nrow,1,ncol,info)
      if (info /= 0) then
        info=4010
        ch_err='psb_coins'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    else
      nrow = psb_cd_get_local_rows(desc_a)
      ncol = psb_cd_get_local_cols(desc_a)
      call a%csput(nz,ia,ja,val,1,nrow,1,ncol,&
           & info,gtl=desc_a%idxmap%glob_to_loc)
      if (info /= 0) then
        info=4010
        ch_err='psb_coins'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if
  else
    info = 1122
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_cspins


subroutine psb_cspins_2desc(nz,ia,ja,val,a,desc_ar,desc_ac,info)
  use psb_sparse_mod, psb_protect_name => psb_cspins_2desc
  use psi_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(in)      :: desc_ar
  type(psb_desc_type), intent(inout)   :: desc_ac
  type(psb_c_sparse_mat), intent(inout) :: a
  integer, intent(in)                  :: nz,ia(:),ja(:)
  complex(psb_spk_), intent(in)        :: val(:)
  integer, intent(out)                 :: info
  !locals.....

  integer :: nrow, err_act, ncol, spstate
  integer                :: ictxt,np,me
  logical, parameter     :: debug=.false.
  integer, parameter     :: relocsz=200
  integer, allocatable   :: ila(:),jla(:)
  character(len=20)  :: name, ch_err

  info = 0
  name = 'psb_cspins'
  call psb_erractionsave(err_act)


  ictxt = psb_cd_get_context(desc_ar)

  call psb_info(ictxt, me, np)

  if (.not.psb_is_ok_desc(desc_ar)) then 
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.psb_is_ok_desc(desc_ac)) then 
    info = 3110
    call psb_errpush(info,name)
    goto 9999
  endif

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
  if (nz==0) return

  if (psb_is_bld_desc(desc_ac)) then 

    allocate(ila(nz),jla(nz),stat=info)
    if (info /= 0) then
      ch_err='allocate'
      call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
      goto 9999
    end if
        ila(1:nz) = ia(1:nz)

    call psb_glob_to_loc(ia(1:nz),ila(1:nz),desc_ar,info,iact='I',owned=.true.)

    call psb_cdins(nz,ja,desc_ac,info,jla=jla, mask=(ila(1:nz)>0))

    if (info /= 0) then
      ch_err='psb_cdins'
      call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
      goto 9999
    end if

    nrow = psb_cd_get_local_rows(desc_ar)
    ncol = psb_cd_get_local_cols(desc_ac)

    call a%csput(nz,ila,jla,val,1,nrow,1,ncol,info)
    if (info /= 0) then
      info=4010
      ch_err='psb_coins'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  else if (psb_is_asb_desc(desc_ac)) then 

    write(0,*) 'Why are you calling me on an assembled desc_ac?'
!!$    if (psb_is_large_desc(desc_a)) then 
!!$
!!$      allocate(ila(nz),jla(nz),stat=info)
!!$      if (info /= 0) then
!!$        ch_err='allocate'
!!$        call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
!!$        goto 9999
!!$      end if
!!$
!!$      ila(1:nz) = ia(1:nz)
!!$      jla(1:nz) = ja(1:nz)
!!$      call psb_glob_to_loc(ila(1:nz),desc_a,info,iact='I')
!!$      call psb_glob_to_loc(jla(1:nz),desc_a,info,iact='I')
!!$      nrow = psb_cd_get_local_rows(desc_a)
!!$      ncol = psb_cd_get_local_cols(desc_a)
!!$
!!$      call psb_coins(nz,ila,jla,val,a,1,nrow,1,ncol,&
!!$           & info,rebuild=rebuild_)
!!$      if (info /= 0) then
!!$        info=4010
!!$        ch_err='psb_coins'
!!$        call psb_errpush(info,name,a_err=ch_err)
!!$        goto 9999
!!$      end if
!!$
!!$    else
!!$      nrow = psb_cd_get_local_rows(desc_a)
!!$      ncol = psb_cd_get_local_cols(desc_a)
!!$      call psb_coins(nz,ia,ja,val,a,1,nrow,1,ncol,&
!!$           & info,gtl=desc_a%idxmap%glob_to_loc,rebuild=rebuild_)
!!$      if (info /= 0) then
!!$        info=4010
!!$        ch_err='psb_coins'
!!$        call psb_errpush(info,name,a_err=ch_err)
!!$        goto 9999
!!$      end if
!!$    end if
  else
    info = 1122
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_cspins_2desc

