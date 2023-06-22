!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
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
! Function: psb_z_getelem
!    Extract entries from a dense  vector. Note: the row indices in index
!    are assumed to be in global numbering and are converted on the fly.
!    Row indices not belonging to the current process have to be in the halo,
!    othewise failure is ensured.
!
! Arguments:
!    x       - type(psb_z_vect_type) The source vector
!    desc_a  - type(psb_desc_type).    The communication descriptor.
!    index   - integer. Row index of x of the value to extract
!    iam     - integer. Index of the process requesting the value
!    info    - integer.                       return code


function psb_z_getelem(x,index,desc_a,info) result(res)
  use psb_base_mod, psb_protect_name => psb_z_getelem
  use psi_mod
  implicit none

  type(psb_z_vect_type), intent(inout) :: x
  type(psb_desc_type), intent(inout)     :: desc_a
  integer(psb_lpk_), intent(in)          :: index
  integer(psb_ipk_), intent(out)         :: info
  complex(psb_dpk_)                        :: res

  !locals
  integer(psb_ipk_) :: localindex(1)
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me, err_act
  integer(psb_lpk_) :: gindex(1)
  integer(psb_lpk_), allocatable :: myidx(:),mylocal(:)
  character(len=20) :: name
  logical, parameter :: debug = .false.

  gindex(1) = index
  res = zzero
  if (psb_errstatus_fatal()) return
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_z_getelem'

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  call desc_a%indxmap%g2l(gindex,localindex,info,owned=.false.)
  if(debug.and.(localindex(1) < 1)) then
    write(*,*)"Process ",me," owns ",desc_a%get_local_rows()," rows"," Global index is ",gindex,"Local index is ",localindex
    myidx = desc_a%get_global_indices(owned=.false.)
    mylocal = desc_a%get_global_indices(owned=.true.)
    write(*,*)"My (local+halo) indexes are: ",myidx
    write(*,*)"My (local) indexes are: ",mylocal
  end if
  if ( localindex(1) < 1) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err="Index not in the HALO")
    goto 9999
  else
    res = x%get_entry(localindex(1))
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end function

function psb_z_getelem_vec(x,index,desc_a,info) result(res)
  use psb_base_mod, psb_protect_name => psb_z_getelem_vec
  use psi_mod
  implicit none

  type(psb_z_vect_type), intent(inout) :: x
  integer(psb_lpk_), intent(in), dimension(:) :: index
  type(psb_desc_type), intent(inout)     :: desc_a
  integer(psb_ipk_), intent(out)         :: info
  complex(psb_dpk_), allocatable, dimension(:) :: res
  ! Local Variables
  integer(psb_ipk_) :: nindex, i
  integer(psb_ipk_) :: np, me, err_act
  integer(psb_ipk_), allocatable, dimension(:) :: localindex
  integer(psb_lpk_) :: gindex(1)
  type(psb_ctxt_type) :: ctxt
  character(len=20) :: name

  if (psb_errstatus_fatal()) return
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_z_getelem_vec'

  nindex = size(index)

  if (.not.allocated(res)) then
    allocate(res(nindex),stat=info)
  end if
  res = zzero

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ctxt = desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  do i = 1,nindex
    gindex(1) = index(i)
    call desc_a%indxmap%g2l(gindex,localindex,info,owned=.true.)
    if ( localindex(1) > 1) then
        res(i) = x%get_entry(localindex(1))
    end if
  end do

  call psb_sum(ctxt,res,root=-1)


  ! Error handling and closing
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end function

