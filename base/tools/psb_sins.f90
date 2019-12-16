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
!    
! Subroutine: psb_sins_vect
!    Insert entries into a dense  vector. Note: the row indices in IRW 
!    are assumed to be in global numbering and are converted on the fly. 
!    Row indices not belonging to the current process are silently discarded.
! 
! Arguments: 
!    m       - integer.        Number of rows of submatrix belonging to 
!                              val to be inserted.
!    irw(:)  - integer(psb_lpk_) Row indices of rows of val (global numbering)
!    val(:)  - real               The source vector
!    x       - type(psb_s_vect_type) The destination vector
!    desc_a  - type(psb_desc_type).         The communication descriptor.
!    info    - integer.                       return code
!    dupl    - integer               What to do with duplicates: 
!                                     psb_dupl_ovwrt_    overwrite
!                                     psb_dupl_add_      add         
subroutine psb_sins_vect(m, irw, val, x, desc_a, info, dupl,local)
  use psb_base_mod, psb_protect_name => psb_sins_vect
  use psi_mod
  implicit none

  !....parameters...
  integer(psb_ipk_), intent(in)                  :: m
  integer(psb_lpk_), intent(in)                  :: irw(:)
  real(psb_spk_), intent(in)        :: val(:)
  type(psb_s_vect_type), intent(inout) :: x
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), optional, intent(in)        :: dupl
  logical, intent(in), optional        :: local

  !locals.....
  integer(psb_ipk_) :: i, loc_rows,loc_cols
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: ictxt, np, me, dupl_,err_act
  integer(psb_ipk_), allocatable   :: irl(:)
  logical :: local_
  character(len=20)      :: name

  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_sinsvi'

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check parameters....
  if (m < 0) then
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/ione,m/))
    goto 9999
  else if (x%get_nrows() < desc_a%get_local_rows()) then
    info = 310
    call psb_errpush(info,name,i_err=(/5_psb_ipk_,4_psb_ipk_/))
    goto 9999
  endif

  if (m == 0) return  
  loc_rows = desc_a%get_local_rows()
  loc_cols = desc_a%get_local_cols()
  mglob    = desc_a%get_global_rows()

  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif



  allocate(irl(m),stat=info) 
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(dupl)) then 
    dupl_ = dupl
  else
    dupl_ = psb_dupl_ovwrt_
  endif
  if (present(local)) then 
    local_ = local
  else
    local_ = .false.
  endif

  if (local_) then 
    irl(1:m) = irw(1:m)
  else
    call desc_a%indxmap%g2l(irw(1:m),irl(1:m),info,owned=.true.)
  end if
  call x%ins(m,irl,val,dupl_,info) 
  if (info /= 0) then 
    call psb_errpush(info,name)
    goto 9999
  end if
  deallocate(irl)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sins_vect

! Subroutine: psb_sins_vect_v
!    Insert entries into a vector from another encapsulated vector.
!    Note: the row indices in IRW 
!    are assumed to be in global numbering and are converted on the fly. 
!    Row indices not belonging to the current process are silently discarded.
! 
! Arguments: 
!    m       - integer.        Number of rows of submatrix belonging to 
!                              val to be inserted.
!    irw     - type(psb_l_vect_type) Row indices of rows of val (global numbering)
!    val     - type(psb_s_vect_type) The source vector
!    x       - type(psb_s_vect_type) The destination vector
!    desc_a  - type(psb_desc_type).         The communication descriptor.
!    info    - integer.                       return code
!    dupl    - integer               What to do with duplicates: 
!                                     psb_dupl_ovwrt_    overwrite
!                                     psb_dupl_add_      add         
subroutine psb_sins_vect_v(m, irw, val, x, desc_a, info, dupl,local)
  use psb_base_mod, psb_protect_name => psb_sins_vect_v
  use psi_mod
  implicit none

  ! m rows number of submatrix belonging to val to be inserted
  ! ix  x global-row corresponding to position at which val submatrix
  !     must be inserted

  !....parameters...
  integer(psb_ipk_), intent(in)                  :: m
  type(psb_l_vect_type), intent(inout)           :: irw
  type(psb_s_vect_type), intent(inout)        :: val
  type(psb_s_vect_type), intent(inout) :: x
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), optional, intent(in)        :: dupl
  logical, intent(in), optional        :: local

  !locals.....
  integer(psb_ipk_) :: i, loc_rows,loc_cols,err_act
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: ictxt, np, me, dupl_
  integer(psb_ipk_), allocatable   :: irl(:)
  real(psb_spk_), allocatable   :: lval(:)
  logical :: local_
  character(len=20)      :: name

  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_sinsvi_vect_v'

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check parameters....
  if (m < 0) then
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/ione,m/))
    goto 9999
  else if (x%get_nrows() < desc_a%get_local_rows()) then
    info = 310
    call psb_errpush(info,name,i_err=(/5_psb_ipk_,4_psb_ipk_/))
    goto 9999
  endif

  if (m == 0) return  
  loc_rows = desc_a%get_local_rows()
  loc_cols = desc_a%get_local_cols()
  mglob    = desc_a%get_global_rows()

  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif



  if (present(dupl)) then 
    dupl_ = dupl
  else
    dupl_ = psb_dupl_ovwrt_
  endif
  if (present(local)) then 
    local_ = local
  else
    local_ = .false.
  endif

  if (irw%is_dev()) call irw%sync()
  if (local_) then 
    irl(1:m) = irw%v%v(1:m)
  else
    call desc_a%indxmap%g2l(irw%v%v(1:m),irl(1:m),info,owned=.true.)
  end if

  call x%ins(m,irl,lval,dupl_,info) 
  if (info /= 0) then 
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sins_vect_v

subroutine psb_sins_vect_r2(m, irw, val, x, desc_a, info, dupl,local)
  use psb_base_mod, psb_protect_name => psb_sins_vect_r2
  use psi_mod
  implicit none

  ! m rows number of submatrix belonging to val to be inserted
  ! ix  x global-row corresponding to position at which val submatrix
  !     must be inserted

  !....parameters...
  integer(psb_ipk_), intent(in)                  :: m
  integer(psb_lpk_), intent(in)                  :: irw(:)
  real(psb_spk_), intent(in)        :: val(:,:)
  type(psb_s_vect_type), intent(inout) :: x(:)
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), optional, intent(in)        :: dupl
  logical, intent(in), optional        :: local

  !locals.....
  integer(psb_ipk_) :: i, loc_rows,loc_cols, n
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: ictxt, np, me, dupl_, err_act
  integer(psb_ipk_), allocatable   :: irl(:)
  logical :: local_
  character(len=20)      :: name

  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_sinsvi'

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x(1)%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check parameters....
  if (m < 0) then
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/ione,m/))
    goto 9999
  else if (x(1)%get_nrows() < desc_a%get_local_rows()) then
    info = 310
    call psb_errpush(info,name,i_err=(/5_psb_ipk_,4_psb_ipk_/))
    goto 9999
  endif

  if (m == 0) return  
  loc_rows = desc_a%get_local_rows()
  loc_cols = desc_a%get_local_cols()
  mglob    = desc_a%get_global_rows()


  n = min(size(x),size(val,2))
  allocate(irl(m),stat=info) 
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(dupl)) then 
    dupl_ = dupl
  else
    dupl_ = psb_dupl_ovwrt_
  endif
  if (present(local)) then 
    local_ = local
  else
    local_ = .false.
  endif

  if (local_) then 
    irl(1:m) = irw(1:m)
  else
    call desc_a%indxmap%g2l(irw(1:m),irl(1:m),info,owned=.true.)
  end if

  do i=1,n
    if (.not.allocated(x(i)%v)) info = psb_err_invalid_vect_state_
    if (info == 0) call x(i)%ins(m,irl,val(:,i),dupl_,info) 
    if (info /= 0) exit
  end do
  if (info /= 0) then 
    call psb_errpush(info,name)
    goto 9999
  end if
  deallocate(irl)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sins_vect_r2

subroutine psb_sins_multivect(m, irw, val, x, desc_a, info, dupl,local)
  use psb_base_mod, psb_protect_name => psb_sins_multivect
  use psi_mod
  implicit none

  ! m rows number of submatrix belonging to val to be inserted
  ! ix  x global-row corresponding to position at which val submatrix
  !     must be inserted

  !....parameters...
  integer(psb_ipk_), intent(in)                  :: m
  integer(psb_lpk_), intent(in)                  :: irw(:)
  real(psb_spk_), intent(in)        :: val(:,:)
  type(psb_s_multivect_type), intent(inout) :: x
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), optional, intent(in)        :: dupl
  logical, intent(in), optional        :: local

  !locals.....
  integer(psb_ipk_) :: i, loc_rows,loc_cols
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: ictxt, np, me, dupl_, err_act
  integer(psb_ipk_), allocatable   :: irl(:)
  logical :: local_
  character(len=20)      :: name

  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_sinsvi'

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check parameters....
  if (m < 0) then
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/ione,m/))
    goto 9999
  else if (x%get_nrows() < desc_a%get_local_rows()) then
    info = 310
    call psb_errpush(info,name,i_err=(/5_psb_ipk_,4_psb_ipk_/))
    goto 9999
  endif

  if (m == 0) return  
  loc_rows = desc_a%get_local_rows()
  loc_cols = desc_a%get_local_cols()
  mglob    = desc_a%get_global_rows()

  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif



  allocate(irl(m),stat=info) 
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(dupl)) then 
    dupl_ = dupl
  else
    dupl_ = psb_dupl_ovwrt_
  endif
  if (present(local)) then 
    local_ = local
  else
    local_ = .false.
  endif

  if (local_) then 
    irl(1:m) = irw(1:m)
  else
    call desc_a%indxmap%g2l(irw(1:m),irl(1:m),info,owned=.true.)
  end if
  call x%ins(m,irl,val,dupl_,info) 
  if (info /= 0) then 
    call psb_errpush(info,name)
    goto 9999
  end if
  deallocate(irl)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_sins_multivect


