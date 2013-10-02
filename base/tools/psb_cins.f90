!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Subroutine: psb_cinsvi
!    Insert dense submatrix to dense matrix. Note: the row indices in IRW 
!    are assumed to be in global numbering and are converted on the fly. 
!    Row indices not belonging to the current process are silently discarded.
! 
! Arguments: 
!    m       - integer.        Number of rows of submatrix belonging to 
!                              val to be inserted.
!    irw(:)  - integer          Row indices of rows of val (global numbering)
!    val(:)  - complex               The source dense submatrix.  
!    x(:)    - complex               The destination dense matrix.  
!    desc_a  - type(psb_desc_type).         The communication descriptor.
!    info    - integer.                       return code
!    dupl    - integer               What to do with duplicates: 
!                                     psb_dupl_ovwrt_    overwrite
!                                     psb_dupl_add_      add         
subroutine psb_cinsvi(m, irw, val, x, desc_a, info, dupl,local)
  use psb_base_mod, psb_protect_name => psb_cinsvi
  use psi_mod
  implicit none

  ! m rows number of submatrix belonging to val to be inserted

  ! ix  x global-row corresponding to position at which val submatrix
  !     must be inserted

  !....parameters...
  integer(psb_ipk_), intent(in)              ::  m
  integer(psb_ipk_), intent(in)              ::  irw(:)
  complex(psb_spk_), intent(in)  ::  val(:)
  complex(psb_spk_),intent(inout)      ::  x(:)
  type(psb_desc_type), intent(in)  ::  desc_a
  integer(psb_ipk_), intent(out)             ::  info
  integer(psb_ipk_), optional, intent(in)    ::  dupl
  logical, intent(in), optional        :: local

  !locals.....
  integer(psb_ipk_) :: ictxt,i,&
       & loc_rows,loc_cols,mglob,err_act, int_err(5)
  integer(psb_ipk_) :: np, me, dupl_
  integer(psb_ipk_), allocatable   :: irl(:)
  logical :: local_
  character(len=20)      :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cinsvi'

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    return
  end if

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check parameters....
  if (m < 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (size(x, dim=1) < desc_a%get_local_rows()) then
    info = 310
    int_err(1) = 5
    int_err(2) = 4
    call psb_errpush(info,name,int_err)
    goto 9999
  endif

  if (m == 0) return
  loc_rows = desc_a%get_local_rows()
  loc_cols = desc_a%get_local_cols()
  mglob    = desc_a%get_global_rows()

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
  select case(dupl_) 
  case(psb_dupl_ovwrt_) 
    do i = 1, m
      !loop over all val's rows

      ! row actual block row 
      if (irl(i) > 0) then
        ! this row belongs to me
        ! copy i-th row of block val in x
        x(irl(i)) = val(i)
      end if
    enddo

  case(psb_dupl_add_) 

    do i = 1, m
      !loop over all val's rows

      if (irl(i) > 0) then
        ! this row belongs to me
        ! copy i-th row of block val in x
        x(irl(i)) = x(irl(i)) +  val(i)
      end if
    enddo

  case default
    info = 321
    call psb_errpush(info,name)
    goto 9999
  end select
  deallocate(irl)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_cinsvi


subroutine psb_cins_vect(m, irw, val, x, desc_a, info, dupl,local)
  use psb_base_mod, psb_protect_name => psb_cins_vect
  use psi_mod
  implicit none

  ! m rows number of submatrix belonging to val to be inserted
  ! ix  x global-row corresponding to position at which val submatrix
  !     must be inserted

  !....parameters...
  integer(psb_ipk_), intent(in)                  :: m
  integer(psb_ipk_), intent(in)                  :: irw(:)
  complex(psb_spk_), intent(in)        :: val(:)
  type(psb_c_vect_type), intent(inout) :: x
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), optional, intent(in)        :: dupl
  logical, intent(in), optional        :: local

  !locals.....
  integer(psb_ipk_) :: ictxt,i,&
       & loc_rows,loc_cols,mglob,err_act, int_err(5)
  integer(psb_ipk_) :: np, me, dupl_
  integer(psb_ipk_), allocatable   :: irl(:)
  logical :: local_
  character(len=20)      :: name

  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cinsvi'

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
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (x%get_nrows() < desc_a%get_local_rows()) then
    info = 310
    int_err(1) = 5
    int_err(2) = 4
    call psb_errpush(info,name,int_err)
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

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_cins_vect

subroutine psb_cins_vect_r2(m, irw, val, x, desc_a, info, dupl,local)
  use psb_base_mod, psb_protect_name => psb_cins_vect_r2
  use psi_mod
  implicit none

  ! m rows number of submatrix belonging to val to be inserted
  ! ix  x global-row corresponding to position at which val submatrix
  !     must be inserted

  !....parameters...
  integer(psb_ipk_), intent(in)                  :: m
  integer(psb_ipk_), intent(in)                  :: irw(:)
  complex(psb_spk_), intent(in)        :: val(:,:)
  type(psb_c_vect_type), intent(inout) :: x(:)
  type(psb_desc_type), intent(in)      :: desc_a
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), optional, intent(in)        :: dupl
  logical, intent(in), optional        :: local

  !locals.....
  integer(psb_ipk_) :: ictxt,i,&
       & loc_rows,loc_cols,mglob,err_act, int_err(5), n
  integer(psb_ipk_) :: np, me, dupl_
  integer(psb_ipk_), allocatable   :: irl(:)
  logical :: local_
  character(len=20)      :: name

  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cinsvi'

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
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (x(1)%get_nrows() < desc_a%get_local_rows()) then
    info = 310
    int_err(1) = 5
    int_err(2) = 4
    call psb_errpush(info,name,int_err)
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

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_cins_vect_r2



!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Subroutine: psb_cinsi
!    Insert dense submatrix to dense matrix. Note: the row indices in IRW 
!    are assumed to be in global numbering and are converted on the fly. 
!    Row indices not belonging to the current process are silently discarded.
! 
! Arguments: 
!    m       - integer.        Number of rows of submatrix belonging to 
!                              val to be inserted.
!    irw(:)   - integer          Row indices of rows of val (global numbering)
!    val(:,:) - complex                 The source dense submatrix.  
!    x(:,:)   - complex                 The destination dense matrix.  
!    desc_a   - type(psb_desc_type).         The communication descriptor.
!    info     - integer.                       return code
!    dupl    - integer               What to do with duplicates: 
!                                     psb_dupl_ovwrt_    overwrite
!                                     psb_dupl_add_      add         
subroutine psb_cinsi(m, irw, val, x, desc_a, info, dupl,local)
  use psb_base_mod, psb_protect_name => psb_cinsi
  use psi_mod
  implicit none

  ! m rows number of submatrix belonging to val to be inserted

  ! ix  x global-row corresponding to position at which val submatrix
  !     must be inserted

  !....parameters...
  integer(psb_ipk_), intent(in)             ::  m
  integer(psb_ipk_), intent(in)             ::  irw(:)
  complex(psb_spk_), intent(in) ::  val(:,:)
  complex(psb_spk_),intent(inout)     ::  x(:,:)
  type(psb_desc_type), intent(in) ::  desc_a
  integer(psb_ipk_), intent(out)            ::  info
  integer(psb_ipk_), optional, intent(in)   ::  dupl
  logical, intent(in), optional        :: local

  !locals.....
  integer(psb_ipk_) :: ictxt,i,loc_row,j,n,&
       & loc_rows,loc_cols,mglob,err_act, int_err(5)
  integer(psb_ipk_) :: np,me,dupl_
  integer(psb_ipk_), allocatable   :: irl(:)
  logical :: local_
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cinsi'

  if (.not.desc_a%is_ok()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    return
  end if

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  !... check parameters....
  if (m < 0) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
    call psb_errpush(info,name,int_err)
    goto 9999
  else if (size(x, dim=1) < desc_a%get_local_rows()) then
    info = 310
    int_err(1) = 5
    int_err(2) = 4
    call psb_errpush(info,name,int_err)
    goto 9999
  endif
  if (m == 0) return 

  loc_rows = desc_a%get_local_rows()
  loc_cols = desc_a%get_local_cols()
  mglob    = desc_a%get_global_rows()

  n = min(size(val,2),size(x,2))

  if (present(dupl)) then 
    dupl_ = dupl
  else
    dupl_ = psb_dupl_ovwrt_
  endif

  allocate(irl(m),stat=info) 
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
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
  
  select case(dupl_) 
  case(psb_dupl_ovwrt_) 
    do i = 1, m
      !loop over all val's rows

      ! row actual block row 
      loc_row = irl(i)
      if (loc_row > 0) then
        ! this row belongs to me
        ! copy i-th row of block val in x
        do j=1,n
          x(loc_row,j) = val(i,j)
        end do
      end if
    enddo

  case(psb_dupl_add_) 

    do i = 1, m
      !loop over all val's rows

      ! row actual block row 
      loc_row = irl(i)
      if (loc_row > 0) then
        ! this row belongs to me
        ! copy i-th row of block val in x
        do j=1,n
          x(loc_row,j) = x(loc_row,j) +  val(i,j)
        end do
      end if
    enddo

  case default
    info = 321
    call psb_errpush(info,name)
    goto 9999
  end select
  deallocate(irl)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_cinsi


