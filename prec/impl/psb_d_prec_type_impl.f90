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
!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine psb_d_apply2_vect(prec,x,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_d_apply2_vect
  implicit none
  type(psb_desc_type),intent(in)       :: desc_data
  class(psb_dprec_type), intent(inout) :: prec
  type(psb_d_vect_type),intent(inout)  :: x
  type(psb_d_vect_type),intent(inout)  :: y
  integer(psb_ipk_), intent(out)                 :: info
  character(len=1), optional           :: trans
  real(psb_dpk_),intent(inout), optional, target :: work(:)

  character     :: trans_
  real(psb_dpk_), pointer :: work_(:)
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: err_act
  character(len=20)   :: name

  name = 'psb_d_apply2v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt = desc_data%get_context()
  call psb_info(ictxt, me, np)

  if (present(trans)) then
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (present(work)) then
    work_ => work
  else
    allocate(work_(4*desc_data%get_local_cols()),stat=info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if

  end if

  if (.not.allocated(prec%prec)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if

  call prec%prec%apply(done,x,dzero,y,desc_data,info,&
       & trans=trans_,work=work_)

  if (present(work)) then
  else
    deallocate(work_,stat=info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999
    end if
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_d_apply2_vect

subroutine psb_d_apply1_vect(prec,x,desc_data,info,trans,work)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_d_apply1_vect
  implicit none
  type(psb_desc_type),intent(in)       :: desc_data
  class(psb_dprec_type), intent(inout) :: prec
  type(psb_d_vect_type),intent(inout)  :: x
  integer(psb_ipk_), intent(out)                 :: info
  character(len=1), optional           :: trans
  real(psb_dpk_),intent(inout), optional, target :: work(:)

  type(psb_d_vect_type)       :: ww
  character     :: trans_
  real(psb_dpk_), pointer :: work_(:)
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: err_act
  character(len=20)   :: name

  name = 'psb_d_apply1v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt = desc_data%get_context()
  call psb_info(ictxt, me, np)

  if (present(trans)) then
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (present(work)) then
    work_ => work
  else
    allocate(work_(4*desc_data%get_local_cols()),stat=info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if

  end if

  if (.not.allocated(prec%prec)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if

  call psb_geasb(ww,desc_data,info,mold=x%v,scratch=.true.)
  if (info == 0) call prec%prec%apply(done,x,dzero,ww,desc_data,info,&
       & trans=trans_,work=work_)
  if (info == 0) call psb_geaxpby(done,ww,dzero,x,desc_data,info)
  call psb_gefree(ww,desc_data,info)
  if (present(work)) then
  else
    deallocate(work_,stat=info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999
    end if
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_d_apply1_vect

subroutine psb_d_apply2v(prec,x,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_d_apply2v
  implicit none
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_dprec_type), intent(inout) :: prec
  real(psb_dpk_),intent(inout)   :: x(:)
  real(psb_dpk_),intent(inout)   :: y(:)
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans
  real(psb_dpk_),intent(inout), optional, target :: work(:)

  character     :: trans_
  real(psb_dpk_), pointer :: work_(:)
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: err_act
  character(len=20)   :: name

  name='psb_d_apply2v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ictxt = desc_data%get_context()
  call psb_info(ictxt, me, np)

  if (present(trans)) then
    trans_=trans
  else
    trans_='N'
  end if

  if (present(work)) then
    work_ => work
  else
    allocate(work_(4*desc_data%get_local_cols()),stat=info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999
    end if

  end if

  if (.not.allocated(prec%prec)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if
  call prec%prec%apply(done,x,dzero,y,desc_data,info,trans_,work=work_)
  if (present(work)) then
  else
    deallocate(work_,stat=info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999
    end if
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_d_apply2v

subroutine psb_d_apply1v(prec,x,desc_data,info,trans)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_d_apply1v
  implicit none
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_dprec_type), intent(inout) :: prec
  real(psb_dpk_),intent(inout)   :: x(:)
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans

  character     :: trans_
  integer(psb_ipk_) :: ictxt,np,me
  integer(psb_ipk_) :: err_act
  real(psb_dpk_), pointer :: WW(:), w1(:)
  character(len=20)   :: name
  name='psb_d_apply1v'
  info = psb_success_
  call psb_erractionsave(err_act)


  ictxt=desc_data%get_context()
  call psb_info(ictxt, me, np)
  if (present(trans)) then
    trans_=psb_toupper(trans)
  else
    trans_='N'
  end if

  if (.not.allocated(prec%prec)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if
  allocate(ww(size(x)),w1(size(x)),stat=info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='Allocate')
    goto 9999
  end if
  call prec%prec%apply(done,x,dzero,ww,desc_data,info,&
       & trans_,work=w1)
  if(info /= psb_success_) goto 9999
  x(:) = ww(:)
  deallocate(ww,W1,stat=info)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='DeAllocate')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_error_handler(err_act)
  return

end subroutine psb_d_apply1v

subroutine psb_dcprecseti(prec,what,val,info,ilev,ilmax,pos,idx)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_dcprecseti
  implicit none

  class(psb_dprec_type), intent(inout)   :: prec
  character(len=*), intent(in)             :: what
  integer(psb_ipk_), intent(in)            :: val
  integer(psb_ipk_), intent(out)           :: info
  ! This optional inputs are backport from the inputs available in AMG4PSBLAS,
  ! they are of no actual use here a part from compatibility reasons.
  integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
  character(len=*), optional, intent(in)   :: pos

  ! Local variables
  character(len=*), parameter            :: name='psb_precseti'

  info = psb_success_

  ! We need to convert from the 'what' string to the corresponding integer
  ! value befor passing the call to the set of the inner method.
  select case (psb_toupper(what))
    case ("SUB_FILLIN")
      call prec%prec%precset(psb_ilu_fill_in_,val,info)
    case('INV_FILLIN')
      call prec%prec%precset(psb_inv_fillin_,val,info)
    case default
      info = psb_err_invalid_args_combination_
      write(psb_err_unit,*) name,&
           & ': Error: uninitialized preconditioner,',&
           &' should call prec%init'
      return
  end select

end subroutine psb_dcprecseti

subroutine psb_dcprecsetr(prec,what,val,info,ilev,ilmax,pos,idx)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_dcprecsetr
  implicit none

  class(psb_dprec_type), intent(inout)   :: prec
  character(len=*), intent(in)             :: what
  real(psb_dpk_), intent(in)             :: val
  integer(psb_ipk_), intent(out)           :: info
  ! This optional inputs are backport from the inputs available in AMG4PSBLAS,
  ! they are of no actual use here a part from compatibility reasons.
  integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
  character(len=*), optional, intent(in)   :: pos

  ! Local variables
  character(len=*), parameter            :: name='amg_precsetr'

  info = psb_success_

  ! We need to convert from the 'what' string to the corresponding integer
  ! value befor passing the call to the set of the inner method.
  select case (psb_toupper(what))
  case('SUB_ILUTHRS')
    call prec%prec%precset(psb_fact_eps_,val,info)
  case('INV_THRESH')
    call prec%prec%precset(psb_inv_thresh_,val,info)
  case default
    info = psb_err_invalid_args_combination_
    write(psb_err_unit,*) name,&
         & ': Error: uninitialized preconditioner,',&
         &' should call prec%init'
    return
  end select

end subroutine psb_dcprecsetr

subroutine psb_dcprecsetc(prec,what,string,info,ilev,ilmax,pos,idx)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_dcprecsetc
  implicit none

  class(psb_dprec_type), intent(inout)   :: prec
  character(len=*), intent(in)             :: what
  character(len=*), intent(in)             :: string
  integer(psb_ipk_), intent(out)           :: info
  ! This optional inputs are backport from the inputs available in AMG4PSBLAS,
  ! they are of no actual use here a part from compatibility reasons.
  integer(psb_ipk_), optional, intent(in)  :: ilev,ilmax,idx
  character(len=*), optional, intent(in)   :: pos

  ! Local variables
  character(len=*), parameter            :: name='amg_precsetc'

  info = psb_success_

  ! We need to convert from the 'what' string to the corresponding integer
  ! value befor passing the call to the set of the inner method.
  select case (psb_toupper(what))
    case ('SUB_SOLVE')
      ! We select here the type of solver on the block
      select case (psb_toupper(string))
        case("ILU")
            call prec%prec%precset(psb_f_type_,psb_f_ilu_k_,info)
            call prec%prec%precset(psb_ilu_ialg_,psb_ilu_n_,info)
        case("ILUT")
            call prec%prec%precset(psb_f_type_,psb_f_ilu_t_,info)
            call prec%prec%precset(psb_ilu_ialg_,psb_ilu_t_,info)
        case("AINV")
            call prec%prec%precset(psb_f_type_,psb_f_ainv_,info)
        case default
          ! Default to ILU(0) factorization
          call prec%prec%precset(psb_f_type_,psb_f_ilu_n_,info)
          call prec%prec%precset(psb_ilu_ialg_,psb_ilu_n_,info)
      end select
    case ("ILU_ALG")
      select case (psb_toupper(string))
        case ("MILU")
          call prec%prec%precset(psb_ilu_ialg_,psb_milu_n_,info)
        case default
          ! Do nothing
      end select
    case ("ILUT_SCALE")
      select case (psb_toupper(string))
      case ("MAXVAL")
        call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_maxval_,info)
      case ("DIAG")
        call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_diag_,info)
      case ("ARWSUM")
        call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_arwsum_,info)
      case ("ARCSUM")
        call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_arcsum_,info)
      case ("ACLSUM")
        call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_aclsum_,info)
      case default
        call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_none_,info)
      end select
    case ("AINV_ALG")
      select case (psb_toupper(string))
      case("LLK")
        call prec%prec%precset(psb_ainv_alg_,psb_ainv_llk_,info)
      case("SYM-LLK")
        call prec%prec%precset(psb_ainv_alg_,psb_ainv_s_llk_,info)
      case("STAB-LLK")
        call prec%prec%precset(psb_ainv_alg_,psb_ainv_s_ft_llk_,info)
      case("MLK","LMX")
        call prec%prec%precset(psb_ainv_alg_,psb_ainv_mlk_,info)
      case default
        call prec%prec%precset(psb_ainv_alg_,psb_ainv_llk_,info)
      end select
    case default

  end select

end subroutine psb_dcprecsetc
