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
!         software without specific prior written permission.
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
!!$       software without specific prior written permission.
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

subroutine psb_d_apply2_vect(prec,x,y,desc_data,info,trans)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_d_apply2_vect
  implicit none
  type(psb_desc_type),intent(in)          :: desc_data
  class(psb_dprec_type), intent(inout)  :: prec
  type(psb_d_vect_type),intent(inout)   :: x
  type(psb_d_vect_type),intent(inout)   :: y
  integer(psb_ipk_), intent(out)          :: info
  character(len=1), optional              :: trans

  character           :: trans_
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np,me
  integer(psb_ipk_)   :: err_act
  character(len=20)   :: name

  name = 'psb_d_apply2v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)

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

  call prec%prec%apply(done,x,dzero,y,desc_data,info,trans=trans_)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_d_apply2_vect

subroutine psb_d_apply1_vect(prec,x,desc_data,info,trans)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_d_apply1_vect
  implicit none
  type(psb_desc_type),intent(in)          :: desc_data
  class(psb_dprec_type), intent(inout)  :: prec
  type(psb_d_vect_type),intent(inout)   :: x
  integer(psb_ipk_), intent(out)          :: info
  character(len=1), optional              :: trans

  type(psb_d_vect_type) :: ww
  character               :: trans_
  type(psb_ctxt_type)     :: ctxt
  integer(psb_ipk_)       :: np,me
  integer(psb_ipk_)       :: err_act
  character(len=20)       :: name

  name = 'psb_d_apply1v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)

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

  call psb_geasb(ww,desc_data,info,mold=x%v,scratch=.true.)
  if (info == 0) call prec%prec%apply(done,x,dzero,ww,desc_data,info,&
       & trans=trans_)
  if (info == 0) call psb_geaxpby(done,ww,dzero,x,desc_data,info)
  call psb_gefree(ww,desc_data,info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_d_apply1_vect

subroutine psb_d_apply2v(prec,x,y,desc_data,info,trans, work)
  use psb_base_mod
  use psb_d_prec_type, psb_protect_name => psb_d_apply2v
  implicit none
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_dprec_type), intent(inout) :: prec
  real(psb_dpk_),intent(inout)   :: x(:)
  real(psb_dpk_),intent(inout)   :: y(:)
  integer(psb_ipk_), intent(out)              :: info
  real(psb_dpk_),intent(inout), optional, target :: work(:)
  character(len=1), optional        :: trans

  character     :: trans_
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: err_act
  character(len=20)   :: name

  name='psb_d_apply2v'
  info = psb_success_
  call psb_erractionsave(err_act)

  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)

  if (present(trans)) then
    trans_=trans
  else
    trans_='N'
  end if

  if (.not.allocated(prec%prec)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if
  call prec%prec%apply(done,x,dzero,y,desc_data,info,trans_, work)

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
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: err_act
  real(psb_dpk_), pointer :: WW(:), w1(:)
  character(len=20)   :: name
  name='psb_d_apply1v'
  info = psb_success_
  call psb_erractionsave(err_act)
  ctxt = desc_data%get_context()
  call psb_info(ctxt, me, np)
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
       & trans_)
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
  use psb_d_nestedprec, only : psb_d_nested_prec_type, &
       & psb_d_nested_schur_maxit_, psb_d_nested_inner_maxit_, &
       & psb_d_nested_inner_itrace_, psb_d_nested_inner_istop_, &
       & psb_d_nested_set_field_precseti, psb_d_nested_set_field_krylov_i
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
  integer(psb_ipk_) :: field

  info = psb_success_

  ! We need to convert from the 'what' string to the corresponding integer
  ! value befor passing the call to the set of the inner method.
  if (.not. allocated(prec%prec)) then
    info = psb_err_invalid_preca_
    return
  end if

  select case (psb_toupper(trim(what)))
    case ('SCHUR_MAXIT','NEST_SCHUR_MAXIT')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        call p%precset(psb_d_nested_schur_maxit_, val, info)
      class default
        info = psb_err_invalid_preca_
      end select
    case ('INNER_MAXIT','INNER_ITMAX','KRYLOV_MAXIT','KRYLOV_ITMAX')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        field = 0
        if (present(idx)) field = idx
        call psb_d_nested_set_field_krylov_i(p, field, psb_d_nested_inner_maxit_, val, info)
      class default
        info = psb_err_invalid_preca_
      end select
    case ('INNER_ITRACE','KRYLOV_ITRACE')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        field = 0
        if (present(idx)) field = idx
        call psb_d_nested_set_field_krylov_i(p, field, psb_d_nested_inner_itrace_, val, info)
      class default
        info = psb_err_invalid_preca_
      end select
    case ('INNER_ISTOP','KRYLOV_ISTOP')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        field = 0
        if (present(idx)) field = idx
        call psb_d_nested_set_field_krylov_i(p, field, psb_d_nested_inner_istop_, val, info)
      class default
        info = psb_err_invalid_preca_
      end select
    case ('SUB_FILLIN')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        field = 0
        if (present(idx)) field = idx
        call psb_d_nested_set_field_precseti(p, field, psb_ilu_fill_in_, val, info)
      class default
        call prec%prec%precset(psb_ilu_fill_in_,val,info)
      end select
    case ('INV_FILLIN')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        field = 0
        if (present(idx)) field = idx
        call psb_d_nested_set_field_precseti(p, field, psb_inv_fillin_, val, info)
      class default
        call prec%prec%precset(psb_inv_fillin_,val,info)
      end select
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
  use psb_d_nestedprec, only : psb_d_nested_prec_type, &
       & psb_d_nested_schur_tol_, psb_d_nested_inner_tol_, &
       & psb_d_nested_set_field_precsetr, psb_d_nested_set_field_krylov_r
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
  character(len=*), parameter            :: name='psb_precsetr'
  integer(psb_ipk_) :: field

  info = psb_success_

  ! We need to convert from the 'what' string to the corresponding integer
  ! value befor passing the call to the set of the inner method.
  if (.not. allocated(prec%prec)) then
    info = psb_err_invalid_preca_
    return
  end if

  select case (psb_toupper(trim(what)))
  case('SCHUR_TOL','NEST_SCHUR_TOL')
    select type (p => prec%prec)
    type is (psb_d_nested_prec_type)
      call p%precset(psb_d_nested_schur_tol_, val, info)
    class default
      info = psb_err_invalid_preca_
    end select
  case('INNER_TOL','KRYLOV_TOL')
    select type (p => prec%prec)
    type is (psb_d_nested_prec_type)
      field = 0
      if (present(idx)) field = idx
      call psb_d_nested_set_field_krylov_r(p, field, psb_d_nested_inner_tol_, val, info)
    class default
      info = psb_err_invalid_preca_
    end select
  case('SUB_ILUTHRS')
    select type (p => prec%prec)
    type is (psb_d_nested_prec_type)
      field = 0
      if (present(idx)) field = idx
      call psb_d_nested_set_field_precsetr(p, field, psb_fact_eps_, val, info)
    class default
      call prec%prec%precset(psb_fact_eps_,val,info)
    end select
  case('INV_THRESH')
    select type (p => prec%prec)
    type is (psb_d_nested_prec_type)
      field = 0
      if (present(idx)) field = idx
      call psb_d_nested_set_field_precsetr(p, field, psb_inv_thresh_, val, info)
    class default
      call prec%prec%precset(psb_inv_thresh_,val,info)
    end select
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
  use psb_d_nestedprec, only : psb_d_nested_prec_type, &
       & psb_d_nested_composition_, psb_d_nested_block_solve_, &
       & psb_d_nested_schur_solve_, psb_d_nested_set_block_solve_field, &
       & psb_d_nested_inner_solve_, psb_d_nested_set_field_precseti, &
       & psb_d_nested_set_field_krylov_c
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
  character(len=*), parameter            :: name='psb_precsetc'
  integer(psb_ipk_) :: field

  info = psb_success_

  ! We need to convert from the 'what' string to the corresponding integer
  ! value befor passing the call to the set of the inner method.
  if (.not. allocated(prec%prec)) then
    info = psb_err_invalid_preca_
    return
  end if

  select case (psb_toupper(trim(what)))
    case ('COMPOSITION','NEST_COMPOSITION')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        call p%precset(psb_d_nested_composition_, string, info)
      class default
        info = psb_err_invalid_preca_
      end select
    case ('SCHUR_SOLVE','NEST_SCHUR_SOLVE')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        call p%precset(psb_d_nested_schur_solve_, string, info)
      class default
        info = psb_err_invalid_preca_
      end select
    case ('BLOCK_SOLVE','NEST_BLOCK_SOLVE')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        if (present(idx)) then
          call psb_d_nested_set_block_solve_field(p, idx, string, info)
        else
          call p%precset(psb_d_nested_block_solve_, string, info)
        end if
      class default
        info = psb_err_invalid_preca_
      end select
    case ('INNER_SOLVE','KRYLOV_SOLVE','FIELD_SOLVE')
      select type (p => prec%prec)
      type is (psb_d_nested_prec_type)
        field = 0
        if (present(idx)) field = idx
        call psb_d_nested_set_field_krylov_c(p, field, psb_d_nested_inner_solve_, string, info)
      class default
        info = psb_err_invalid_preca_
      end select
    case ('SUB_SOLVE')
      ! We select here the type of solver on the block
      field = 0
      if (present(idx)) field = idx
      select case (psb_toupper(trim(string)))
        case("ILU")
          select type (p => prec%prec)
          type is (psb_d_nested_prec_type)
            call psb_d_nested_set_field_precseti(p, field, psb_f_type_, psb_f_ilu_k_, info)
            if (info == psb_success_) &
                 & call psb_d_nested_set_field_precseti(p, field, psb_ilu_ialg_, psb_ilu_n_, info)
          class default
            call prec%prec%precset(psb_f_type_,psb_f_ilu_k_,info)
            if (info == psb_success_) call prec%prec%precset(psb_ilu_ialg_,psb_ilu_n_,info)
          end select
        case("ILUT")
          select type (p => prec%prec)
          type is (psb_d_nested_prec_type)
            call psb_d_nested_set_field_precseti(p, field, psb_f_type_, psb_f_ilu_t_, info)
            if (info == psb_success_) &
                 & call psb_d_nested_set_field_precseti(p, field, psb_ilu_ialg_, psb_ilu_t_, info)
          class default
            call prec%prec%precset(psb_f_type_,psb_f_ilu_t_,info)
            if (info == psb_success_) call prec%prec%precset(psb_ilu_ialg_,psb_ilu_t_,info)
          end select
        case("AINV")
          select type (p => prec%prec)
          type is (psb_d_nested_prec_type)
            call psb_d_nested_set_field_precseti(p, field, psb_f_type_, psb_f_ainv_, info)
          class default
            call prec%prec%precset(psb_f_type_,psb_f_ainv_,info)
          end select
        case("INVK")
          select type (p => prec%prec)
          type is (psb_d_nested_prec_type)
            call psb_d_nested_set_field_precseti(p, field, psb_f_type_, psb_f_invk_, info)
          class default
            call prec%prec%precset(psb_f_type_,psb_f_invk_,info)
          end select
        case("INVT")
          select type (p => prec%prec)
          type is (psb_d_nested_prec_type)
            call psb_d_nested_set_field_precseti(p, field, psb_f_type_, psb_f_invt_, info)
          class default
            call prec%prec%precset(psb_f_type_,psb_f_invt_,info)
          end select
        case default
          ! Default to ILU(0) factorization
          select type (p => prec%prec)
          type is (psb_d_nested_prec_type)
            call psb_d_nested_set_field_precseti(p, field, psb_f_type_, psb_f_ilu_n_, info)
            if (info == psb_success_) &
                 & call psb_d_nested_set_field_precseti(p, field, psb_ilu_ialg_, psb_ilu_n_, info)
          class default
            call prec%prec%precset(psb_f_type_,psb_f_ilu_n_,info)
            if (info == psb_success_) call prec%prec%precset(psb_ilu_ialg_,psb_ilu_n_,info)
          end select
      end select
    case ("ILU_ALG")
      field = 0
      if (present(idx)) field = idx
      select case (psb_toupper(trim(string)))
        case ("MILU")
          select type (p => prec%prec)
          type is (psb_d_nested_prec_type)
            call psb_d_nested_set_field_precseti(p, field, psb_ilu_ialg_, psb_milu_n_, info)
          class default
            call prec%prec%precset(psb_ilu_ialg_,psb_milu_n_,info)
          end select
        case default
          ! Do nothing
      end select
    case ("ILUT_SCALE")
      field = 0
      if (present(idx)) field = idx
      select case (psb_toupper(trim(string)))
      case ("MAXVAL")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ilu_scale_, psb_ilu_scale_maxval_, info)
        class default
          call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_maxval_,info)
        end select
      case ("DIAG")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ilu_scale_, psb_ilu_scale_diag_, info)
        class default
          call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_diag_,info)
        end select
      case ("ARWSUM")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ilu_scale_, psb_ilu_scale_arwsum_, info)
        class default
          call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_arwsum_,info)
        end select
      case ("ARCSUM")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ilu_scale_, psb_ilu_scale_arcsum_, info)
        class default
          call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_arcsum_,info)
        end select
      case ("ACLSUM")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ilu_scale_, psb_ilu_scale_aclsum_, info)
        class default
          call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_aclsum_,info)
        end select
      case ("NONE")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ilu_scale_, psb_ilu_scale_none_, info)
        class default
          call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_none_,info)
        end select
      case default
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ilu_scale_, psb_ilu_scale_none_, info)
        class default
          call prec%prec%precset(psb_ilu_scale_,psb_ilu_scale_none_,info)
        end select
      end select
    case ("AINV_ALG")
      field = 0
      if (present(idx)) field = idx
      select case (psb_toupper(trim(string)))
      case("LLK")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ainv_alg_, psb_ainv_llk_, info)
        class default
          call prec%prec%precset(psb_ainv_alg_,psb_ainv_llk_,info)
        end select
      case("SYM-LLK")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ainv_alg_, psb_ainv_s_llk_, info)
        class default
          call prec%prec%precset(psb_ainv_alg_,psb_ainv_s_llk_,info)
        end select
      case("STAB-LLK")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ainv_alg_, psb_ainv_s_ft_llk_, info)
        class default
          call prec%prec%precset(psb_ainv_alg_,psb_ainv_s_ft_llk_,info)
        end select
      case("MLK","LMX")
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ainv_alg_, psb_ainv_mlk_, info)
        class default
          call prec%prec%precset(psb_ainv_alg_,psb_ainv_mlk_,info)
        end select
      case default
        select type (p => prec%prec)
        type is (psb_d_nested_prec_type)
          call psb_d_nested_set_field_precseti(p, field, psb_ainv_alg_, psb_ainv_llk_, info)
        class default
          call prec%prec%precset(psb_ainv_alg_,psb_ainv_llk_,info)
        end select
      end select
    case default

  end select

end subroutine psb_dcprecsetc
