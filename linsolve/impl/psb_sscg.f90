subroutine psb_sscg_vect(a, prec, b, x, s, eps, desc_a, info, &
                      & itmax, iter, err, itrace, istop, &
                      & base_type, eigext, Gram_solver, FGS_sweeps)
  use psb_base_mod
  use psb_prec_mod
  use psb_s_linsolve_conv_mod
  use psb_linsolve_mod
  use psb_pMPK_mod
  implicit none
  type(psb_sspmat_type), intent(in)     :: a
  class(psb_sprec_type), intent(inout)  :: prec
  type(psb_s_vect_type), intent(inout)  :: b, x
  integer(psb_ipk_), intent(in)         :: s
  real(psb_spk_), intent(in)            :: eps
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
  integer(psb_ipk_), optional, intent(out)  :: iter
  real(psb_spk_), optional, intent(out)     :: err
  character, optional, intent(in)           :: base_type
  real(psb_spk_), optional, intent(in)      :: eigext(2)
  character(len=3), optional, intent(in)    :: Gram_solver
  integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps

  ! Local vars
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: istop_, itmax_, itrace_, FGS_sweeps_
  character(len=3)    :: base_type_, Gram_solver_
  integer(psb_ipk_)   :: err_act, np, me, debug_level, debug_unit, n_col, n_row
  integer(psb_lpk_)   :: mglob
  character(len=20)           :: name = 'psb_sscg'
  character(len=*), parameter :: methdbasename = 'sStepCG'
  character(len=20)           :: methdfullname

  real(psb_spk_), allocatable :: alpha(:), beta(:, :), W(:, :), pW(:), temp_fa(:, :)
  type(psb_s_vect_type)       :: r  
  type(psb_s_multivect_type)  :: Z, Q, P, V, temp_mv
  real(psb_spk_)              :: cheb_coeff(3)
  integer(psb_ipk_)           :: itidx
  
  type(psb_itconv_type)         :: stopdat
  real(psb_dpk_)                :: derr

  type(psb_s_multivect_type), target  :: aux_mv
  real(psb_spk_), allocatable, target :: aux_fa(:)

  character(len=3), parameter   :: forwardGS = "FGS"
  character(len=3), parameter   :: lapackLU = "LLU"
  character(len=3), parameter   :: lapackCC = "LCC"

  info = psb_success_
  call psb_erractionsave(err_act)
  
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if(debug_level >= psb_debug_ext_) &
       & write(debug_unit, *) me, ' ', trim(name), ': from psb_info ', np
  
  if(s < 1) then
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info, name)
    goto 9999
  endif

  write(methdfullname, '(A, "(", I0, ")")') methdbasename, s

  if((.not. allocated(b%v)) .or. (.not.allocated(x%v))) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  istop_ = 2
  if(present(istop)) istop_ = istop
  
  !  ISTOP_ = 1:  Normwise backward error, infinity norm 
  !  ISTOP_ = 2:  ||r||/||b||, 2-norm 
  if((istop_ < 1 ) .or. (istop_ > 2 )) then
    info = psb_err_invalid_istop_
    err = info
    call psb_errpush(info, name, i_err = (/istop_/))
    goto 9999
  endif

  itmax_ = 1000
  if(present(itmax)) itmax_ = itmax

  itrace_ = 0
  if(present(itrace)) itrace_ = itrace

  base_type_ = "C"
  if(present(base_type)) base_type_ = base_type

  Gram_solver_ = lapackCC
  if(present(Gram_solver)) Gram_solver_ = Gram_solver
  
  FGS_sweeps_ = 30
  if(present(FGS_sweeps)) FGS_sweeps_ = FGS_sweeps

  mglob = desc_a%get_global_rows()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()

  call psb_chkvect(mglob, lone, x%get_nrows(), lone, lone, desc_a, info)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info, name,  a_err = 'psb_chkvect on x')
    goto 9999
  end if

  call psb_chkvect(mglob, lone, b%get_nrows(), lone, lone, desc_a, info)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_    
    call psb_errpush(info, name, a_err='psb_chkvect on b')
    goto 9999
  end if

  !Allocate and assembly data structure
  allocate(alpha(s), beta(s, s), W(s, s), pW(s), temp_fa(s, s + 1), aux_fa(4*n_col), stat = info)
  if(info == psb_success_) call psb_geall(r, desc_a, info)
  if(info == psb_success_) call psb_geall(Z, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(Q, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(P, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(V, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(temp_mv, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(aux_mv, desc_a, info, n = 3)
  if(info == psb_success_) call psb_geasb(r, desc_a, info)
  if(info == psb_success_) call psb_geasb(Z, desc_a, info)
  if(info == psb_success_) call psb_geasb(Q, desc_a, info)
  if(info == psb_success_) call psb_geasb(P, desc_a, info)
  if(info == psb_success_) call psb_geasb(V, desc_a, info)
  if(info == psb_success_) call psb_geasb(temp_mv, desc_a, info)
  if(info == psb_success_) call psb_geasb(aux_mv, desc_a, info)

  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! First residual calculation
  call psb_geaxpby(sone, b, szero, r, desc_a, info)
  if(info == psb_success_) call psb_spmm(-sone, a, x, sone, r, desc_a, info)
  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! Init converence
  call psb_init_conv(methdfullname, istop_, itrace_, itmax_, a, x, b, eps, desc_a, stopdat, info)
  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! check convergence here?

  ! Chebyshev coefficient calculation
  select case (base_type_)
    case("M")
      cheb_coeff = szero
    case("C")
      cheb_coeff = psb_s_chebyshev_coefficients(a, prec, desc_a, info, eigext)
    case default
      info = psb_err_invalid_input_ 
      call psb_errpush(info, name)
      goto 9999
  end select

  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  
  ! First matrix power kernel
  call psb_pMPK(a, prec, r, P, V, s, desc_a, info, base_type = base_type_, &
                  & alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3), &
                  & mvec_temp = aux_mv, farr_temp = aux_fa)
  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! Loop until convergence (or maxiter)
  do itidx = 1, itmax_
    ! Compute matrix W and rhs for alpha
    call psb_gedots(P, V, temp_fa(:, 1 : s), desc_a, info, global = .false.)
    call psb_gedots(P, r, temp_fa(:, s + 1), desc_a, info, global = .false.)
    call psb_sum(desc_a%get_context(), temp_fa)
    W = temp_fa(:, 1 : s)
    alpha = temp_fa(:, s + 1)

    ! Factor matrix W (if soving with LU or Cholesky factorization)
    if(Gram_solver_ == lapackLU) call sgetrf(s, s, W, s, pW, info)
    if(Gram_solver_ == lapackCC) call spotrf('L', s, W, s, info)

    ! Solve for alpha
    select case(Gram_solver_)
      case(forwardGS);  call inner_solver_fgs_1D(W, alpha, FGS_sweeps_)
      case(lapackLU);   call sgetrs('N', s, 1, W, s, pW, alpha, s, info)
      case(lapackCC);   call spotrs('L', s, 1, W, s, alpha, s, info)
      case default
        info = psb_err_invalid_input_ 
        call psb_errpush(info, name)
        goto 9999
    end select

    ! Update solution and residual
    call psb_geaxpby(P, alpha, x, desc_a, info, upd_flag = .true.)
    call psb_geaxpby(V, -alpha, r, desc_a, info, upd_flag = .true.)

    ! Check convergence. 
    if(psb_check_conv(methdfullname, itidx, x, r, desc_a, stopdat, info)) exit
    
    ! Matrix power kernel
    call psb_pMPK(a, prec, r, Z, Q, s, desc_a, info, base_type = base_type_, &
                  & alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3), &
                  & mvec_temp = aux_mv, farr_temp = aux_fa)

    ! Compute rhs for beta
    call psb_gedots(P, Q, beta, desc_a, info, .true.)
    beta = -beta;

    ! Solve for beta
    select case(Gram_solver_)
      case(forwardGS);  call inner_solver_fgs_2D(W, beta, FGS_sweeps_)
      case(lapackLU);   call sgetrs('N', s, s, W, s, pW, beta, s, info)
      case(lapackCC);   call spotrs('L', s, s, W, s, beta, s, info)
      case default
        info = psb_err_invalid_input_ 
        call psb_errpush(info, name)
        goto 9999
    end select

    ! Update P and V. Use of temp_mv in needed because internal dgemm constraint
    call psb_geaxpby(P, beta, temp_mv, desc_a, info, .false.)
    call psb_geaxpby(sone, Z, sone, temp_mv, P, desc_a, info)
    call psb_geaxpby(V, beta, temp_mv, desc_a, info, .false.)
    call psb_geaxpby(sone, Q, sone, temp_mv, V, desc_a, info)
  end do

  call psb_end_conv(methdfullname, itidx, desc_a, stopdat, info, derr, iter)
  if(present(err)) err = derr
  if(present(iter)) iter = iter * s

  if(info == psb_success_) call psb_gefree(r, desc_a, info)
  if(info == psb_success_) call psb_gefree(Z, desc_a, info)
  if(info == psb_success_) call psb_gefree(Q, desc_a, info)
  if(info == psb_success_) call psb_gefree(P, desc_a, info)
  if(info == psb_success_) call psb_gefree(V, desc_a, info)
  if(info == psb_success_) call psb_gefree(temp_mv, desc_a, info)
  if(info == psb_success_) call psb_gefree(aux_mv, desc_a, info)

  if(info == psb_success_) deallocate(alpha, beta, W, pW, temp_fa, aux_fa, stat = info)
  if(info /= psb_success_) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains
  function psb_s_chebyshev_coefficients(a, prec, desc, info, eigext) result(coeff)
    use psb_eigsolve_mod
    type(psb_sspmat_type), intent(in)     :: a
    class(psb_sprec_type), intent(inout)  :: prec
    type(psb_desc_type), intent(in)       :: desc
    integer(psb_ipk_), intent(out)        :: info
    real(psb_spk_), optional, intent(in)  :: eigext(2)
    
    real(psb_spk_) :: lambda_max, lambda_min
    real(psb_spk_) :: coeff(3)
 
    info = psb_success_

    if(present(eigext)) then
      lambda_min = eigext(1)
      lambda_max = eigext(2)
    else 
      call psb_powermethod(a, prec, lambda_max, desc, info)
      lambda_min = szero
    end if

    if(info /= psb_success_) then 
      info = psb_err_from_subroutine_
      return
    end if

    if(lambda_min < 0) then 
      info = psb_err_fatal_ ! TODO: set the correct error
      return
    end if

    if(lambda_max < lambda_min) then 
      info = psb_err_fatal_ ! TODO: set the correct error
      return
    end if

    coeff(1) = 2_psb_spk_ / (lambda_max - lambda_min)
    coeff(2) = (lambda_max + lambda_min) / (lambda_max - lambda_min)
    coeff(3) = sone
  end function psb_s_chebyshev_coefficients

  subroutine inner_solver_fgs_1D(M, rhs, num_iter)
    real(psb_spk_), intent(in)    :: M(:, :)
    real(psb_spk_), intent(inout) :: rhs(:)
    integer(psb_ipk_), intent(in) :: num_iter

    integer(psb_ipk_) :: iter_idx, i, j, n
    real(psb_spk_)    :: sol(size(rhs))

    sol = szero
    n = size(sol)
    do iter_idx = 1, num_iter
      do i = 1, n
        sol(i) = rhs(i)
        do j = 1, n
          if(j /= i) sol(i) = sol(i) - M(i, j) * sol(j)
        end do
        sol(i) = sol(i) / M(i, i)
      end do
    end do
    rhs = sol
  end subroutine inner_solver_fgs_1D

  subroutine inner_solver_fgs_2D(M, rhs, num_iter)
    real(psb_spk_), intent(in)    :: M(:, :)
    real(psb_spk_), intent(inout) :: rhs(:, :)
    integer(psb_ipk_), intent(in) :: num_iter

    integer(psb_ipk_) :: iter_idx, i, j, n
    real(psb_spk_)    :: sol(size(rhs, 1), size(rhs, 2))

    sol = szero
    n = size(sol, 1)
    do iter_idx = 1, num_iter
      do i = 1, n
        sol(i, :) = rhs(i, :)
        do j = 1, n
          if(j /= i) sol(i, :) = sol(i, :) - M(i, j) * sol(j, :)
        end do
        sol(i, :) = sol(i, :) / M(i, i)
      end do
    end do
    rhs = sol
  end subroutine inner_solver_fgs_2D
end subroutine psb_sscg_vect

subroutine psb_sscg2_vect(a, prec, b, x, s, eps, desc_a, info, &
                        & itmax, iter, err, itrace, istop, &
                        & base_type, eigext, Gram_solver, FGS_sweeps)
  use psb_base_mod
  use psb_prec_mod
  use psb_s_linsolve_conv_mod
  use psb_linsolve_mod
  use psb_pMPK_mod

  implicit none
  type(psb_sspmat_type), intent(in)     :: a
  class(psb_sprec_type), intent(inout)  :: prec
  type(psb_s_vect_type), intent(inout)  :: b, x
  integer(psb_ipk_), intent(in)         :: s
  real(psb_spk_), intent(in)            :: eps
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
  integer(psb_ipk_), optional, intent(out)  :: iter
  real(psb_spk_), optional, intent(out)     :: err
  character, optional, intent(in)           :: base_type
  real(psb_spk_), optional, intent(in)      :: eigext(2)
  character(len=3), optional, intent(in)    :: Gram_solver
  integer(psb_ipk_), optional, intent(in)   :: FGS_sweeps

  ! Local vars
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: istop_, itmax_, itrace_, FGS_sweeps_
  character(len=3)    :: base_type_, Gram_solver_
  integer(psb_ipk_)   :: err_act, np, me, debug_level, debug_unit, &
                          & n_col, n_row
  integer(psb_lpk_)   :: mglob
  character(len=20)           :: name = 'psb_sscg'
  character(len=*), parameter :: methdbasename = 'sStepCGv2'
  character(len=20)           :: methdfullname

  real(psb_spk_), allocatable :: alpha(:), beta(:, :), W(:, :), pW(:), temp_fa(:, :), B2(:, :), c0(:)
  type(psb_s_vect_type)       :: r  
  type(psb_s_multivect_type)  :: Z, Q, P, V, temp_mv
  real(psb_spk_)              :: cheb_coeff(3)
  integer(psb_ipk_)           :: itidx
  
  type(psb_itconv_type)         :: stopdat
  real(psb_dpk_)                :: derr 

  type(psb_s_multivect_type), target  :: aux_mv
  real(psb_spk_), allocatable, target :: aux_fa(:)

  character(len=3), parameter   :: forwardGS = "FGS"
  character(len=3), parameter   :: lapackLU = "LLU"
  character(len=3), parameter   :: lapackCC = "LCC"

  info = psb_success_
  call psb_erractionsave(err_act)
  
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if(debug_level >= psb_debug_ext_) &
       & write(debug_unit, *) me, ' ', trim(name), ': from psb_info ', np
  
  if(s < 1) then
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info, name)
    goto 9999
  endif

  write(methdfullname, '(A, "(", I0, ")")') methdbasename, s

  if((.not. allocated(b%v)) .or. (.not.allocated(x%v))) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  istop_ = 2
  if(present(istop)) istop_ = istop
  
  !  ISTOP_ = 1:  Normwise backward error, infinity norm 
  !  ISTOP_ = 2:  ||r||/||b||, 2-norm 
  if((istop_ < 1 ) .or. (istop_ > 2 )) then
    info = psb_err_invalid_istop_
    err = info
    call psb_errpush(info, name, i_err = (/istop_/))
    goto 9999
  endif

  itmax_ = 1000
  if(present(itmax)) itmax_ = itmax

  itrace_ = 0
  if(present(itrace)) itrace_ = itrace

  base_type_ = "C"
  if(present(base_type)) base_type_ = base_type

  Gram_solver_ = lapackCC
  if(present(Gram_solver)) Gram_solver_ = Gram_solver
  
  FGS_sweeps_ = 30
  if(present(FGS_sweeps)) FGS_sweeps_ = FGS_sweeps

  mglob = desc_a%get_global_rows()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()

  call psb_chkvect(mglob, lone, x%get_nrows(), lone, lone, desc_a, info)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info, name,  a_err = 'psb_chkvect on x')
    goto 9999
  end if

  call psb_chkvect(mglob, lone, b%get_nrows(), lone, lone, desc_a, info)
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_    
    call psb_errpush(info, name, a_err='psb_chkvect on b')
    goto 9999
  end if

  !Allocate and assembly data structure
  allocate(alpha(s), beta(s, s), W(s, s), pW(s), temp_fa(s, 2*s + 1), B2(s, s), c0(s), aux_fa(4*n_col), stat = info)
  if(info == psb_success_) call psb_geall(r, desc_a, info)
  if(info == psb_success_) call psb_geall(Z, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(Q, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(P, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(V, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(temp_mv, desc_a, info, n = s)
  if(info == psb_success_) call psb_geall(aux_mv, desc_a, info, n = 3)
  if(info == psb_success_) call psb_geasb(r, desc_a, info)
  if(info == psb_success_) call psb_geasb(Z, desc_a, info)
  if(info == psb_success_) call psb_geasb(Q, desc_a, info)
  if(info == psb_success_) call psb_geasb(P, desc_a, info)
  if(info == psb_success_) call psb_geasb(V, desc_a, info)
  if(info == psb_success_) call psb_geasb(temp_mv, desc_a, info)
  if(info == psb_success_) call psb_geasb(aux_mv, desc_a, info)

  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! First residual calculation
  call psb_geaxpby(sone, b, szero, r, desc_a, info)
  if(info == psb_success_) call psb_spmm(-sone, a, x, sone, r, desc_a, info)
  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! Init converence
  call psb_init_conv(methdfullname, istop_, itrace_, itmax_, a, x, b, eps, desc_a, stopdat, info)
  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! check convergence here?

  ! Chebyshev coefficient calculation
  select case (base_type_)
    case("M")
      cheb_coeff = szero
    case("C")
      cheb_coeff = psb_s_chebyshev_coefficients(a, prec, desc_a, info, eigext)
    case default
      info = psb_err_invalid_input_ 
      call psb_errpush(info, name)
      goto 9999
  end select

  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  
  ! First matrix power kernel
  call psb_pMPK(a, prec, r, P, V, s, desc_a, info, base_type = base_type_, &
                  & alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3), &
                  & mvec_temp = aux_mv, farr_temp = aux_fa)
  if(info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! Compute first Gram system components
  call psb_gedots(P, V, temp_fa(:, 1 : s), desc_a, info, global = .false.)
  call psb_gedots(P, r, temp_fa(:, 2*s + 1), desc_a, info, global = .false.)
  call psb_sum(desc_a%get_context(), temp_fa)

  W = temp_fa(:, 1 : s)
  alpha = temp_fa(:, 2*s + 1)

  ! Loop until convergence (or maxiter)
  do itidx = 1, itmax_
    ! Factor matrix W (if soving with LU or Cholesky factorization)
    if(Gram_solver_ == lapackLU) call sgetrf(s, s, W, s, pW, info)
    if(Gram_solver_ == lapackCC) call spotrf('L', s, W, s, info)

    ! Solve for alpha
    select case(Gram_solver_)
      case(forwardGS);  call inner_solver_fgs_1D(W, alpha, FGS_sweeps_)
      case(lapackLU);   call sgetrs('N', s, 1, W, s, pW, alpha, s, info)
      case(lapackCC);   call spotrs('L', s, 1, W, s, alpha, s, info)
      case default
        info = psb_err_invalid_input_ 
        call psb_errpush(info, name)
        goto 9999
    end select

    ! Update solution and residual
    call psb_geaxpby(P, alpha, x, desc_a, info, upd_flag = .true.)
    call psb_geaxpby(V, -alpha, r, desc_a, info, upd_flag = .true.)

    ! Check convergence
    if(psb_check_conv(methdfullname, itidx, x, r, desc_a, stopdat, info)) exit

    ! Matrix power kernel
    call psb_pMPK(a, prec, r, Z, Q, s, desc_a, info, base_type = base_type_, &
                  & alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3), &
                  & mvec_temp = aux_mv, farr_temp = aux_fa)

    ! Compute dot products
    call psb_gedots(P, Q, temp_fa(:, 1 : s), desc_a, info, global = .false.)
    call psb_gedots(Z, Q, temp_fa(:, s+1 : 2*s), desc_a, info, global = .false.)
    call psb_gedots(Z, r, temp_fa(:, 2*s + 1), desc_a, info, global = .false.)
    call psb_sum(desc_a%get_context(), temp_fa)

    !Compute new rhs for beta
    B2 = -temp_fa(:, 1 : s)

    ! Solve for beta
    beta = B2
    select case(Gram_solver_)
      case(forwardGS);  call inner_solver_fgs_2D(W, beta, FGS_sweeps_)
      case(lapackLU);   call sgetrs('N', s, s, W, s, pW, beta, s, info)
      case(lapackCC);   call spotrs('L', s, s, W, s, beta, s, info)
      case default
        info = psb_err_invalid_input_ 
        call psb_errpush(info, name)
        goto 9999
    end select

    ! Update P and V. Use of temp_mv in needed because internal dgemm constraint
    call psb_geaxpby(P, beta, temp_mv, desc_a, info, upd_flag = .false.)
    call psb_geaxpby(sone, Z, sone, temp_mv, P, desc_a, info)
    call psb_geaxpby(V, beta, temp_mv, desc_a, info, upd_flag = .false.)
    call psb_geaxpby(sone, Q, sone, temp_mv, V, desc_a, info)

    !Compute new Gram matrix
    W = temp_fa(:, s+1 : 2*s)
    call sgemm('T', 'N', s, s, s, -sone, beta, s, B2, s, sone, W, s)
    
    !Compute new rhs for alpha
    alpha = temp_fa(:, 2*s + 1)
  end do

  call psb_end_conv(methdfullname, itidx, desc_a, stopdat, info, derr, iter)
  if(present(err)) err = derr
  if(present(iter)) iter = iter * s

  if(info == psb_success_) call psb_gefree(r, desc_a, info)
  if(info == psb_success_) call psb_gefree(Z, desc_a, info)
  if(info == psb_success_) call psb_gefree(Q, desc_a, info)
  if(info == psb_success_) call psb_gefree(P, desc_a, info)
  if(info == psb_success_) call psb_gefree(V, desc_a, info)
  if(info == psb_success_) call psb_gefree(temp_mv, desc_a, info)
  if(info == psb_success_) call psb_gefree(aux_mv, desc_a, info)

  if(info == psb_success_) deallocate(alpha, beta, W, pW, temp_fa, B2, c0, aux_fa, stat = info)
  if(info /= psb_success_) then
    call psb_errpush(info, name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains 
  function psb_s_chebyshev_coefficients(a, prec, desc, info, eigext) result(coeff)
    use psb_eigsolve_mod
    type(psb_sspmat_type), intent(in)     :: a
    class(psb_sprec_type), intent(inout)  :: prec
    type(psb_desc_type), intent(in)       :: desc
    integer(psb_ipk_), intent(out)        :: info
    real(psb_spk_), optional, intent(in)  :: eigext(2)
    
    real(psb_spk_) :: lambda_max, lambda_min
    real(psb_spk_) :: coeff(3)
 
    info = psb_success_

    if(present(eigext)) then
      lambda_min = eigext(1)
      lambda_max = eigext(2)
    else 
      call psb_powermethod(a, prec, lambda_max, desc, info)
      lambda_min = szero
    end if

    if(info /= psb_success_) then 
      info = psb_err_from_subroutine_
      return
    end if

    if(lambda_min < 0) then 
      info = psb_err_fatal_ ! TODO: set the correct error
      return
    end if

    if(lambda_max < lambda_min) then 
      info = psb_err_fatal_ ! TODO: set the correct error
      return
    end if

    coeff(1) = 2_psb_spk_ / (lambda_max - lambda_min)
    coeff(2) = (lambda_max + lambda_min) / (lambda_max - lambda_min)
    coeff(3) = sone
  end function psb_s_chebyshev_coefficients

  subroutine inner_solver_fgs_1D(M, rhs, num_iter)
    real(psb_spk_), intent(in)    :: M(:, :)
    real(psb_spk_), intent(inout) :: rhs(:)
    integer(psb_ipk_), intent(in) :: num_iter

    integer(psb_ipk_) :: iter_idx, i, j, n
    real(psb_spk_)    :: sol(size(rhs))

    sol = szero
    n = size(sol)
    do iter_idx = 1, num_iter
      do i = 1, n
        sol(i) = rhs(i)
        do j = 1, n
          if(j /= i) sol(i) = sol(i) - M(i, j) * sol(j)
        end do
        sol(i) = sol(i) / M(i, i)
      end do
    end do
    rhs = sol
  end subroutine inner_solver_fgs_1D

  subroutine inner_solver_fgs_2D(M, rhs, num_iter)
    real(psb_spk_), intent(in)    :: M(:, :)
    real(psb_spk_), intent(inout) :: rhs(:, :)
    integer(psb_ipk_), intent(in) :: num_iter

    integer(psb_ipk_) :: iter_idx, i, j, n
    real(psb_spk_)    :: sol(size(rhs, 1), size(rhs, 2))

    sol = szero
    n = size(sol, 1)
    do iter_idx = 1, num_iter
      do i = 1, n
        sol(i, :) = rhs(i, :)
        do j = 1, n
          if(j /= i) sol(i, :) = sol(i, :) - M(i, j) * sol(j, :)
        end do
        sol(i, :) = sol(i, :) / M(i, i)
      end do
    end do
    rhs = sol
  end subroutine inner_solver_fgs_2D
end subroutine psb_sscg2_vect