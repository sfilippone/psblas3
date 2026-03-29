subroutine psb_dscg_vect(a, prec, b, x, s, eps, base_type, desc_a, info, itmax, iter, err, itrace, istop, eigext)
  use psb_base_mod
  use psb_prec_mod
  use psb_d_linsolve_conv_mod
  use psb_linsolve_mod
  use psb_pMPK_mod

  implicit none
  type(psb_dspmat_type), intent(in)     :: a
  class(psb_dprec_type), intent(inout)  :: prec
  type(psb_d_vect_type), intent(inout)  :: b, x
  integer(psb_ipk_), intent(in)         :: s
  real(psb_dpk_), intent(in)            :: eps
  character, intent(in)                 :: base_type
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
  integer(psb_ipk_), optional, intent(out)  :: iter
  real(psb_dpk_), optional, intent(out)     :: err
  real(psb_dpk_), optional, intent(in)      :: eigext(2)

  ! Local vars
  integer(psb_ipk_)   :: istop_, itmax_, itrace_
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: err_act, np, me, debug_level, debug_unit, &
                          & n_col, n_row
  integer(psb_lpk_)   :: mglob
  character(len=20)           :: name = 'psb_dscg'
  character(len=*), parameter :: methdname = 'sStepCG'

  real(psb_dpk_), allocatable :: alpha(:), beta(:, :), W(:, :), pW(:), temp_a(:, :)
  type(psb_d_vect_type)       :: r  
  type(psb_d_multivect_type)  :: Z, Q, P, V, temp
  integer(psb_ipk_)           :: itidx
  
  type(psb_itconv_type)         :: stopdat
  integer(psb_ipk_), parameter  :: Gram_solver_type = ione
  real(psb_dpk_)                :: derr 
  real(psb_dpk_)                :: cheb_coeff(3)


  info = psb_success_
  call psb_erractionsave(err_act)
  
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit, *) me, ' ', trim(name), ': from psb_info ', np
  
  if(s < 1) then
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info, name)
    goto 9999
  endif

  if ((.not. allocated(b%v)) .or. (.not.allocated(x%v))) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info, name)
    goto 9999
  endif

  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 2
  endif
  
  !  ISTOP_ = 1:  Normwise backward error, infinity norm 
  !  ISTOP_ = 2:  ||r||/||b||, 2-norm 
  if ((istop_ < 1 ) .or. (istop_ > 2 )) then
    info = psb_err_invalid_istop_
    err = info
    call psb_errpush(info, name, i_err = (/istop_/))
    goto 9999
  endif

  if (present(itmax)) then 
    itmax_ = itmax
  else
    itmax_ = 1000
  endif

  if (present(itrace)) then
    itrace_ = itrace
  else
    itrace_ = 0
  end if

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
  allocate(alpha(s), beta(s, s), W(s, s), pW(s), temp_a(s, s + 1), stat = info)
  if (info == psb_success_) call psb_geall(r, desc_a, info)
  if (info == psb_success_) call psb_geall(Z, desc_a, info, n = s)
  if (info == psb_success_) call psb_geall(Q, desc_a, info, n = s)
  if (info == psb_success_) call psb_geall(P, desc_a, info, n = s)
  if (info == psb_success_) call psb_geall(V, desc_a, info, n = s)
  if (info == psb_success_) call psb_geall(temp, desc_a, info, n = s)
  if (info == psb_success_) call psb_geasb(r, desc_a, info)
  if (info == psb_success_) call psb_geasb(Z, desc_a, info)
  if (info == psb_success_) call psb_geasb(Q, desc_a, info)
  if (info == psb_success_) call psb_geasb(P, desc_a, info)
  if (info == psb_success_) call psb_geasb(V, desc_a, info)
  if (info == psb_success_) call psb_geasb(temp, desc_a, info)

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! First residual calculation
  call psb_geaxpby(done, b, dzero, r, desc_a, info)
  if (info == psb_success_) call psb_spmm(-done, a, x, done, r, desc_a, info)
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! Init converence
  call psb_init_conv(methdname, istop_, itrace_, itmax_, a, x, b, eps, desc_a, stopdat, info)
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! check convergence here?

  ! Chebyshev coefficient calculation
  select case (base_type)
    case("M")
      cheb_coeff = dzero
    case("C")
      cheb_coeff = psb_d_chebyshev_coefficients(a, prec, desc_a, info, eigext)
    case default
      info = psb_err_invalid_input_ 
      call psb_errpush(info, name)
      goto 9999
  end select
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  
  ! First matrix power kernel
  call psb_pMPK(a, prec, r, Z, Q, s, desc_a, info, base_type = base_type, &
                  & alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3))
  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! Inizialization of P and V
  call psb_geaxpby(done, Z, dzero, P, desc_a, info)
  call psb_geaxpby(done, Q, dzero, V, desc_a, info)

  ! Loop until convergence (or maxiter)
  do itidx = 1, itmax_
    ! Compute matrix W and rhs for alpha 
    ! Compute matrix W and rhs for alpha 
    call psb_gedots(P, V, W, desc_a, info, global = .true.)
    call psb_gedots(P, r, alpha, desc_a, info, global = .true.)
    ! call psb_gedots(P, V, temp_a(:, 1 : s), desc_a, info, global = .false.)
    ! call psb_gedots(P, r, temp_a(:, s + 1), desc_a, info, global = .false.)
    ! call psb_sum(desc_a%get_ctxt(), temp_a, info)
    ! W = temp_a(:, 1 : s)
    ! alpha = temp_a(:, s + 1)

    ! Factor matrix W
    call dgetrf(s, s, W, s, pW, info)

    ! Solve for alpha
    call dgetrs('N', s, 1, W, s, pW, alpha, s, info)

    ! Update solution and residual
    call psb_geaxpby(P, alpha, x, desc_a, info, .true.)
    call psb_geaxpby(V, -alpha, r, desc_a, info, .true.)

    ! Check convergence. 
    if(psb_check_conv(methdname, itidx, x, r, desc_a, stopdat, info)) exit
    
    ! Matrix power kernel
    call psb_pMPK(a, prec, r, Z, Q, s, desc_a, info, base_type = base_type, &
                  & alpha = cheb_coeff(1), beta = cheb_coeff(2), gamma = cheb_coeff(3))

    ! Compute rhs for beta
    call psb_gedots(P, Q, beta, desc_a, info, .true.)

    ! Solve for beta
    beta = -beta;
    call dgetrs('N', s, s, W, s, pW, beta, s, info)

    ! Update P and V. Use of temp in needed because internal dgemm constraint
    call psb_geaxpby(P, beta, temp, desc_a, info, .false.)
    call psb_geaxpby(done, Z, done, temp, P, desc_a, info)
    call psb_geaxpby(V, beta, temp, desc_a, info, .false.)
    call psb_geaxpby(done, Q, done, temp, V, desc_a, info)
  end do

  call psb_end_conv(methdname, itidx, desc_a, stopdat, info, derr, iter)
  if (present(err)) err = derr
  if (present(iter)) iter = iter * s

  if (info == psb_success_) call psb_gefree(r, desc_a, info)
  if (info == psb_success_) call psb_gefree(Z, desc_a, info)
  if (info == psb_success_) call psb_gefree(Q, desc_a, info)
  if (info == psb_success_) call psb_gefree(P, desc_a, info)
  if (info == psb_success_) call psb_gefree(V, desc_a, info)
  if (info == psb_success_) call psb_gefree(temp, desc_a, info)

  if (info == psb_success_) deallocate(alpha, beta, W, pW, stat = info)
  if (info /= psb_success_) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

contains
  function psb_d_chebyshev_coefficients(a, prec, desc, info, eigext) result(coeff)
    use psb_eigsolve_mod
    type(psb_dspmat_type), intent(in)     :: a
    class(psb_dprec_type), intent(inout)  :: prec
    type(psb_desc_type), intent(in)       :: desc
    integer(psb_ipk_), intent(out)        :: info
    real(psb_dpk_), optional, intent(in)  :: eigext(2)
    
    real(psb_dpk_) :: lambda_max, lambda_min
    real(psb_dpk_) :: coeff(3)
 
    info = psb_success_

    if(present(eigext)) then
      lambda_min = eigext(1)
      lambda_max = eigext(2)
    else 
      call psb_powermethod(a, prec, lambda_max, desc, info)
      lambda_min = dzero
    end if

    if (info /= psb_success_) then 
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

    coeff(1) = 2_psb_dpk_ / (lambda_max - lambda_min)
    coeff(2) = (lambda_max + lambda_min) / (lambda_max - lambda_min)
    coeff(3) = done
  end function psb_d_chebyshev_coefficients
end subroutine psb_dscg_vect