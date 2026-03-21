subroutine psb_dscg_vect(a, prec, b, x, s, eps, desc_a, info, itmax, iter, err, itrace, istop)
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
  type(psb_desc_type), intent(in)       :: desc_a
  integer(psb_ipk_), intent(out)        :: info
  integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
  integer(psb_ipk_), optional, intent(out)  :: iter
  real(psb_dpk_), optional, intent(out)     :: err

  ! Local vars
  integer(psb_ipk_)   :: istop_, itmax_, itrace_
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: err_act, np, me, debug_level, debug_unit, &
                          & n_col, n_row
  integer(psb_lpk_)   :: mglob
  character(len=20)           :: name = 'psb_dscg'
  character(len=*), parameter :: methdname = 'sStepCG'

  real(psb_dpk_), allocatable :: alpha(:), beta(:, :), W(:, :), pW(:)
  type(psb_d_vect_type)       :: r  
  type(psb_d_multivect_type)  :: Z, Q, P, V, temp
  integer(psb_ipk_)           :: itidx
  
  type(psb_itconv_type)       :: stopdat
  real(psb_dpk_)              :: derr 


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

  !TO DO: allocate work arrays

  !Allocate and assembly data structure
  allocate(alpha(s), beta(s, s), W(s, s), pW(s), stat = info)
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
  if (info /= psb_success_) Then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info, name)
    goto 9999
  end if

  ! check convergence here?

  ! TODO: chebyshev coefficient calculation

  ! First matrix power kernel (now monomial)
  call psb_pMPK(a, prec, r, Z, Q, s, desc_a, info)
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
    ! Compute and factor matrix W
    call psb_gedots(P, V, W, desc_a, info, .true.)
    call dgetrf(s, s, W, s, pW, info)

    ! Compute rhs for alpha
    call psb_gedots(P, r, alpha, desc_a, info, .true.)

    ! Solve for alpha
    call dgetrs('N', s, 1, W, s, pW, alpha, s, info)

    ! Update solution and residual
    call psb_geaxpby(P, alpha, x, desc_a, info, .true.)
    call psb_geaxpby(V, -alpha, r, desc_a, info, .true.)

    ! Check convergence. 
    if(psb_check_conv(methdname, itidx, x, r, desc_a, stopdat, info)) exit
    
    ! Matrix power kernel
    call psb_pMPK(a, prec, r, Z, Q, s, desc_a, info)

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
end subroutine psb_dscg_vect