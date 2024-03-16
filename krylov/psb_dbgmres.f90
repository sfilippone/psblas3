! File:  psb_dbgmres.f90
!
! Subroutine: psb_dbgmres
!    This subroutine implements the BGMRES method with right preconditioning.
!
! Arguments:
!
!    a      -  type(psb_dspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_dprec_type)       Input: preconditioner
!    b      -  real,dimension(:,:)       Input: vector containing the
!                                         right hand side B
!    x      -  real,dimension(:,:)       Input/Output: vector containing the
!                                         initial guess and final solution X.
!    eps    -  real                       Input: Stopping tolerance; the iteration is
!                                         stopped when the error estimate |err| <= eps
!    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
!    info   -  integer.                   Output: Return code
!
!    itmax  -  integer(optional)          Input: maximum number of iterations to be
!                                         performed.
!
!    iter   -  integer(optional)          Output: how many iterations have been
!                                         performed.
!                                         performed.
!    err    -  real   (optional)          Output: error estimate on exit. If the
!                                         denominator of the estimate is exactly
!                                         0, it is changed into 1.
!    itrace -  integer(optional)          Input: print an informational message
!                                         with the error estimate every itrace
!                                         iterations
!    itrs   -  integer(optional)          Input: iteration number parameter
!    istop  -  integer(optional)          Input: stopping criterion, or how
!                                         to estimate the error.
!                                         1: err =  |r|/(|a||x|+|b|);  here the iteration is
!                                            stopped when  |r| <= eps * (|a||x|+|b|)
!                                         2: err =  |r|/|b|; here the iteration is
!                                            stopped when  |r| <= eps * |b|
!                                         where r is the (preconditioned, recursive
!                                         estimate of) residual.
!

subroutine psb_dbgmres_multivect(a, prec, b, x, eps, desc_a, info, itmax, iter, err, itrace, itrs, istop)

   use psb_base_mod
   use psb_prec_mod
   use psb_d_krylov_conv_mod
   use psb_krylov_mod

   implicit none

   type(psb_dspmat_type), intent(in)         :: a
   type(psb_desc_type), Intent(in)           :: desc_a
   class(psb_dprec_type), intent(inout)      :: prec

   type(psb_d_multivect_type), Intent(inout) :: b
   type(psb_d_multivect_type), Intent(inout) :: x

   real(psb_dpk_), Intent(in)                :: eps
   integer(psb_ipk_), intent(out)            :: info
   integer(psb_ipk_), Optional, Intent(in)   :: itmax, itrace, itrs, istop
   integer(psb_ipk_), Optional, Intent(out)  :: iter
   real(psb_dpk_), Optional, Intent(out)     :: err

   real(psb_dpk_), allocatable               :: aux(:), h(:,:), beta(:,:), beta_e1(:,:)

   type(psb_d_multivect_type), allocatable   :: v(:)
   type(psb_d_multivect_type)                :: v_tot, w

   real(psb_dpk_)                            :: t1, t2

   real(psb_dpk_)                            :: rti, rti1
   integer(psb_ipk_)                         :: litmax, naux, itrace_, n_row, n_col, nrhs, nrep
   integer(psb_lpk_)                         :: mglob, n_add

   integer(psb_ipk_)                         :: i, j, k, istop_, err_act, idx_i, idx_j, idx
   integer(psb_ipk_)                         :: debug_level, debug_unit

   type(psb_ctxt_type)                       :: ctxt
   integer(psb_ipk_)                         :: np, me, itx
   real(psb_dpk_)                            :: rni, xni, bni, ani, bn2, r0n2
   real(psb_dpk_)                            :: errnum, errden, deps, derr
   character(len=20)                         :: name
   character(len=*), parameter               :: methdname='BGMRES'

   info = psb_success_
   name = 'psb_dbgmres'
   call psb_erractionsave(err_act)
   debug_unit  = psb_get_debug_unit()
   debug_level = psb_get_debug_level()

   ctxt = desc_a%get_context()
   call psb_info(ctxt, me, np)

   if (debug_level >= psb_debug_ext_) &
   & write(debug_unit,*) me,' ',trim(name),': from psb_info',np

   if (.not.allocated(b%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
   endif
   if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
   endif

   mglob = desc_a%get_global_rows()
   n_row = desc_a%get_local_rows()
   n_col = desc_a%get_local_cols()

   if (present(istop)) then
      istop_ = istop
   else
      istop_ = 2
   endif

   if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
      info=psb_err_invalid_istop_
      err=info
      call psb_errpush(info,name,i_err=(/istop_/))
      goto 9999
   endif

   if (present(itmax)) then
      litmax = itmax
   else
      litmax = 1000
   endif

   if (present(itrace)) then
      itrace_ = itrace
   else
      itrace_ = 0
   end if

   if (present(itrs)) then
      nrep = itrs
      if (debug_level >= psb_debug_ext_) &
      & write(debug_unit,*) me,' ',trim(name),&
      & ' present: itrs: ',itrs,nrep
   else
      nrep = 10
      if (debug_level >= psb_debug_ext_) &
      & write(debug_unit,*) me,' ',trim(name),&
      & ' not present: itrs: ',itrs,nrep
   endif
   if (nrep <=0 ) then
      info=psb_err_invalid_irst_
      err=info
      call psb_errpush(info,name,i_err=(/nrep/))
      goto 9999
   endif

   call psb_chkvect(mglob,x%get_ncols(),x%get_nrows(),lone,lone,desc_a,info)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_chkvect on X')
      goto 9999
   end if
   call psb_chkvect(mglob,b%get_ncols(),b%get_nrows(),lone,lone,desc_a,info)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_chkvect on B')
      goto 9999
   end if

   naux = 4*n_col
   nrhs = x%get_ncols()
   allocate(aux(naux),h((nrep+1)*nrhs,nrep*nrhs),stat=info)
   if (info == psb_success_) call psb_geall(v,desc_a,info,m=nrep+1,n=nrhs)
   if (info == psb_success_) call psb_geall(v_tot,desc_a,info,n=nrep*nrhs)
   if (info == psb_success_) call psb_geall(w,desc_a,info,n=nrhs)
   if (info == psb_success_) call psb_geasb(v,desc_a,info,mold=x%v,n=nrhs)
   if (info == psb_success_) call psb_geasb(v_tot,desc_a,info,mold=x%v,n=nrep*nrhs)
   if (info == psb_success_) call psb_geasb(w,desc_a,info,mold=x%v,n=nrhs)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   if (debug_level >= psb_debug_ext_) &
   & write(debug_unit,*) me,' ',trim(name),&
   & ' Size of V,W ',v(1)%get_nrows(),size(v),&
   & w%get_nrows()

   if (istop_ == 1) then
      ani = psb_spnrmi(a,desc_a,info)
      bni = psb_geamax(b,desc_a,info)
   else if (istop_ == 2) then
      bn2 = psb_genrm2(b,desc_a,info)
   endif
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   h      = dzero
   errnum = dzero
   errden = done
   deps   = eps
   itx    = 0
   n_add  = nrhs-1

   if ((itrace_ > 0).and.(me == psb_root_)) call log_header(methdname)

   ! BGMRES algorithm

   ! TODO QR fact seriale per ora

   ! STEP 1: Compute R(0) = B - A*X(0)

   ! Store B in V(1)
   call psb_geaxpby(done,b,dzero,v(1),desc_a,info)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   ! Store R(0) in V(1)
   call psb_spmm(-done,a,x,done,v(1),desc_a,info,work=aux)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   ! STEP 2: Compute QR_fact(R(0))
   beta = psb_geqrfact(v(1),desc_a,info)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   ! STEP 3: Outer loop
   outer: do j=1,nrep

      ! TODO Check convergence
      ! if (istop_ == 1) then
      !    rni = psb_geamax(v(1),desc_a,info)
      !    xni = psb_geamax(x,desc_a,info)
      !    errnum = rni
      !    errden = (ani*xni+bni)
      ! else if (istop_ == 2) then
      !    rni = psb_genrm2(v(1),desc_a,info)
      !    errnum = rni
      !    errden = bn2
      ! endif
      ! if (info /= psb_success_) then
      !    info=psb_err_from_subroutine_non_
      !    call psb_errpush(info,name)
      !    goto 9999
      ! end if

      ! if (errnum <= eps*errden) exit outer

      ! if (itrace_ > 0) call log_conv(methdname,me,itx,itrace_,errnum,errden,deps)

      itx = itx + 1

      ! Compute j index for H operations
      idx_j = (j-1)*nrhs+1

      ! STEP 4: Compute W = AV(j)
      call psb_spmm(done,a,v(j),dzero,w,desc_a,info,work=aux)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      if (itx >= litmax) exit outer

      ! STEP 5: Inner loop
      inner: do i=1,j

         ! Compute i index for H operations
         idx_i = (i-1)*nrhs+1

         ! STEP 6: Compute H(i,j) = V(i)_T*W
         h(idx_i:idx_i+n_add,idx_j:idx_j+n_add) = psb_geprod(v(i),w,desc_a,info,trans=.true.)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         ! STEP 7: Compute W = W - V(i)*H(i,j)

         ! TODO si blocca con NRHS grandi?
         !temp = psb_geprod(v(i),h(idx_i:idx_i+n_add,idx_j:idx_j+n_add),desc_a,info,global=.false.)
         call psb_geaxpby(-done,psb_geprod(v(i),h(idx_i:idx_i+n_add,idx_j:idx_j+n_add),desc_a,info,global=.false.),done,w,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

      end do inner

      ! STEP 8: Compute QR_fact(W)

      ! Store R in H(j+1,j)
      h(idx_j+nrhs:idx_j+nrhs+n_add,idx_j:idx_j+n_add) = psb_geqrfact(w,desc_a,info)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      ! Store Q in V(j+1)
      call psb_geaxpby(done,w,dzero,v(j+1),desc_a,info)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

   end do outer

   ! STEP 9: Compute Y(m)
   call frobenius_norm_min()
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   ! TODO V_tot comprende V(nrep+1)?
   ! STEP 10: Compute V = {V(1),...,V(m)}
   do i=1,nrep
      idx = (i-1)*nrhs+1
      v_tot%v%v(1:n_row,idx:idx+n_add) = v(i)%v%v(1:n_row,1:nrhs)
   enddo

   ! STEP 11: X(m) = X(0) + V*Y(m)
   call psb_geaxpby(done,psb_geprod(v_tot,beta_e1,desc_a,info,global=.false.),done,x,desc_a,info)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   ! END algorithm

   if (itrace_ > 0) call log_conv(methdname,me,itx,ione,errnum,errden,deps)

   call log_end(methdname,me,itx,itrace_,errnum,errden,deps,err=derr,iter=iter)
   if (present(err)) err = derr

   if (info == psb_success_) call psb_gefree(v,desc_a,info)
   if (info == psb_success_) call psb_gefree(v_tot,desc_a,info)
   if (info == psb_success_) call psb_gefree(w,desc_a,info)
   if (info == psb_success_) deallocate(aux,h,stat=info)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   call psb_erractionrestore(err_act)
   return

9999 call psb_error_handler(err_act)
   return

contains

   ! Minimize Frobenius norm
   subroutine frobenius_norm_min()

      implicit none

      integer(psb_ipk_)           :: lwork
      real(psb_dpk_), allocatable :: work(:), beta_temp(:,:)

      integer(psb_ipk_) :: m_h, n_h, mn

      ! Initialize params
      m_h   = (nrep+1)*nrhs
      n_h   = nrep*nrhs
      mn    = min(m_h,n_h)
      lwork = max(1,mn+max(mn,nrhs))
      allocate(work(lwork))

      ! Compute E1*beta
      allocate(beta_temp(m_h,nrhs))
      beta_temp = dzero
      beta_temp(1:nrhs,1:nrhs) = beta

      ! Compute min Frobenius norm
      call dgels('N',m_h,n_h,nrhs,h,m_h,beta_temp,m_h,work,lwork,info)

      ! Set solution
      allocate(beta_e1(n_h,nrhs))
      beta_e1 = beta_temp(1:n_h,1:nrhs)

      ! Deallocate
      deallocate(work,beta,beta_temp)

      return

   end subroutine frobenius_norm_min

end subroutine psb_dbgmres_multivect
