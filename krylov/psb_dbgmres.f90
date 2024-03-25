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

   real(psb_dpk_), allocatable               :: aux(:), h(:,:), vt(:,:), beta(:,:), y(:,:)

   type(psb_d_multivect_type), allocatable   :: v(:)
   type(psb_d_multivect_type)                :: w, xt, r

   real(psb_dpk_)                            :: t1, t2

   real(psb_dpk_)                            :: rti, rti1
   integer(psb_ipk_)                         :: litmax, naux, itrace_, n_row, n_col, nrhs, nrep
   integer(psb_lpk_)                         :: mglob, n_add

   integer(psb_ipk_)                         :: i, j, k, col, istop_, err_act, idx_i, idx_j, idx
   integer(psb_ipk_)                         :: debug_level, debug_unit

   type(psb_ctxt_type)                       :: ctxt
   integer(psb_ipk_)                         :: np, me, itx
   real(psb_dpk_), allocatable               :: r0n2(:), rmn2(:)
   real(psb_dpk_), allocatable               :: errnum(:), errden(:)
   real(psb_dpk_)                            :: deps, derr
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
      istop_ = 1
   endif

   if (istop_ /= 1) then
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
   allocate(aux(naux),h((nrep+1)*nrhs,nrep*nrhs),y(nrep*nrhs,nrhs),&
         & vt(n_row,(nrep+1)*nrhs),r0n2(nrhs),rmn2(nrhs),&
         & errnum(nrhs),errden(nrhs),stat=info)
   if (info == psb_success_) call psb_geall(v,desc_a,info,m=nrep+1,n=nrhs)
   if (info == psb_success_) call psb_geall(w,desc_a,info,n=nrhs)
   if (info == psb_success_) call psb_geall(xt,desc_a,info,n=nrhs)
   if (info == psb_success_) call psb_geall(r,desc_a,info,n=nrhs)
   if (info == psb_success_) call psb_geasb(v,desc_a,info,mold=x%v,n=nrhs)
   if (info == psb_success_) call psb_geasb(w,desc_a,info,mold=x%v,n=nrhs)
   if (info == psb_success_) call psb_geasb(xt,desc_a,info,mold=x%v,n=nrhs)
   if (info == psb_success_) call psb_geasb(r,desc_a,info,mold=x%v,n=nrhs)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   if (debug_level >= psb_debug_ext_) &
   & write(debug_unit,*) me,' ',trim(name),&
   & ' Size of V,W ',v(1)%get_nrows(),size(v),&
   & w%get_nrows()

   ! Compute norm2 of R(0)
   if (istop_ == 1) then
      call psb_geaxpby(done,b,dzero,r,desc_a,info)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if
      call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if
      r0n2 = psb_genrm2(r,desc_a,info)
   endif
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   h      = dzero
   y      = dzero
   errnum = dzero
   errden = done
   deps   = eps
   itx    = 0
   n_add  = nrhs-1

   if ((itrace_ > 0).and.(me == psb_root_)) call log_header(methdname)

   ! BGMRES algorithm

   ! TODO Con tanti ITRS e tanti NRHS si ottengono NaN
   ! TODO Deflazione e restart dopo aver trovato una colonna, difficile...

   ! TODO L'algo converge abbastanza bene. Capire come fare check residui

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

   ! Add V(1) to VT
   vt(:,1:nrhs) = v(1)%get_vect()

   ! STEP 3: Outer loop
   outer: do j=1,nrep

      ! Update itx counter
      itx = itx + 1
      if (itx >= litmax) exit outer

      ! Compute j index for H operations
      idx_j = (j-1)*nrhs+1

      ! STEP 4: Compute W = AV(j)
      call psb_spmm(done,a,v(j),dzero,w,desc_a,info,work=aux)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      ! STEP 5: Inner loop
      inner: do i=1,j

         ! Compute i index for H operations
         idx_i = (i-1)*nrhs+1

         ! STEP 6: Compute H(i,j) = (V(i)**T)*W
         h(idx_i:idx_i+n_add,idx_j:idx_j+n_add) = psb_geprod(v(i),w,desc_a,info,trans=.true.)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         ! STEP 7: Compute W = W - V(i)*H(i,j)
         call psb_geaxpby(-done,&
         & psb_geprod(v(i),h(idx_i:idx_i+n_add,idx_j:idx_j+n_add),desc_a,info,global=.false.),&
         & done,w,desc_a,info)
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

      ! Add V(j+1) to VT
      idx = j*nrhs+1
      vt(:,idx:idx+n_add) = v(j+1)%get_vect()

      ! STEP 9: Compute Y(j)
      call frobenius_norm_min(j)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      ! Compute residues
      if (istop_ == 1) then

         ! TODO Compute R(j) = R(0) - VT(j+1)*H(j)*Y(j)
         call psb_geaxpby(-done,psb_geprod(psb_geprod(vt(:,1:(j+1)*nrhs),h(1:(j+1)*nrhs,1:j*nrhs),&
                           & desc_a,info,global=.false.),&
                           & y(1:j*nrhs,1:nrhs),desc_a,info,global=.false.),done,r,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         write(*,*)
         do col=1,r%get_nrows()
            write(*,*) r%v%v(col,:)
         end do
         write(*,*)

         ! TODO Calcolo soluzione al passo J e vedo i residui (se minore esco dal ciclo)
         ! TODO Compute R(j) = B - A*X(j)

         ! Copy X in XT
         call psb_geaxpby(done,x,dzero,xt,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         ! Compute current solution X(j) = X(0) + VT(j)*Y(j)
         call psb_geaxpby(done,psb_geprod(vt(:,1:j*nrhs),y(1:j*nrhs,:),desc_a,info,global=.false.),done,xt,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         ! Copy B in R
         call psb_geaxpby(done,b,dzero,r,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         ! Compute R(j) = B - A*X(j)
         call psb_spmm(-done,a,xt,done,r,desc_a,info,work=aux)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         write(*,*)
         do col=1,r%get_nrows()
            write(*,*) r%v%v(col,:)
         end do
         write(*,*)

         ! Compute nrm2 of each column of R(j)
         rmn2 = psb_genrm2(r,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         ! Set error num and den
         errnum = rmn2
         errden = r0n2

         ! TODO Ogni entrata della norma2 di R(m) deve essere più piccola di tolleranza*nrm2(r0)
         do col=1,nrhs
            write(*,*) rmn2(col), r0n2(col)
         end do
      end if

      ! TODO Norma dei residui con Xm devono essere minori di tolleranza * nrm2(R0)?
      ! Check convergence
      if (maxval(errnum) <= eps*maxval(errden)) then

         ! Compute result and exit
         if (istop_ == 1) then

            ! Compute X(j) = X(0) + VT(j)*Y(j)
            ! call psb_geaxpby(done,psb_geprod(vt(:,1:j*nrhs),y(1:j*nrhs,:),desc_a,info,global=.false.),done,x,desc_a,info)
            ! if (info /= psb_success_) then
            !    info=psb_err_from_subroutine_non_
            !    call psb_errpush(info,name)
            !    goto 9999
            ! end if

            ! Copy current solution XT in X
            call psb_geaxpby(done,xt,dzero,x,desc_a,info)
            if (info /= psb_success_) then
               info=psb_err_from_subroutine_non_
               call psb_errpush(info,name)
               goto 9999
            end if

         end if

         ! Exit algorithm
         exit outer

      end if

      ! Log update
      if (itrace_ > 0) call log_conv(methdname,me,itx,ione,maxval(errnum),maxval(errden),deps)

   end do outer

   ! STEP 10: X(m) = X(0) + VT(m)*Y(m)
   ! call psb_geaxpby(done,psb_geprod(vt(:,1:nrep*nrhs),y,desc_a,info,global=.false.),done,x,desc_a,info)
   ! if (info /= psb_success_) then
   !    info=psb_err_from_subroutine_non_
   !    call psb_errpush(info,name)
   !    goto 9999
   ! end if

   ! END algorithm

   ! TODO log_conv passa scalari errnum,errden,deps (servono vettori)
   ! TODO Inizialmente versione verbosa che stampa errore per tutte le colonne
   ! TODO Versione finale che stampa errore massimo (si può usare log_conv con questo)
   if (itrace_ > 0) call log_conv(methdname,me,itx,ione,maxval(errnum),maxval(errden),deps)

   call log_end(methdname,me,itx,itrace_,maxval(errnum),maxval(errden),deps,err=derr,iter=iter)
   if (present(err)) err = derr

   if (info == psb_success_) call psb_gefree(v,desc_a,info)
   if (info == psb_success_) call psb_gefree(w,desc_a,info)
   if (info == psb_success_) call psb_gefree(xt,desc_a,info)
   if (info == psb_success_) call psb_gefree(r,desc_a,info)
   if (info == psb_success_) deallocate(aux,h,y,vt,stat=info)
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
   subroutine frobenius_norm_min(rep)

      implicit none

      integer(psb_ipk_), intent(in) :: rep

      integer(psb_ipk_)           :: lwork
      real(psb_dpk_), allocatable :: work(:), beta_e1(:,:)

      real(psb_dpk_), allocatable :: h_temp(:,:)
      integer(psb_ipk_)           :: m_h, n_h, mn

      ! Initialize params
      h_temp = h
      m_h    = (rep+1)*nrhs
      n_h    = rep*nrhs
      mn     = min(m_h,n_h)
      lwork  = max(1,mn+max(mn,nrhs))
      allocate(work(lwork))

      ! Compute E1*beta
      allocate(beta_e1(m_h,nrhs))
      beta_e1 = dzero
      beta_e1(1:nrhs,1:nrhs) = beta

      ! Compute min Frobenius norm
      call dgels('N',m_h,n_h,nrhs,h_temp(1:m_h,1:n_h),m_h,beta_e1,m_h,work,lwork,info)

      ! Set solution
      y = beta_e1(1:n_h,1:nrhs)

      ! Deallocate
      deallocate(work,h_temp,beta_e1)

      return

   end subroutine frobenius_norm_min

end subroutine psb_dbgmres_multivect
