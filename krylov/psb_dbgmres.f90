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
!    istop  -  integer(optional)          Input: stopping criterion, or how
!                                         to estimate the error.
!                                         1: err =  |r|/(|a||x|+|b|);  here the iteration is
!                                            stopped when  |r| <= eps * (|a||x|+|b|)
!                                         2: err =  |r|/|b|; here the iteration is
!                                            stopped when  |r| <= eps * |b|
!                                         where r is the (preconditioned, recursive
!                                         estimate of) residual.
!

subroutine psb_dbgmres_multivect(a, prec, b, x, eps, desc_a, info, itmax, iter, err, itrace, istop)

   use psb_base_mod
   use psb_prec_mod
   use psb_d_krylov_conv_mod
   use psb_krylov_mod

   implicit none

   type(psb_dspmat_type), intent(in)         :: a
   type(psb_desc_type), intent(in)           :: desc_a
   class(psb_dprec_type), intent(inout)      :: prec

   type(psb_d_multivect_type), intent(inout) :: b
   type(psb_d_multivect_type), intent(inout) :: x

   real(psb_dpk_), intent(in)                :: eps
   integer(psb_ipk_), intent(out)            :: info
   integer(psb_ipk_), optional, intent(in)   :: itmax, itrace, istop
   integer(psb_ipk_), optional, intent(out)  :: iter
   real(psb_dpk_), optional, intent(out)     :: err

   real(psb_dpk_), allocatable               :: aux(:), h(:,:), beta(:,:), y(:,:)

   type(psb_d_multivect_type), allocatable   :: v(:)
   type(psb_d_multivect_type)                :: w, r0, rm, pd

   integer(psb_ipk_)                         :: naux, itrace_, n_row, n_col, nrhs, nrep
   integer(psb_lpk_)                         :: mglob, n_add, ncv

   integer(psb_ipk_)                         :: i, j, istop_, err_act, idx_i, idx_j, idx
   integer(psb_ipk_)                         :: debug_level, debug_unit

   type(psb_ctxt_type)                       :: ctxt
   integer(psb_ipk_)                         :: np, me, itx
   real(psb_dpk_), allocatable               :: r0n2(:), rmn2(:)
   real(psb_dpk_)                            :: errnum, errden
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
      nrep = itmax
   else
      nrep = 10
   endif

   if (present(itrace)) then
      itrace_ = itrace
   else
      itrace_ = 0
   end if

   ncv = x%get_ncols()
   call psb_chkvect(mglob,ncv,x%get_nrows(),lone,lone,desc_a,info)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_chkvect on X')
      goto 9999
   end if
   ncv = b%get_ncols()
   call psb_chkvect(mglob,ncv,b%get_nrows(),lone,lone,desc_a,info)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_chkvect on B')
      goto 9999
   end if

   naux = 4*n_col
   nrhs = x%get_ncols()
   allocate(aux(naux),h((nrep+1)*nrhs,nrep*nrhs),y(nrep*nrhs,nrhs),r0n2(nrhs),rmn2(nrhs),stat=info)
   if (info == psb_success_) call psb_geall(v,desc_a,info,m=nrep+1,n=nrhs)
   if (info == psb_success_) call psb_geall(w,desc_a,info,n=nrhs)
   if (info == psb_success_) call psb_geall(r0,desc_a,info,n=nrhs)
   if (info == psb_success_) call psb_geall(rm,desc_a,info,n=nrhs)
   if (info == psb_success_) call psb_geall(pd,desc_a,info,n=nrhs)
   if (info == psb_success_) call psb_geasb(v,desc_a,info,mold=x%v,n=nrhs)
   if (info == psb_success_) call psb_geasb(w,desc_a,info,mold=x%v,n=nrhs)
   if (info == psb_success_) call psb_geasb(r0,desc_a,info,mold=x%v,n=nrhs)
   if (info == psb_success_) call psb_geasb(rm,desc_a,info,mold=x%v,n=nrhs)
   if (info == psb_success_) call psb_geasb(pd,desc_a,info,mold=x%v,n=nrhs)
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
      call psb_geaxpby(done,b,dzero,r0,desc_a,info)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if
      call psb_spmm(-done,a,x,done,r0,desc_a,info,work=aux)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if
      r0n2 = psb_genrm2(r0,desc_a,info)
   endif
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   h     = dzero
   y     = dzero
   deps  = eps
   itx   = 0
   n_add = nrhs-1

   if ((itrace_ > 0).and.(me == psb_root_)) call log_header(methdname)

   ! BGMRES algorithm

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
   beta = qr_fact(v(1))
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   ! STEP 3: Outer loop
   outer: do j=1,nrep

      ! Update itx counter
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

      ! STEP 5: Inner loop
      inner: do i=1,j

         ! Compute i index for H operations
         idx_i = (i-1)*nrhs+1

         ! STEP 6: Compute H(i,j) = (V(i)**T)*W
         h(idx_i:idx_i+n_add,idx_j:idx_j+n_add) = psb_gedot(v(i),w,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         ! STEP 7: Compute W = W - V(i)*H(i,j)

         ! Compute product V(i)*H(i,j)
         call psb_geprod(v(i),h(idx_i:idx_i+n_add,idx_j:idx_j+n_add),pd,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

         ! Compute W
         call psb_geaxpby(-done,pd,done,w,desc_a,info)
         if (info /= psb_success_) then
            info=psb_err_from_subroutine_non_
            call psb_errpush(info,name)
            goto 9999
         end if

      end do inner

      ! STEP 8: Compute QR_fact(W)

      ! Store R in H(j+1,j)
      h(idx_j+nrhs:idx_j+nrhs+n_add,idx_j:idx_j+n_add) = qr_fact(w)
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

      ! STEP 9: Compute Y(j)
      rmn2 = frobenius_norm_min(j)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      ! Compute residues
      if (istop_ == 1) then
         errnum = sum(rmn2)/size(rmn2)
         errden = sum(r0n2)/size(r0n2)
      end if

      ! Check convergence
      if (errnum <= eps*errden) then
         exit outer
      end if

      ! Log update
      if (itrace_ > 0) call log_conv(methdname,me,itx,ione,errnum,errden,deps)

   end do outer

   ! STEP 10: X(m) = X(0) + VT(m)*Y(m)
   do i=1,j

      ! Compute index for Y products
      idx = (i-1)*nrhs+1

      ! Compute product V(i)*Y(i)
      call psb_geprod(v(i),y(idx:idx+n_add,:),pd,desc_a,info)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      ! Add product to X(m)
      call psb_geaxpby(done,pd,done,x,desc_a,info)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

   end do

   ! END algorithm

   if (itrace_ > 0) call log_conv(methdname,me,itx,ione,errnum,errden,deps)

   call log_end(methdname,me,itx,itrace_,errnum,errden,deps,err=derr,iter=iter)
   if (present(err)) err = derr

   if (info == psb_success_) call psb_gefree(v,desc_a,info)
   if (info == psb_success_) call psb_gefree(w,desc_a,info)
   if (info == psb_success_) call psb_gefree(r0,desc_a,info)
   if (info == psb_success_) call psb_gefree(rm,desc_a,info)
   if (info == psb_success_) call psb_gefree(pd,desc_a,info)
   if (info == psb_success_) deallocate(aux,h,y,r0n2,rmn2,stat=info)
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

   ! QR factorization
   function qr_fact(mv) result(res)
      implicit none

      ! I/O parameters
      type(psb_d_multivect_type), intent(inout)  :: mv
      real(psb_dpk_), allocatable                :: res(:,:)

      ! Utils
      real(psb_dpk_), allocatable :: vec(:,:)
      real(psb_dpk_), allocatable :: tau(:), work(:)
      integer(psb_ipk_)           :: ii, jj, m, lda, lwork

      ! Allocate output
      allocate(res(nrhs,nrhs))

      ! Gather multivector to factorize
      call psb_gather(vec,mv,desc_a,info,root=psb_root_)

      ! If root
      if (me == psb_root_) then

         ! Initialize params
         m = size(vec,1)
         lda = m
         lwork = nrhs
         allocate(tau(nrhs),work(lwork))

         ! Perform QR factorization
         call dgeqrf(m,nrhs,vec,lda,tau,work,lwork,info)

         ! Set R
         res = vec(1:nrhs,1:nrhs)
         do ii = 2, nrhs
            do jj = 1, ii-1
               res(ii,jj) = dzero
            enddo
         enddo

         ! Generate Q matrix
         call dorgqr(m,nrhs,nrhs,vec,lda,tau,work,lwork,info)

         ! Deallocate
         deallocate(tau,work)

      end if

      ! Share R
      call psb_bcast(ctxt,res)

      ! Scatter Q
      call psb_scatter(vec,mv,desc_a,info,root=psb_root_,mold=mv%v)

   end function qr_fact

   ! TODO Loop con Givens rotation su ogni colonna
   ! Minimize Frobenius norm
   function frobenius_norm_min(rep) result(res)
      implicit none

      ! I/O parameters
      real(psb_dpk_), allocatable   :: res(:)
      integer(psb_ipk_), intent(in) :: rep

      ! Utils
      integer(psb_ipk_)           :: lwork, m_h, n_h, mn
      real(psb_dpk_), allocatable :: work(:), beta_e1(:,:), h_temp(:,:)

      ! Initialize params
      m_h    = (rep+1)*nrhs
      n_h    = rep*nrhs
      h_temp = h(1:m_h,1:n_h)
      mn     = min(m_h,n_h)
      lwork  = max(1,mn+max(mn,nrhs))
      allocate(work(lwork),beta_e1(m_h,nrhs),res(nrhs))

      ! Compute E1*beta
      beta_e1 = dzero
      beta_e1(1:nrhs,1:nrhs) = beta

      ! Compute min Frobenius norm
      call dgels('N',m_h,n_h,nrhs,h_temp,m_h,beta_e1,m_h,work,lwork,info)

      ! Set solution
      y = beta_e1(1:n_h,1:nrhs)

      ! Set residues
      res = sum(beta_e1(n_h+1:m_h,:)**2,dim=1)

      ! Deallocate
      deallocate(work,h_temp,beta_e1)

   end function frobenius_norm_min

end subroutine psb_dbgmres_multivect
