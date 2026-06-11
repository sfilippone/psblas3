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
! File:  psb_dminres.f90
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   C                                                                      C
!   C  References:                                                         C
!   C          [1] Duff, I., Marrone, M., Radicati, G., and Vittoli, C.    C
!   C              Level 3 basic linear algebra subprograms for sparse     C
!   C              matrices: a user level interface                        C
!   C              ACM Trans. Math. Softw., 23(3), 379-401, 1997.          C
!   C                                                                      C
!   C                                                                      C
!   C         [2]  S. Filippone, M. Colajanni                              C
!   C              PSBLAS: A library for parallel linear algebra           C
!   C              computation on sparse matrices                          C
!   C              ACM Trans. on Math. Softw., 26(4), 527-550, Dec. 2000.  C
!   C                                                                      C
!   C         [3] M. Arioli, I. Duff, M. Ruiz                              C
!   C             Stopping criteria for iterative solvers                  C
!   C             SIAM J. Matrix Anal. Appl., Vol. 13, pp. 138-144, 1992   C
!   C                                                                      C
!   C                                                                      C
!   C         [4] R. Barrett et al                                         C
!   C             Templates for the solution of linear systems             C
!   C             SIAM, 1993
!   C                                                                      C
!   C                                                                      C
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_dminres.f90
!
! Subroutine: psb_dminres
!    This subroutine implements the MINRES method.
!
!
! Arguments:
!
!    a      -  type(psb_dspmat_type)    Input: sparse matrix containing A.
!    prec   -  class(psb_dprec_type)    Input: preconditioner
!    b(:)   -  real                     Input: vector containing the
!                                         right hand side B
!    x(:)   -  real                     Input/Output: vector containing the
!                                         initial guess and final solution X.
!    eps    -  real                       Input: Stopping tolerance; the iteration is
!                                         stopped when the error estimate |err| <= eps
!    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
!    info   -  integer.                   Output: Return code
!
!    itmax  -  integer(optional)          Input: maximum number of iterations to be
!                                         performed.
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
!
subroutine psb_dminres_vect(a,prec,b,x,eps,desc_a,info,&
& itmax,iter,err,itrace,istop,cond)
   use psb_base_mod
   use psb_prec_mod
   use psb_d_linsolve_conv_mod
   use psb_linsolve_mod
   implicit none
   type(psb_dspmat_type), intent(in)    :: a
   Type(psb_desc_type), Intent(in)      :: desc_a
   class(psb_dprec_type), intent(inout) :: prec
   type(psb_d_vect_type), Intent(inout) :: b
   type(psb_d_vect_type), Intent(inout) :: x
   Real(psb_dpk_), Intent(in)           :: eps
   integer(psb_ipk_), intent(out)                 :: info
   integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, istop
   integer(psb_ipk_), Optional, Intent(out)       :: iter
   Real(psb_dpk_), Optional, Intent(out) :: err, cond
! =   Local data
   real(psb_dpk_), allocatable, target   :: aux(:)
   type(psb_d_vect_type), allocatable, target :: wwrk(:)
   type(psb_d_vect_type), pointer  :: y, r1, r2, v, w, w1, w2, wtmp, res
   real(psb_dpk_)   :: alfa, beta, beta1, dbar, delta, denom
   real(psb_dpk_)   :: epsln, gamma, gbar, oldb, oldeps, phi, phibar
   real(psb_dpk_)   :: rhs1, rhs2, s, zeta, betatmp
   real(psb_dpk_)   :: rotx(1), roty(1)
   real(psb_dpk_)  :: cs
   real(psb_dpk_)   :: sn
   integer(psb_ipk_) :: itmax_, istop_, naux, itx, itrace_, n_col, n_row, err_act
   integer(psb_lpk_) :: mglob
   integer(psb_ipk_) :: debug_level, debug_unit
   type(psb_ctxt_type) :: ctxt
   integer(psb_ipk_) :: np, me
   real(psb_dpk_)     :: derr, errnum, errden, deps
   logical                     :: do_cond
   logical                     :: converged
   real(psb_dpk_)               :: epsmach
   real(psb_dpk_)               :: gmax, gmin, tnorm2, ynorm2
   real(psb_dpk_)               :: rni, xni, bni, ani, bn2, r0n2
   real(psb_dpk_)               :: dotv
   character(len=20)           :: name
   character(len=*), parameter :: methdname='MINRES'

   info = psb_success_
   name = 'psb_dminres'
   call psb_erractionsave(err_act)
   debug_unit  = psb_get_debug_unit()
   debug_level = psb_get_debug_level()

   ctxt = desc_a%get_context()

   call psb_info(ctxt, me, np)
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
      istop_ = psb_get_istop_default()
   endif
   if (.not.psb_is_valid_istop(istop_)) then
      info=psb_err_invalid_istop_
      err=info
      call psb_errpush(info,name,i_err=(/istop_/))
      goto 9999
   end if
   !
   !  istop_ = 1:  normwise backward error, infinity norm
   !  istop_ = 2:  ||r||/||b||   norm 2
   !
   select case(istop_)
    case(psb_istop_ani_,psb_istop_bn2_,&
    & psb_istop_rn2_abs_, psb_istop_rrn2_)
      ! nothing needed
    case default
      ! should never get here
      info=psb_err_internal_error_
      err=info
      call psb_errpush(info,name,a_err="invalid istop_")
      goto 9999
   end select

   call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
   if (info == psb_success_)&
   & call psb_chkvect(mglob,lone,b%get_nrows(),lone,lone,desc_a,info)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_chkvect on X/B')
      goto 9999
   end if

   naux=4*n_col
   allocate(aux(naux), stat=info)
   if (info == psb_success_) call psb_geall(wwrk,desc_a,info,n=9_psb_ipk_)
   if (info == psb_success_) call psb_geasb(wwrk,desc_a,info,mold=x%v,scratch=.true.)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

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

   y    => wwrk(1)
   r1   => wwrk(2)
   r2   => wwrk(3)
   v    => wwrk(4)
   w    => wwrk(5)
   w1   => wwrk(6)
   w2   => wwrk(7)
   wtmp => wwrk(8)
   res  => wwrk(9)

   do_cond = present(cond)
   epsmach = epsilon(done)

   ! res = b - A*x
   call psb_geaxpby(done,b,dzero,res,desc_a,info)
   if (info == psb_success_) call psb_spmm(-done,a,x,done,res,desc_a,info,work=aux)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   select case(istop_)
   case(psb_istop_ani_)
      ani = psb_spnrmi(a,desc_a,info)
      bni = psb_geamax(b,desc_a,info)
   case(psb_istop_bn2_)
      bn2 = psb_genrm2(b,desc_a,info)
   case(psb_istop_rn2_abs_)
      ! nothing needed
   case(psb_istop_rrn2_)
      r0n2 = psb_genrm2(res,desc_a,info)
   end select
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   errnum = dzero
   errden = done
   deps   = eps
   converged = .false.
   if ((itrace_ > 0).and.(me == 0)) call log_header(methdname)

   ! y  = beta1 * P' * v1, with v1 the first Lanczos vector.
   call psb_geaxpby(done,res,dzero,y,desc_a,info)
   if (info == psb_success_) call psb_geaxpby(done,res,dzero,r1,desc_a,info)
   if (info == psb_success_) call prec%apply(res,y,desc_a,info,work=aux)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   dotv = psb_gedot(res,y,desc_a,info)
   if (dotv < -epsmach*max(done,abs(dotv))) then
      info=psb_err_pivot_too_small_
      call psb_errpush(info,name,a_err='MINRES breakdown: negative (r,z)')
      goto 9999
   end if
   if (dotv < dzero) dotv = dzero
   beta1 = sqrt(dotv)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   ! If beta1 is numerically zero, x is already acceptable.
   if (abs(beta1) <= epsmach) then
      itx = 0
      converged = .true.
      if (do_cond) cond = dzero
      select case(istop_)
      case(psb_istop_ani_)
         rni = psb_geamax(res,desc_a,info)
         xni = psb_geamax(x,desc_a,info)
         errnum = rni
         errden = ani*xni + bni
      case(psb_istop_bn2_)
         errnum = abs(beta1)
         errden = bn2
      case(psb_istop_rn2_abs_)
         errnum = abs(beta1)
         errden = done
      case(psb_istop_rrn2_)
         errnum = abs(beta1)
         errden = r0n2
      end select
      if (errden == dzero) errden = done
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if
   end if

   beta   = beta1
   oldb   = dzero
   dbar   = dzero
   epsln  = dzero
   phibar = beta1
   rhs1   = beta1
   rhs2   = dzero
   cs     = -done
   sn     = dzero
   itx    = 0

   if (do_cond) then
      tnorm2 = 0.0_psb_dpk_
      ynorm2 = 0.0_psb_dpk_
      gmax   = 0.0_psb_dpk_
      gmin   = huge(done)
   end if

   call psb_geaxpby(done,r1,dzero,r2,desc_a,info)
   if (info == psb_success_) call psb_geaxpby(dzero,x,dzero,w,desc_a,info)
   if (info == psb_success_) call psb_geaxpby(dzero,x,dzero,w1,desc_a,info)
   if (info == psb_success_) call psb_geaxpby(dzero,x,dzero,w2,desc_a,info)
   if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
   end if

   if (.not.converged) then
   do while (itx < itmax_)
      if (abs(beta) <= epsmach) exit
      itx = itx + 1

      s = done/beta
      call psb_geaxpby(s,y,dzero,v,desc_a,info)
      if (info == psb_success_) call psb_spmm(done,a,v,dzero,y,desc_a,info,work=aux)
      if (itx >= 2 .and. info == psb_success_) then
         call psb_geaxpby((-beta/oldb),r1,done,y,desc_a,info)
      end if
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      alfa = psb_gedot(v,y,desc_a,info)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      call psb_geaxpby((-alfa/beta),r2,done,y,desc_a,info)
      if (info == psb_success_) call psb_geaxpby(done,r2,dzero,r1,desc_a,info)
      if (info == psb_success_) call psb_geaxpby(done,y,dzero,r2,desc_a,info)
      if (info == psb_success_) call prec%apply(r2,y,desc_a,info,work=aux)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      oldb = beta
      dotv = psb_gedot(r2,y,desc_a,info)
      if (dotv < -epsmach*max(done,abs(dotv))) then
         info=psb_err_pivot_too_small_
         call psb_errpush(info,name,a_err='MINRES breakdown: negative (r2,z2)')
         goto 9999
      end if
      if (dotv < dzero) dotv = dzero
      beta = sqrt(dotv)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      if (do_cond) then
        tnorm2 = tnorm2 + abs(alfa)**2 + abs(oldb)**2 + abs(beta)**2
        if (itx == 1) then
          gmax = abs(alfa)
          gmin = gmax
        end if
      end if

      oldeps = epsln
      delta  = dbar
      gbar   = alfa
      rotx(1) = delta
      roty(1) = gbar
      call drot(1,rotx,1,roty,1,cs,sn)
      delta = rotx(1)
      gbar  = -roty(1)

      epsln  = dzero
      dbar   = beta
      rotx(1) = epsln
      roty(1) = dbar
      call drot(1,rotx,1,roty,1,cs,sn)
      epsln = rotx(1)
      dbar  = -roty(1)

      gamma   = gbar
      betatmp = beta
      call drotg(gamma,betatmp,cs,sn)
      if (abs(gamma) <= epsmach) exit

      phi    = phibar
      phibar = dzero
      rotx(1) = phi
      roty(1) = phibar
      call drot(1,rotx,1,roty,1,cs,-sn)
      phi    = rotx(1)
      phibar = roty(1)

      denom = done/gamma
      call psb_geaxpby(done,w2,dzero,w1,desc_a,info)
      if (info == psb_success_) call psb_geaxpby(done,w,dzero,w2,desc_a,info)
      if (info == psb_success_) call psb_geaxpby(done,v,dzero,wtmp,desc_a,info)
      if (info == psb_success_) call psb_geaxpby(-oldeps,w1,done,wtmp,desc_a,info)
      if (info == psb_success_) call psb_geaxpby(-delta,w2,done,wtmp,desc_a,info)
      if (info == psb_success_) call psb_geaxpby(denom,wtmp,dzero,w,desc_a,info)
      if (info == psb_success_) call psb_geaxpby(phi,w,done,x,desc_a,info)
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      if (do_cond) then
        gmax = max(gmax,abs(gamma))
        gmin = min(gmin,abs(gamma))
        zeta = rhs1/gamma
        ynorm2 = ynorm2 + abs(zeta)**2
        rhs1 = rhs2 - delta*zeta
        rhs2 = -epsln*zeta
      end if

      select case(istop_)
      case(psb_istop_ani_)
         ! Compute true residual only for the ANI stopping criterion.
         call psb_geaxpby(done,b,dzero,res,desc_a,info)
         if (info == psb_success_) call psb_spmm(-done,a,x,done,res,desc_a,info,work=aux)
         if (info == psb_success_) rni = psb_geamax(res,desc_a,info)
         if (info == psb_success_) xni = psb_geamax(x,desc_a,info)
         errnum = rni
         errden = ani*xni + bni
      case(psb_istop_bn2_)
         errnum = abs(phibar)
         errden = bn2
      case(psb_istop_rn2_abs_)
         errnum = abs(phibar)
         errden = done
      case(psb_istop_rrn2_)
         errnum = abs(phibar)
         errden = r0n2
      end select
      if (errden == dzero) errden = done
      if (info /= psb_success_) then
         info=psb_err_from_subroutine_non_
         call psb_errpush(info,name)
         goto 9999
      end if

      if (itrace_ > 0) &
           & call log_conv(methdname,me,itx,itrace_,errnum,errden,deps)

      if (errnum <= eps*errden) then
         converged = .true.
         exit
      end if
   end do
   end if

   if ((.not.converged) .and. (itx >= itmax_)) then
      if (debug_level >= psb_debug_ext_)&
           & write(debug_unit,*) me,' ',trim(name),': MINRES reached iteration limit'
   end if

   if (do_cond) then
      if (gmin > 0.0_psb_dpk_) then
         cond = gmax/gmin
      else
         cond = huge(done)
      end if
   end if

   call log_end(methdname,me,itx,itrace_,errnum,errden,deps,err=derr,iter=iter)
   if (present(err)) err = derr

   if (info == psb_success_) call psb_gefree(wwrk,desc_a,info)
   if (info == psb_success_) deallocate(aux,stat=info)
   if (info /= psb_success_) then
      call psb_errpush(info,name)
      goto 9999
   end if


   call psb_erractionrestore(err_act)
   return

9999 call psb_error_handler(err_act)
   return

   end subroutine psb_dminres_vect
