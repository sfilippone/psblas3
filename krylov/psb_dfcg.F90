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
! File:  psb_dfcg.f90
!!
!! Contributors: Ambra Abdullahi (UNITOV) and Pasqua D’Ambra (IAC-CNR)
!!
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
!   C         [4] Notay, Yvan                                              C
!   C             Flexible Conjugate gradients                             C
!   C             SIAM Journal on Scientific Computing 22(4),              C 
!   C             pp. 1444-1460, 2000                                      C    
!   C                                                                      C
!   C                                                                      C
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_dfcg.f90
!
! Subroutine: psb_dfcg
!    This subroutine implements the Flexible Conjugate Gradient method.
!
!
! Arguments:
!
!    a      -  type(psb_dspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_dprec_type)       Input: preconditioner
!    b(:)   -  real                    Input: vector containing the
!                                         right hand side B
!    x(:)   -  real                    Input/Output: vector containing the
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
subroutine psb_dfcg_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,istop,cond)
  use psb_base_mod
  use psb_prec_mod
  use psb_d_krylov_conv_mod
  use psb_krylov_mod
  implicit none
  type(psb_dspmat_type), intent(in)    :: a
  Type(psb_desc_type), Intent(in)      :: desc_a
  class(psb_dprec_type), intent(inout) :: prec
  type(psb_d_vect_type), Intent(inout) :: b
  type(psb_d_vect_type), Intent(inout) :: x
  real(psb_dpk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  real(psb_dpk_), Optional, Intent(out) :: err,cond
! =   Local data
  type(psb_d_vect_type)  :: v, w, d , q, r
  real(psb_dpk_) :: alpha, beta, delta, gamma, theta
  real(psb_dpk_) :: derr
  integer(psb_ipk_) ::  i, idx, nc2l, it, itx, istop_, itmax_, itrace_
  integer(psb_ipk_) :: n_col, naux, err_act
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: np, me, ictxt
  real(psb_dpk_), allocatable, target   :: aux(:)
  real(psb_dpk_)   :: vres(3)
  character(len=20)           :: name
  type(psb_itconv_type)       :: stopdat
  character(len=*), parameter :: methdname='FCG'


  info = psb_success_
  name = 'psb_dfcg'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)
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
  n_col = desc_a%get_local_cols()

  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 2
  endif


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


  !Assemble w, v, d, q, r, u 
  call psb_geasb(w, desc_a,info,&
       & scratch=.true.,mold=x%v)
  call psb_geasb(v, desc_a,info,&
       & scratch=.true.,mold=x%v)
  call psb_geasb(d, desc_a,info,&
       & scratch=.true.,mold=x%v)
  call psb_geasb(q, desc_a,info,&
       & scratch=.true.,mold=x%v)
  call psb_geasb(r, desc_a,info,&
       & scratch=.true.,mold=x%v)

  call psb_init_conv(methdname,istop_,itrace_,itmax_,&
       & a,x,b,eps,desc_a,stopdat,info)
  itx = 0 

  restart: do 
    if (itx>= itmax_) exit restart 

    ! r=b -Ax
    call psb_geaxpby(done,b,dzero,r, desc_a,info)
    if (info == psb_success_) call psb_spmm(-done,a,x,done,r,desc_a,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error during residual')
      goto 9999
    end if
    

    if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart


    ! Apply the preconditioner v=Pr
    ! Compute w = Av  
    call prec%apply(r,v,desc_a,info,work=aux)  
    if (info == psb_success_) call psb_spmm(done,a,v,dzero,w,desc_a,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
           & a_err='Error during residual')
      goto 9999
    end if


    vres(1) = psb_gedot(r, v, desc_a, info, global = .false.) 
    vres(2) = psb_gedot(w, v, desc_a, info, global = .false.) 


    call psb_sum(ictxt, vres(1:2))

    alpha = vres(1)
    beta  = vres(2)

    ! d = v
    call psb_geaxpby(done, v, dzero, d, desc_a, info)   
    ! q = w
    call psb_geaxpby(done, w, dzero, q, desc_a, info)   

    ! compute delta=beta
    ! then 
    ! x = x + (alpha/delta)*d
    ! r = r - (alpha/delta)*q

    delta = beta
    theta = alpha/delta
    call psb_geaxpby(theta, d, done, x, desc_a, info)   
    call psb_geaxpby(-theta, q, done, r, desc_a, info)   

    iteration: do 

      itx = itx + 1

      if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart

      ! Apply the preconditioner v = Pr
      ! Compute w = Av  
      call prec%apply(r,v,desc_a,info,work=aux)  
      if (info == psb_success_) call psb_spmm(done,a,v,dzero,w,desc_a,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_internal_error_,name,&
             & a_err='Error during residual'); goto 9999
      end if

      vres(1) = psb_gedot(r, v, desc_a, info, global = .false.) 
      vres(2) = psb_gedot(w, v, desc_a, info, global = .false.) 
      vres(3) = psb_gedot(q, v, desc_a, info, global = .false.) 

      call psb_sum(ictxt, vres(1:3))

      alpha = vres(1)
      beta  = vres(2)
      gamma = vres(3)

      ! Compute d = v-(gamma/delta)*d
      !         q = w-(gamma/delta)*q
      theta= gamma/delta
      call psb_geaxpby(done, v, -theta, d, desc_a, info)   
      call psb_geaxpby(done, w, -theta, q , desc_a, info)   

      ! update delta
      delta = beta - (gamma*gamma)/delta

      ! update u and r      
      ! u = u + (alpha/delta)*d
      ! r = r - (alpha/delta)*q
      theta= alpha/delta
      call psb_geaxpby(theta, d, done, x, desc_a, info)   
      call psb_geaxpby(-theta, q, done, r, desc_a, info)   

    end do iteration
  end do restart


  call psb_end_conv(methdname,itx ,desc_a,stopdat,info,derr,iter)
  if (present(err)) err = derr
  return
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if

  return
end subroutine psb_dfcg_vect
