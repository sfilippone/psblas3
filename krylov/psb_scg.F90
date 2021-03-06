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
! File:  psb_scg.f90
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
! File:  psb_scg.f90
!
! Subroutine: psb_scg
!    This subroutine implements the Conjugate Gradient method.
!
!
! Arguments:
!
!    a      -  type(psb_sspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_sprec_type)       Input: preconditioner
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
subroutine psb_scg_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,istop,cond)
  use psb_base_mod
  use psb_prec_mod
  use psb_s_krylov_conv_mod
  use psb_krylov_mod
  implicit none
  type(psb_sspmat_type), intent(in)    :: a
  Type(psb_desc_type), Intent(in)      :: desc_a
  class(psb_sprec_type), intent(inout) :: prec
  type(psb_s_vect_type), Intent(inout) :: b
  type(psb_s_vect_type), Intent(inout) :: x
  Real(psb_spk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  Real(psb_spk_), Optional, Intent(out) :: err,cond
! =   Local data
  real(psb_spk_), allocatable, target   :: aux(:),td(:),tu(:),eig(:),ewrk(:)
  integer(psb_mpk_), allocatable :: ibl(:), ispl(:), iwrk(:)
  type(psb_s_vect_type), allocatable, target :: wwrk(:)
  type(psb_s_vect_type), pointer  :: q, p, r, z, w
  real(psb_spk_)   :: alpha, beta, rho, rho_old, sigma,alpha_old,beta_old
  integer(psb_ipk_) :: itmax_, istop_, naux, it, itx, itrace_,&
       &  n_col, n_row,err_act, ieg,nspl, istebz
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: debug_level, debug_unit
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me
  real(psb_dpk_)     :: derr  
  type(psb_itconv_type)       :: stopdat
  logical                     :: do_cond
  character(len=20)           :: name
  character(len=*), parameter :: methdname='CG'

  info = psb_success_
  name = 'psb_scg'
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
  if (info == psb_success_) call psb_geall(wwrk,desc_a,info,n=5_psb_ipk_)
  if (info == psb_success_) call psb_geasb(wwrk,desc_a,info,mold=x%v,scratch=.true.)  
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if

  p  => wwrk(1)
  q  => wwrk(2)
  r  => wwrk(3)
  z  => wwrk(4) 
  w  => wwrk(5)


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

  do_cond=present(cond)
  if (do_cond) then 
    istebz = 0
    allocate(td(itmax_),tu(itmax_), eig(itmax_),&
         & ibl(itmax_),ispl(itmax_),iwrk(3*itmax_),ewrk(4*itmax_),&
         & stat=info)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if
  end if
  itx=0
  alpha = szero

  restart: do 
! =   
! =    r0 = b-Ax0
! =   
    if (itx>= itmax_) exit restart 

    it = 0
    call psb_geaxpby(sone,b,szero,r,desc_a,info)
    if (info == psb_success_) call psb_spmm(-sone,a,x,sone,r,desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

    rho = szero
    
    call psb_init_conv(methdname,istop_,itrace_,itmax_,a,x,b,eps,desc_a,stopdat,info)
    if (info /= psb_success_) Then 
      call psb_errpush(psb_err_from_subroutine_non_,name)
      goto 9999
    End If

    iteration:  do 

      it   = it + 1
      itx = itx + 1

      call prec%apply(r,z,desc_a,info,work=aux)
      rho_old = rho
      rho     = psb_gedot(r,z,desc_a,info)

      if (it == 1) then
        call psb_geaxpby(sone,z,szero,p,desc_a,info)
      else
        if (rho_old == szero) then
          if (debug_level >= psb_debug_ext_)&
               & write(debug_unit,*) me,' ',trim(name),&
               & ': CG Iteration breakdown rho'
          exit iteration
        endif
        beta = rho/rho_old
        call psb_geaxpby(sone,z,beta,p,desc_a,info)
      end if

      call psb_spmm(sone,a,p,szero,q,desc_a,info,work=aux)
      sigma = psb_gedot(p,q,desc_a,info)
      if (sigma == szero) then
          if (debug_level >= psb_debug_ext_)&
               & write(debug_unit,*) me,' ',trim(name),&
               & ': CG Iteration breakdown sigma'
        exit iteration
      endif
      alpha_old = alpha
      alpha = rho/sigma

      if (do_cond) then 
        istebz = istebz + 1
        if (istebz == 1) then 
          td(istebz) = sone/alpha
        else 
          td(istebz) = sone/alpha + beta/alpha_old
          tu(istebz-1) = sqrt(beta)/alpha_old
        end if
      end if


      call psb_geaxpby(alpha,p,sone,x,desc_a,info)
      call psb_geaxpby(-alpha,q,sone,r,desc_a,info)

      if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
      if (info /= psb_success_) Then 
        call psb_errpush(psb_err_from_subroutine_non_,name)
        goto 9999
      End If

    end do iteration
  end do restart
  if (do_cond) then 
    if (me == psb_root_) then 
#if defined(HAVE_LAPACK) 
      call sstebz('A','E',istebz,szero,szero,0,0,-sone,td,tu,&
           & ieg,nspl,eig,ibl,ispl,ewrk,iwrk,info)
      if (info < 0) then 
        call psb_errpush(psb_err_from_subroutine_ai_,name,&
             & a_err='sstebz',i_err=(/info/))
        info=psb_err_from_subroutine_ai_
        goto 9999
      end if
      cond = eig(ieg)/eig(1)
#else 
      cond = szero
#endif
      info=psb_success_
    end if
    call psb_bcast(ctxt,cond)
  end if


  call psb_end_conv(methdname,itx,desc_a,stopdat,info,derr,iter)
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

end subroutine psb_scg_vect

