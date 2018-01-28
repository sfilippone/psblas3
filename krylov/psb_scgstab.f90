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
! File:  psb_scgstab.f90
!
! Subroutine: psb_scgstab
!    This subroutine implements the BiCG Stabilized method.
!
!    
! Arguments:
!
!    a      -  type(psb_sspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_sprec_type)       Input: preconditioner
!    b      -  real,dimension(:)       Input: vector containing the
!                                         right hand side B
!    x      -  real,dimension(:)       Input/Output: vector containing the
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
Subroutine psb_scgstab_vect(a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_s_krylov_conv_mod
  use psb_krylov_mod
  implicit none
  type(psb_sspmat_type), intent(in)  :: a
  class(psb_sprec_type), Intent(inout)  :: prec 
  Type(psb_desc_type), Intent(in)    :: desc_a
  type(psb_s_vect_type), Intent(inout)   :: b
  type(psb_s_vect_type), Intent(inout)   :: x
  Real(psb_spk_), Intent(in)         :: eps
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_), Optional, Intent(in)      :: itmax, itrace, istop
  integer(psb_ipk_), Optional, Intent(out)     :: iter
  Real(psb_spk_), Optional, Intent(out) :: err
! =   Local data
  real(psb_spk_), allocatable, target   :: aux(:),wwrk(:,:)
  type(psb_s_vect_type) :: q, r, p, v, s, t, z, f

  integer(psb_ipk_) :: itmax_, naux, mglob, it,itrace_,&
       & n_row, n_col
  integer(psb_ipk_) :: debug_level, debug_unit
  Logical, Parameter :: exchange=.True., noexchange=.False., debug1 = .False.
  integer(psb_ipk_), Parameter :: irmax = 8
  integer(psb_ipk_) :: itx, err_act, i
  integer(psb_ipk_) :: istop_
  integer(psb_ipk_) :: ictxt, np, me
  real(psb_dpk_)     :: derr
  real(psb_spk_)     :: alpha, beta, rho, rho_old, sigma, omega, tau
  type(psb_itconv_type) :: stopdat

  character(len=20)           :: name
  character(len=*), parameter :: methdname='BiCGStab'

  info = psb_success_
  name = 'psb_scgstab'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit,*) me,' ',trim(name),': from psb_info',np
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(b%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  mglob = desc_a%get_global_rows()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()

  If (Present(istop)) Then 
    istop_ = istop 
  Else
    istop_ = 2
  Endif
  !
  !  ISTOP_ = 1:  Normwise backward error, infinity norm 
  !  ISTOP_ = 2:  ||r||/||b||   norm 2 
  !
  ! =  if (.not.same_type_as(x,b)) then 
  ! =    write(0,*) 'Warning: different dynamic types for X and B '
  ! =  end if

  call psb_chkvect(mglob,ione,x%get_nrows(),ione,ione,desc_a,info)
  if (psb_errstatus_fatal()) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_chkvect on X')
    goto 9999
  end if
  call psb_chkvect(mglob,ione,b%get_nrows(),ione,ione,desc_a,info)
  if (psb_errstatus_fatal()) then
    info=psb_err_from_subroutine_     
    call psb_errpush(info,name,a_err='psb_chkvect on B')
    goto 9999
  end if
  naux=6*n_col 
  if (info == psb_success_) allocate(aux(naux),stat=info)
  if (info /= psb_success_) then  
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  End If


  call psb_geasb(q,desc_a,info,mold=x%v,scratch=.true.) 
  call psb_geasb(r,desc_a,info,mold=x%v,scratch=.true.)
  call psb_geasb(p,desc_a,info,mold=x%v,scratch=.true.)
  call psb_geasb(v,desc_a,info,mold=x%v,scratch=.true.)
  call psb_geasb(s,desc_a,info,mold=x%v,scratch=.true.)
  call psb_geasb(t,desc_a,info,mold=x%v,scratch=.true.)
  call psb_geasb(z,desc_a,info,mold=x%v,scratch=.true.)
  call psb_geasb(f,desc_a,info,mold=x%v,scratch=.true.)


  if (psb_errstatus_fatal()) then 
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  End If

  If (Present(itmax)) Then 
    itmax_ = itmax
  Else
    itmax_ = 1000
  Endif

  If (Present(itrace)) Then
    itrace_ = itrace
  Else
    itrace_ = 0
  End If

  itx   = 0
  call psb_init_conv(methdname,istop_,itrace_,itmax_,a,x,b,eps,desc_a,stopdat,info)
  if (psb_errstatus_fatal()) Then 
    call psb_errpush(psb_err_from_subroutine_non_,name)
    goto 9999
  End If

  restart: Do 

    if (itx >= itmax_) exit restart  

    it = 0      
    call psb_geaxpby(sone,b,szero,r,desc_a,info)

    call psb_spmm(-sone,a,x,sone,r,desc_a,info,work=aux)
    call psb_geaxpby(sone,r,szero,q,desc_a,info)

    ! Perhaps we already satisfy the convergence criterion...
    if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart

    if (psb_errstatus_fatal()) then 
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='Init residual chk')
      goto 9999
    end if


    rho = szero

    iteration:  Do 
      it   = it + 1
      itx = itx + 1

      if (debug_level >= psb_debug_ext_)&
           & write(debug_unit,*) me,' ',trim(name),&
           & ' Iteration: ',itx

      rho_old = rho    
      rho     = psb_gedot(q,r,desc_a,info)

      if (debug_level >= psb_debug_ext_) then 
        call q%sync()
        call r%sync()
      end if
      if (debug_level >= psb_debug_ext_)&
           & write(debug_unit,*) me,' ',trim(name),&
           & ' Rho: ',rho

      if (rho == szero) then
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' Iteration breakdown R',rho
        exit iteration
      endif

      if (it == 1) then
        call psb_geaxpby(sone,r,szero,p,desc_a,info)
      else
        beta = (rho/rho_old)*(alpha/omega)
        call psb_geaxpby(-omega,v,sone,p,desc_a,info)
        call psb_geaxpby(sone,r,beta,p,desc_a,info)
      End If

      call prec%apply(p,f,desc_a,info,work=aux)

      call psb_spmm(sone,a,f,szero,v,desc_a,info,&
           & work=aux)


      sigma = psb_gedot(q,v,desc_a,info)

      if (debug_level >= psb_debug_ext_) then 
        call q%sync()
        call v%sync()
      end if
      if (debug_level >= psb_debug_ext_)&
           & write(debug_unit,*) me,' ',trim(name),&
           & ' Sigma: ',sigma

      if (sigma == szero) then
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' Iteration breakdown S1', sigma
        exit iteration
      endif

      alpha = rho/sigma
      call psb_geaxpby(sone,r,szero,s,desc_a,info)
      call psb_geaxpby(-alpha,v,sone,s,desc_a,info)
      if (debug_level >= psb_debug_ext_)&
           & write(debug_unit,*) me,' ',trim(name),&
           & ' alpha: ',alpha


      if (psb_errstatus_fatal()) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_geaxpby')
        goto 9999
      end if


      call prec%apply(s,z,desc_a,info,work=aux)
      Call psb_spmm(sone,a,z,szero,t,desc_a,info,work=aux)

      if(psb_errstatus_fatal()) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='precaply/spmm')
        goto 9999
      end if

      sigma = psb_gedot(t,t,desc_a,info)
      if (sigma == szero) then
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' Iteration breakdown S2', sigma
        exit iteration
      endif

      tau   = psb_gedot(t,s,desc_a,info)
      omega = tau/sigma

      if (debug_level >= psb_debug_ext_) then 
        call t%sync()
        call s%sync()
      end if
      if (debug_level >= psb_debug_ext_)&
           & write(debug_unit,*) me,' ',trim(name),&
           & ' sigma, tau, omega: ',sigma, tau, omega

      if (omega == szero) then
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' Iteration breakdown O',omega
        exit iteration
      endif

      call psb_geaxpby(alpha,f,sone,x,desc_a,info)
      call psb_geaxpby(omega,z,sone,x,desc_a,info)
      call psb_geaxpby(sone,s,szero,r,desc_a,info)
      call psb_geaxpby(-omega,t,sone,r,desc_a,info)
      if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart

      if (psb_errstatus_fatal()) Then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='X/R update ')
        goto 9999
      End If

    end do iteration
  end do restart

  call psb_end_conv(methdname,itx,desc_a,stopdat,info,derr,iter)
  if (present(err)) err = derr

  deallocate(aux,stat=info)

  call x%sync()
  call psb_gefree(q,desc_a,info) 
  call psb_gefree(r,desc_a,info)
  call psb_gefree(p,desc_a,info)
  call psb_gefree(v,desc_a,info)
  call psb_gefree(s,desc_a,info)
  call psb_gefree(t,desc_a,info)
  call psb_gefree(z,desc_a,info)
  call psb_gefree(f,desc_a,info)

  if(psb_errstatus_fatal()) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

End Subroutine psb_scgstab_vect
