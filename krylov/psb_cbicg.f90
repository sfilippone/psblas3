!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006, 2010, 2015, 2017
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
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
!   C             SIAM, 1993                                               C
!   C                                                                      C
!   C                                                                      C
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! File:  psb_cbicg.f90
!
! Subroutine: psb_cbicg
!    This subroutine implements the BiCG method.
!
! Arguments:
!
!    a      -  type(psb_cspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_cprec_type)       Input: preconditioner
!    b(:)   -  complex                    Input: vector containing the
!                                         right hand side B
!    x(:)   -  complex                    Input/Output: vector containing the
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

subroutine psb_cbicg_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_c_krylov_conv_mod
  use psb_krylov_mod
  implicit none
  type(psb_cspmat_type), intent(in)    :: a
  type(psb_desc_type), intent(in)      :: desc_a
  class(psb_cprec_type), intent(inout) :: prec
  type(psb_c_vect_type), Intent(inout) :: b
  type(psb_c_vect_type), Intent(inout) :: x
  real(psb_spk_), intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), optional, intent(in)        :: itmax, itrace, istop
  integer(psb_ipk_), optional, intent(out)       :: iter
  real(psb_spk_), optional, intent(out) :: err
! !$   local data
  complex(psb_spk_), allocatable, target   :: aux(:)
  type(psb_c_vect_type), allocatable, target :: wwrk(:)
  type(psb_c_vect_type), pointer  :: ww, q, r, p,&
       & zt, pt, z, rt, qt
  integer(psb_ipk_) :: int_err(5)
  integer(psb_ipk_) :: itmax_, naux, mglob, it, itrace_,&
       & n_row, n_col, istop_, err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  logical, parameter :: exchange=.true., noexchange=.false.  
  integer(psb_ipk_), parameter :: irmax = 8
  integer(psb_ipk_) :: itx
  integer(psb_ipk_) :: ictxt, np, me
  complex(psb_spk_)     :: alpha, beta, rho, rho_old, sigma
  real(psb_dpk_)     :: derr  
  type(psb_itconv_type) :: stopdat
  character(len=20)           :: name,ch_err
  character(len=*), parameter :: methdname='BiCG'

  info = psb_success_
  name = 'psb_cbicg'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit,*) me,' ',trim(name),': from psb_info',np

  mglob = desc_a%get_global_rows()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()

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

  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 2
  endif
  !
  !  istop_ = 1:  normwise backward error, infinity norm 
  !  istop_ = 2:  ||r||/||b||   norm 2 
  !

  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
    info=psb_err_invalid_istop_
    int_err=istop_
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  call psb_chkvect(mglob,ione,x%get_nrows(),ione,ione,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_chkvect on X')
    goto 9999
  end if
  call psb_chkvect(mglob,ione,b%get_nrows(),ione,ione,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_    
    call psb_errpush(info,name,a_err='psb_chkvect on B')
    goto 9999
  end if


  naux=4*n_col 

  allocate(aux(naux),stat=info)
  if (info == psb_success_) call psb_geall(wwrk,desc_a,info,n=9_psb_ipk_)
  if (info == psb_success_) call psb_geasb(wwrk,desc_a,info,mold=x%v)  
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_non_
    ch_err='psb_asb'
    err=info
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  q  => wwrk(1)
  qt => wwrk(2)
  r  => wwrk(3)
  rt => wwrk(4)
  p  => wwrk(5)
  pt => wwrk(6)
  z  => wwrk(7)
  zt => wwrk(8)
  ww => wwrk(9)

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

  itx   = 0


  call psb_init_conv(methdname,istop_,itrace_,itmax_,a,x,b,eps,desc_a,stopdat,info)
  if (info /= psb_success_) Then 
     call psb_errpush(psb_err_from_subroutine_non_,name)
     goto 9999
  End If

  restart: do 
    ! !$   
    ! !$   r0 = b-ax0
    ! !$ 
    if (itx >= itmax_) exit restart  
    it = 0      
    call psb_geaxpby(cone,b,czero,r,desc_a,info)
    if (info == psb_success_) call psb_spmm(-cone,a,x,cone,r,desc_a,info,work=aux)
    if (debug_level >= psb_debug_ext_)&
         & write(debug_unit,*) me,' ',trim(name),' Done spmm',info
    if (info == psb_success_) call psb_geaxpby(cone,r,czero,rt,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

    rho = czero
    
    ! Perhaps we already satisfy the convergence criterion...
    if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
    if (info /= psb_success_) Then 
      call psb_errpush(psb_err_from_subroutine_non_,name)
      goto 9999
    End If

    iteration:  do 
      it   = it + 1
      itx = itx + 1

      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),'iteration: ',itx

      call prec%apply(r,z,desc_a,info,work=aux)
      if (info == psb_success_) call prec%apply(rt,zt,desc_a,info,trans='c',work=aux)

      rho_old = rho    
      rho = psb_gedot(rt,z,desc_a,info)
      if (rho == czero) then
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' iteration breakdown r',rho
        exit iteration
      endif

      if (it == 1) then
        call psb_geaxpby(cone,z,czero,p,desc_a,info)
        call psb_geaxpby(cone,zt,czero,pt,desc_a,info)
      else
        beta = (rho/rho_old)
        call psb_geaxpby(cone,z,beta,p,desc_a,info)
        call psb_geaxpby(cone,zt,beta,pt,desc_a,info)
      end if

      call psb_spmm(cone,a,p,czero,q,desc_a,info,&
           & work=aux)
      call psb_spmm(cone,a,pt,czero,qt,desc_a,info,&
           & work=aux,trans='c')

      sigma = psb_gedot(pt,q,desc_a,info)
      if (sigma == czero) then
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' iteration breakdown s1', sigma
        exit iteration
      endif

      alpha = rho/sigma


      call psb_geaxpby(alpha,p,cone,x,desc_a,info)
      call psb_geaxpby(-alpha,q,cone,r,desc_a,info)
      call psb_geaxpby(-alpha,qt,cone,rt,desc_a,info)

      if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
      if (info /= psb_success_) Then 
        call psb_errpush(psb_err_from_subroutine_non_,name)
        goto 9999
      End If

    end do iteration
  end do restart

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

end subroutine psb_cbicg_vect

