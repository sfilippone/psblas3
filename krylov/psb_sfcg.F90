!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File:  psb_sfcg.f90
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
! File:  psb_sfcg.f90
!
! Subroutine: psb_sfcg
!    This subroutine implements the Flexible Conjugate Gradient method.
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
subroutine psb_sfcg_vect(a,prec,b,x,eps,desc_a,info,&
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
  real(psb_spk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  real(psb_spk_), Optional, Intent(out) :: err,cond
! =   Local data
  type(psb_s_vect_type)  :: v, w
  type(psb_s_vect_type), dimension(0:1) ::  d
  real(psb_spk_) :: alpha, tau, tau1, beta, delta
  real(psb_dpk_) :: derr
  integer(psb_ipk_) ::  i, idx, nc2l, it, itx, istop_, itmax_, itrace_
  integer(psb_ipk_) :: n_col, mglob, naux, err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: np, me, ictxt
  real(psb_spk_), allocatable, target   :: aux(:)
  character(len=20)           :: name
  type(psb_itconv_type)       :: stopdat
  character(len=*), parameter :: methdname='FCG'


  info = psb_success_
  name = 'psb_sfcg'
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


  call psb_chkvect(mglob,ione,x%get_nrows(),ione,ione,desc_a,info)
  if (info == psb_success_)&
       & call psb_chkvect(mglob,ione,b%get_nrows(),ione,ione,desc_a,info)
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


  !Assemble w, v

  call psb_geasb(w,&
       & desc_a,info,&
       & scratch=.true.,mold=b%v)
  call psb_geasb(v,&
       & desc_a,info,&
       & scratch=.true.,mold=b%v)

  !Assemble d(0) and d(1)
  call psb_geasb(d(0),&
       & desc_a,info,&
       & scratch=.true.,mold=x%v)
  call psb_geasb(d(1),&
       & desc_a,info,&
       & scratch=.true.,mold=x%v)


    call psb_init_conv(methdname,istop_,itrace_,itmax_,a,b,eps,desc_a,stopdat,info)
  itx=0

  restart: do 
    if (itx>= itmax_) exit restart 

    ! w=b
    call psb_geaxpby(sone,b,szero,w,&
         &   desc_a,info)

      if (psb_errstatus_fatal()) then 
        nc2l = desc_a%get_local_cols()
        info=psb_err_alloc_request_
        call psb_errpush(info,name,i_err=(/2*nc2l,izero,izero,izero,izero/),&
             & a_err='real(psb_spk_)')
        goto 9999      
      end if

   !Compute v = Ax  

    call psb_spmm(sone,a,x,szero,v,desc_a,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_internal_error_,name,&
        & a_err='Error during residue')
      goto 9999
    end if

    !Compute w = -Ax + b

    call psb_geaxpby(-sone, v, sone, w, desc_a, info)   

    !Apply the preconditioner

    idx=0

    call prec%apply(w,d(idx),desc_a,info,work=aux)  

    delta = psb_gedot(d(idx), w, desc_a, info)


    !Loop

    if (psb_check_conv(methdname,itx ,x,w,desc_a,stopdat,info)) exit restart

    if (info /= psb_success_) Then 
      call psb_errpush(psb_err_from_subroutine_non_,name)
      goto 9999
    End If

    iteration: do 

      call psb_spmm(sone,a,d(idx),szero,v,desc_a,info)
      if (info /= psb_success_) then
       call psb_errpush(psb_err_internal_error_,name,&
          & a_err='Error during residue')
       goto 9999
      end if
      tau = psb_gedot(d(idx), v, desc_a, info) 


      alpha = delta/tau
      !Update solution x
      call psb_geaxpby(alpha, d(idx), sone, x, desc_a, info)   
      !Update residual w
      call psb_geaxpby(-alpha, v, sone, w, desc_a, info) 

      itx = itx + 1
      idx=mod(itx ,2)

      call d(idx)%set(szero)  
      call prec%apply(w,d(idx),desc_a,info,work=aux)    

      tau1= psb_gedot(d(idx), v, desc_a, info)
      beta=tau1/tau

      if (idx == 1) then
        call psb_geaxpby(-beta, d(idx - 1), sone, d(idx), desc_a, info)   
      else
        call psb_geaxpby(-beta, d(idx + 1), sone, d(idx), desc_a, info)         
      endif
    
      delta = psb_gedot(w, d(idx), desc_a, info)

      if (psb_check_conv(methdname,itx ,x,w,desc_a,stopdat,info)) exit restart
      if (info /= psb_success_) Then 
        call psb_errpush(psb_err_from_subroutine_non_,name)
        goto 9999
      End If

    end do iteration
  end do restart


  call psb_end_conv(methdname,itx ,desc_a,stopdat,info,derr,iter)
  if (present(err)) err = derr

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if

  return
end subroutine psb_sfcg_vect
