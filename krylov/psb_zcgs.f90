!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$ C                                                                      C
!!$ C  References:                                                         C
!!$ C          [1] Duff, I., Marrone, M., Radicati, G., and Vittoli, C.    C
!!$ C              Level 3 basic linear algebra subprograms for sparse     C
!!$ C              matrices: a user level interface                        C
!!$ C              ACM Trans. Math. Softw., 23(3), 379-401, 1997.          C
!!$ C                                                                      C
!!$ C                                                                      C
!!$ C         [2]  S. Filippone, M. Colajanni                              C
!!$ C              PSBLAS: A library for parallel linear algebra           C
!!$ C              computation on sparse matrices                          C
!!$ C              ACM Trans. on Math. Softw., 26(4), 527-550, Dec. 2000.  C
!!$ C                                                                      C
!!$ C         [3] M. Arioli, I. Duff, M. Ruiz                              C
!!$ C             Stopping criteria for iterative solvers                  C
!!$ C             SIAM J. Matrix Anal. Appl., Vol. 13, pp. 138-144, 1992   C
!!$ C                                                                      C
!!$ C                                                                      C
!!$ C         [4] R. Barrett et al                                         C
!!$ C             Templates for the solution of linear systems             C
!!$ C             SIAM, 1993                                          
!!$ C                                                                      C
!!$ C                                                                      C
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_zcgs.f90
!
! Subroutine: psb_zcgs
!    Implements the Conjugate Gradient Squared method.
!    
! Arguments:
!
!    a      -  type(psb_zspmat_type)      Input: sparse matrix containing A.
!    prec   -  type(psb_zprec_type)       Input: preconditioner
!    b      -  complex,dimension(:)       Input: vector containing the
!                                         right hand side B
!    x      -  complex,dimension(:)       Input/Output: vector containing the
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
!                                         1: err =  |r|/|b|; here the iteration is
!                                            stopped when  |r| <= eps * |b|
!                                         2: err =  |r|/(|a||x|+|b|);  here the iteration is
!                                            stopped when  |r| <= eps * (|a||x|+|b|)
!                                         where r is the (preconditioned, recursive
!                                         estimate of) residual. 
!
Subroutine psb_zcgs(a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod, psb_protect_name => psb_zcgs
  implicit none

!!$  parameters 
  Type(psb_zspmat_type), Intent(in)  :: a
  Type(psb_desc_type), Intent(in)    :: desc_a 
  Type(psb_zprec_type), Intent(in)   :: prec 
  Complex(Kind(1.d0)), Intent(in)       :: b(:)
  Complex(Kind(1.d0)), Intent(inout)    :: x(:)
  Real(Kind(1.d0)), Intent(in)       :: eps
  integer, intent(out)               :: info
  Integer, Optional, Intent(in)      :: itmax, itrace,istop
  Integer, Optional, Intent(out)     :: iter
  Real(Kind(1.d0)), Optional, Intent(out) :: err
!!$   local data
  Complex(Kind(1.d0)), allocatable, target   :: aux(:),wwrk(:,:)
  Complex(Kind(1.d0)), Pointer  :: ww(:), q(:),&
       & r(:), p(:), v(:), s(:), z(:), f(:), rt(:),qt(:),uv(:)
  Integer       :: itmax_, naux, mglob, it, itrace_,int_err(5),&
       & np,me, n_row, n_col,istop_, err_act
  Integer            :: itx, isvch, ictxt
  integer            :: debug_level, debug_unit
  complex(Kind(1.d0))   :: alpha, beta, rho, rho_old, sigma 
  type(psb_itconv_type) :: stopdat
  character(len=20)           :: name
  character(len=*), parameter :: methdname='CGS'

  info = 0
  name = 'psb_zcgs'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = psb_cd_get_context(desc_a)
  Call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit,*) me,' ',trim(name),': from psb_info',np

  mglob = psb_cd_get_global_rows(desc_a)
  n_row = psb_cd_get_local_rows(desc_a)
  n_col = psb_cd_get_local_cols(desc_a)

  If (Present(istop)) Then 
    istop_ = istop 
  Else
    istop_ = 1
  Endif

  call psb_chkvect(mglob,1,size(x,1),1,1,desc_a,info)
  if (info == 0) call psb_chkvect(mglob,1,size(b,1),1,1,desc_a,info)
  if(info /= 0) then
    info=4010    
    call psb_errpush(info,name,a_err='psb_chkvect on X/B')
    goto 9999
  end if

  naux=4*n_col 
  Allocate(aux(naux),stat=info)
  if (info == 0) Call psb_geall(wwrk,desc_a,info,n=11)
  if (info == 0) Call psb_geasb(wwrk,desc_a,info)  
  if (info /= 0) Then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  End If

  q  => wwrk(:,1)
  qt => wwrk(:,2)
  r  => wwrk(:,3)
  rt => wwrk(:,4)
  p  => wwrk(:,5)
  v  => wwrk(:,6)
  uv => wwrk(:,7)
  z  => wwrk(:,8)
  f  => wwrk(:,9)
  s  => wwrk(:,10)
  ww => wwrk(:,11)


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

  ! Ensure global coherence for convergence checks.
  call psb_set_coher(ictxt,isvch)
  
  itx   = 0

  call psb_init_conv(methdname,istop_,itrace_,itmax_,a,b,eps,desc_a,stopdat,info)
  if (info /= 0) Then 
     call psb_errpush(4011,name)
     goto 9999
  End If

  restart: Do 
!!$
!!$   r0 = b-ax0
!!$ 
    if (itx >= itmax_) exit restart  
    it = 0      
    call psb_geaxpby(zone,b,zzero,r,desc_a,info)
    if (info == 0) call psb_spmm(-zone,a,x,zone,r,desc_a,info,work=aux)
    if (info == 0) call psb_geaxpby(zone,r,zzero,rt,desc_a,info)
    if (info/=0) then
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if
    

    ! Perhaps we already satisfy the convergence criterion...
    if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
    if (info /= 0) Then 
      call psb_errpush(4011,name)
      goto 9999
    End If

    rho = zzero

    iteration:  do 
      it   = it + 1
      itx = itx + 1
      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),'iteration: ',itx

      rho_old = rho    
      rho = psb_gedot(rt,r,desc_a,info)

      if (rho==zzero) then
         if (debug_level >= psb_debug_ext_) &
              & write(debug_unit,*) me,' ',trim(name),&
              & ' iteration breakdown r',rho
        exit iteration
      endif

      if (it==1) then
        call psb_geaxpby(zone,r,zzero,uv,desc_a,info)
        if (info == 0) call psb_geaxpby(zone,r,zzero,p,desc_a,info)
      else
        beta = (rho/rho_old)
        call psb_geaxpby(zone,r,zzero,uv,desc_a,info)
        if (info == 0) call psb_geaxpby(beta,q,zone,uv,desc_a,info)
        if (info == 0) call psb_geaxpby(zone,q,beta,p,desc_a,info)
        if (info == 0) call psb_geaxpby(zone,uv,beta,p,desc_a,info)
      end if

      if (info == 0) call psb_precaply(prec,p,f,desc_a,info,work=aux)

      if (info == 0) call psb_spmm(zone,a,f,zzero,v,desc_a,info,&
           & work=aux)
      
      if (info /= 0) then
         call psb_errpush(4010,name,a_err='First loop part ')
         goto 9999
      end if
      
      sigma = psb_gedot(rt,v,desc_a,info)
      if (sigma==zzero) then
         if (debug_level >= psb_debug_ext_) &
              & write(debug_unit,*) me,' ',trim(name),&
              & ' iteration breakdown s1', sigma
         exit iteration
      endif
      
      alpha = rho/sigma

      if (info == 0) call psb_geaxpby(zone,uv,zzero,q,desc_a,info)
      if (info == 0) call psb_geaxpby(-alpha,v,zone,q,desc_a,info)
      if (info == 0) call psb_geaxpby(zone,uv,zzero,s,desc_a,info)
      if (info == 0) call psb_geaxpby(zone,q,zone,s,desc_a,info)
      
      if (info == 0) call psb_precaply(prec,s,z,desc_a,info,work=aux)

      if (info == 0) call psb_geaxpby(alpha,z,zone,x,desc_a,info)

      if (info == 0) call psb_spmm(zone,a,z,zzero,qt,desc_a,info,&
           & work=aux)
      
      if (info == 0) call psb_geaxpby(-alpha,qt,zone,r,desc_a,info)
      
      if (info /= 0) then
         call psb_errpush(4010,name,a_err='X update ')
         goto 9999
      end if

      if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
      if (info /= 0) Then 
        call psb_errpush(4011,name)
        goto 9999
      End If

    end do iteration
  end do restart

  call psb_end_conv(methdname,itx,desc_a,stopdat,info,err,iter)

  deallocate(aux,stat=info)
  if (info == 0) call psb_gefree(wwrk,desc_a,info)
  if (info /= 0) then
    call psb_errpush(info,name)
    goto 9999
  end if

  ! restore external global coherence behaviour
  call psb_restore_coher(ictxt,isvch)
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_zcgs
