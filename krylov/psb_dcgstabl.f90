!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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
!!$ C             SIAM, 1993                                               C
!!$ C                                                                      C
!!$ C                                                                      C
!!$ C         [5] G. Sleijpen, D. Fokkema                                  C
!!$ C             BICGSTAB(L) for linear equations involving unsymmetric   C
!!$ C             matrices with complex spectrum                           C
!!$ C             Electronic Trans. on Numer. Analysis, Vol. 1, pp. 11-32, C
!!$ C             Sep. 1993                                                C
!!$ C                                                                      C
!!$ C                                                                      C
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_dcgstabl.f90
!
! Subroutine: psb_dcgstabl
!   Implements the BICGSTAB(L) method
!
!
! Arguments:
!
!    a      -  type(psb_dspmat_type)      Input: sparse matrix containing A.
!    prec   -  type(psb_dprec_type)       Input: preconditioner
!    b      -  real,dimension(:)          Input: vector containing the
!                                         right hand side B
!    x      -  real,dimension(:)          Input/Output: vector containing the
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
!    irst   -  integer(optional)          Input: restart parameter L 
!
!
!
Subroutine psb_dcgstabl(a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod, psb_protect_name => psb_dcgstabl
  implicit none

!!$  parameters 
  Type(psb_dspmat_type), Intent(in)  :: a
  Type(psb_dprec_type), Intent(in)   :: prec 
  Type(psb_desc_type), Intent(in)    :: desc_a
  Real(Kind(1.d0)), Intent(in)       :: b(:)
  Real(Kind(1.d0)), Intent(inout)    :: x(:)
  Real(Kind(1.d0)), Intent(in)       :: eps
  integer, intent(out)               :: info
  Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
  Integer, Optional, Intent(out)     :: iter
  Real(Kind(1.d0)), Optional, Intent(out) :: err
!!$   local data
  Real(Kind(1.d0)), allocatable, target   :: aux(:),wwrk(:,:),uh(:,:), rh(:,:)
  Real(Kind(1.d0)), Pointer  :: ww(:), q(:), r(:), rt0(:), p(:), v(:), &
       & s(:), t(:), z(:), f(:), gamma(:), gamma1(:), gamma2(:), taum(:,:), sigma(:)

  Integer       :: litmax, naux, mglob, it, itrace_,&
       & np,me, n_row, n_col, nl, err_act
  Logical, Parameter :: exchange=.True., noexchange=.False.  
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ictxt,istop_,j, int_err(5)
  integer            :: debug_level, debug_unit
  Real(Kind(1.d0)) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,bn2,& 
       & omega
  real(kind(1.d0))   :: errnum, errden
  character(len=20)             :: name
  character(len=*), parameter   :: methdname='BiCGStab(L)'

  info = 0
  name = 'psb_dcgstabl'
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

  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 1
  endif
!
!  ISTOP_ = 1:  Normwise backward error, infinity norm 
!  ISTOP_ = 2:  ||r||/||b||   norm 2 
!

  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
    info=5001
    int_err=istop_
    err=info
    call psb_errpush(info,name,i_err=int_err)
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
  
  if (present(irst)) then
    nl = irst
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & 'present: irst: ',irst,nl
  else
    nl = 1 
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' not present: irst: ',irst,nl
  endif
  if (nl <=0 ) then 
    info=5001
    int_err(1)=nl
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  call psb_chkvect(mglob,1,size(x,1),1,1,desc_a,info)
  if(info /= 0) then
    info=4010
    call psb_errpush(info,name,a_err='psb_chkvect on X')
    goto 9999
  end if
  call psb_chkvect(mglob,1,size(b,1),1,1,desc_a,info)
  if(info /= 0) then
    info=4010    
    call psb_errpush(info,name,a_err='psb_chkvect on B')
    goto 9999
  end if

  naux=4*n_col 
  allocate(aux(naux),gamma(0:nl),gamma1(nl),&
       &gamma2(nl),taum(nl,nl),sigma(nl), stat=info)

  if (info /= 0) then 
     info=4000
     call psb_errpush(info,name)
     goto 9999
  end if
  if (info == 0) Call psb_geall(wwrk,desc_a,info,n=10)
  if (info == 0) Call psb_geall(uh,desc_a,info,n=nl+1)
  if (info == 0) Call psb_geall(rh,desc_a,info,n=nl+1)
  if (info == 0) Call psb_geasb(wwrk,desc_a,info)  
  if (info == 0) Call psb_geasb(uh,desc_a,info)  
  if (info == 0) Call psb_geasb(rh,desc_a,info)  
  if (info /= 0) then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  end if

  q   => wwrk(:,1)
  r   => wwrk(:,2)
  p   => wwrk(:,3)
  v   => wwrk(:,4)
  f   => wwrk(:,5)
  s   => wwrk(:,6)
  t   => wwrk(:,7)
  z   => wwrk(:,8)
  ww  => wwrk(:,9)
  rt0 => wwrk(:,10)
  
  ! Ensure global coherence for convergence checks.
  call psb_set_coher(ictxt,isvch)

  if (istop_ == 1) then 
    ani = psb_spnrmi(a,desc_a,info)
    bni = psb_geamax(b,desc_a,info)
  else if (istop_ == 2) then 
    bn2 = psb_genrm2(b,desc_a,info)
  endif
  errnum = dzero
  errden = done
  if (info /= 0) then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  end if

  itx   = 0
  restart: do 
!!$   
!!$   r0 = b-ax0
!!$ 
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),' restart: ',itx,it
    if (itx >= litmax) exit restart  
    it = 0      
    call psb_geaxpby(done,b,dzero,r,desc_a,info)
    call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
    
    call psb_precaply(prec,r,desc_a,info)

    call psb_geaxpby(done,r,dzero,rt0,desc_a,info)
    call psb_geaxpby(done,r,dzero,rh(:,0),desc_a,info)
    call psb_geaxpby(dzero,r,dzero,uh(:,0),desc_a,info)
    if (info /= 0) then 
       info=4011 
       call psb_errpush(info,name)
       goto 9999
    end if
   
    rho   = done
    alpha = dzero
    omega = done 

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' on entry to amax: b: ',Size(b)

    if (istop_ == 1) then 
      rni = psb_geamax(r,desc_a,info)
      xni = psb_geamax(x,desc_a,info)
      errnum = rni
      errden = (ani*xni+bni)
    else if (istop_ == 2) then 
      rni = psb_genrm2(r,desc_a,info)
      errnum = rni
      errden = bn2
    endif
    if (info /= 0) Then 
       info=4011 
       call psb_errpush(info,name)
       goto 9999
    end if

    if (errnum <= eps*errden) exit restart
    if (itrace_ > 0) &
         & call log_conv(methdname,me,itx,itrace_,errnum,errden,eps)
     
    iteration:  do 
      it   = it + nl
      itx  = itx + nl
      rho = -omega*rho 

      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' iteration: ',itx, rho,rh(1,0)

      do j = 0, nl -1 
        If (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),'bicg part:  ',j, nl

        rho_old = rho
        rho = psb_gedot(rh(:,j),rt0,desc_a,info)
        if (rho==dzero) then
          if (debug_level >= psb_debug_ext_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' bi-cgstab iteration breakdown r',rho
          exit iteration
        endif

        beta = alpha*rho/rho_old 
        rho_old = rho
        call psb_geaxpby(done,rh(:,0:j),-beta,uh(:,0:j),desc_a,info)
        call psb_spmm(done,a,uh(:,j),dzero,uh(:,j+1),desc_a,info,work=aux)

        call psb_precaply(prec,uh(:,j+1),desc_a,info)

        gamma(j) = psb_gedot(uh(:,j+1),rt0,desc_a,info)

        if (gamma(j)==dzero) then
          if (debug_level >= psb_debug_ext_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' bi-cgstab iteration breakdown s2',gamma(j)
          exit iteration
        endif
        alpha = rho/gamma(j)
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' bicg part: alpha=r/g ',alpha,rho,gamma(j)

        call psb_geaxpby(-alpha,uh(:,1:j+1),done,rh(:,0:j),desc_a,info)        
        call psb_geaxpby(alpha,uh(:,0),done,x,desc_a,info)
        call psb_spmm(done,a,rh(:,j),dzero,rh(:,j+1),desc_a,info,work=aux)

        call psb_precaply(prec,rh(:,j+1),desc_a,info)
                
      enddo
      
      do j=1, nl 
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' mod g-s part:  ',j, nl,rh(1,0)

        do i=1, j-1 
          taum(i,j) = psb_gedot(rh(:,i),rh(:,j),desc_a,info)
          taum(i,j) = taum(i,j)/sigma(i) 
          call psb_geaxpby(-taum(i,j),rh(:,i),done,rh(:,j),desc_a,info)        
        enddo        
        sigma(j)  = psb_gedot(rh(:,j),rh(:,j),desc_a,info)
        gamma1(j) = psb_gedot(rh(:,0),rh(:,j),desc_a,info)
        gamma1(j) = gamma1(j)/sigma(j)
      enddo
      
      gamma(nl) = gamma1(nl) 
      omega     = gamma(nl) 

      do j=nl-1,1,-1
        gamma(j) = gamma1(j)
        do i=j+1,nl
          gamma(j) = gamma(j) - taum(j,i) * gamma(i) 
        enddo
      enddo

      do j=1,nl-1
        gamma2(j) = gamma(j+1)
        do i=j+1,nl-1
          gamma2(j) = gamma2(j) + taum(j,i) * gamma(i+1) 
        enddo
      enddo
      
      call psb_geaxpby(gamma(1),rh(:,0),done,x,desc_a,info)        
      call psb_geaxpby(-gamma1(nl),rh(:,nl),done,rh(:,0),desc_a,info)        
      call psb_geaxpby(-gamma(nl),uh(:,nl),done,uh(:,0),desc_a,info)        

      do j=1, nl-1
        call psb_geaxpby(-gamma(j),uh(:,j),done,uh(:,0),desc_a,info)        
        call psb_geaxpby(gamma2(j),rh(:,j),done,x,desc_a,info)        
        call psb_geaxpby(-gamma1(j),rh(:,j),done,rh(:,0),desc_a,info)        
      enddo
      
      if (istop_ == 1) then 
        rni = psb_geamax(rh(:,0),desc_a,info)
        xni = psb_geamax(x,desc_a,info)
        errnum = rni
        errden = (ani*xni+bni)
      else  if (istop_ == 2) then 
        rni = psb_genrm2(rh(:,0),desc_a,info)
        errnum = rni
        errden = bn2
      endif

      if (errnum <= eps*errden) exit restart
      if (itx >= litmax) exit restart
      if (itrace_ > 0) &
           & call log_conv(methdname,me,itx,itrace_,errnum,errden,eps)
      
    end do iteration
  end do restart

  if (itrace_ > 0) &
       & call log_conv(methdname,me,itx,1,errnum,errden,eps)

  if (present(err)) then 
    if (errden /= dzero) then 
      err = errnum/errden
    else
      err = errnum
    end if
  end if

  if (present(iter)) iter = itx
  if (errnum > eps*errden) &
       & call end_log(methdname,me,itx,errnum,errden,eps)

  deallocate(aux,stat=info)
  if (info == 0) call psb_gefree(wwrk,desc_a,info)
  if (info == 0) call psb_gefree(uh,desc_a,info)
  if (info == 0) call psb_gefree(rh,desc_a,info)
  if(info/=0) then
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

End Subroutine psb_dcgstabl


