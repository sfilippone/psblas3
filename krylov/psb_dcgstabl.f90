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
!
! Parameters:
!    a       -  type(<psb_dspmat_type>).     The sparse matrix containing A.
!    prec    -  type(<psb_prec_type>).       The data structure containing the preconditioner.
!    b       -  real,dimension(:).           The right hand side.
!    x       -  real,dimension(:).           The vector of unknowns.
!    eps     -  real.                        The error tolerance.
!    desc_a  -  type(<psb_desc_type>).       The communication descriptor.
!    info    -  integer.                     Eventually returns an error code.
!    itmax   -  integer(optional).           The maximum number of iterations.
!    iter    -  integer(optional).           The number of iterations performed.
!    err     -  real(optional).              The error on return.
!    itrace  -  integer(optional).           The unit to write messages onto.
!    istop   -  integer(optional).           The stopping criterium.
!
Subroutine psb_dcgstabl(a,prec,b,x,eps,desc_a,info,&
     &itmax,iter,err,itrace,irst,istop)
  use psb_base_mod
  use psb_prec_mod
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

  Real(Kind(1.d0)) :: rerr
  Integer       :: litmax, naux, m, mglob, it, itrace_,&
       & np,me, n_row, n_col, nl, err_act
  Character     ::diagl, diagu
  Logical, Parameter :: exchange=.True., noexchange=.False.  
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ich, ictxt,istop_,j, int_err(5)
  Logical, Parameter :: debug = .False.
  Real(Kind(1.d0)) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,bn2,& 
       & omega, tau 
  character(len=20)             :: name

  info = 0
  name = 'psb_dcgstabl'
  call psb_erractionsave(err_act)

  If (debug) Write(0,*) 'entering psb_dbicgstabl'
  ictxt = psb_cd_get_context(desc_a)
  Call psb_info(ictxt, me, np)

  If (debug) Write(0,*) 'psb_dbicgstabl: from gridinfo',np,me

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
    write(0,*) 'psb_bicgstabl: invalid istop',istop_ 
    info=5001
    int_err=istop_
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  If (Present(itmax)) Then 
    litmax = itmax
  Else
    litmax = 1000
  Endif

  If (Present(itrace)) Then
     itrace_ = itrace
  Else
     itrace_ = 0
  End If
  
  If (Present(irst)) Then
    nl = irst
    If (debug) Write(0,*) 'present: irst: ',irst,nl
  Else
    nl = 1 
    If (debug) Write(0,*) 'not present: irst: ',irst,nl
  Endif

  naux=4*n_col 
  Allocate(aux(naux),gamma(0:nl),gamma1(nl),&
       &gamma2(nl),taum(nl,nl),sigma(nl), stat=info)

  If (info.Ne.0) Then 
     info=4000
     call psb_errpush(info,name)
     goto 9999
  End If
  if (info == 0) Call psb_geall(wwrk,desc_a,info,n=10)
  if (info == 0) Call psb_geall(uh,desc_a,info,n=nl+1)
  if (info == 0) Call psb_geall(rh,desc_a,info,n=nl+1)
  if (info == 0) Call psb_geasb(wwrk,desc_a,info)  
  if (info == 0) Call psb_geasb(uh,desc_a,info)  
  if (info == 0) Call psb_geasb(rh,desc_a,info)  
  if (info.ne.0) Then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  End If

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
  if (info.ne.0) Then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  End If

  diagl  = 'u'
  diagu  = 'u'
  itx   = 0
  restart: Do 
!!$   
!!$   r0 = b-ax0
!!$ 
    If (debug) Write(0,*) 'restart: ',itx,it
    If (itx.Ge.litmax) Exit restart  
    it = 0      
    Call psb_geaxpby(done,b,dzero,r,desc_a,info)
    Call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
    
    call psb_precaply(prec,r,desc_a,info)

    Call psb_geaxpby(done,r,dzero,rt0,desc_a,info)
    Call psb_geaxpby(done,r,dzero,rh(:,0),desc_a,info)
    Call psb_geaxpby(dzero,r,dzero,uh(:,0),desc_a,info)
    if (info.ne.0) Then 
       info=4011 
       call psb_errpush(info,name)
       goto 9999
    End If
   
    rho   = done
    alpha = dzero
    omega = done 

    If (debug) Write(0,*) 'on entry to amax: b: ',Size(b)

    if (istop_ == 1) then 
      rni = psb_geamax(r,desc_a,info)
      xni = psb_geamax(x,desc_a,info)
      rerr =  rni/(ani*xni+bni)
    else if (istop_ == 2) then 
      rni = psb_genrm2(r,desc_a,info)
      rerr = rni/bn2
    endif
    if (info.ne.0) Then 
       info=4011 
       call psb_errpush(info,name)
       goto 9999
    End If

    If (rerr<=eps) Then 
      Exit restart
    End If
    If (itrace_ > 0) then 
      if ((mod(itx,itrace_)==0).and.(me == 0))&
           & write(*,'(a,i4,3(2x,es10.4))') 'bicgstab(l): ',itx,rerr
    end If
     
    iteration:  Do 
      it   = it + nl
      itx  = itx + nl
      rho = -omega*rho 
      If (debug) Write(0,*) 'iteration: ',itx, rho,rh(1,0)

      Do j = 0, nl -1 
        If (debug) Write(0,*) 'bicg part:  ',j, nl
        rho_old = rho
        rho = psb_gedot(rh(:,j),rt0,desc_a,info)
        If (rho==dzero) Then
          If (debug) Write(0,*) 'bi-cgstab iteration breakdown r',rho
          Exit iteration
        Endif
        beta = alpha*rho/rho_old 
        If (debug) Write(0,*) 'bicg part:  ',alpha,beta,rho,rho_old
        rho_old = rho
        Call psb_geaxpby(done,rh(:,0:j),-beta,uh(:,0:j),desc_a,info)
        If (debug) Write(0,*) 'bicg part:  ',rh(1,0),beta
        Call psb_spmm(done,a,uh(:,j),dzero,uh(:,j+1),desc_a,info,work=aux)

        call psb_precaply(prec,uh(:,j+1),desc_a,info)

        gamma(j) = psb_gedot(uh(:,j+1),rt0,desc_a,info)
        If (gamma(j)==dzero) Then
          If (debug) Write(0,*) 'bi-cgstab iteration breakdown s2',gamma(j)
          Exit iteration
        Endif
        alpha = rho/gamma(j)
        If (debug) Write(0,*) 'bicg part: alpha=r/g ',alpha,rho,gamma(j)

        Call psb_geaxpby(-alpha,uh(:,1:j+1),done,rh(:,0:j),desc_a,info)        
        Call psb_geaxpby(alpha,uh(:,0),done,x,desc_a,info)
        Call psb_spmm(done,a,rh(:,j),dzero,rh(:,j+1),desc_a,info,work=aux)

        call psb_precaply(prec,rh(:,j+1),desc_a,info)
                
      Enddo
      
      Do j=1, nl 
        If (debug) Write(0,*) 'mod g-s part:  ',j, nl,rh(1,0)
        Do i=1, j-1 
          taum(i,j) = psb_gedot(rh(:,i),rh(:,j),desc_a,info)
          taum(i,j) = taum(i,j)/sigma(i) 
          Call psb_geaxpby(-taum(i,j),rh(:,i),done,rh(:,j),desc_a,info)        
        Enddo        
        If (debug) Write(0,*) 'mod g-s part:  dot prod '
        sigma(j)  = psb_gedot(rh(:,j),rh(:,j),desc_a,info)
        gamma1(j) = psb_gedot(rh(:,0),rh(:,j),desc_a,info)
        If (debug) Write(0,*) 'mod g-s part: gamma1 ', &
             &gamma1(j), sigma(j)
        gamma1(j) = gamma1(j)/sigma(j)
      Enddo
      
      gamma(nl) = gamma1(nl) 
      omega     = gamma(nl) 

      Do j=nl-1,1,-1
        gamma(j) = gamma1(j)
        Do i=j+1,nl
          gamma(j) = gamma(j) - taum(j,i) * gamma(i) 
        Enddo
      Enddo
      If (debug) Write(0,*) 'first solve: ', gamma(:)
      
      Do j=1,nl-1
        gamma2(j) = gamma(j+1)
        Do i=j+1,nl-1
          gamma2(j) = gamma2(j) + taum(j,i) * gamma(i+1) 
        Enddo
      Enddo
      If (debug) Write(0,*) 'second solve: ', gamma(:)
      
      Call psb_geaxpby(gamma(1),rh(:,0),done,x,desc_a,info)        
      Call psb_geaxpby(-gamma1(nl),rh(:,nl),done,rh(:,0),desc_a,info)        
      Call psb_geaxpby(-gamma(nl),uh(:,nl),done,uh(:,0),desc_a,info)        

      Do j=1, nl-1
        Call psb_geaxpby(-gamma(j),uh(:,j),done,uh(:,0),desc_a,info)        
        Call psb_geaxpby(gamma2(j),rh(:,j),done,x,desc_a,info)        
        Call psb_geaxpby(-gamma1(j),rh(:,j),done,rh(:,0),desc_a,info)        
      Enddo
      
      if (istop_ == 1) then 
        rni = psb_geamax(rh(:,0),desc_a,info)
        xni = psb_geamax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)
      else  if (istop_ == 2) then 
        rni = psb_genrm2(rh(:,0),desc_a,info)
        rerr = rni/bn2
      endif

      If (rerr<=eps) Then 
        Exit restart
      End If
      If (itx.Ge.litmax) Exit restart

      If (itrace_ > 0) then 
        if ((mod(itx,itrace_)==0).and.(me == 0))&
             & write(*,'(a,i4,3(2x,es10.4))') 'bicgstab(l): ',itx,rerr
      end If
      
    End Do iteration
  End Do restart

  If (itrace_ > 0) then 
    if (me == 0) write(*,'(a,i4,3(2x,es10.4))') 'bicgstab(l): ',itx,rerr
  end If
  If (Present(err)) err=rerr
  If (Present(iter)) iter = itx
  If (rerr>eps) Then
    Write(0,*) 'bi-cgstabl failed to converge to ',eps,&
         & ' in ',itx,' iterations  '
  End If

  Deallocate(aux)
  Call psb_gefree(wwrk,desc_a,info)
  Call psb_gefree(uh,desc_a,info)
  Call psb_gefree(rh,desc_a,info)

  ! restore external global coherence behaviour
  call psb_restore_coher(ictxt,isvch)

  if(info/=0) then
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

End Subroutine psb_dcgstabl


