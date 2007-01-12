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
!!$ C             SIAM, 1993                                          
!!$ C                                                                      C
!!$ C                                                                      C
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_dcgstab.f90
!
! Subroutine: psb_dcgstab
!    This subroutine implements the CG Stabilized method.
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
Subroutine psb_dcgstab(a,prec,b,x,eps,desc_a,info,&
     &itmax,iter,err,itrace, istop)
  use psb_base_mod
  use psb_prec_mod
  Implicit None
!!$  parameters 
  Type(psb_dspmat_type), Intent(in)  :: a
  Type(psb_dprec_type), Intent(in)   :: prec 
  Type(psb_desc_type), Intent(in)    :: desc_a
  Real(Kind(1.d0)), Intent(in)       :: b(:)
  Real(Kind(1.d0)), Intent(inout)    :: x(:)
  Real(Kind(1.d0)), Intent(in)       :: eps
  integer, intent(out)               :: info
  Integer, Optional, Intent(in)      :: itmax, itrace, istop
  Integer, Optional, Intent(out)     :: iter
  Real(Kind(1.d0)), Optional, Intent(out) :: err
!!$   Local data
  Real(Kind(1.d0)), allocatable, target   :: aux(:),wwrk(:,:)
  Real(Kind(1.d0)), Pointer  :: q(:),&
       & r(:), p(:), v(:), s(:), t(:), z(:), f(:)
  Real(Kind(1.d0)) :: rerr
  Integer       :: litmax, naux, m, mglob, it,itrace_,&
       & np,me, n_row, n_col
  Character     ::diagl, diagu
  Logical, Parameter :: debug = .false.
  Logical, Parameter :: exchange=.True., noexchange=.False., debug1 = .False.
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ictxt, err_act, int_err(5),ii
  Integer            :: istop_
  Real(Kind(1.d0)) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,& 
       & sigma, omega, tau, rn0, bn2
!!$  Integer   istpb, istpe, ifctb, ifcte, imerr, irank, icomm,immb,imme
!!$  Integer mpe_log_get_event_number,mpe_Describe_state,mpe_log_event
  character(len=20)             :: name

  info = 0
  name = 'psb_dcgstab'
  call psb_erractionsave(err_act)

  If (debug) Write(*,*) 'Entering PSB_DCGSTAB',present(istop)
  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (debug) write(*,*) 'PSB_DCGSTAB: From GRIDINFO',np,me

  mglob = psb_cd_get_global_rows(desc_a)
  n_row = psb_cd_get_local_rows(desc_a)
  n_col = psb_cd_get_local_cols(desc_a)

  If (Present(istop)) Then 
    istop_ = istop 
  Else
    istop_ = 1
  Endif
!
!  ISTOP_ = 1:  Normwise backward error, infinity norm 
!  ISTOP_ = 2:  ||r||/||b||   norm 2 
!

!!$  If ((prec%prec < min_prec_).Or.(prec%prec > max_prec_) ) Then
!!$     Write(0,*) 'PSB_CGSTAB: Invalid IPREC',prec%prec
!!$     info=5002
!!$     int_err(1)=prec%prec
!!$     err=info
!!$     call psb_errpush(info,name,i_err=int_err)
!!$     goto 9999
!!$  Endif
  
  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
    write(0,*) 'psb_bicgstab: invalid istop',istop_ 
    info=5001
    int_err(1)=istop_
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  naux=6*n_col 
  allocate(aux(naux),stat=info)
  if (info == 0) call psb_geall(wwrk,desc_a,info,n=8)
  if (info == 0) call psb_geasb(wwrk,desc_a,info)  
  if (info /= 0) then 
     info=4011
     call psb_errpush(info,name)
     goto 9999
  End If

  Q => WWRK(:,1)
  R => WWRK(:,2)
  P => WWRK(:,3)
  V => WWRK(:,4)
  F => WWRK(:,5)
  S => WWRK(:,6)
  T => WWRK(:,7)
  Z => WWRK(:,8)

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
  
  diagl = 'U'
  diagu = 'U'

  ! Ensure global coherence for convergence checks.
  call psb_set_coher(ictxt,isvch)

  itx   = 0

  If (istop_ == 1) Then 
    ani = psb_spnrmi(a,desc_a,info)
    bni = psb_geamax(b,desc_a,info)
  Else If (istop_ == 2) Then 
    bn2 = psb_genrm2(b,desc_a,info)
  Endif
  if (info /= 0) Then 
     info=4011
     call psb_errpush(info,name)
     goto 9999
  End If

  restart: Do 
!!$   
!!$   r0 = b-Ax0
!!$ 
    If (itx >= litmax) Exit restart  
    it = 0      
    Call psb_geaxpby(done,b,dzero,r,desc_a,info)
    Call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
    Call psb_geaxpby(done,r,dzero,q,desc_a,info)
    if (info /= 0) Then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    End If
    
    rho = dzero
    If (debug) Write(*,*) 'On entry to AMAX: B: ',Size(b)
    
!
!   Must always provide norm of R into RNI below for first check on 
!   residual
!
    If (istop_ == 1) Then 
      rni = psb_geamax(r,desc_a,info)
      xni = psb_geamax(x,desc_a,info)
    Else If (istop_ == 2) Then 
      rni = psb_genrm2(r,desc_a,info)
    Endif
    if (info /= 0) Then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    End If

    If (itx == 0) Then 
      rn0 = rni
    End If
    If (rn0 == 0.d0 ) Then 
      If (itrace_ > 0 ) Then 
        If (me == 0) Write(*,*) 'BiCGSTAB: ',itx,rn0
      Endif
      Exit restart
    End If
    
    If (istop_ == 1) Then 
      xni  = psb_geamax(x,desc_a,info)
      rerr =  rni/(ani*xni+bni)
    Else  If (istop_ == 2) Then 
      rerr = rni/bn2
    Endif
    if (info /= 0) Then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    End If

    If (rerr<=eps) Then 
      Exit restart
    End If
    
    If (itrace_ > 0) then 
      if (((itx==0).or.(mod(itx,itrace_)==0)).and.(me == 0)) &
           & write(*,'(a,i4,3(2x,es10.4))') 'bicgstab: ',itx,rerr
    Endif

    iteration:  Do 
      it   = it + 1
      itx = itx + 1
      If (debug) Write(*,*) 'Iteration: ',itx
      rho_old = rho    
      rho = psb_gedot(q,r,desc_a,info)
!!$      write(0,'(i2," rho old ",2(f,2x))')me,rho,rho_old
      If (debug) Write(0,*) 'Bi-CGSTAB RHO:',rho
      If (rho==dzero) Then
         If (debug) Write(0,*) 'Bi-CGSTAB Itxation breakdown R',rho
        Exit iteration
      Endif

      If (it==1) Then
        Call psb_geaxpby(done,r,dzero,p,desc_a,info)
      Else
        beta = (rho/rho_old)*(alpha/omega)
        Call psb_geaxpby(-omega,v,done,p,desc_a,info)
        Call psb_geaxpby(done,r,beta,p,desc_a,info)
      End If

      Call psb_prc_aply(prec,p,f,desc_a,info,work=aux)

      Call psb_spmm(done,a,f,dzero,v,desc_a,info,&
           & work=aux)

      sigma = psb_gedot(q,v,desc_a,info)
      If (sigma==dzero) Then
         If (debug) Write(0,*) 'Bi-CGSTAB Iteration breakdown S1', sigma
         Exit iteration
      Endif
      If (debug) Write(0,*) 'Bi-CGSTAB SIGMA:',sigma
      alpha = rho/sigma
      Call psb_geaxpby(done,r,dzero,s,desc_a,info)
      if(info.ne.0) then
         call psb_errpush(4010,name,a_err='psb_geaxpby')
         goto 9999
      end if
      Call psb_geaxpby(-alpha,v,done,s,desc_a,info)
      if(info.ne.0) then
         call psb_errpush(4010,name,a_err='psb_geaxpby')
         goto 9999
      end if
      
      Call psb_prc_aply(prec,s,z,desc_a,info,work=aux)
      if(info.ne.0) then
         call psb_errpush(4010,name,a_err='psb_prc_aply')
         goto 9999
      end if

      Call psb_spmm(done,a,z,dzero,t,desc_a,info,&
           & work=aux)

      if(info.ne.0) then
         call psb_errpush(4010,name,a_err='psb_spmm')
         goto 9999
      end if
      
      sigma = psb_gedot(t,t,desc_a,info)
      If (sigma==dzero) Then
         If (debug) Write(0,*) 'BI-CGSTAB ITERATION BREAKDOWN S2', sigma
        Exit iteration
      Endif
      
      tau  = psb_gedot(t,s,desc_a,info)
      omega = tau/sigma
      If (debug) Write(0,*) 'Bi-CGSTAB OMEGA:',omega      
      If (omega==dzero) Then
         If (debug) Write(0,*) 'BI-CGSTAB ITERATION BREAKDOWN O',omega
        Exit iteration
      Endif

      Call psb_geaxpby(alpha,f,done,x,desc_a,info)
      Call psb_geaxpby(omega,z,done,x,desc_a,info)
      Call psb_geaxpby(done,s,dzero,r,desc_a,info)
      Call psb_geaxpby(-omega,t,done,r,desc_a,info)
      
      If (istop_ == 1) Then 
        rni = psb_geamax(r,desc_a,info)
        xni = psb_geamax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)

      Else  If (istop_ == 2) Then 
        rni = psb_genrm2(r,desc_a,info)
        rerr = rni/bn2
      Endif
      If (rerr<=eps) Then 
        Exit restart
      End If
      
      If (itx.Ge.litmax) Exit restart

      If (itrace_ > 0) then 
        if ((mod(itx,itrace_)==0).and.(me == 0)) &
             & write(*,'(a,i4,3(2x,es10.4))') &
             & 'bicgstab: ',itx,rerr
      Endif
      
    End Do iteration
  End Do restart

  If (itrace_ > 0) then 
    if (me == 0) write(*,'(a,i4,3(2x,es10.4))') 'bicgstab: ',itx,rerr
  Endif
  
  If (Present(err)) err=rerr
  If (Present(iter)) iter = itx
  If (rerr>eps) Then
    Write(0,*) 'BI-CGSTAB failed to converge to ',EPS,&
         & ' in ',ITX,' iterations. '
  End If

  Deallocate(aux)
  Call psb_gefree(wwrk,desc_a,info)
  ! restore external global coherence behaviour
  call psb_restore_coher(ictxt,isvch)
!!$  imerr = MPE_Log_event( istpe, 0, "ed CGSTAB" )
  if(info/=0) then
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(ictxt)
     return
  end if
  return

End Subroutine psb_dcgstab

