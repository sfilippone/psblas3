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
! File:  psb_zcgstab.f90
!
! Subroutine: psb_zcgstab
!    This subroutine implements the CG Stabilized method.
!
! Parameters:
!    a       -  type(<psb_zspmat_type>).     The sparse matrix containing A.
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
Subroutine psb_zcgstab(a,prec,b,x,eps,desc_a,info,&
     &itmax,iter,err,itrace, istop)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_tools_mod
  use psb_const_mod
  use psb_prec_mod
  use psb_error_mod
  Implicit None
!!$  parameters 
  Type(psb_zspmat_type), Intent(in)  :: a
  Type(psb_zprec_type), Intent(in)   :: prec 
  Type(psb_desc_type), Intent(in)    :: desc_a
  Complex(Kind(1.d0)), Intent(in)       :: b(:)
  Complex(Kind(1.d0)), Intent(inout)    :: x(:)
  Real(Kind(1.d0)), Intent(in)       :: eps
  integer, intent(out)               :: info
  Integer, Optional, Intent(in)      :: itmax, itrace, istop
  Integer, Optional, Intent(out)     :: iter
  Real(Kind(1.d0)), Optional, Intent(out) :: err
!!$   Local data
  Complex(Kind(1.d0)), Pointer  :: aux(:),wwrk(:,:)
  Complex(Kind(1.d0)), Pointer  :: q(:),&
       & r(:), p(:), v(:), s(:), t(:), z(:), f(:)
  Integer, Pointer           :: iperm(:), ipnull(:), ipsave(:)
  Real(Kind(1.d0)) :: rerr
  Integer          :: litmax, liter, naux, m, mglob, it,itrac,&
       & nprows,npcols,myrow,mycol, n_row, n_col
  Character     ::diagl, diagu
  Logical, Parameter :: debug = .false.
  Logical, Parameter :: exchange=.True., noexchange=.False., debug1 = .False.
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ich, icontxt, err_act, int_err(5),ii
  Integer            :: listop
  Logical            :: do_renum_left
  complex(Kind(1.d0)) :: alpha, beta, rho, rho_old, sigma, omega, tau
  Real(Kind(1.d0)) :: rni, xni, bni, ani, rn0, bn2
!!$  Integer   istpb, istpe, ifctb, ifcte, imerr, irank, icomm,immb,imme
!!$  Integer mpe_log_get_event_number,mpe_Describe_state,mpe_log_event
  character(len=20)             :: name,ch_err

  info = 0
  name = 'psb_zcgstab'
  call psb_erractionsave(err_act)

  If (debug) Write(*,*) 'Entering PSB_ZCGSTAB',present(istop)
  icontxt = desc_a%matrix_data(psb_ctxt_)
  CALL blacs_gridinfo(icontxt,nprows,npcols,myrow,mycol)
  if (debug) write(*,*) 'PSB_ZCGSTAB: From GRIDINFO',nprows,npcols,myrow

  mglob = desc_a%matrix_data(psb_m_)
  n_row = desc_a%matrix_data(psb_n_row_)
  n_col = desc_a%matrix_data(psb_n_col_)

  If (Present(istop)) Then 
    listop = istop 
  Else
    listop = 1
  Endif
!
!  LISTOP = 1:  Normwise backward error, infinity norm 
!  LISTOP = 2:  ||r||/||b||   norm 2 
!

  If ((prec%prec < min_prec_).Or.(prec%prec > max_prec_) ) Then
     Write(0,*) 'PSB_CGSTAB: Invalid IPREC',prec%prec
     info=5002
     int_err(1)=prec%prec
     err=info
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  Endif
  
  if ((listop < 1 ).or.(listop > 2 ) ) then
    write(0,*) 'psb_bicgstab: invalid istop',listop 
    info=5001
    int_err(1)=listop
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  naux=6*n_col 
  allocate(aux(naux),stat=info)
  call psb_geall(mglob,8,wwrk,desc_a,info)
  call psb_geasb(wwrk,desc_a,info)  
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
     itrac = itrace
  Else
     itrac = -1
  End If
  
  diagl = 'U'
  diagu = 'U'

  ! Ensure global coherence for convergence checks.
  Call blacs_get(icontxt,16,isvch)
  ich = 1 
  Call blacs_set(icontxt,16,ich)

  itx   = 0

  If (listop == 1) Then 
    ani = psb_spnrmi(a,desc_a,info)
    bni = psb_geamax(b,desc_a,info)
  Else If (listop == 2) Then 
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
    Call psb_geaxpby(zone,b,zzero,r,desc_a,info)
    Call psb_spmm(-zone,a,x,zone,r,desc_a,info,work=aux)
    Call psb_geaxpby(zone,r,zzero,q,desc_a,info)
    if (info /= 0) Then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    End If
    
    rho = zzero
    If (debug) Write(*,*) 'On entry to AMAX: B: ',Size(b)
    
!
!   Must always provide norm of R into RNI below for first check on 
!   residual
!
    If (listop == 1) Then 
      rni = psb_geamax(r,desc_a,info)
      xni = psb_geamax(x,desc_a,info)
    Else If (listop == 2) Then 
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
      If (itrac /= -1) Then 
        If (myrow == 0) Write(itrac,*) 'BiCGSTAB: ',itx,rn0
      Endif
      Exit restart
    End If
    
    If (listop == 1) Then 
      xni  = psb_geamax(x,desc_a,info)
      rerr =  rni/(ani*xni+bni)
      If (itrac /= -1) Then 
        If (myrow == 0) Write(itrac,'(a,i4,5(2x,es10.4))') 'bicgstab: ',itx,rerr,rni,bni,&
             &xni,ani
      Endif
    Else  If (listop == 2) Then 
      rerr = rni/bn2
      If (itrac /= -1) Then 
        If (myrow == 0) Write(itrac,'(a,i4,3(2x,es10.4))') 'bicgstab: ',itx,rerr,rni,bn2
      Endif
    Endif
    if (info /= 0) Then 
       info=4011
       call psb_errpush(info,name)
       goto 9999
    End If

    
    If (rerr<=eps) Then 
      Exit restart
    End If

    iteration:  Do 
      it   = it + 1
      itx = itx + 1
      If (debug) Write(*,*) 'Iteration: ',itx

      rho_old = rho    
      rho = psb_gedot(q,r,desc_a,info)

      If (rho==zzero) Then
         If (debug) Write(0,*) 'Bi-CGSTAB Itxation breakdown R',rho
        Exit iteration
      Endif

      If (it==1) Then
        Call psb_geaxpby(zone,r,zzero,p,desc_a,info)
      Else
        beta = (rho/rho_old)*(alpha/omega)
        Call psb_geaxpby(-omega,v,zone,p,desc_a,info)
        Call psb_geaxpby(zone,r,beta,p,desc_a,info)
      End If

      Call psb_prc_aply(prec,p,f,desc_a,info,work=aux)

      Call psb_spmm(zone,a,f,zzero,v,desc_a,info,&
           & work=aux)

      sigma = psb_gedot(q,v,desc_a,info)
      If (sigma==zzero) Then
         If (debug) Write(0,*) 'Bi-CGSTAB Iteration breakdown S1', sigma
         Exit iteration
      Endif
      
      alpha = rho/sigma
      Call psb_geaxpby(zone,r,zzero,s,desc_a,info)
      if(info.ne.0) then
         call psb_errpush(4010,name,a_err='psb_geaxpby')
         goto 9999
      end if
      Call psb_geaxpby(-alpha,v,zone,s,desc_a,info)
      if(info.ne.0) then
         call psb_errpush(4010,name,a_err='psb_geaxpby')
         goto 9999
      end if
      
      Call psb_prc_aply(prec,s,z,desc_a,info,work=aux)
      if(info.ne.0) then
         call psb_errpush(4010,name,a_err='psb_prc_aply')
         goto 9999
      end if

      Call psb_spmm(zone,a,z,zzero,t,desc_a,info,&
           & work=aux)

      if(info.ne.0) then
         call psb_errpush(4010,name,a_err='psb_spmm')
         goto 9999
      end if
      
      sigma = psb_gedot(t,t,desc_a,info)
      If (sigma==zzero) Then
         If (debug) Write(0,*) 'BI-CGSTAB ITERATION BREAKDOWN S2', sigma
        Exit iteration
      Endif
      
      tau  = psb_gedot(t,s,desc_a,info)
      omega = tau/sigma
      
      If (omega==zzero) Then
         If (debug) Write(0,*) 'BI-CGSTAB ITERATION BREAKDOWN O',omega
        Exit iteration
      Endif

      Call psb_geaxpby(alpha,f,zone,x,desc_a,info)
      Call psb_geaxpby(omega,z,zone,x,desc_a,info)
      Call psb_geaxpby(zone,s,zzero,r,desc_a,info)
      Call psb_geaxpby(-omega,t,zone,r,desc_a,info)
      
      If (listop == 1) Then 
        rni = psb_geamax(r,desc_a,info)
        xni = psb_geamax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)
        If (itrac /= -1) Then 
          If (myrow == 0) Write(itrac,'(a,i4,5(2x,es10.4))') &
               & 'bicgstab: ',itx,rerr,rni,bni,xni,ani
        Endif

      Else  If (listop == 2) Then 
        rni = psb_genrm2(r,desc_a,info)
        rerr = rni/bn2
        If (itrac /= -1) Then 
          If (myrow == 0) Write(itrac,'(a,i4,3(2x,es10.4)))') &
               & 'bicgstab: ',itx,rerr,rni,bn2
        Endif
      Endif
      
      If (rerr<=eps) Then 
        Exit restart
      End If
      
      If (itx.Ge.litmax) Exit restart
    End Do iteration
  End Do restart

  If (Present(err)) err=rerr
  If (Present(iter)) iter = itx
  If (rerr>eps) Then
    Write(0,*) 'BI-CGSTAB FAILED TO CONVERGE TO ',EPS,&
         & ' IN ',ITX,' ITERATIONS  '
  End If

  Deallocate(aux)
  Call psb_gefree(wwrk,desc_a,info)
  ! restore external global coherence behaviour
  Call blacs_set(icontxt,16,isvch)
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
     call psb_error(icontxt)
     return
  end if
  return

End Subroutine psb_zcgstab

