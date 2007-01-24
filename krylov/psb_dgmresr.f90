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
! File:  psb_dgmresr.f90
!
! Subroutine: psb_dgmres
!    This subroutine implements the restarted GMRES method.
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
!    irst    -  integer(optional).           The restart value.
!    istop   -  integer(optional).           The stopping criterium.
!
Subroutine psb_dgmresr(a,prec,b,x,eps,desc_a,info,&
     &itmax,iter,err,itrace,irst,istop)
  use psb_base_mod
  use psb_prec_mod
  implicit none

!!$  Parameters 
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
  Real(Kind(1.d0)), allocatable, target   :: aux(:),w(:), v(:,:)
  Real(Kind(1.d0)), allocatable   ::  c(:),s(:), h(:,:), rs(:)
  Real(Kind(1.d0)) :: rerr, scal, gm 
  Integer       ::litmax, liter, naux, m, mglob, it,k, itrace_,&
       & np,me, n_row, n_col, nl, int_err(5)
  Character     :: diagl, diagu
  Logical, Parameter :: exchange=.True., noexchange=.False.  
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ich, ictxt,istop_, err_act
  Logical, Parameter :: debug = .false.
  Real(Kind(1.d0)) :: rni, xni, bni, ani,bn2, dt
  real(kind(1.d0)), external :: dnrm2
  character(len=20)          :: name

  info = 0
  name = 'psb_dgmres'
  call psb_erractionsave(err_act)

  If (debug) Write(0,*) 'entering psb_dgmres'
  ictxt = psb_cd_get_context(desc_a)
  Call psb_info(ictxt, me, np)

  If (debug) Write(0,*) 'psb_dgmres: from gridinfo',np,me

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
    write(0,*) 'psb_dgmres: invalid istop',istop_ 
    info=5001
    int_err(1)=istop_
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
    nl = 10 
    If (debug) Write(0,*) 'not present: irst: ',irst,nl
  Endif


  naux=4*n_col 
  Allocate(aux(naux),h(nl+1,nl+1),&
       &c(nl+1),s(nl+1),rs(nl+1), stat=info)

  if (info == 0) Call psb_geall(v,desc_a,info,n=nl+1)
  if (info == 0) Call psb_geall(w,desc_a,info)
  if (info == 0) Call psb_geasb(v,desc_a,info)  
  if (info == 0) Call psb_geasb(w,desc_a,info)  
  if (info.ne.0) Then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  End If
  if (debug) write(0,*) 'Size of V,W ',size(v),size(v,1),&
       &size(w),size(w,1), size(v(:,1))

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
    it = 0      
    Call psb_geaxpby(done,b,dzero,v(:,1),desc_a,info)
    if (info.ne.0) Then 
       info=4011 
       call psb_errpush(info,name)
       goto 9999
    End If
    Call psb_spmm(-done,a,x,done,v(:,1),desc_a,info,work=aux)
    
    call psb_precaply(prec,v(:,1),desc_a,info)
    rs(1) = psb_genrm2(v(:,1),desc_a,info)
    if (info.ne.0) Then 
       info=4011 
       call psb_errpush(info,name)
       goto 9999
    End If

    scal=done/rs(1)
    If (debug) Write(0,*) 'on entry to amax: b: ',Size(b),rs(1),scal

    if (istop_ == 1) then 
      rni = psb_geamax(v(:,1),desc_a,info)
      xni = psb_geamax(x,desc_a,info)
      rerr =  rni/(ani*xni+bni)
    else if (istop_ == 2) then 
      rni = psb_genrm2(v(:,1),desc_a,info)
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
           & write(*,'(a,i4,3(2x,es10.4))') 'gmres(l): ',itx,rerr
    end If
     
    If (itx.Ge.litmax) Exit restart  

    v(:,1) = v(:,1) * scal

    inner:  Do i=1,nl
      itx  = itx + 1

      Call psb_spmm(done,a,v(:,i),dzero,w,desc_a,info,work=aux)
      call psb_precaply(prec,w,desc_a,info)

      do k = 1, i
        h(k,i) = psb_gedot(v(:,k),w,desc_a,info)
        call psb_geaxpby(-h(k,i),v(:,k),done,w,desc_a,info)
      end do
      h(i+1,i) = psb_genrm2(w,desc_a,info)
      scal=done/h(i+1,i)
      call psb_geaxpby(scal,w,dzero,v(:,i+1),desc_a,info)
      do k=2,i
        dt       = h(k-1,i)
        h(k-1,i) =  c(k-1)*dt + s(k-1)*h(k,i)
        h(k,i)   = -s(k-1)*dt + c(k-1)*h(k,i)
      enddo
      gm =  safe_dn2(h(i,i),h(i+1,i))
      if (debug) write(0,*) 'GM : ',gm
      gm = max(gm,epstol)
      
      c(i) = h(i,i)/gm
      s(i) = h(i+1,i)/gm
      rs(i+1) = -s(i)*rs(i)
      rs(i)   = c(i)*rs(i)
      h(i,i)  = c(i)*h(i,i)+s(i)*h(i+1,i)
            
      if (istop_ == 1) then 
!da modificare, la norma infinito del residuo va calcolata a parte
        rni = abs(rs(i+1))
        xni = psb_geamax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)
      else if (istop_ == 2) then 
        rni = abs(rs(i+1))
        rerr = rni/bn2
      endif

      if (rerr < eps ) then 
        call dtrsm('l','u','n','n',i,1,done,h,size(h,1),rs,nl)
        if (debug) write(0,*) 'Rebuild x-> RS:',rs(21:nl)
        do k=1, i
          call psb_geaxpby(rs(k),v(:,k),done,x,desc_a,info)
        end do
        exit restart
      end if
      If (itrace_ > 0) then 
        if ((mod(itx,itrace_)==0).and.(me == 0))&
             & write(*,'(a,i4,3(2x,es10.4))') 'gmres(l): ',itx,rerr
      end If

    end Do inner
    if (debug) write(0,*) 'Before DTRSM :',rs(1:nl)
    call dtrsm('l','u','n','n',nl,1,done,h,size(h,1),rs,nl)
    if (debug) write(0,*) 'Rebuild x-> RS:',rs(21:nl)
    do k=1, nl
      call psb_geaxpby(rs(k),v(:,k),done,x,desc_a,info)
    end do
     
  End Do restart
  If (itrace_ > 0) then 
    if (me == 0) write(*,'(a,i4,3(2x,es10.4))') 'gmres(l): ',itx,rerr
  end If

  If (Present(err)) err=rerr
  If (Present(iter)) iter = itx
  If ((rerr>eps).and. (me == 0))  Then
    Write(0,*) 'gmresr(l) failed to converge to ',eps,&
         & ' in ',itx,' iterations  '
  End If


  Deallocate(aux,h,c,s,rs, stat=info)
  Call psb_gefree(v,desc_a,info)
  Call psb_gefree(w,desc_a,info)

  ! restore external global coherence behaviour
  call psb_restore_coher(ictxt,isvch)

  if (info /= 0) then
     info=4011
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


contains
  function safe_dn2(a,b)
    real(kind(1.d0)), intent(in) :: a, b
    real(kind(1.d0))  :: safe_dn2
    real(kind(1.d0))  :: t
    
    t = max(abs(a),abs(b))
    if (t==0.d0) then 
      safe_dn2 = 0.d0
    else
      safe_dn2 = t * sqrt(abs(a/t)**2 + abs(b/t)**2)
    endif
    return
  end function safe_dn2
    

End Subroutine psb_dgmresr


