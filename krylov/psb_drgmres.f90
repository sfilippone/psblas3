!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$
!!$     Contributions to this routine:
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR
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
! File:  psb_drgmres.f90
!
! Subroutine: psb_drgmres
!    This subroutine implements the restarted GMRES method with right
!    preconditioning.
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
!    irst   -  integer(optional)          Input: restart parameter
!                                         
! 
Subroutine psb_drgmres(a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop)
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
  Real(Kind(1.d0)), allocatable, target   :: aux(:),w(:),w1(:), v(:,:)
  Real(Kind(1.d0)), allocatable   ::  c(:),s(:), h(:,:), rs(:),rst(:),xt(:)
  Real(Kind(1.d0)) :: scal, gm, rti, rti1
  Integer       ::litmax, naux, mglob, it,k, itrace_,&
       & np,me, n_row, n_col, nl, int_err(5)
  Logical, Parameter :: exchange=.True., noexchange=.False., use_drot=.true.
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ictxt,istop_, err_act
  integer            :: debug_level, debug_unit
  Real(Kind(1.d0)) :: rni, xni, bni, ani,bn2, dt
  real(kind(1.d0)), external :: dnrm2
  real(kind(1.d0))   :: errnum, errden
  character(len=20)          :: name

  info = 0
  name = 'psb_dgmres'
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
!  ISTOP_ = 2:  ||r||/||b||, 2-norm 
!

  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
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
    If (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' present: irst: ',irst,nl
  Else
    nl = 10 
    If (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' not present: irst: ',irst,nl
  Endif
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
  Allocate(aux(naux),h(nl+1,nl+1),&
       &c(nl+1),s(nl+1),rs(nl+1), rst(nl+1),stat=info)

  if (info == 0) Call psb_geall(v,desc_a,info,n=nl+1)
  if (info == 0) Call psb_geall(w,desc_a,info)
  if (info == 0) Call psb_geall(w1,desc_a,info)
  if (info == 0) Call psb_geall(xt,desc_a,info)
  if (info == 0) Call psb_geasb(v,desc_a,info)  
  if (info == 0) Call psb_geasb(w,desc_a,info)  
  if (info == 0) Call psb_geasb(w1,desc_a,info)
  if (info == 0) Call psb_geasb(xt,desc_a,info)
  if (info.ne.0) Then 
    info=4011 
    call psb_errpush(info,name)
    goto 9999
  End If
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Size of V,W,W1 ',size(v),size(v,1),&
       & size(w),size(w,1),size(w1),size(w1,1), size(v(:,1))

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
  if (info.ne.0) Then 
    info=4011 
    call psb_errpush(info,name)
    goto 9999
  End If

  itx   = 0
  restart: Do 
  
    ! compute r0 = b-ax0
    ! check convergence
    ! compute v1 = r0/||r0||_2

    If (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' restart: ',itx,it
    it = 0      
    Call psb_geaxpby(done,b,dzero,v(:,1),desc_a,info)
    if (info.ne.0) Then 
      info=4011 
      call psb_errpush(info,name)
      goto 9999
    End If

    Call psb_spmm(-done,a,x,done,v(:,1),desc_a,info,work=aux)
    if (info.ne.0) Then 
      info=4011 
      call psb_errpush(info,name)
      goto 9999
    End If

    rs(1) = psb_genrm2(v(:,1),desc_a,info)
    rs(2:) = dzero
    if (info.ne.0) Then 
      info=4011 
      call psb_errpush(info,name)
      goto 9999
    End If
    scal=done/rs(1)  ! rs(1) MIGHT BE VERY SMALL - USE DSCAL TO DEAL WITH IT?

    If (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' on entry to amax: b: ',Size(b),rs(1),scal

    !
    ! check convergence
    !
    if (istop_ == 1) then 
      rni = psb_geamax(v(:,1),desc_a,info)
      xni = psb_geamax(x,desc_a,info)
      errnum = rni
      errden = (ani*xni+bni)
    else if (istop_ == 2) then 
      rni = psb_genrm2(v(:,1),desc_a,info)
      errnum = rni
      errden = bn2
    endif
    if (info.ne.0) Then 
      info=4011 
      call psb_errpush(info,name)
      goto 9999
    End If
    
    If (errnum <= eps*errden) Then 
      Exit restart
    End If

    If (itrace_ > 0) then 
      if ((mod(itx,itrace_)==0).and.(me == 0))&
           & write(*,'(a,i4,3(2x,es10.4))') 'gmres(l): ',itx,errnum,eps*errden
    end If
     
    v(:,1) = v(:,1) * scal

    If (itx.Ge.litmax) Exit restart  

    !
    ! inner iterations
    !

    inner:  Do i=1,nl
      itx  = itx + 1

      call psb_precaply(prec,v(:,i),w1,desc_a,info)
      Call psb_spmm(done,a,w1,dzero,w,desc_a,info,work=aux)
      !

      do k = 1, i
        h(k,i) = psb_gedot(v(:,k),w,desc_a,info)
        call psb_geaxpby(-h(k,i),v(:,k),done,w,desc_a,info)
      end do
      h(i+1,i) = psb_genrm2(w,desc_a,info)
      scal=done/h(i+1,i)
      call psb_geaxpby(scal,w,dzero,v(:,i+1),desc_a,info)
      if (use_drot) then 
        do k=2,i
          call drot(1,h(k-1,i),1,h(k,i),1,c(k-1),s(k-1))
        enddo
        
        rti  = h(i,i)
        rti1 = h(i+1,i) 
        call drotg(rti,rti1,c(i),s(i))
        call drot(1,h(i,i),1,h(i+1,i),1,c(i),s(i))
        h(i+1,i) = dzero
        call drot(1,rs(i),1,rs(i+1),1,c(i),s(i))

      else
        do k=2,i
          dt       = h(k-1,i)
          h(k-1,i) =  c(k-1)*dt + s(k-1)*h(k,i)
          h(k,i)   = -s(k-1)*dt + c(k-1)*h(k,i)
        enddo
        gm =  safe_dn2(h(i,i),h(i+1,i))
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),' GM : ',gm
        gm = max(gm,epstol)

        c(i) = h(i,i)/gm
        s(i) = h(i+1,i)/gm
        rs(i+1) = -s(i)*rs(i)
        rs(i)   = c(i)*rs(i)
        h(i,i)  = c(i)*h(i,i)+s(i)*h(i+1,i)
      endif
      if (istop_ == 1) then 
        !
        ! build x and then compute the residual and its infinity norm
        !
        rst = rs
        xt = dzero
        call dtrsm('l','u','n','n',i,1,done,h,size(h,1),rst,size(rst,1))
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' Rebuild x-> RS:',rst(1:nl)
        do k=1, i
          call psb_geaxpby(rst(k),v(:,k),done,xt,desc_a,info)
        end do
        call psb_precaply(prec,xt,desc_a,info)
        call psb_geaxpby(done,x,done,xt,desc_a,info)
        call psb_geaxpby(done,b,dzero,w1,desc_a,info)
        call psb_spmm(-done,a,xt,done,w1,desc_a,info,work=aux)
        rni = psb_geamax(w1,desc_a,info)
        xni = psb_geamax(xt,desc_a,info)
        errnum = rni
        errden = (ani*xni+bni)
        !

      else if (istop_ == 2) then 
        !
        ! compute the residual 2-norm as byproduct of the solution
        ! procedure of the least-squares problem
        !
        rni = abs(rs(i+1))
        errnum = rni
        errden = bn2
      endif

      If (errnum <= eps*errden) Then 

        if (istop_ == 1) then 
          x = xt 
        else if (istop_ == 2) then
          !
          ! build x
          !
          call dtrsm('l','u','n','n',i,1,done,h,size(h,1),rs,size(rs,1))
          if (debug_level >= psb_debug_ext_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' Rebuild x-> RS:',rs(1:nl)
          w1 = dzero 
          do k=1, i
            call psb_geaxpby(rs(k),v(:,k),done,w1,desc_a,info)
          end do
          call psb_precaply(prec,w1,w,desc_a,info)
          call psb_geaxpby(done,w,done,x,desc_a,info)
        end if

        exit restart

      end if

      If (itrace_ > 0) then 
        if ((mod(itx,itrace_)==0).and.(me == 0))&
             & write(*,'(a,i4,3(2x,es10.4))') 'gmres(l): ',itx,errnum,eps*errden
      end If

    end Do inner

    if (istop_ == 1) then 
      x = xt 
    else if (istop_ == 2) then
      !
      ! build x
      !
      call dtrsm('l','u','n','n',nl,1,done,h,size(h,1),rs,size(rs,1))
      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' Rebuild x-> RS:',rs(1:nl)
      w1 = dzero 
      do k=1, nl
        call psb_geaxpby(rs(k),v(:,k),done,w1,desc_a,info)
      end do
      call psb_precaply(prec,w1,w,desc_a,info)
      call psb_geaxpby(done,w,done,x,desc_a,info)
    end if
     
  End Do restart
  If (itrace_ > 0) then 
    if (me == 0) write(*,'(a,i4,3(2x,es10.4))') 'gmres(l): ',itx,errnum,eps*errden
  end If

  If (Present(err)) then 
    if (errden /= dzero) then 
      err = errnum/errden
    else
      err = errnum
    end if
  end If

  If (Present(iter)) iter = itx
  If ((errnum > eps*errden).and.(me==0)) Then
    write(debug_unit,*) 'gmresr(l) failed to converge to ',eps*errden,&
         & ' in ',itx,' iterations  '
  End If


  Deallocate(aux,h,c,s,rs,rst, stat=info)
  Call psb_gefree(v,desc_a,info)
  Call psb_gefree(w,desc_a,info)
  Call psb_gefree(w1,desc_a,info)
  Call psb_gefree(xt,desc_a,info)

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
    

End Subroutine psb_drgmres


