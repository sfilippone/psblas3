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
! File:  psb_zrgmres.f90
!
! Subroutine: psb_zrgmres
!    This subroutine implements the restarted GMRES method with right
!    preconditioning.
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
!    irst   -  integer(optional)          Input: restart parameter 
!
Subroutine psb_zrgmres(a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod, psb_protect_name => psb_zrgmres
  implicit none

!!$  Parameters 
  Type(psb_zspmat_type), Intent(in)  :: a
  Type(psb_zprec_type), Intent(in)   :: prec 
  Type(psb_desc_type), Intent(in)    :: desc_a
  complex(Kind(1.d0)), Intent(in)       :: b(:)
  complex(Kind(1.d0)), Intent(inout)    :: x(:)
  Real(Kind(1.d0)), Intent(in)       :: eps
  integer, intent(out)               :: info
  Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
  Integer, Optional, Intent(out)     :: iter
  Real(Kind(1.d0)), Optional, Intent(out) :: err
!!$   local data
  complex(Kind(1.d0)), allocatable, target   :: aux(:),w(:),w1(:), v(:,:)
  complex(Kind(1.d0)), allocatable   ::  c(:),s(:), h(:,:), rs(:),rst(:),xt(:)
  Real(Kind(1.d0)) :: tmp
  complex(kind(1.d0)) :: rti, rti1, scal
  Integer       ::litmax, naux, mglob, it,k, itrace_,&
       & np,me, n_row, n_col, nl, int_err(5)
  Logical, Parameter :: exchange=.True., noexchange=.False.
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ictxt,istop_, err_act
  integer            :: debug_level, debug_unit
  Real(Kind(1.d0)) :: rni, xni, bni, ani,bn2
  real(kind(1.d0)), external :: dznrm2
  real(kind(1.d0))   :: errnum, errden
  character(len=20)           :: name
  character(len=*), parameter :: methdname='RGMRES'

  info = 0
  name = 'psb_zgmres'
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
         & ' present: irst: ',irst,nl
  else
    nl = 10 
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
  allocate(aux(naux),h(nl+1,nl+1),&
       &c(nl+1),s(nl+1),rs(nl+1), rst(nl+1),stat=info)

  if (info == 0) Call psb_geall(v,desc_a,info,n=nl+1)
  if (info == 0) Call psb_geall(w,desc_a,info)
  if (info == 0) Call psb_geall(w1,desc_a,info)
  if (info == 0) Call psb_geall(xt,desc_a,info)
  if (info == 0) Call psb_geasb(v,desc_a,info)  
  if (info == 0) Call psb_geasb(w,desc_a,info)  
  if (info == 0) Call psb_geasb(w1,desc_a,info)
  if (info == 0) Call psb_geasb(xt,desc_a,info)
  if (info /= 0) then 
    info=4011 
    call psb_errpush(info,name)
    goto 9999
  end if
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
  if (info /= 0) then 
    info=4011 
    call psb_errpush(info,name)
    goto 9999
  end if

  itx   = 0
  restart: Do 

    ! compute r0 = b-ax0
    ! check convergence
    ! compute v1 = r0/||r0||_2

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' restart: ',itx,it
    it = 0      
    call psb_geaxpby(zone,b,zzero,v(:,1),desc_a,info)
    if (info /= 0) then 
      info=4011 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_spmm(-zone,a,x,zone,v(:,1),desc_a,info,work=aux)
    if (info /= 0) Then 
      info=4011 
      call psb_errpush(info,name)
      goto 9999
    End If

    rs(1) = psb_genrm2(v(:,1),desc_a,info)
    rs(2:) = zzero
    if (info /= 0) Then 
      info=4011 
      call psb_errpush(info,name)
      goto 9999
    end if
    scal=done/rs(1)  ! rs(1) MIGHT BE VERY SMALL - USE DSCAL TO DEAL WITH IT?

    if (debug_level >= psb_debug_ext_) &
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
    if (info /= 0) then 
      info=4011 
      call psb_errpush(info,name)
      goto 9999
    end if

    if (errnum <= eps*errden) exit restart

    if (itrace_ > 0) &
         & call log_conv(methdname,me,itx,itrace_,errnum,errden,eps)

    v(:,1) = v(:,1) * scal

    if (itx >= litmax) exit restart  

    !
    ! inner iterations
    !

    inner:  Do i=1,nl
      itx  = itx + 1

      call psb_precaply(prec,v(:,i),w1,desc_a,info)
      Call psb_spmm(zone,a,w1,zzero,w,desc_a,info,work=aux)
      !

      do k = 1, i
        h(k,i) = psb_gedot(v(:,k),w,desc_a,info)
        call psb_geaxpby(-h(k,i),v(:,k),zone,w,desc_a,info)
      end do
      h(i+1,i) = psb_genrm2(w,desc_a,info)
      scal=done/h(i+1,i)
      call psb_geaxpby(scal,w,zzero,v(:,i+1),desc_a,info)
      do k=2,i
        call zrot(1,h(k-1,i),1,h(k,i),1,real(c(k-1)),s(k-1))
      enddo

      rti  = h(i,i)
      rti1 = h(i+1,i) 
      call zrotg(rti,rti1,tmp,s(i))
      c(i) = cmplx(tmp,dzero)
      call zrot(1,h(i,i),1,h(i+1,i),1,real(c(i)),s(i))
      h(i+1,i) = zzero
      call zrot(1,rs(i),1,rs(i+1),1,real(c(i)),s(i))

      if (istop_ == 1) then 
        !
        ! build x and then compute the residual and its infinity norm
        !
        rst = rs
        xt = zzero
        call ztrsm('l','u','n','n',i,1,zone,h,size(h,1),rst,size(rst,1))
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' Rebuild x-> RS:',rst(1:nl)
        do k=1, i
          call psb_geaxpby(rst(k),v(:,k),zone,xt,desc_a,info)
        end do
        call psb_precaply(prec,xt,desc_a,info)
        call psb_geaxpby(zone,x,zone,xt,desc_a,info)
        call psb_geaxpby(zone,b,zzero,w1,desc_a,info)
        call psb_spmm(-zone,a,xt,zone,w1,desc_a,info,work=aux)
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
          call ztrsm('l','u','n','n',i,1,zone,h,size(h,1),rs,size(rs,1))
          if (debug_level >= psb_debug_ext_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' Rebuild x-> RS:',rs(1:nl)
          w1 = zzero 
          do k=1, i
            call psb_geaxpby(rs(k),v(:,k),zone,w1,desc_a,info)
          end do
          call psb_precaply(prec,w1,w,desc_a,info)
          call psb_geaxpby(zone,w,zone,x,desc_a,info)
        end if

        exit restart

      end if

      if (itrace_ > 0) &
           & call log_conv(methdname,me,itx,itrace_,errnum,errden,eps)

    end do inner

    if (istop_ == 1) then 
      x = xt 
    else if (istop_ == 2) then
      !
      ! build x
      !
      call ztrsm('l','u','n','n',nl,1,zone,h,size(h,1),rs,size(rs,1))
      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' Rebuild x-> RS:',rs(1:nl)
      w1 = zzero 
      do k=1, nl
        call psb_geaxpby(rs(k),v(:,k),zone,w1,desc_a,info)
      end do
      call psb_precaply(prec,w1,w,desc_a,info)
      call psb_geaxpby(zone,w,zone,x,desc_a,info)
    end if

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

  deallocate(aux,h,c,s,rs,rst, stat=info)
  if (info == 0) call psb_gefree(v,desc_a,info)
  if (info == 0) call psb_gefree(w,desc_a,info)
  if (info == 0) call psb_gefree(w1,desc_a,info)
  if (info == 0) call psb_gefree(xt,desc_a,info)
  if (info /= 0) then
    info=4011
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


contains

  subroutine zrot( n, cx, incx, cy, incy, c, s )
    !
    !  -- lapack auxiliary routine (version 3.0) --
    !     univ. of tennessee, univ. of california berkeley, nag ltd.,
    !     courant institute, argonne national lab, and rice university
    !     october 31, 1992
    !
    !     .. scalar arguments ..
    integer            incx, incy, n
    real(kind(1.d0))    c
    complex(kind(1.d0))   s
    !     ..
    !     .. array arguments ..
    complex(kind(1.d0)) cx( * ), cy( * )
    !     ..
    !
    !  purpose
    !  =======
    !
    !  zrot   applies a plane rotation, where the cos (c) is real and the
    !  sin (s) is complex, and the vectors cx and cy are complex.
    !
    !  arguments
    !  =========
    !
    !  n       (input) integer
    !          the number of elements in the vectors cx and cy.
    !
    !  cx      (input/output) complex*16 array, dimension (n)
    !          on input, the vector x.
    !          on output, cx is overwritten with c*x + s*y.
    !
    !  incx    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  cy      (input/output) complex*16 array, dimension (n)
    !          on input, the vector y.
    !          on output, cy is overwritten with -conjg(s)*x + c*y.
    !
    !  incy    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  c       (input) double precision
    !  s       (input) complex*16
    !          c and s define a rotation
    !             [  c          s  ]
    !             [ -conjg(s)   c  ]
    !          where c*c + s*conjg(s) = 1.0.
    !
    ! =====================================================================
    !
    !     .. local scalars ..
    integer            i, ix, iy
    complex(kind(1.d0))         stemp
    !     ..
    !     .. intrinsic functions ..
    intrinsic          dconjg
    !     ..
    !     .. executable statements ..
    !
    if( n <= 0 ) return
    if( incx == 1 .and. incy == 1 ) then 
      !
      !     code for both increments equal to 1
      !
      do  i = 1, n
        stemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - dconjg(s)*cx(i)
        cx(i) = stemp
      end do
    else
      !
      !     code for unequal increments or equal increments not equal to 1
      !
      ix = 1
      iy = 1
      if( incx < 0 )ix = ( -n+1 )*incx + 1
      if( incy < 0 )iy = ( -n+1 )*incy + 1
      do  i = 1, n
        stemp  = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - dconjg(s)*cx(ix)
        cx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    return
    return
  end subroutine zrot
  !
  !
  subroutine zrotg(ca,cb,c,s)
    complex(kind(1.d0)) ca,cb,s
    real(kind(1.d0)) c
    real(kind(1.d0)) norm,scale
    complex(kind(1.d0)) alpha
    !
    if (cdabs(ca) == 0.0d0) then 
      !
      c = 0.0d0
      s = (1.0d0,0.0d0)
      ca = cb
      return
    end if
    !

    scale = cdabs(ca) + cdabs(cb)
    norm = scale*dsqrt((cdabs(ca/dcmplx(scale,0.0d0)))**2 +&
         &   (cdabs(cb/dcmplx(scale,0.0d0)))**2)
    alpha = ca /cdabs(ca)
    c = cdabs(ca) / norm
    s = alpha * dconjg(cb) / norm
    ca = alpha * norm
    !

    return
  end subroutine zrotg
End Subroutine psb_zrgmres
