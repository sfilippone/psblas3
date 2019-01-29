!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!  
!       Contributions to this routine:
!                         Daniela di Serafino    Second University of Naples
!                         Pasqua D'Ambra         ICAR-CNR
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
!   C         [5] G. Sleijpen, D. Fokkema                                  C
!   C             BICGSTAB(L) for linear equations involving unsymmetric   C
!   C             matrices with complex spectrum                           C
!   C             Electronic Trans. on Numer. Analysis, Vol. 1, pp. 11-32, C
!   C             Sep. 1993                                                C
!   C                                                                      C
!   C                                                                      C
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_srgmres.f90
!
! Subroutine: psb_srgmres
!    This subroutine implements the restarted GMRES method with right
!    preconditioning.
!
! Arguments:
!
!    a      -  type(psb_sspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_sprec_type)       Input: preconditioner
!    b      -  real,dimension(:)       Input: vector containing the
!                                         right hand side B
!    x      -  real,dimension(:)       Input/Output: vector containing the
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
!    irst   -  integer(optional)          Input: restart parameter 
!

recursive subroutine psb_srgmres_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,irst,istop)
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
  Real(psb_spk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, irst,istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  Real(psb_spk_), Optional, Intent(out) :: err
! =   local data
  real(psb_spk_), allocatable   :: aux(:)
  real(psb_spk_), allocatable   :: c(:), s(:), h(:,:), rs(:), rst(:)
  type(psb_s_vect_type), allocatable :: v(:)
  type(psb_s_vect_type)              :: w, w1, xt
  real(psb_spk_) :: tmp 
  real(psb_spk_) :: scal, gm, rti, rti1
  integer(psb_ipk_) ::litmax, naux, mglob, it,k, itrace_,&
       & n_row, n_col, nl, int_err(5)
  Logical, Parameter :: exchange=.True., noexchange=.False., use_srot=.true.
  integer(psb_ipk_), Parameter :: irmax = 8
  integer(psb_ipk_) :: itx, i, istop_, err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: ictxt, np, me  
  Real(psb_spk_)     :: rni, xni, bni, ani,bn2, dt, r0n2
  real(psb_dpk_)     :: errnum, errden, deps, derr
  character(len=20)           :: name
  character(len=*), parameter :: methdname='RGMRES'

  info = psb_success_
  name = 'psb_sgmres'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()
  Call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit,*) me,' ',trim(name),': from psb_info',np
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
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()

  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 2
  endif
!
!  ISTOP_ = 1:  Normwise backward error, infinity norm 
!  ISTOP_ = 2:  ||r||/||b||, 2-norm 
!

  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
    info=psb_err_invalid_istop_
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
    info=psb_err_invalid_istop_
    int_err(1)=nl
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
  allocate(aux(naux),h(nl+1,nl+1),&
       &c(nl+1),s(nl+1),rs(nl+1), rst(nl+1),stat=info)

  if (info == psb_success_) call psb_geall(v,desc_a,info,n=nl+1)
  if (info == psb_success_) call psb_geall(w,desc_a,info)
  if (info == psb_success_) call psb_geall(w1,desc_a,info)
  if (info == psb_success_) call psb_geall(xt,desc_a,info)
  if (info == psb_success_) call psb_geasb(v,desc_a,info,mold=x%v)  
  if (info == psb_success_) call psb_geasb(w,desc_a,info,mold=x%v)  
  if (info == psb_success_) call psb_geasb(w1,desc_a,info,mold=x%v)  
  if (info == psb_success_) call psb_geasb(xt,desc_a,info,mold=x%v)  
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_ 
    call psb_errpush(info,name)
    goto 9999
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ' Size of V,W,W1 ',v(1)%get_nrows(),size(v),&
       & w%get_nrows(),w1%get_nrows()


  if (istop_ == 1) then 
    ani = psb_spnrmi(a,desc_a,info)
    bni = psb_geamax(b,desc_a,info)
  else if (istop_ == 2) then 
    bn2 = psb_genrm2(b,desc_a,info)
  else if (istop_ == 3) then
    call psb_geaxpby(sone,b,szero,v(1),desc_a,info)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call psb_spmm(-sone,a,x,sone,v(1),desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    r0n2 = psb_genrm2(v(1),desc_a,info)
  endif
  errnum = szero
  errden = sone
  deps   = eps
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_ 
    call psb_errpush(info,name)
    goto 9999
  end if
  if ((itrace_ > 0).and.(me == 0)) call log_header(methdname)

  itx   = 0
  restart: do 
  
    ! compute r0 = b-ax0
    ! check convergence

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' restart: ',itx,it
    it = 0      
    call psb_geaxpby(sone,b,szero,v(1),desc_a,info)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_spmm(-sone,a,x,sone,v(1),desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if

    rs(1) = psb_genrm2(v(1),desc_a,info)
    rs(2:) = szero
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    scal=sone/rs(1)  ! rs(1) MIGHT BE VERY SMALL - USE DSCAL TO DEAL WITH IT?

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' on entry to amax: b: ',b%get_nrows(),rs(1),scal

    !
    ! check convergence
    !
    if (istop_ == 1) then 
      rni = psb_geamax(v(1),desc_a,info)
      xni = psb_geamax(x,desc_a,info)
      errnum = rni
      errden = (ani*xni+bni)
    else if (istop_ == 2) then 
      rni = psb_genrm2(v(1),desc_a,info)
      errnum = rni
      errden = bn2
    else if (istop_ == 3) then 
      rni = psb_genrm2(v(1),desc_a,info)
      errnum = rni
      errden = r0n2
    endif
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    if (errnum <= eps*errden) exit restart

    if (itrace_ > 0) &
         & call log_conv(methdname,me,itx,itrace_,errnum,errden,deps)
     
    call v(1)%scal(scal) !v(1) = v(1) * scal

    if (itx >= litmax) exit restart  

    !
    ! inner iterations
    !

    inner:  Do i=1,nl
      itx  = itx + 1

      call prec%apply(v(i),w1,desc_a,info)
      call psb_spmm(sone,a,w1,szero,w,desc_a,info,work=aux)
      !

      do k = 1, i
        h(k,i) = psb_gedot(v(k),w,desc_a,info)
        call psb_geaxpby(-h(k,i),v(k),sone,w,desc_a,info)
      end do
      h(i+1,i) = psb_genrm2(w,desc_a,info)
      scal=sone/h(i+1,i)
      call psb_geaxpby(scal,w,szero,v(i+1),desc_a,info)
      do k=2,i
        call srot(1,h(k-1,i),1,h(k,i),1,real(c(k-1),kind=psb_spk_),s(k-1))
      enddo


      rti  = h(i,i)
      rti1 = h(i+1,i) 
      call srotg(rti,rti1,tmp,s(i))
      c(i) = cmplx(tmp,szero)
      call srot(1,h(i,i),1,h(i+1,i),1,real(c(i),kind=psb_spk_),s(i))
      h(i+1,i) = szero
      call srot(1,rs(i),1,rs(i+1),1,real(c(i),kind=psb_spk_),s(i))
      
      if (istop_ == 1) then 
        !
        ! build x and then compute the residual and its infinity norm
        !
        rst = rs
        call w1%set(szero)
        call strsm('l','u','n','n',i,1,sone,h,size(h,1),rst,size(rst,1))
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' Rebuild x-> RS:',rst(1:i)
        do k=1, i
          call psb_geaxpby(rst(k),v(k),sone,xt,desc_a,info)
        end do
        call prec%apply(xt,desc_a,info)
        call psb_geaxpby(sone,x,sone,xt,desc_a,info)
        call psb_geaxpby(sone,b,szero,w1,desc_a,info)
        call psb_spmm(-sone,a,xt,sone,w1,desc_a,info,work=aux)
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
      else if (istop_ == 3) then 
        !
        ! compute the residual 2-norm as byproduct of the solution
        ! procedure of the least-squares problem
        !
        rni = abs(rs(i+1))
        errnum = rni
        errden = r0n2
      endif

      if (errnum <= eps*errden) then 

        if (istop_ == 1) then 
          call psb_geaxpby(sone,xt,szero,x,desc_a,info)
! =          x = xt 
        else if (istop_ == 2) then
          !
          ! build x
          !
          call strsm('l','u','n','n',i,1,sone,h,size(h,1),rs,size(rs,1))
          if (debug_level >= psb_debug_ext_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' Rebuild x-> RS:',rs(1:i)
          call w1%set(szero)
          do k=1, i
            call psb_geaxpby(rs(k),v(k),sone,w1,desc_a,info)
          end do
          call prec%apply(w1,w,desc_a,info)
          call psb_geaxpby(sone,w,sone,x,desc_a,info)
        end if

        exit restart

      end if

      if (itrace_ > 0) &
           & call log_conv(methdname,me,itx,itrace_,errnum,errden,deps)

    end do inner

    if (istop_ == 1) then 
      call psb_geaxpby(sone,xt,szero,x,desc_a,info)!      x = xt 
    else if (istop_ == 2) then
      !
      ! build x
      !
      call strsm('l','u','n','n',nl,1,sone,h,size(h,1),rs,size(rs,1))
      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' Rebuild x-> RS:',rs(1:nl)
      call w1%set(szero)
      do k=1, nl
        call psb_geaxpby(rs(k),v(k),sone,w1,desc_a,info)
      end do
      call prec%apply(w1,w,desc_a,info)
      call psb_geaxpby(sone,w,sone,x,desc_a,info)
    end if
     
  end do restart
  if (itrace_ > 0) &
       & call log_conv(methdname,me,itx,ione,errnum,errden,deps)

  call log_end(methdname,me,itx,errnum,errden,deps,err=derr,iter=iter)
  if (present(err)) err = derr

  
  if (info == psb_success_) call psb_gefree(v,desc_a,info)
  if (info == psb_success_) call psb_gefree(w,desc_a,info)
  if (info == psb_success_) call psb_gefree(w1,desc_a,info)
  if (info == psb_success_) call psb_gefree(xt,desc_a,info)
  if (info == psb_success_) deallocate(aux,h,c,s,rs,rst, stat=info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_srgmres_vect

