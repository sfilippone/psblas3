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
!!$ CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! File:  psb_dbicg.f90
!
! Subroutine: psb_dbicg
!    This subroutine implements the BiCG method.
!
! Arguments:
!
!    a      -  type(psb_dspmat_type)      Input: sparse matrix containing A.
!    prec   -  type(psb_dprec_type)       Input: preconditioner
!    b      -  real,dimension(:)            Input: vector containing the
!                                           right hand side B
!    x      -  real,dimension(:)            Input/Output: vector containing the
!                                           initial guess and final solution X.
!    eps    -  real                         Input: Stopping tolerance; the iteration is
!                                           stopped when the error estimate |err| <= eps
!    desc_a -  type(psb_desc_type).       Input: The communication descriptor.
!    info   -  integer.                     Output: Return code
!
!    itmax  -  integer(optional)            Input: maximum number of iterations to be
!                                           performed.
!    iter   -  integer(optional)            Output: how many iterations have been
!                                           performed.
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
!
subroutine psb_dbicg(a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,istop)
  use psb_base_mod
  use psb_prec_mod
  implicit none

!!$  parameters 
  type(psb_dspmat_type), intent(in)  :: a
  type(psb_dprec_type), intent(in)   :: prec 
  type(psb_desc_type), intent(in)    :: desc_a
  real(kind(1.d0)), intent(in)       :: b(:)
  real(kind(1.d0)), intent(inout)    :: x(:)
  real(kind(1.d0)), intent(in)       :: eps
  integer, intent(out)               :: info
  integer, optional, intent(in)      :: itmax, itrace, istop
  integer, optional, intent(out)     :: iter
  real(kind(1.d0)), optional, intent(out) :: err
!!$   local data
  real(kind(1.d0)), allocatable, target  :: aux(:),wwrk(:,:)
  real(kind(1.d0)), pointer  :: ww(:), q(:),&
       & r(:), p(:), zt(:), pt(:), z(:), rt(:),qt(:)
  integer           :: int_err(5)
  integer       ::litmax, naux, mglob, it, itrace_,&
       & np,me, n_row, n_col, istop_, err_act
  integer            :: debug_level, debug_unit
  logical, parameter :: exchange=.true., noexchange=.false.  
  integer, parameter :: irmax = 8
  integer            :: itx, isvch, ictxt
  real(kind(1.d0)) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,& 
       & sigma,bn2
  real(kind(1.d0))   :: errnum, errden
  character(len=20)             :: name,ch_err

  info = 0
  name = 'psb_dbicg'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit,*) me,' ',trim(name),': from psb_info',np

  mglob = psb_cd_get_global_rows(desc_a)
  n_row = psb_cd_get_local_rows(desc_a)
  n_col = psb_cd_get_local_cols(desc_a)

  ! Ensure global coherence for convergence checks.
  call psb_set_coher(ictxt,isvch)


  if (present(istop)) then 
    istop_ = istop 
  else
    istop_ = 1
  endif
  !
  !  istop_ = 1:  normwise backward error, infinity norm 
  !  istop_ = 2:  ||r||/||b||   norm 2 
  !

  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
    info=5001
    int_err=istop_
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

  allocate(aux(naux),stat=info)
  if (info == 0) call psb_geall(wwrk,desc_a,info,n=9)
  if (info == 0) call psb_geasb(wwrk,desc_a,info)  
  if(info.ne.0) then
    info=4011
    ch_err='psb_asb'
    err=info
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  q  => wwrk(:,1)
  qt => wwrk(:,2)
  r  => wwrk(:,3)
  rt => wwrk(:,4)
  p  => wwrk(:,5)
  pt => wwrk(:,6)
  z  => wwrk(:,7)
  zt => wwrk(:,8)
  ww => wwrk(:,9)

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

  itx   = 0

  if (istop_ == 1) then 
    ani = psb_spnrmi(a,desc_a,info)
    bni = psb_geamax(b,desc_a,info)
  else if (istop_ == 2) then 
    bn2 = psb_genrm2(b,desc_a,info)
  endif
  errnum = dzero
  errden = done

  if(info.ne.0) then
    info=4011
    err=info
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  restart: do 
!!$   
!!$   r0 = b-ax0
!!$ 
    if (itx.ge.litmax) exit restart  
    it = 0      
    call psb_geaxpby(done,b,dzero,r,desc_a,info)
    if (info == 0) call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
    if (debug_level >= psb_debug_ext_)&
         & write(debug_unit,*) me,' ',trim(name),' Done spmm',info
    if (info == 0) call psb_geaxpby(done,r,dzero,rt,desc_a,info)
    if(info.ne.0) then
      info=4011
      call psb_errpush(info,name)
      goto 9999
    end if

    rho = dzero
    if (debug_level >= psb_debug_ext_)&
         & write(debug_unit,*) me,' ',trim(name),'on entry to amax: b: ',size(b)
    if (istop_ == 1) then 
      rni = psb_geamax(r,desc_a,info)
      xni = psb_geamax(x,desc_a,info)
    else if (istop_ == 2) then 
      rni = psb_genrm2(r,desc_a,info)
    endif
    if(info.ne.0) then
      info=4011
      call psb_errpush(info,name)
      goto 9999
    end if

    if (istop_ == 1) then 
      xni  = psb_geamax(x,desc_a,info)
      errnum = rni
      errden = (ani*xni+bni)
    else  if (istop_ == 2) then 
      errnum = rni
      errden = bn2
    endif

    if(info.ne.0) then
      info=4011
      call psb_errpush(info,name)
      goto 9999
    end if

    if (errnum <= eps*errden) Then 
      exit restart
    end if
    If (itrace_ > 0) then 
      if ((mod(itx,itrace_)==0).and.(me == 0))&
           & write(*,'(a,i4,3(2x,es10.4))') 'bicg: ',itx,errnum,eps*errden
    end If


    iteration:  do 
      it   = it + 1
      itx = itx + 1
      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),'iteration: ',itx

      call psb_precaply(prec,r,z,desc_a,info,work=aux)
      call psb_precaply(prec,rt,zt,desc_a,info,trans='t',work=aux)

      rho_old = rho    
      rho = psb_gedot(rt,z,desc_a,info)
      if (rho==dzero) then
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' iteration breakdown r',rho
        exit iteration
      endif

      if (it==1) then
        call psb_geaxpby(done,z,dzero,p,desc_a,info)
        call psb_geaxpby(done,zt,dzero,pt,desc_a,info)
      else
        beta = (rho/rho_old)
        call psb_geaxpby(done,z,beta,p,desc_a,info)
        call psb_geaxpby(done,zt,beta,pt,desc_a,info)
      end if

      call psb_spmm(done,a,p,dzero,q,desc_a,info,&
           & work=aux)
      call psb_spmm(done,a,pt,dzero,qt,desc_a,info,&
           & work=aux,trans='t')

      sigma = psb_gedot(pt,q,desc_a,info)
      if (sigma==dzero) then
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' iteration breakdown s1', sigma
        exit iteration
      endif

      alpha = rho/sigma


      call psb_geaxpby(alpha,p,done,x,desc_a,info)
      call psb_geaxpby(-alpha,q,done,r,desc_a,info)
      call psb_geaxpby(-alpha,qt,done,rt,desc_a,info)


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
      If (errnum <= eps*errden) Then 
        exit restart
      end if

      if (itx.ge.litmax) exit restart

      If (itrace_ > 0) then 
        if ((mod(itx,itrace_)==0).and.(me == 0))&
             & write(*,'(a,i4,3(2x,es10.4))') 'bicg: ',itx,errnum,eps*errden
      end If
    end do iteration
  end do restart
  If (itrace_ > 0) then 
    if (me == 0) write(*,'(a,i4,3(2x,es10.4))') 'bicg: ',itx,errnum,eps*errden
  end If

  if (present(err)) then 
    if (errden /= dzero) then 
      err = errnum/errden
    else
      err = errnum
    end if
  end if

  if (present(iter)) iter = itx
  If ((errnum > eps*errden).and.(me==0)) Then
    write(debug_unit,*) 'bicg failed to converge to ',eps*errden,&
         & ' in ',itx,' iterations  '
  end if


  deallocate(aux)
  call psb_gefree(wwrk,desc_a,info)

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

end subroutine psb_dbicg


