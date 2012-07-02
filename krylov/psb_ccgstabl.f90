!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File:  psb_ccgstabl.f90
!
! Subroutine: psb_ccgstabl
!   Implements the BICGSTAB(L) method
!
!
! Arguments:
!
!    a      -  type(psb_cspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_cprec_type)       Input: preconditioner
!    b(:)   -  complex                    Input: vector containing the
!                                         right hand side B
!    x(:)   -  complex                    Input/Output: vector containing the
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
!                                         1: err =  |r|/(|a||x|+|b|);  here the iteration
!                                            is stopped when  |r| <= eps * (|a||x|+|b|)
!                                         2: err =  |r|/|b|; here the iteration is
!                                            stopped when  |r| <= eps * |b|
!                                         where r is the (preconditioned, recursive
!                                         estimate of) residual. 
!    irst   -  integer(optional)          Input: restart parameter L 
!
!
!
!!$Subroutine psb_ccgstabl(a,prec,b,x,eps,desc_a,info,itmax,iter,err,itrace,irst,istop)
!!$  use psb_base_mod
!!$  use psb_prec_mod
!!$  use psb_c_krylov_conv_mod
!!$  use psb_krylov_mod
!!$  implicit none
!!$
!!$! =  parameters 
!!$  Type(psb_cspmat_type), Intent(in)  :: a
!!$  class(psb_cprec_type), Intent(in)   :: prec 
!!$  Type(psb_desc_type), Intent(in)    :: desc_a
!!$  complex(psb_spk_), Intent(in)    :: b(:)
!!$  complex(psb_spk_), Intent(inout) :: x(:)
!!$  Real(psb_spk_), Intent(in)       :: eps
!!$  integer(psb_ipk_), intent(out)               :: info
!!$  integer(psb_ipk_), Optional, Intent(in)      :: itmax, itrace, irst,istop
!!$  integer(psb_ipk_), Optional, Intent(out)     :: iter
!!$  Real(psb_spk_), Optional, Intent(out) :: err
!!$! =   local data
!!$  complex(psb_spk_), allocatable, target   :: aux(:),wwrk(:,:),uh(:,:), rh(:,:)
!!$  complex(psb_spk_), Pointer  :: ww(:), q(:), r(:), rt0(:), p(:), v(:), &
!!$       & s(:), t(:), z(:), f(:), gamma(:), gamma1(:), gamma2(:), taum(:,:), sigma(:)
!!$
!!$  integer(psb_ipk_) :: itmax_, naux, mglob, it, itrace_,&
!!$       & np,me, n_row, n_col, nl, err_act
!!$  Logical, Parameter :: exchange=.True., noexchange=.False.  
!!$  integer(psb_ipk_), Parameter :: irmax = 8
!!$  integer(psb_ipk_) :: itx, i, isvch, ictxt,istop_,j, int_err(5)
!!$  integer(psb_ipk_) :: debug_level, debug_unit
!!$  complex(psb_spk_) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,bn2,& 
!!$       & omega
!!$  type(psb_itconv_type)        :: stopdat
!!$  real(psb_dpk_)               :: derr
!!$  character(len=20)            :: name
!!$  character(len=*), parameter  :: methdname='BiCGStab(L)'
!!$
!!$  info = psb_success_
!!$  name = 'psb_ccgstabl'
!!$  call psb_erractionsave(err_act)
!!$  debug_unit  = psb_get_debug_unit()
!!$  debug_level = psb_get_debug_level()
!!$
!!$  ictxt = desc_a%get_context()
!!$  Call psb_info(ictxt, me, np)
!!$  if (debug_level >= psb_debug_ext_)&
!!$       & write(debug_unit,*) me,' ',trim(name),': from psb_info',np
!!$
!!$
!!$  mglob = desc_a%get_global_rows()
!!$  n_row = desc_a%get_local_rows()
!!$  n_col = desc_a%get_local_cols()
!!$
!!$  if (present(istop)) then 
!!$    istop_ = istop 
!!$  else
!!$    istop_ = 2
!!$  endif
!!$
!!$  if (present(itmax)) then 
!!$    itmax_ = itmax
!!$  else
!!$    itmax_ = 1000
!!$  endif
!!$
!!$  if (present(itrace)) then
!!$     itrace_ = itrace
!!$  else
!!$     itrace_ = 0
!!$  end if
!!$  
!!$  if (present(irst)) then
!!$    nl = irst
!!$    if (debug_level >= psb_debug_ext_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & 'present: irst: ',irst,nl
!!$  else
!!$    nl = 1 
!!$    if (debug_level >= psb_debug_ext_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & ' not present: irst: ',irst,nl
!!$  endif
!!$  if (nl <=0 ) then 
!!$    info=psb_err_invalid_istop_
!!$    int_err(1)=nl
!!$    err=info
!!$    call psb_errpush(info,name,i_err=int_err)
!!$    goto 9999
!!$  endif
!!$
!!$  call psb_chkvect(mglob,ione,size(x,ione),ione,ione,desc_a,info)
!!$  if (info == psb_success_) call psb_chkvect(mglob,ione,size(b,ione),ione,ione,desc_a,info)
!!$  if (info /= psb_success_) then
!!$    info=psb_err_from_subroutine_    
!!$    call psb_errpush(info,name,a_err='psb_chkvect on X/B')
!!$    goto 9999
!!$  end if
!!$
!!$  naux=4*n_col 
!!$  allocate(aux(naux),gamma(0:nl),gamma1(nl),&
!!$       &gamma2(nl),taum(nl,nl),sigma(nl), stat=info)
!!$
!!$  if (info /= psb_success_) then 
!!$     info=psb_err_alloc_dealloc_
!!$     call psb_errpush(info,name)
!!$     goto 9999
!!$  end if
!!$  if (info == psb_success_) Call psb_geall(wwrk,desc_a,info,n=psb_err_iarg_neg_)
!!$  if (info == psb_success_) Call psb_geall(uh,desc_a,info,n=nl+1,lb=0)
!!$  if (info == psb_success_) Call psb_geall(rh,desc_a,info,n=nl+1,lb=0)
!!$  if (info == psb_success_) Call psb_geasb(wwrk,desc_a,info)  
!!$  if (info == psb_success_) Call psb_geasb(uh,desc_a,info)  
!!$  if (info == psb_success_) Call psb_geasb(rh,desc_a,info)  
!!$  if (info /= psb_success_) then 
!!$     info=psb_err_from_subroutine_non_ 
!!$     call psb_errpush(info,name)
!!$     goto 9999
!!$  end if
!!$
!!$  q   => wwrk(:,1)
!!$  r   => wwrk(:,2)
!!$  p   => wwrk(:,3)
!!$  v   => wwrk(:,4)
!!$  f   => wwrk(:,5)
!!$  s   => wwrk(:,6)
!!$  t   => wwrk(:,7)
!!$  z   => wwrk(:,8)
!!$  ww  => wwrk(:,9)
!!$  rt0 => wwrk(:,10)
!!$  
!!$
!!$  call psb_init_conv(methdname,istop_,itrace_,itmax_,a,b,eps,desc_a,stopdat,info)
!!$  if (info /= psb_success_) Then 
!!$     call psb_errpush(psb_err_from_subroutine_non_,name)
!!$     goto 9999
!!$  End If
!!$
!!$  itx   = 0
!!$  restart: do 
!!$! =   
!!$! =   r0 = b-ax0
!!$! = 
!!$    if (debug_level >= psb_debug_ext_) &
!!$         & write(debug_unit,*) me,' ',trim(name),' restart: ',itx,it
!!$    if (itx >= itmax_) exit restart  
!!$
!!$    it = 0      
!!$    call psb_geaxpby(cone,b,czero,r,desc_a,info)
!!$    if (info == psb_success_) call psb_spmm(-cone,a,x,cone,r,desc_a,info,work=aux)
!!$    
!!$    if (info == psb_success_) call prec%apply(r,desc_a,info)
!!$
!!$    if (info == psb_success_) call psb_geaxpby(cone,r,czero,rt0,desc_a,info)
!!$    if (info == psb_success_) call psb_geaxpby(cone,r,czero,rh(:,0),desc_a,info)
!!$    if (info == psb_success_) call psb_geaxpby(czero,r,czero,uh(:,0),desc_a,info)
!!$    if (info /= psb_success_) then 
!!$       info=psb_err_from_subroutine_non_ 
!!$       call psb_errpush(info,name)
!!$       goto 9999
!!$    end if
!!$   
!!$    rho   = cone
!!$    alpha = czero
!!$    omega = cone 
!!$
!!$    if (debug_level >= psb_debug_ext_) &
!!$         & write(debug_unit,*) me,' ',trim(name),&
!!$         & ' on entry to amax: b: ',Size(b)
!!$
!!$    if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
!!$    if (info /= psb_success_) Then 
!!$      call psb_errpush(psb_err_from_subroutine_non_,name)
!!$      goto 9999
!!$    End If
!!$    
!!$    iteration:  do 
!!$      it   = it + nl
!!$      itx  = itx + nl
!!$      rho = -omega*rho 
!!$
!!$      if (debug_level >= psb_debug_ext_) &
!!$           & write(debug_unit,*) me,' ',trim(name),&
!!$           & ' iteration: ',itx, rho,rh(1,0)
!!$
!!$      do j = 0, nl -1 
!!$        If (debug_level >= psb_debug_ext_) &
!!$             & write(debug_unit,*) me,' ',trim(name),'bicg part:  ',j, nl
!!$
!!$        rho_old = rho
!!$        rho = psb_gedot(rh(:,j),rt0,desc_a,info)
!!$        if (rho == czero) then
!!$          if (debug_level >= psb_debug_ext_) &
!!$               & write(debug_unit,*) me,' ',trim(name),&
!!$               & ' bi-cgstab iteration breakdown r',rho
!!$          exit iteration
!!$        endif
!!$
!!$        beta = alpha*rho/rho_old 
!!$        rho_old = rho
!!$        call psb_geaxpby(cone,rh(:,0:j),-beta,uh(:,0:j),desc_a,info)
!!$        call psb_spmm(cone,a,uh(:,j),czero,uh(:,j+1),desc_a,info,work=aux)
!!$
!!$        call prec%apply(uh(:,j+1),desc_a,info)
!!$
!!$        gamma(j) = psb_gedot(uh(:,j+1),rt0,desc_a,info)
!!$
!!$        if (gamma(j) == czero) then
!!$          if (debug_level >= psb_debug_ext_) &
!!$               & write(debug_unit,*) me,' ',trim(name),&
!!$               & ' bi-cgstab iteration breakdown s2',gamma(j)
!!$          exit iteration
!!$        endif
!!$        alpha = rho/gamma(j)
!!$        if (debug_level >= psb_debug_ext_) &
!!$             & write(debug_unit,*) me,' ',trim(name),&
!!$             & ' bicg part: alpha=r/g ',alpha,rho,gamma(j)
!!$
!!$        call psb_geaxpby(-alpha,uh(:,1:j+1),cone,rh(:,0:j),desc_a,info)        
!!$        call psb_geaxpby(alpha,uh(:,0),cone,x,desc_a,info)
!!$        call psb_spmm(cone,a,rh(:,j),czero,rh(:,j+1),desc_a,info,work=aux)
!!$
!!$        call prec%apply(rh(:,j+1),desc_a,info)
!!$                
!!$      enddo
!!$      
!!$      do j=1, nl 
!!$        if (debug_level >= psb_debug_ext_) &
!!$             & write(debug_unit,*) me,' ',trim(name),&
!!$             & ' mod g-s part:  ',j, nl,rh(1,0)
!!$
!!$        do i=1, j-1 
!!$          taum(i,j) = psb_gedot(rh(:,i),rh(:,j),desc_a,info)
!!$          taum(i,j) = taum(i,j)/sigma(i) 
!!$          call psb_geaxpby(-taum(i,j),rh(:,i),cone,rh(:,j),desc_a,info)        
!!$        enddo        
!!$        sigma(j)  = psb_gedot(rh(:,j),rh(:,j),desc_a,info)
!!$        gamma1(j) = psb_gedot(rh(:,0),rh(:,j),desc_a,info)
!!$        gamma1(j) = gamma1(j)/sigma(j)
!!$      enddo
!!$      
!!$      gamma(nl) = gamma1(nl) 
!!$      omega     = gamma(nl) 
!!$
!!$      do j=nl-1,1,-1
!!$        gamma(j) = gamma1(j)
!!$        do i=j+1,nl
!!$          gamma(j) = gamma(j) - taum(j,i) * gamma(i) 
!!$        enddo
!!$      enddo
!!$
!!$      do j=1,nl-1
!!$        gamma2(j) = gamma(j+1)
!!$        do i=j+1,nl-1
!!$          gamma2(j) = gamma2(j) + taum(j,i) * gamma(i+1) 
!!$        enddo
!!$      enddo
!!$      
!!$      call psb_geaxpby(gamma(1),rh(:,0),cone,x,desc_a,info)        
!!$      call psb_geaxpby(-gamma1(nl),rh(:,nl),cone,rh(:,0),desc_a,info)        
!!$      call psb_geaxpby(-gamma(nl),uh(:,nl),cone,uh(:,0),desc_a,info)        
!!$
!!$      do j=1, nl-1
!!$        call psb_geaxpby(-gamma(j),uh(:,j),cone,uh(:,0),desc_a,info)        
!!$        call psb_geaxpby(gamma2(j),rh(:,j),cone,x,desc_a,info)        
!!$        call psb_geaxpby(-gamma1(j),rh(:,j),cone,rh(:,0),desc_a,info)        
!!$      enddo
!!$      
!!$      if (psb_check_conv(methdname,itx,x,rh(:,0),desc_a,stopdat,info)) exit restart
!!$      if (info /= psb_success_) Then 
!!$        call psb_errpush(psb_err_from_subroutine_non_,name)
!!$        goto 9999
!!$      End If
!!$      
!!$    end do iteration
!!$  end do restart
!!$
!!$  call psb_end_conv(methdname,itx,desc_a,stopdat,info,derr,iter)
!!$
!!$  if (present(err)) then 
!!$    err = derr
!!$  end if
!!$
!!$  deallocate(aux,stat=info)
!!$  if (info == psb_success_) call psb_gefree(wwrk,desc_a,info)
!!$  if (info == psb_success_) call psb_gefree(uh,desc_a,info)
!!$  if (info == psb_success_) call psb_gefree(rh,desc_a,info)
!!$  if (info /= psb_success_) then
!!$     call psb_errpush(info,name)
!!$     goto 9999
!!$  end if
!!$
!!$  call psb_erractionrestore(err_act)
!!$  return
!!$
!!$9999 continue
!!$  call psb_erractionrestore(err_act)
!!$  if (err_act == psb_act_abort_) then
!!$     call psb_error()
!!$     return
!!$  end if
!!$  return
!!$
!!$End Subroutine psb_ccgstabl
!!$


Subroutine psb_ccgstabl_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,irst,istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_c_krylov_conv_mod
  use psb_krylov_mod
  implicit none
  type(psb_cspmat_type), intent(in)    :: a
  class(psb_cprec_type), Intent(inout) :: prec 
  Type(psb_desc_type), Intent(in)      :: desc_a
  type(psb_c_vect_type), Intent(inout) :: b
  type(psb_c_vect_type), Intent(inout) :: x
  Real(psb_spk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, irst,istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  Real(psb_spk_), Optional, Intent(out) :: err
! =   local data
  complex(psb_spk_), allocatable, target   :: aux(:), gamma(:),&
       & gamma1(:), gamma2(:), taum(:,:), sigma(:)
  type(psb_c_vect_type), allocatable, target :: wwrk(:),uh(:), rh(:)
  type(psb_c_vect_type), Pointer  :: ww, q, r, rt0, p, v, &
       & s, t, z, f

  integer(psb_ipk_) :: itmax_, naux, mglob, it, itrace_,&
       & n_row, n_col, nl, err_act
  Logical, Parameter :: exchange=.True., noexchange=.False.  
  integer(psb_ipk_), Parameter :: irmax = 8
  integer(psb_ipk_) :: itx, i, istop_,j, k, int_err(5)
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: ictxt, np, me
  complex(psb_spk_) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,bn2,& 
       & omega
  real(psb_dpk_)     :: derr  
  type(psb_itconv_type)        :: stopdat
  character(len=20)            :: name
  character(len=*), parameter  :: methdname='BiCGStab(L)'

  info = psb_success_
  name = 'psb_ccgstabl'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()
  Call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit,*) me,' ',trim(name),': from psb_info',np
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(b%v)) then 
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

  if (present(itmax)) then 
    itmax_ = itmax
  else
    itmax_ = 1000
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
    info=psb_err_invalid_istop_
    int_err(1)=nl
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  call psb_chkvect(mglob,ione,x%get_nrows(),ione,ione,desc_a,info)
  if (info == psb_success_) call psb_chkvect(mglob,ione,b%get_nrows(),ione,ione,desc_a,info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_    
    call psb_errpush(info,name,a_err='psb_chkvect on X/B')
    goto 9999
  end if

  naux=4*n_col 
  allocate(aux(naux),gamma(0:nl),gamma1(nl),&
       &gamma2(nl),taum(nl,nl),sigma(nl), stat=info)

  if (info /= psb_success_) then 
     info=psb_err_alloc_dealloc_
     call psb_errpush(info,name)
     goto 9999
  end if
  if (info == psb_success_) Call psb_geall(wwrk,desc_a,info,n=10_psb_ipk_)
  if (info == psb_success_) Call psb_geall(uh,desc_a,info,n=nl+1,lb=izero)
  if (info == psb_success_) Call psb_geall(rh,desc_a,info,n=nl+1,lb=izero)
  if (info == psb_success_) Call psb_geasb(wwrk,desc_a,info,mold=x%v)  
  if (info == psb_success_) Call psb_geasb(uh,desc_a,info,mold=x%v)  
  if (info == psb_success_) Call psb_geasb(rh,desc_a,info,mold=x%v)    
  if (info /= psb_success_) then 
     info=psb_err_from_subroutine_non_ 
     call psb_errpush(info,name)
     goto 9999
  end if

  q   => wwrk(1)
  r   => wwrk(2)
  p   => wwrk(3)
  v   => wwrk(4)
  f   => wwrk(5)
  s   => wwrk(6)
  t   => wwrk(7)
  z   => wwrk(8)
  ww  => wwrk(9)
  rt0 => wwrk(10)
  

  call psb_init_conv(methdname,istop_,itrace_,itmax_,a,b,eps,desc_a,stopdat,info)
  if (info /= psb_success_) Then 
     call psb_errpush(psb_err_from_subroutine_non_,name)
     goto 9999
  End If

  itx   = 0
  restart: do 
! =   
! =   r0 = b-ax0
! = 
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),' restart: ',itx,it
    if (itx >= itmax_) exit restart  

    it = 0      
    call psb_geaxpby(cone,b,czero,r,desc_a,info)
    if (info == psb_success_) call psb_spmm(-cone,a,x,cone,r,desc_a,info,work=aux)
    
    if (info == psb_success_) call prec%apply(r,desc_a,info)

    if (info == psb_success_) call psb_geaxpby(cone,r,czero,rt0,desc_a,info)
    if (info == psb_success_) call psb_geaxpby(cone,r,czero,rh(0),desc_a,info)
    if (info == psb_success_) call psb_geaxpby(czero,r,czero,uh(0),desc_a,info)
    if (info /= psb_success_) then 
       info=psb_err_from_subroutine_non_ 
       call psb_errpush(info,name)
       goto 9999
    end if
   
    rho   = cone
    alpha = czero
    omega = cone 

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' on entry to amax: b: ',b%get_nrows()

    if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
    if (info /= psb_success_) Then 
      call psb_errpush(psb_err_from_subroutine_non_,name)
      goto 9999
    End If
    
    iteration:  do 
      it   = it  + nl
      itx  = itx + nl
      rho = -omega*rho 

      if (debug_level >= psb_debug_ext_) &
           & write(debug_unit,*) me,' ',trim(name),&
           & ' iteration: ',itx, rho

      do j = 0, nl -1 
        If (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),'bicg part:  ',j, nl

        rho_old = rho
        rho = psb_gedot(rh(j),rt0,desc_a,info)
        if (rho == czero) then
          if (debug_level >= psb_debug_ext_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' bi-cgstab iteration breakdown r',rho
          exit iteration
        endif

        beta = alpha*rho/rho_old 
        rho_old = rho
        do k=0, j
! =          call psb_geaxpby(cone,rh(:,0:j),-beta,uh(:,0:j),desc_a,info)
          call psb_geaxpby(cone,rh(k),-beta,uh(k),desc_a,info)
        end do
        call psb_spmm(cone,a,uh(j),czero,uh(j+1),desc_a,info,work=aux)

        call prec%apply(uh(j+1),desc_a,info)

        gamma(j) = psb_gedot(uh(j+1),rt0,desc_a,info)

        if (gamma(j) == czero) then
          if (debug_level >= psb_debug_ext_) &
               & write(debug_unit,*) me,' ',trim(name),&
               & ' bi-cgstab iteration breakdown s2',gamma(j)
          exit iteration
        endif
        alpha = rho/gamma(j)
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' bicg part: alpha=r/g ',alpha,rho,gamma(j)

        do k=0,j
! =        call psb_geaxpby(-alpha,uh(:,1:j+1),cone,rh(:,0:j),desc_a,info)        
          call psb_geaxpby(-alpha,uh(k+1),cone,rh(k),desc_a,info)        
        end do
        call psb_geaxpby(alpha,uh(0),cone,x,desc_a,info)
        call psb_spmm(cone,a,rh(j),czero,rh(j+1),desc_a,info,work=aux)

        call prec%apply(rh(j+1),desc_a,info)
                
      enddo
      
      do j=1, nl 
        if (debug_level >= psb_debug_ext_) &
             & write(debug_unit,*) me,' ',trim(name),&
             & ' mod g-s part:  ',j, nl

        do i=1, j-1 
          taum(i,j) = psb_gedot(rh(i),rh(j),desc_a,info)
          taum(i,j) = taum(i,j)/sigma(i) 
          call psb_geaxpby(-taum(i,j),rh(i),cone,rh(j),desc_a,info)        
        enddo        
        sigma(j)  = psb_gedot(rh(j),rh(j),desc_a,info)
        gamma1(j) = psb_gedot(rh(0),rh(j),desc_a,info)
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
      
      call psb_geaxpby(gamma(1),rh(0),cone,x,desc_a,info)        
      call psb_geaxpby(-gamma1(nl),rh(nl),cone,rh(0),desc_a,info)        
      call psb_geaxpby(-gamma(nl),uh(nl),cone,uh(0),desc_a,info)        

      do j=1, nl-1
        call psb_geaxpby(-gamma(j),uh(j),cone,uh(0),desc_a,info)        
        call psb_geaxpby(gamma2(j),rh(j),cone,x,desc_a,info)        
        call psb_geaxpby(-gamma1(j),rh(j),cone,rh(0),desc_a,info)        
      enddo
      
      if (psb_check_conv(methdname,itx,x,rh(0),desc_a,stopdat,info)) exit restart
      if (info /= psb_success_) Then 
        call psb_errpush(psb_err_from_subroutine_non_,name)
        goto 9999
      End If
      
    end do iteration
  end do restart

  call psb_end_conv(methdname,itx,desc_a,stopdat,info,derr,iter)
  if (present(err)) err = derr

  if (info == psb_success_) call psb_gefree(uh,desc_a,info)
  if (info == psb_success_) call psb_gefree(rh,desc_a,info)
  if (info == psb_success_) call psb_gefree(wwrk,desc_a,info)
  if (info == psb_success_) deallocate(aux,stat=info)
  if (info /= psb_success_) then
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
     call psb_error()
     return
  end if
  return

End Subroutine psb_ccgstabl_vect



