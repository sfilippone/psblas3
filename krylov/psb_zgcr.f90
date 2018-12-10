!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
! File:  psb_zgcr.f90
!!
!! Contributors: Ambra Abdullahi (UNITOV) and Pasqua D’Ambra (IAC-CNR)
!!
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
!   C             SIAM, 1993                                        
!   C                                                                      C  
!   C         [4] Notay, Yvan                                              C
!   C             Aggregation-based algebraic multigrid method             C
!   C             SIAM Journal on Scientific Computing 34,                 C 
!   C             pp. A2288-A2316, 2012                                    C    
!   C                                                                      C
!   C                                                                      C
!   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! File:  psb_zgcr.f90
!
! Subroutine: psb_zgcr
!    This subroutine implements the GCR method.
!
!
! Arguments:
!
!    a      -  type(psb_zspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_zprec_type)       Input: preconditioner
!    b(:)   -  real                    Input: vector containing the
!                                         right hand side B
!    x(:)   -  real                    Input/Output: vector containing the
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
! 
!    irst   -  integer(optional)          Input: restart parameter 
!

subroutine psb_zgcr_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace, irst, istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_z_krylov_conv_mod
  use psb_krylov_mod
  implicit none
  
  
  type(psb_zspmat_type), intent(in)    :: a
  Type(psb_desc_type), Intent(in)      :: desc_a
  class(psb_zprec_type), intent(inout) :: prec
  type(psb_z_vect_type), Intent(inout) :: b
  type(psb_z_vect_type), Intent(inout) :: x
  real(psb_dpk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, irst, istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  real(psb_dpk_), Optional, Intent(out) :: err
  ! =   local data
  complex(psb_dpk_), allocatable   :: alpha(:), h(:,:)
  type(psb_z_vect_type), allocatable :: z(:), c(:), c_scale(:)
  type(psb_z_vect_type)   ::  r
  
  real(psb_dpk_) :: r_norm, b_norm, a_norm, derr
  integer(psb_ipk_) :: n_col, naux, err_act
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: np, me, ictxt
  integer(psb_ipk_) ::  i, j, it, itx, istop_, itmax_, itrace_, nl, m, nrst
  complex(psb_dpk_) :: hjj
  complex(psb_dpk_), allocatable, target   :: aux(:)
  character(len=20)           :: name
  type(psb_itconv_type)       :: stopdat
  character(len=*), parameter :: methdname='GCR'
  info = psb_success_
  name = 'psb_zgcr'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  
  ictxt = desc_a%get_context()
  
  call psb_info(ictxt, me, np)
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
    err=info
    call psb_errpush(info,name,i_err=(/istop_/))
    goto 9999
  endif
  
  
  call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
  if (info == psb_success_)&
       & call psb_chkvect(mglob,lone,b%get_nrows(),lone,lone,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_    
    call psb_errpush(info,name,a_err='psb_chkvect on X/B')
    goto 9999
  end if
  
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
         & ' present: irst: ',irst,nl
  else
    nl = 10 
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' not present: irst: ',irst,nl
  endif
  
  if (nl <=0 ) then 
    info=psb_err_invalid_istop_
    err=info
    call psb_errpush(info,name,i_err=(/nl/))
    goto 9999
  endif
  
  naux=4*n_col 
  allocate(aux(naux),h(nl+1,nl+1),&
       &c_scale(nl+1),c(nl+1),z(nl+1), alpha(nl+1), stat=info)
  
  h = zzero
  if (info /= psb_success_) then 
    info=psb_err_from_subroutine_non_ 
    call psb_errpush(info,name)
    goto 9999
  end if
  
  call psb_geasb(r, desc_a,info, scratch=.true.,mold=x%v)
  
  do i =1,nl+1
    call psb_geasb(c(i), desc_a,info, scratch=.true.,mold=x%v)
    call psb_geasb(z(i), desc_a,info, scratch=.true.,mold=x%v)
    call psb_geasb(c_scale(i), desc_a,info, scratch=.true.,mold=x%v)
  end do
  
  itx = 0
  
  nrst = -1
  call psb_init_conv(methdname,istop_,itrace_,itmax_,a,x,b,eps,desc_a,stopdat,info)
  restart: do 
    if (itx>= itmax_) exit restart 
    h = zzero
    
    it = 0
    ! compute r0 = b-ax0
    
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    
    call psb_geaxpby(zone, b, zzero, r, desc_a, info) 
    call psb_spmm(-zone,a,x,zone,r,desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
        
    nrst = nrst + 1 
    
    iteration: do 
      
      itx = itx + 1
      it = it + 1
      j = it    
      !Apply preconditioner
      call prec%apply(r,z(j),desc_a,info,work=aux)  
      
      call psb_spmm(zone,a,z(j),zzero,c(1),desc_a,info,work=aux)
      do i =1, j - 1
        
        h(i,j) = psb_gedot(c_scale(i), c(i), desc_a, info)   
        
        call psb_geaxpby(zone, c(i), zzero, c(i+1), desc_a, info)
        call psb_geaxpby(-h(i,j), c_scale(i), zone, c(i+1), desc_a, info)   
      end do
      
      h(j,j) = psb_norm2(c(j), desc_a, info)
      hjj = zone/h(j,j)
      call psb_geaxpby(hjj, c(j), zzero, c_scale(j), desc_a, info)   
      
      alpha(j) = psb_gedot(c_scale(j), r, desc_a, info) 
      
      !Update residual
      call psb_geaxpby(zone, r, zzero, r, desc_a, info)   
      call psb_geaxpby(-alpha(j), c_scale(j), zone, r, desc_a, info)   
      
      if (psb_check_conv(methdname,itx,x,r,desc_a,stopdat,info)) exit restart
      
      if (j >= irst) exit iteration
      
      
    end do iteration
    
    m = j
    
    !Compute solution
    
    call ztrsm('l','u','n','n',m,1,zone,h,size(h,1),alpha,size(alpha,1))
    
    if (nrst == 0 ) then    
      call x%set(zzero)
    endif
    do i=1,m
      call psb_geaxpby(alpha(i), z(i), zone, x, desc_a, info)   
    enddo
    
    
    
    
  end do restart
  m = j
  !Compute solution
  call ztrsm('l','u','n','n',m,1,zone,h,size(h,1),alpha,size(alpha,1))
  call x%set(zzero)
  do i=1,m
    call psb_geaxpby(alpha(i), z(i), zone, x, desc_a, info)   
  enddo
  
  iter = j
  
  call psb_end_conv(methdname,itx,desc_a,stopdat,info,derr,iter)
  if (present(err)) err = derr
  
  if (info == psb_success_) call psb_gefree(r,desc_a,info)
  
  do j = 1,m
    if (info == psb_success_) call psb_gefree(z(j),desc_a,info)
    if (info == psb_success_) call psb_gefree(c_scale(j),desc_a,info)
  enddo
  
  do i =1,nl+1
    if (info == psb_success_) call psb_gefree(c(i),desc_a,info)   
  end do
  
  if (info == psb_success_) deallocate(aux,h,c_scale,z,c,alpha,stat=info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_non_
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
end subroutine psb_zgcr_vect


