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
!         software without specific prior written permission.
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
! File:  psb_drgmres.f90
!
! Subroutine: psb_drgmres
!    This subroutine implements the restarted GMRES method with right
!    preconditioning.
!
! Arguments:
!
!    a      -  type(psb_dspmat_type)      Input: sparse matrix containing A.
!    prec   -  class(psb_dprec_type)       Input: preconditioner
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

subroutine psb_drgmres_vect(a,prec,b,x,eps,desc_a,info,&
     & itmax,iter,err,itrace,irst,istop)
  use psb_base_mod
  use psb_prec_mod
  use psb_d_linsolve_conv_mod
  use psb_linsolve_mod
  implicit none
  type(psb_dspmat_type), intent(in)    :: a
  Type(psb_desc_type), Intent(in)      :: desc_a
  class(psb_dprec_type), intent(inout) :: prec
  type(psb_d_vect_type), Intent(inout) :: b
  type(psb_d_vect_type), Intent(inout) :: x
  Real(psb_dpk_), Intent(in)           :: eps
  integer(psb_ipk_), intent(out)                 :: info
  integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace, irst,istop
  integer(psb_ipk_), Optional, Intent(out)       :: iter
  Real(psb_dpk_), Optional, Intent(out) :: err
! =   local data
  real(psb_dpk_), allocatable   :: aux(:)
  real(psb_dpk_), allocatable   :: c(:), s(:), h(:,:), rs(:), rst(:)
  type(psb_d_vect_type), allocatable :: v(:)
  type(psb_d_vect_type)              :: w, w1, xt
  real(psb_dpk_) :: tmp 
  real(psb_dpk_) :: scal, gm, rti, rti1
  integer(psb_ipk_) ::itmax_, naux, it, k, itrace_,&
       & n_row, n_col, nl
  integer(psb_lpk_) :: mglob
  Logical, Parameter :: exchange=.True., noexchange=.False., use_srot=.true.
  integer(psb_ipk_), Parameter :: irmax = 8
  integer(psb_ipk_) :: itx, i, istop_, err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me  
  Real(psb_dpk_)     :: rni, xni, bni, ani,bn2, dt, r0n2
  real(psb_dpk_)     :: errnum, errden, deps, derr
  character(len=20)           :: name
  character(len=*), parameter :: methdname='RGMRES'

  info = psb_success_
  name = 'psb_dgmres'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt = desc_a%get_context()
  Call psb_info(ctxt, me, np)
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
    istop_ = psb_get_istop_default()
  endif
  
  if (.not.psb_is_valid_istop(istop_)) then
    info=psb_err_invalid_istop_
    err=info
    call psb_errpush(info,name,i_err=(/istop_/))
    goto 9999
  end if
  !
  !  istop_ = 1:  normwise backward error, infinity norm 
  !  istop_ = 2:  ||r||/||b||   norm 2
  !
  select case(istop_)
  case(psb_istop_ani_,psb_istop_bn2_,&
       & psb_istop_rn2_abs_,psb_istop_rrn2_)
    ! nothing needed
  case default
    ! should never get here
    info=psb_err_internal_error_
    err=info
    call psb_errpush(info,name,a_err="invalid istop_")
    goto 9999
  end select

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
    info=psb_err_invalid_irst_
    err=info
    call psb_errpush(info,name,i_err=(/nl/))
    goto 9999
  endif

  call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_chkvect on X')
    goto 9999
  end if
  call psb_chkvect(mglob,lone,b%get_nrows(),lone,lone,desc_a,info)
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


  select case(istop_)
  case(psb_istop_ani_)
    ani = psb_spnrmi(a,desc_a,info)
    bni = psb_geamax(b,desc_a,info)
  case(psb_istop_bn2_)
    bn2 = psb_genrm2(b,desc_a,info)    
  case(psb_istop_rn2_abs_)
    ! do nothing
  case(psb_istop_rrn2_)
    call psb_geaxpby(done,b,dzero,v(1),desc_a,info)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    call psb_spmm(-done,a,x,done,v(1),desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    r0n2 = psb_genrm2(v(1),desc_a,info)
  end select
  
  errnum = dzero
  errden = done
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
    call psb_geaxpby(done,b,dzero,v(1),desc_a,info)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_spmm(-done,a,x,done,v(1),desc_a,info,work=aux)
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if

    rs(1) = psb_genrm2(v(1),desc_a,info)
    rs(2:) = dzero
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    scal=done/rs(1)  ! rs(1) MIGHT BE VERY SMALL - USE DSCAL TO DEAL WITH IT?

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' on entry to amax: b: ',b%get_nrows(),rs(1),scal

    !
    ! check convergence
    !
    select case(istop_)
    case(psb_istop_ani_)
      rni = psb_geamax(v(1),desc_a,info)
      xni = psb_geamax(x,desc_a,info)
      errnum = rni
      errden = (ani*xni+bni)
    case(psb_istop_bn2_)
      rni = psb_genrm2(v(1),desc_a,info)
      errnum = rni
      errden = bn2
    case(psb_istop_rn2_abs_)
      rni = psb_genrm2(v(1),desc_a,info)
      errnum = rni
      errden = done
    case(psb_istop_rrn2_) 
      rni = psb_genrm2(v(1),desc_a,info)
      errnum = rni
      errden = r0n2
    end select
    if (info /= psb_success_) then 
      info=psb_err_from_subroutine_non_ 
      call psb_errpush(info,name)
      goto 9999
    end if
    
    if ((errnum <= eps*errden).or.(itx >= itmax_)) exit restart  

    if (itrace_ > 0) &
         & call log_conv(methdname,me,itx,itrace_,errnum,errden,deps)
     
    call v(1)%scal(scal) !v(1) = v(1) * scal

    !
    ! inner iterations
    !
    inner:  Do i=1,nl
      itx  = itx + 1

      call prec%apply(v(i),w1,desc_a,info)
      call psb_spmm(done,a,w1,dzero,w,desc_a,info,work=aux)
      !
      call mgs(i,h,v,w,rs,c,s,desc_a,info)

      
      select case(istop_)
      case(psb_istop_ani_)
        !
        ! build x and then compute the residual and its infinity norm
        !
        rst = rs
        call psb_geaxpby(done,x,dzero,xt,desc_a,info)
        call rebuildx(i,h,v,w,w1,xt,rst,c,s,prec,desc_a,info)

        call psb_geaxpby(done,b,dzero,w1,desc_a,info)
        call psb_spmm(-done,a,xt,done,w1,desc_a,info,work=aux)
        rni = psb_geamax(w1,desc_a,info)
        xni = psb_geamax(xt,desc_a,info)
        errnum = rni
        errden = (ani*xni+bni)
        !
      case(psb_istop_bn2_)
        !
        ! compute the residual 2-norm as byproduct of the solution
        ! procedure of the least-squares problem
        !
        rni = abs(rs(i+1))
        errnum = rni
        errden = bn2
        
      case(psb_istop_rn2_abs_)
        rni = abs(rs(i+1))
        errnum = rni
        errden = done

      case(psb_istop_rrn2_)
        
        !
        ! compute the residual 2-norm as byproduct of the solution
        ! procedure of the least-squares problem
        !
        rni = abs(rs(i+1))
        errnum = rni
        errden = r0n2
      end select

      if (errnum <= eps*errden) then 

        select case(istop_)
        case(psb_istop_ani_)
          call psb_geaxpby(done,xt,dzero,x,desc_a,info)
! =          x = xt 
        case(psb_istop_bn2_, psb_istop_rn2_abs_,psb_istop_rrn2_)
          call  rebuildx(i,h,v,w,w1,x,rs,c,s,prec,desc_a,info)          !

        end select

        if (itrace_ > 0) &
             & call log_conv(methdname,me,itx,ione,errnum,errden,deps)
        exit restart

      end if

      if (itrace_ > 0) &
           & call log_conv(methdname,me,itx,itrace_,errnum,errden,deps)

    end do inner

    select case(istop_)
    case(psb_istop_ani_)
      call psb_geaxpby(done,xt,dzero,x,desc_a,info)!      x = xt 

    case(psb_istop_bn2_, psb_istop_rn2_abs_,psb_istop_rrn2_)
      call  rebuildx(nl,h,v,w,w1,x,rs,c,s,prec,desc_a,info)          !
    end select
    
    if (itx >= itmax_) then 
      if (itrace_ > 0) then 
        if (mod(itx,itrace_)/=0) &
             &  call log_conv(methdname,me,itx,ione,errnum,errden,deps)
      end if
      exit restart
    end if
    
  end do restart

  call log_end(methdname,me,itx,itrace_,errnum,errden,deps,err=derr,iter=iter)
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
contains
  !
  ! Advance one step the Gram-Schmidt process, and apply
  ! Givens rotations to keep H triangular
  !
  subroutine mgs(n,h,v,w,rs,c,s,desc_a,info)
    real(psb_dpk_)   :: c(:), s(:), h(:,:), rs(:)
    type(psb_d_vect_type) :: v(:), w
    type(psb_desc_type) :: desc_a
    real(psb_dpk_) :: scal, gm, rti, rti1
    real(psb_dpk_) :: tmp
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: k,n
    !
    
    do k = 1, n
      h(k,n) = psb_gedot(v(k),w,desc_a,info)
      call psb_geaxpby(-h(k,n),v(k),done,w,desc_a,info)
    end do
    h(n+1,n) = psb_genrm2(w,desc_a,info)
    scal=done/h(n+1,n)
    call psb_geaxpby(scal,w,dzero,v(n+1),desc_a,info)
    do k=2,n
      call drot(1,h(k-1:k-1,n),1,h(k:k,n),1,real(c(k-1),kind=psb_dpk_),s(k-1))
    enddo
    
    rti  = h(n,n)
    rti1 = h(n+1,n) 
    call drotg(rti,rti1,tmp,s(n))
    c(n) = cmplx(tmp,dzero)
    call drot(1,h(n:n,n),1,h(n+1:n+1,n),1,real(c(n),kind=psb_dpk_),s(n))
    call drot(1,rs(n:n),1,rs(n+1:n+1),1,real(c(n),kind=psb_dpk_),s(n))
  end subroutine mgs

  !
  ! Rebuild solution X from the space V using the factor
  ! stored in R
  !
  subroutine rebuildx(n,h,v,w,w1,x,rs,c,s,prec,desc_a,info)
    real(psb_dpk_)   :: c(:), s(:), rs(:), h(:,:)
    type(psb_d_vect_type) :: v(:), w, w1, x
    type(psb_desc_type) :: desc_a
    class(psb_dprec_type) :: prec
    integer(psb_ipk_) :: info
    integer(psb_ipk_) :: k,n

    !
    !
    ! build x
    !
    call dtrsm('l','u','n','n',n,1,done,h,size(h,1),rs,size(rs,1))
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),&
         & ' Rebuild x-> RS:',rs(1:n)
    call w1%zero()
    do k=1, n
      call psb_geaxpby(rs(k),v(k),done,w1,desc_a,info)
    end do
    call prec%apply(w1,w,desc_a,info)
    call psb_geaxpby(done,w,done,x,desc_a,info)
  end subroutine rebuildx

end subroutine psb_drgmres_vect

