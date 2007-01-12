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
! File:  psb_dcg.f90
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
! File:  psb_dcg.f90
!
! Subroutine: psb_dcg
!    This subroutine implements the Conjugate Gradient method.
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
Subroutine psb_dcg(a,prec,b,x,eps,desc_a,info,&
     &itmax,iter,err, itrace, istop)
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
  Integer, Optional, Intent(in)      :: itmax, itrace, istop
  Integer, Optional, Intent(out)     :: iter
  Real(Kind(1.d0)), Optional, Intent(out) :: err
!!$   Local data
  real(kind(1.d0)), allocatable, target   :: aux(:), wwrk(:,:)
  real(kind(1.d0)), pointer  :: q(:), p(:), r(:), z(:), w(:)
  real(kind(1.d0))    ::rerr
  real(kind(1.d0))    ::alpha, beta, rho, rho_old, rni, xni, bni, ani,bn2,& 
       & sigma
  integer         :: litmax, liter, istop_, naux, m, mglob, it, itx, itrace_,&
       & np,me, n_col, isvch, ich, ictxt, n_row,err_act, int_err(5)
  character          :: diagl, diagu
  logical, parameter :: exchange=.true., noexchange=.false.  
  character(len=20)             :: name

  info = 0
  name = 'psb_dcg'
  call psb_erractionsave(err_act)


  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)

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

!!$  If ((prec%prec < min_prec_).Or.(prec%prec > max_prec_) ) Then
!!$    Write(0,*) 'F90_CG: Invalid IPREC',prec%prec
!!$    If (Present(ierr)) ierr=-1
!!$    Return
!!$  Endif

  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
    write(0,*) 'psb_cg: invalid istop',istop_ 
    info=5001
    int_err(1)=istop_
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  naux=4*n_col
  allocate(aux(naux), stat=info)
  if (info == 0) call psb_geall(wwrk,desc_a,info,n=5)
  if (info == 0) call psb_geasb(wwrk,desc_a,info)  
  if (info.ne.0) then 
    info=4011
    call psb_errpush(info,name)
    goto 9999
  end if

  p  => wwrk(:,1)
  q  => wwrk(:,2)
  r  => wwrk(:,3)
  z  => wwrk(:,4) 
  w  => wwrk(:,5)


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

  itx=0

  ! Ensure global coherence for convergence checks.
  call psb_set_coher(ictxt,isvch)

  restart: do 
!!$   
!!$    r0 = b-Ax0
!!$   
    if (itx>= litmax) exit restart 
    it = 0
    call psb_geaxpby(done,b,dzero,r,desc_a,info)
    call psb_spmm(-done,a,x,done,r,desc_a,info,work=aux)
    if (info.ne.0) then 
      info=4011
      call psb_errpush(info,name)
      goto 9999
    end if

    rho = dzero
    if (istop_ == 1) then 
      ani = psb_spnrmi(a,desc_a,info)
      bni = psb_geamax(b,desc_a,info)
    else if (istop_ == 2) then 
      bn2 = psb_genrm2(b,desc_a,info)
    endif
    if (info.ne.0) then 
      info=4011
      call psb_errpush(info,name)
      goto 9999
    end if


    iteration:  do 
      it   = it + 1
      itx = itx + 1

      Call psb_prc_aply(prec,r,z,desc_a,info,work=aux)
      rho_old = rho
      rho     = psb_gedot(r,z,desc_a,info)

      if (it==1) then
        call psb_geaxpby(done,z,dzero,p,desc_a,info)
      else
        if (rho_old==dzero) then
          write(0,*) 'CG Iteration breakdown'
          exit iteration
        endif
        beta = rho/rho_old
        call psb_geaxpby(done,z,beta,p,desc_a,info)
      end if

      call psb_spmm(done,a,p,dzero,q,desc_a,info,work=aux)
      sigma = psb_gedot(p,q,desc_a,info)
      if (sigma==dzero) then
        write(0,*) 'CG Iteration breakdown'
        exit iteration
      endif

      alpha = rho/sigma
      call psb_geaxpby(alpha,p,done,x,desc_a,info)
      call psb_geaxpby(-alpha,q,done,r,desc_a,info)


      if (istop_ == 1) Then 
        rni = psb_geamax(r,desc_a,info)
        xni = psb_geamax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)
      Else  If (istop_ == 2) Then 

        rni = psb_genrm2(r,desc_a,info)
        rerr = rni/bn2
      Endif
      if (rerr<=eps) exit restart

      if (itx>= litmax) exit restart 

      If (itrace_ > 0) then 
        if ((mod(itx,itrace_)==0).and.(me == 0))&
             & write(*,'(a,i4,3(2x,es10.4))') 'cg: ',itx,rerr
      end If
    end do iteration
  end do restart
  If (itrace_ > 0) then 
    if (me == 0) write(*,'(a,i4,3(2x,es10.4))') 'cg: ',itx,rerr
  end If

  if (present(err)) err=rerr
  if (present(iter)) iter = itx
  if (rerr>eps) then
    write(0,*) 'CG Failed to converge to ',eps,&
         & ' in ',litmax,' iterations '
    info=itx
  end if

  deallocate(aux)
  call psb_gefree(wwrk,desc_a,info)

  ! restore external global coherence behaviour
  call psb_restore_coher(ictxt,isvch)

  if (info.ne.0) then 
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dcg


