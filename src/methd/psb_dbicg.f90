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
subroutine psb_dbicg(a,prec,b,x,eps,desc_a,info,&
     &itmax,iter,err, itrace,istop)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_tools_mod
  use psb_const_mod
  use psb_prec_mod
  use psb_error_mod
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
  real(kind(1.d0)), pointer  :: aux(:),wwrk(:,:)
  real(kind(1.d0)), pointer  :: ww(:), q(:),&
       & r(:), p(:), zt(:), pt(:), z(:), rt(:),qt(:)
  integer, pointer           :: iperm(:), ipnull(:), ipsave(:), int_err(:)
  real(kind(1.d0)) ::rerr
  integer       ::litmax, liter, naux, m, mglob, it, itrac,&
       & nprows,npcols,me,mecol, n_row, n_col, listop, err_act
  character     ::diagl, diagu
  logical, parameter :: debug = .false.
  logical, parameter :: exchange=.true., noexchange=.false.  
  integer, parameter :: irmax = 8
  integer            :: itx, i, isvch, ich, icontxt
  logical            :: do_renum_left
  real(kind(1.d0)), parameter :: one=1.d0, zero=0.d0, epstol=1.d-35
  real(kind(1.d0)) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,& 
       & sigma, omega, tau,bn2
  character(len=20)             :: name,ch_err

  info = 0
  name = 'psb_dbicg'
  call psb_erractionsave(err_act)

  if (debug) write(*,*) 'entering psb_dbicg'
  icontxt = desc_a%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprows,npcols,me,mecol)
  if (debug) write(*,*) 'psb_dbicg: from gridinfo',nprows,npcols,me

  mglob = desc_a%matrix_data(psb_m_)
  n_row = desc_a%matrix_data(psb_n_row_)
  n_col = desc_a%matrix_data(psb_n_col_)

  ! ensure global coherence for convergence checks.
  call blacs_get(icontxt,16,isvch)
  ich = 1 
  call blacs_set(icontxt,16,ich)


  if (present(istop)) then 
    listop = istop 
  else
    listop = 1
  endif
  !
  !  listop = 1:  normwise backward error, infinity norm 
  !  listop = 2:  ||r||/||b||   norm 2 
  !
!!$
!!$  if ((prec%prec < min_prec_).or.(prec%prec > max_prec_) ) then
!!$    write(0,*) 'f90_bicg: invalid iprec',prec%prec
!!$    if (present(ierr)) ierr=-1
!!$    return
!!$  endif

  if ((listop < 1 ).or.(listop > 2 ) ) then
    write(0,*) 'psb_bicg: invalid istop',listop 
    info=5001
    int_err=listop
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  naux=4*n_col 

  allocate(aux(naux),stat=info)
  call psb_dalloc(mglob,9,wwrk,desc_a,info)
  call psb_asb(wwrk,desc_a,info)  
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
    itrac = itrace
  else
    itrac = -1
  end if

  diagl  = 'u'
  diagu  = 'u'
  itx   = 0
  
  if (listop == 1) then 
     ani = psb_nrmi(a,desc_a,info)
     bni = psb_amax(b,desc_a,info)
  else if (listop == 2) then 
     bn2 = psb_nrm2(b,desc_a,info)
  endif
 
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
    call psb_axpby(one,b,zero,r,desc_a,info)
    call psb_spmm(-one,a,x,one,r,desc_a,info,work=aux)
    call psb_axpby(one,r,zero,rt,desc_a,info)
    if(info.ne.0) then
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if

    rho = zero
    if (debug) write(*,*) 'on entry to amax: b: ',size(b)
    if (listop == 1) then 
      rni = psb_amax(r,desc_a,info)
      xni = psb_amax(x,desc_a,info)
    else if (listop == 2) then 
      rni = psb_nrm2(r,desc_a,info)
    endif
    if(info.ne.0) then
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if

    if (listop == 1) then 
      xni  = psb_amax(x,desc_a,info)
      rerr =  rni/(ani*xni+bni)
      if (itrac /= -1) then 
        if (me.eq.0) write(itrac,'(a,i4,5(2x,es10.4))') 'bicg: ',itx,rerr,rni,bni,&
             &xni,ani
      endif
    else  if (listop == 2) then 
      rerr = rni/bn2
      if (itrac /= -1) then 
        if (me.eq.0) write(itrac,'(a,i4,3(2x,es10.4))') 'bicg: ',itx,rerr,rni,bn2
      endif
    endif

    if(info.ne.0) then
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if
    
    if (rerr<=eps) then 
      exit restart
    end if

    iteration:  do 
      it   = it + 1
      itx = itx + 1
      if (debug) write(*,*) 'iteration: ',itx

      call psb_prc_aply(prec,r,z,desc_a,info,work=aux)
      call psb_prc_aply(prec,rt,zt,desc_a,info,trans='t',work=aux)

      rho_old = rho    
      rho = psb_dot(rt,z,desc_a,info)
      if (rho==zero) then
        if (debug) write(0,*) 'bicg itxation breakdown r',rho
        exit iteration
      endif

      if (it==1) then
        call psb_axpby(one,z,zero,p,desc_a,info)
        call psb_axpby(one,zt,zero,pt,desc_a,info)
      else
        beta = (rho/rho_old)
        call psb_axpby(one,z,beta,p,desc_a,info)
        call psb_axpby(one,zt,beta,pt,desc_a,info)
      end if

      call psb_spmm(one,a,p,zero,q,desc_a,info,&
           & work=aux)
      call psb_spmm(one,a,pt,zero,qt,desc_a,info,&
           & work=aux,trans='t')

      sigma = psb_dot(pt,q,desc_a,info)
      if (sigma==zero) then
        if (debug) write(0,*) 'cgs iteration breakdown s1', sigma
        exit iteration
      endif

      alpha = rho/sigma


      call psb_axpby(alpha,p,one,x,desc_a,info)
      call psb_axpby(-alpha,q,one,r,desc_a,info)
      call psb_axpby(-alpha,qt,one,rt,desc_a,info)


      if (listop == 1) then 
        rni = psb_amax(r,desc_a,info)
        xni = psb_amax(x,desc_a,info)
      else if (listop == 2) then 
        rni = psb_nrm2(r,desc_a,info)
      endif

      if (listop == 1) then 
        xni  = psb_amax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)
        if (itrac /= -1) then 
          if (me.eq.0) write(itrac,'(a,i4,5(2x,es10.4))') 'bicg: ',itx,rerr,rni,bni,&
               &xni,ani
        endif
      else  if (listop == 2) then 
        rerr = rni/bn2
        if (itrac /= -1) then 
          if (me.eq.0) write(itrac,'(a,i4,3(2x,es10.4))') 'bicg: ',itx,rerr,rni,bn2
        endif
      endif
      if (rerr<=eps) then 
        exit restart
      end if
      if (itx.ge.litmax) exit restart
    end do iteration
  end do restart

  if (present(err)) err=rerr
  if (present(iter)) iter = itx
  if (rerr>eps) then
    write(0,*) 'bicg failed to converge to ',eps,&
         & ' in ',itx,' iterations  '
  end if


  deallocate(aux)
  call psb_free(wwrk,desc_a,info)
  ! restore external global coherence behaviour
  call blacs_set(icontxt,16,isvch)

  if(info/=0) then
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

end subroutine psb_dbicg


