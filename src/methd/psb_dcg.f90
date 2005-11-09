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
  use psb_serial_mod
  use psb_descriptor_type
  use psb_prec_type
  use psb_psblas_mod
  use psb_tools_mod
  use psb_const_mod
  use psb_prec_mod
  use psb_error_mod
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
  real(kind(1.d0)), pointer  :: aux(:), q(:), p(:),&
       & r(:), z(:), w(:), wwrk(:,:)
  real(kind(1.d0))    ::rerr
  real(kind(1.d0))    ::alpha, beta, rho, rho_old, rni, xni, bni, ani,bn2,& 
       & sigma
  integer         :: litmax, liter, listop, naux, m, mglob, it, itx, itrac,&
       & nprows,npcols,me,mecol, n_col, isvch, ich, icontxt, n_row,err_act, int_err(5)
  character          ::diagl, diagu
  logical, parameter :: exchange=.true., noexchange=.false.  
  real(kind(1.d0)), parameter :: one=1.d0, zero=0.d0, epstol=1.d-35
  character(len=20)             :: name,ch_err

  info = 0
  name = 'psb_dcg'
  call psb_erractionsave(err_act)


  icontxt = desc_a%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprows,npcols,me,mecol)

  mglob = desc_a%matrix_data(psb_m_)
  n_row = desc_a%matrix_data(psb_n_row_)
  n_col = desc_a%matrix_data(psb_n_col_)


  if (present(istop)) then 
    listop = istop 
  else
    listop = 1
  endif
  !
  !  LISTOP = 1:  Normwise backward error, infinity norm 
  !  LISTOP = 2:  ||r||/||b||   norm 2 
  !

!!$  If ((prec%prec < min_prec_).Or.(prec%prec > max_prec_) ) Then
!!$    Write(0,*) 'F90_CG: Invalid IPREC',prec%prec
!!$    If (Present(ierr)) ierr=-1
!!$    Return
!!$  Endif

  if ((listop < 1 ).or.(listop > 2 ) ) then
    write(0,*) 'psb_cg: invalid istop',listop 
    info=5001
    int_err(1)=listop
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  naux=4*n_col
  allocate(aux(naux), stat=info)
  call psb_dalloc(mglob,5,wwrk,desc_a,info)
  call psb_asb(wwrk,desc_a,info)  
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
    itrac = itrace
  else
    itrac = -1
  end if

!!$  DIAGL  = 'U'
!!$  DIAGU  = 'R'

  ! Ensure global coherence for convergence checks.
  call blacs_get(icontxt,16,isvch)
  ich = 1 
  call blacs_set(icontxt,16,ich)

  restart: do 
!!$   
!!$    r0 = b-Ax0
!!$   
    if (itx>= litmax) exit restart 
    it = 0
    call psb_axpby(one,b,zero,r,desc_a,info)
    call psb_spmm(-one,a,x,one,r,desc_a,info,work=aux)
    if (info.ne.0) then 
      info=4011
      call psb_errpush(info,name)
      goto 9999
    end if

    rho = zero
    if (listop == 1) then 
      ani = psb_nrmi(a,desc_a,info)
      bni = psb_amax(b,desc_a,info)
    else if (listop == 2) then 
      bn2 = psb_nrm2(b,desc_a,info)
    endif
    if (info.ne.0) then 
      info=4011
      call psb_errpush(info,name)
      goto 9999
    end if


    iteration:  do 
      it   = it + 1
      itx = itx + 1

      Call psb_prcaply(prec,r,z,desc_a,info,work=aux)
      rho_old = rho
      rho     = psb_dot(r,z,desc_a,info)

      if (it==1) then
        call psb_axpby(one,z,zero,p,desc_a,info)
      else
        if (rho_old==zero) then
          write(0,*) 'CG Iteration breakdown'
          exit iteration
        endif
        beta = rho/rho_old
        call psb_axpby(one,z,beta,p,desc_a,info)
      end if

      call psb_spmm(one,a,p,zero,q,desc_a,info,work=aux)
      sigma = psb_dot(p,q,desc_a,info)
      if (sigma==zero) then
        write(0,*) 'CG Iteration breakdown'
        exit iteration
      endif

      alpha = rho/sigma
      call psb_axpby(alpha,p,one,x,desc_a,info)
      call psb_axpby(-alpha,q,one,r,desc_a,info)


      if (listop == 1) Then 
        rni = psb_amax(r,desc_a,info)
        xni = psb_amax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)
        If (itrac /= -1) Then 
          If (me.Eq.0) Write(itrac,'(a,i4,5(2x,es10.4))') 'cg: ',itx,rerr,rni,bni,&
               &xni,ani
        Endif

      Else  If (listop == 2) Then 

        rni = psb_nrm2(r,desc_a,info)
        rerr = rni/bn2
        If (itrac /= -1) Then 
          If (me.Eq.0) Write(itrac,'(a,i4,3(2x,es10.4)))') 'cg: ',itx,rerr,rni,bn2
        Endif
      Endif
      if (rerr<=eps) exit restart
      if (itx>= litmax) exit restart 
    end do iteration
  end do restart

  if (present(err)) err=rerr
  if (present(iter)) iter = itx
  if (rerr>eps) then
    write(0,*) 'CG Failed to converge to ',eps,&
         & ' in ',litmax,' iterations '
    info=itx
  end if

  deallocate(aux)
  call psb_free(wwrk,desc_a,info)
  ! restore external global coherence behaviour
  call blacs_set(icontxt,16,isvch)

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


