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
! File:  psb_dcgs.f90
!
! Subroutine: psb_dcgs
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
Subroutine psb_dcgs(a,prec,b,x,eps,desc_a,info,&
     &itmax,iter,err,itrace,istop)
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
  Type(psb_dspmat_type), Intent(in)  :: a
  Type(psb_desc_type), Intent(in)    :: desc_a 
  Type(psb_dprec_type), Intent(in)   :: prec 
  Real(Kind(1.d0)), Intent(in)       :: b(:)
  Real(Kind(1.d0)), Intent(inout)    :: x(:)
  Real(Kind(1.d0)), Intent(in)       :: eps
  integer, intent(out)               :: info
  Integer, Optional, Intent(in)      :: itmax, itrace,istop
  Integer, Optional, Intent(out)     :: iter
  Real(Kind(1.d0)), Optional, Intent(out) :: err
!!$   local data
  Real(Kind(1.d0)), Pointer  :: aux(:),wwrk(:,:)
  Real(Kind(1.d0)), Pointer  :: ww(:), q(:),&
       & r(:), p(:), v(:), s(:), t(:), z(:), f(:), rt(:),qt(:),uv(:)
  Integer, Pointer           :: iperm(:), ipnull(:), ipsave(:)
  Real(Kind(1.d0)) ::rerr
  Integer       ::litmax, liter, naux, m, mglob, it, itrac,int_err(5),&
       & nprows,npcols,me,mecol, n_row, n_col,listop, err_act
  Character     ::diagl, diagu
  Logical, Parameter :: exchange=.True., noexchange=.False.  
  Integer, Parameter :: ione=1
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ich, icontxt
  Logical            :: do_renum_left
  Logical, Parameter :: debug = .false.
  Real(Kind(1.d0)), Parameter :: one=1.d0, zero=0.d0, epstol=1.d-35
  Real(Kind(1.d0)) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,bn2,& 
       & sigma, omega, tau 
  character(len=20)             :: name,ch_err

  info = 0
  name = 'psb_dcgs'
  call psb_erractionsave(err_act)

  If (debug) Write(*,*) 'entering psb_dcgs'
  icontxt = desc_a%matrix_data(psb_ctxt_)
  Call blacs_gridinfo(icontxt,nprows,npcols,me,mecol)
  If (debug) Write(*,*) 'psb_dcgs: from gridinfo',nprows,npcols,me

  mglob = desc_a%matrix_data(psb_m_)
  n_row = desc_a%matrix_data(psb_n_row_)
  n_col  = desc_a%matrix_data(psb_n_col_)

  If (Present(istop)) Then 
    listop = istop 
  Else
    listop = 1
  Endif
!
!  listop = 1:  normwise backward error, infinity norm 
!  listop = 2:  ||r||/||b||   norm 2 
!
!!$
!!$  If ((prec%prec < 0).Or.(prec%prec > 6) ) Then
!!$     Write(0,*) 'f90_cgstab: invalid iprec',prec%prec
!!$     If (Present(ierr)) ierr=-1
!!$     Return
!!$  Endif
  
  if ((listop < 1 ).or.(listop > 2 ) ) then
    write(0,*) 'psb_cgs: invalid istop',listop 
    info=5001
    int_err=listop
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  naux=4*n_col 
  Allocate(aux(naux),stat=info)

  Call psb_alloc(mglob,11,wwrk,desc_a,info)
  Call psb_asb(wwrk,desc_a,info)  
  if (info.ne.0) Then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  End If

  q  => wwrk(:,1)
  qt => wwrk(:,2)
  r  => wwrk(:,3)
  rt => wwrk(:,4)
  p  => wwrk(:,5)
  v  => wwrk(:,6)
  uv => wwrk(:,7)
  z  => wwrk(:,8)
  f  => wwrk(:,9)
  s  => wwrk(:,10)
  ww => wwrk(:,11)


  If (Present(itmax)) Then 
    litmax = itmax
  Else
    litmax = 1000
  Endif

  If (Present(itrace)) Then
     itrac = itrace
  Else
     itrac = -1
  End If

  ! ensure global coherence for convergence checks.
  Call blacs_get(icontxt,16,isvch)
  ich = 1 
  Call blacs_set(icontxt,16,ich)
  
  diagl  = 'u'
  diagu  = 'u'
  itx   = 0

  if (listop == 1) then 
    ani = psb_nrmi(a,desc_a,info)
    bni = psb_amax(b,desc_a,info)
  else if (listop == 2) then 
    bn2 = psb_nrm2(b,desc_a,info)
  endif
  if(info/=0)then
     info=4011
     call psb_errpush(info,name)
     goto 9999
  end if

  restart: Do 
!!$
!!$   r0 = b-ax0
!!$ 
    If (itx.Ge.itmax) Exit restart  
    it = 0      
    Call psb_axpby(one,b,zero,r,desc_a,info)
    Call psb_spmm(-one,a,x,one,r,desc_a,info,work=aux)
    Call psb_axpby(one,r,zero,rt,desc_a,info)
    if(info/=0)then
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if
    
    rho = zero
    If (debug) Write(*,*) 'on entry to amax: b: ',Size(b)

    if (listop == 1) then 
      rni = psb_amax(r,desc_a,info)
      xni = psb_amax(x,desc_a,info)
      rerr =  rni/(ani*xni+bni)
      if (itrac /= -1) then 
        If (me == 0) Write(itrac,'(a,i4,5(2x,es10.4))') 'cgs: ',&
             & itx,rerr,rni,bni,xni,ani
      endif
    else if (listop == 2) then 
      rni = psb_nrm2(r,desc_a,info)
      rerr = rni/bn2
      if (itrac /= -1) then 
        If (me == 0) Write(itrac,'(a,i4,3(2x,es10.4))') 'cgs: ',itx,rerr,rni,bn2
      endif
    endif
    if(info/=0)then
       info=4011
       call psb_errpush(info,name)
       goto 9999
    end if
    
    If (rerr<=eps) Then 
      Exit restart
    End If

    iteration:  Do 
      it   = it + 1
      itx = itx + 1
      If (debug) Write(*,*) 'iteration: ',itx
      rho_old = rho    
      rho = psb_dot(rt,r,desc_a,info)
      If (rho==zero) Then
         If (debug) Write(0,*) 'cgs iteration breakdown r',rho
        Exit iteration
      Endif

      If (it==1) Then
        Call psb_axpby(one,r,zero,uv,desc_a,info)
        Call psb_axpby(one,r,zero,p,desc_a,info)
      Else
        beta = (rho/rho_old)
        Call psb_axpby(one,r,zero,uv,desc_a,info)
        Call psb_axpby(beta,q,one,uv,desc_a,info)
        Call psb_axpby(one,q,beta,p,desc_a,info)
        Call psb_axpby(one,uv,beta,p,desc_a,info)

      End If

      Call psb_prcaply(prec,p,f,desc_a,info,work=aux)

      Call psb_spmm(one,a,f,zero,v,desc_a,info,&
           & work=aux)

      sigma = psb_dot(rt,v,desc_a,info)
      If (sigma==zero) Then
         If (debug) Write(0,*) 'cgs iteration breakdown s1', sigma
         Exit iteration
      Endif
      
      alpha = rho/sigma

      Call psb_axpby(one,uv,zero,q,desc_a,info)
      Call psb_axpby(-alpha,v,one,q,desc_a,info)
      Call psb_axpby(one,uv,zero,s,desc_a,info)
      Call psb_axpby(one,q,one,s,desc_a,info)
      
      Call psb_prcaply(prec,s,z,desc_a,info,work=aux)

      Call psb_axpby(alpha,z,one,x,desc_a,info)

      Call psb_spmm(one,a,z,zero,qt,desc_a,info,&
           & work=aux)
      
      Call psb_axpby(-alpha,qt,one,r,desc_a,info)
      
     
      if (listop == 1) then 
        rni = psb_amax(r,desc_a,info)
        xni = psb_amax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)
        if (itrac /= -1) then 
        If (me == 0) Write(itrac,'(a,i4,5(2x,es10.4))') 'cgs: ',&
             & itx,rerr,rni,bni,xni,ani
        endif

      else  if (listop == 2) then 

        rni = psb_nrm2(r,desc_a,info)
        rerr = rni/bn2
        if (itrac /= -1) then 
        If (me == 0) Write(itrac,'(a,i4,3(2x,es10.4))') 'cgs: ',&
             & itx,rerr,rni,bn2
        endif
      endif

      If (rerr<=eps) Then 
        Exit restart
      End If
      If (itx.Ge.itmax) Exit restart
    End Do iteration
  End Do restart

  If (Present(err)) err=rerr
  If (Present(iter)) iter = itx
  If (rerr>eps) Then
    Write(0,*) 'cgs failed to converge to ',eps,&
         & ' in ',itx,' iterations  '
  End If

  Deallocate(aux)
  Call psb_free(wwrk,desc_a,info)
  ! restore external global coherence behaviour
  Call blacs_set(icontxt,16,isvch)

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

End Subroutine psb_dcgs


