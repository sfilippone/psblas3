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
! File:  psb_dgmresr.f90
!
! Subroutine: psb_dgmres
!    This subroutine implements the restarted GMRES method.
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
!    irst    -  integer(optional).           The restart value.
!    istop   -  integer(optional).           The stopping criterium.
!
Subroutine psb_dgmresr(a,prec,b,x,eps,desc_a,info,&
     &itmax,iter,err,itrace,irst,istop)
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
  Integer, Optional, Intent(in)      :: itmax, itrace, irst,istop
  Integer, Optional, Intent(out)     :: iter
  Real(Kind(1.d0)), Optional, Intent(out) :: err
!!$   local data
  Real(Kind(1.d0)), Pointer  :: aux(:),wwrk(:,:)
  Real(Kind(1.d0)), Pointer  :: w(:), q(:), r(:), rt0(:), p(:), v(:,:), &
       & c(:),s(:), t(:), z(:), f(:), uh(:,:), h(:,:), rs(:),&
       & gamma(:), gamma1(:), gamma2(:), taum(:,:), sigma(:),&
       &pv1(:),  pv2(:), pm1(:,:), rr(:,:)
  Integer, Pointer           :: iperm(:), ipnull(:), ipsave(:), ierrv(:)
  Real(Kind(1.d0)) :: rerr, scal, gm 
  Integer       ::litmax, liter, naux, m, mglob, it,k, itrac,&
       & nprows,npcols,me,mecol, n_row, n_col, nl, int_err(5)
  Character     ::diagl, diagu
  Logical, Parameter :: exchange=.True., noexchange=.False.  
  Integer, Parameter :: irmax = 8
  Integer            :: itx, i, isvch, ich, icontxt,listop, err_act
  Logical            :: do_renum_left,inner_stop
  Real(Kind(1.d0)), Parameter :: one=1.d0, zero=0.d0, epstol=1.d-35
  Logical, Parameter :: debug = .false.
  Real(Kind(1.d0)) :: alpha, beta, rho, rho_old, rni, xni, bni, ani,bn2,& 
       & omega, tau 
  real(kind(1.d0)), external :: dnrm2
  character(len=20)             :: name,ch_err

  info = 0
  name = 'psb_dgmres'
  call psb_erractionsave(err_act)

  If (debug) Write(0,*) 'entering psb_dgmres'
  icontxt = desc_a%matrix_data(psb_ctxt_)
  Call blacs_gridinfo(icontxt,nprows,npcols,me,mecol)

  If (debug) Write(0,*) 'psb_dgmres: from gridinfo',nprows,npcols,me

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

  if ((listop < 1 ).or.(listop > 2 ) ) then
    write(0,*) 'psb_dgmres: invalid istop',listop 
    info=5001
    int_err(1)=listop
    err=info
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

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
  
  If (Present(irst)) Then
    nl = irst
    If (debug) Write(0,*) 'present: irst: ',irst,nl
  Else
    nl = 10 
    If (debug) Write(0,*) 'not present: irst: ',irst,nl
  Endif


  naux=4*n_col 
  Allocate(aux(naux),h(nl+1,nl+1),rr(nl+1,nl+1),&
       &c(nl+1),s(nl+1),rs(nl+1), stat=info)

  If (info.Ne.0) Then 
     info = 4000
     call psb_errpush(info,name)
     goto 9999
  End If

  Call psb_dsall(mglob,nl+1,v,desc_a,info)
  Call psb_dsall(mglob,w,desc_a,info)
  Call psb_dsasb(v,desc_a,info)  
  Call psb_dsasb(w,desc_a,info)  
  if (info.ne.0) Then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  End If

  ! ensure global coherence for convergence checks.
  Call blacs_get(icontxt,16,isvch)
  ich = 1 
  Call blacs_set(icontxt,16,ich)

  if (listop == 1) then 
    ani = psb_nrmi(a,desc_a,info)
    bni = psb_amax(b,desc_a,info)
  else if (listop == 2) then 
    bn2 = psb_nrm2(b,desc_a,info)
  endif
  if (info.ne.0) Then 
     info=4011 
     call psb_errpush(info,name)
     goto 9999
  End If

  diagl  = 'u'
  diagu  = 'u'
  itx   = 0
  restart: Do 
!!$   
!!$   r0 = b-ax0
!!$ 
    If (debug) Write(0,*) 'restart: ',itx,it
    it = 0      
    Call psb_axpby(one,b,zero,v(:,1),desc_a,info)
    Call psb_spmm(-one,a,x,one,v(:,1),desc_a,info,work=aux)
    
    call psb_prcaply(prec,v(:,1),desc_a,info)
    rs(1) = psb_nrm2(v(:,1),desc_a,info)
    if (info.ne.0) Then 
       info=4011 
       call psb_errpush(info,name)
       goto 9999
    End If

    scal=one/rs(1)
    If (debug) Write(0,*) 'on entry to amax: b: ',Size(b),rs(1),scal

    if (listop == 1) then 
      rni = psb_amax(v(:,1),desc_a,info)
      xni = psb_amax(x,desc_a,info)
      rerr =  rni/(ani*xni+bni)
      if (itrac /= -1) then 
          If (me == 0) Write(itrac,'(a,i4,5(2x,es10.4))') 'gmresr(l): ',&
               & itx,rerr,rni,bni,xni,ani
      endif
    else if (listop == 2) then 
      rni = psb_nrm2(v(:,1),desc_a,info)
      rerr = rni/bn2
      if (itrac /= -1) then  
        If (me == 0) Write(itrac,'(a,i4,3(2x,es10.4))') 'gmresr(l): ',&
             & itx,rerr,rni,bn2
      endif
    endif
    if (info.ne.0) Then 
       info=4011 
       call psb_errpush(info,name)
       goto 9999
    End If
    
    If (rerr<=eps) Then 
      Exit restart
    End If
     
    If (itx.Ge.litmax) Exit restart  

    v(:,1) = v(:,1) * scal

    inner:  Do i=1,nl
      itx  = itx + 1

      Call psb_spmm(one,a,v(:,i),zero,w,desc_a,info,work=aux)
      call psb_prcaply(prec,w,desc_a,info)

      do k = 1, i
        h(k,i) = psb_dot(v(:,k),w,desc_a,info)
        call psb_axpby(-h(k,i),v(:,k),one,w,desc_a,info)
      end do
      h(i+1,i) = psb_nrm2(w,desc_a,info)
      scal=one/h(i+1,i)
      call psb_axpby(scal,w,zero,v(:,i+1),desc_a,info)
      do k=2,i
        rr(k-1,i) =  c(k-1)*h(k-1,i) + s(k-1)*h(k,i)
        rr(k,i)   = -s(k-1)*h(k-1,i) + c(k-1)*h(k,i)
      enddo
      gm =  safe_dn2(h(i,i),h(i+1,i))
      if (debug) write(0,*) 'GM : ',gm
      gm = max(gm,epstol)
      
      c(i) = h(i,i)/gm
      s(i) = h(i+1,i)/gm
      rs(i+1) = -s(i)*rs(i)
      rs(i)   = c(i)*rs(i)
      rr(i,i)  = c(i)*h(i,i)+s(i)*h(i+1,i)
            
      if (listop == 1) then 
        rni = abs(rs(i+1))
        xni = psb_amax(x,desc_a,info)
        rerr =  rni/(ani*xni+bni)
        if (itrac /= -1) then 
          If (me == 0) Write(itrac,'(a,i4,5(2x,es10.4))') 'gmresr(l): ',&
               & itx,rerr,rni,bni,xni,ani
        endif
      else if (listop == 2) then 
        rni = abs(rs(i+1))
        rerr = rni/bn2
        if (itrac /= -1) then  
          If (me == 0) Write(itrac,'(a,i4,3(2x,es10.4))') 'gmresr(l): ',&
               & itx,rerr,rni,bn2
        endif
      endif

      if (rerr < eps ) then 
        call dtrsm('l','u','n','n',i,1,one,rr,size(rr,1),rs,nl)
        if (debug) write(0,*) 'Rebuild x-> RS:',rs(21:nl)
        do k=1, i
          call psb_axpby(rs(k),v(:,k),one,x,desc_a,info)
        end do
        exit restart
      end if

    end Do inner
    if (debug) write(0,*) 'Before DTRSM :',rs(1:nl)
    call dtrsm('l','u','n','n',nl,1,one,rr,size(rr,1),rs,nl)
    if (debug) write(0,*) 'Rebuild x-> RS:',rs(21:nl)
    do k=1, nl
      call psb_axpby(rs(k),v(:,k),one,x,desc_a,info)
    end do
     
  End Do restart

  If (Present(err)) err=rerr
  If (Present(iter)) iter = itx
  If (rerr>eps) Then
    Write(0,*) 'gmresr(l) failed to converge to ',eps,&
         & ' in ',itx,' iterations  '
  End If


  Deallocate(aux,h,c,s,rs,rr, stat=info)
  Call psb_free(v,desc_a,info)
  Call psb_free(w,desc_a,info)
  ! restore external global coherence behaviour
  Call blacs_set(icontxt,16,isvch)

  if (info /= 0) then
     info=4011
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


contains
  function safe_dn2(a,b)
    real(kind(1.d0)), intent(in) :: a, b
    real(kind(1.d0))  :: safe_dn2
    real(kind(1.d0))  :: t
    
    t = max(abs(a),abs(b))
    if (t==0.d0) then 
      safe_dn2 = 0.d0
    else
      safe_dn2 = t * sqrt(abs(a/t)**2 + abs(b/t)**2)
    endif
    return
  end function safe_dn2
    

End Subroutine psb_dgmresr


