! File: psb_dnrm2.f90
!
! Function: psb_dnrm2
!    Forms the norm2 of a distributed vector,
!
!    norm2 := sqrt ( sub( X )**T * sub( X ) )
!
!    where sub( X ) denotes X(:,JX).
!
! Parameters:
!    x      -  real,dimension(:,:).       The input vector containing the entries of sub( X ).
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset for sub( X ).
!
function psb_dnrm2(x, desc_a, info, jx)  
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      ::  x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(in), optional     :: jx
  integer, intent(out)              :: info
  real(kind(1.D0))                  :: psb_dnrm2

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, me, mycol,&
       & err_act, n, iix, jjx, temp(2), ndim
  real(kind(1.d0))         :: nrm2
  real(kind(1.d0)),pointer :: tmpx(:)
  external dcombnrm2
  character(len=20)        :: name, ch_err

  name='psb_dnrm2'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (nprow == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= 1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif
  
  ix = 1
  if (present(jx)) then
     ijx = jx
  else
     ijx = 1
  endif

  m = desc_data(m_)

  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_data%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  if(m.ne.0) then
     if (desc_a%matrix_data(psb_n_row_) .gt. 0) then 
        ndim = desc_a%matrix_data(psb_n_row_)
        nrm2 = dnrm2( ndim, x(iix,jjx), ione )
        i=1
        do while (desc_a%ovrlap_elem(i) .ne. -1)
           id = desc_a%ovrlap_elem(i+n_dom_ovr_)
           dd = dble(id-1)/dble(id)
           nrm2 = nrm2 * sqrt(&
                &  one - dd * ( &
                &  x(desc_a%ovrlap_elem(i+ovrlp_elem_), jjx) &
                &  / nrm2 &
                &  ) ** 2 &
                &  ) 
           i = i+2
        end do
     else 	    
        nrm2 = zero
     end if
  else 	    
     nrm2 = zero
  end if
  
  call pdtreecomb(icontxt,'All',1,nrm2,-1,-1,dcombnrm2)
  
  psb_dnrm2 = nrm2  
  
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_dnrm2



! Function: psb_dnrm2
!    Forms the norm2 of a distributed vector,
!
!    norm2 := sqrt ( X**T * X)
!
! Parameters:
!    x      -  real,dimension(:).         The input vector containing the entries of X.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
function psb_dnrm2v(x, desc_a, info)  
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  real(kind(1.D0))                  :: psb_dnrm2v

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, me, mycol,&
       & err_act, n, iix, jjx, temp(2), ndim
  real(kind(1.d0))         :: nrm2
  real(kind(1.d0)),pointer :: tmpx(:)
  external dcombnrm2
  character(len=20)        :: name, ch_err

  name='psb_dnrm2v'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (nprow == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= 1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif
  
  ix = 1
  jx=1

  m = desc_data(m_)

  call psb_chkvect(m,1,size(x),ix,jx,desc_data%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  if(m.ne.0) then
     if (desc_a%matrix_data(psb_n_row_) .gt. 0) then 
        ndim = desc_a%matrix_data(psb_n_row_)
        nrm2 = dnrm2( ndim, x, ione )
        i=1
        do while (desc_a%ovrlap_elem(i) .ne. -1)
           id = desc_a%ovrlap_elem(i+n_dom_ovr_)
           dd = dble(id-1)/dble(id)
           nrm2 = nrm2 * sqrt(&
                &  one - dd * ( &
                &  x(desc_a%ovrlap_elem(i+ovrlp_elem_)) &
                &  / nrm2 &
                &  ) ** 2 &
                &  ) 
           i = i+2
        end do
     else 	    
        nrm2 = zero
     end if
  else 	    
     nrm2 = zero
  end if
  
  call pdtreecomb(icontxt,'All',1,nrm2,-1,-1,dcombnrm2)
  
  psb_dnrm2v = nrm2  
  
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_dnrm2v




! Subroutine: psb_dnrm2
!    Forms the norm2 of a distributed vector,
!
!    norm2 := sqrt ( X**T * X)
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:).         The input vector containing the entries of X.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
subroutine psb_dnrm2vs(res, x, desc_a, info)
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:)
  real(kind(1.d0)), intent(out)     :: res
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, me, mycol,&
       & err_act, n, iix, jjx, temp(2), ndim
  real(kind(1.d0))         :: nrm2
  real(kind(1.d0)),pointer :: tmpx(:)
  external dcombnrm2
  character(len=20)        :: name, ch_err

  name='psb_dnrm2'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (nprow == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= 1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif
  
  ix = 1
  jx = 1
  m = desc_data(m_)

  call psb_chkvect(m,1,size(x),ix,jx,desc_data%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  if(m.ne.0) then
     if (desc_a%matrix_data(psb_n_row_) .gt. 0) then 
        ndim = desc_a%matrix_data(psb_n_row_)
        nrm2 = dnrm2( ndim, x, ione )
        i=1
        do while (desc_a%ovrlap_elem(i) .ne. -1)
           id = desc_a%ovrlap_elem(i+n_dom_ovr_)
           dd = dble(id-1)/dble(id)
           nrm2 = nrm2 * sqrt(&
                &  one - dd * ( &
                &  x(desc_a%ovrlap_elem(i+ovrlp_elem_)) &
                &  / nrm2 &
                &  ) ** 2 &
                &  ) 
           i = i+2
        end do
     else 	    
        nrm2 = zero
     end if
  else 	    
     nrm2 = zero
  end if
  
  call pdtreecomb(icontxt,'All',1,nrm2,-1,-1,dcombnrm2)
  
  res = nrm2  
  
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_dnrm2vs
