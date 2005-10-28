! File: psb_damax.f90
!
! Function: psb_damax
!    Searches the absolute max of X.
!
!    normi := max(abs(sub(X)(i))  
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Parameters:
!    x      -  real,dimension(:,:).       The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset.
!
function psb_damax (x,desc_a, info, jx)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  integer, optional, intent(in)     :: jx
  real(kind(1.d0))                  :: psb_damax

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, temp(2), ix, ijx, m, i, k, imax, idamax
  real(kind(1.d0))         :: locmax(2), amax
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_damax'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  locmax(:)=0.d0
  amax=0.d0

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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

  m = desc_a%matrix_data(psb_m_)

  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  ! compute local max
  if ((desc_a%matrix_data(psb_n_row_).gt.0).and.(m.ne.0)) then
     imax=idamax(desc_a%matrix_data(psb_n_row_)-iix+1,x(iix,jjx),1)
     amax=abs(x(iix+imax-1,jjx))
  end if
  
  ! compute global max
  call dgamx2d(icontxt, 'A', ' ', ione, ione, amax, ione,&
           &temp ,temp,-ione ,-ione,-ione)

  psb_damax=amax

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_damax




! Function: psb_damaxv
!    Searches the absolute max of X.
!
!    normi := max(abs(X(i))  
!
! Parameters:
!    x      -  real,dimension(:).         The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
function psb_damaxv (x,desc_a, info)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  real(kind(1.d0))                  :: psb_damaxv

  ! locals
  integer                  :: int_err(5), err, icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, jx, temp(2), ix, ijx, m, imax, idamax
  real(kind(1.d0))         :: locmax(2), amax
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_damaxv'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  locmax(:)=0.d0
  amax=0.d0

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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

  m = desc_a%matrix_data(psb_m_)

  call psb_chkvect(m,1,size(x,1),ix,jx,desc_a%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  ! compute local max
  if ((desc_a%matrix_data(psb_n_row_).gt.0).and.(m.ne.0)) then
     imax=idamax(desc_a%matrix_data(psb_n_row_)-iix+1,x(iix),1)
     amax=abs(x(iix+imax-1))
  end if
  
  ! compute global max
  call dgamx2d(icontxt, 'A', ' ', ione, ione, amax, ione,&
           &temp ,temp,-ione ,-ione,-ione)

  psb_damaxv=amax

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_damaxv


! Subroutine: psb_damaxvs
!    Searches the absolute max of X.
!
!    normi := max(abs(sub(X)(i))  
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:,:).       The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset.
!
subroutine psb_damaxvs (res,x,desc_a, info)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  real(kind(1.D0)), intent(out)     :: res

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, temp(2), ix, ijx, m, imax, idamax
  real(kind(1.d0))         :: locmax(2), amax
  character(len=20)        :: name, ch_err

  name='psb_damaxvs'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  locmax(:)=0.d0
  amax=0.d0

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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
  ijx=1

  m = desc_a%matrix_data(psb_m_)

  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  ! compute local max
  if ((desc_a%matrix_data(psb_n_row_).gt.0).and.(m.ne.0)) then
     imax=idamax(desc_a%matrix_data(psb_n_row_)-iix+1,x(iix),1)
     amax=abs(x(iix+imax-1))
  end if
  
  ! compute global max
  call dgamx2d(icontxt, 'A', ' ', ione, ione, amax, ione,&
           &temp ,temp,-ione ,-ione,-ione)

  res = amax

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_damaxvs




! Subroutine: psb_dmamaxs
!    Searches the absolute max of X.
!
!    normi := max(abs(X(i))  
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:).         The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
subroutine psb_dmamaxs (res,x,desc_a, info,jx)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  integer, optional, intent(in)     :: jx
  real(kind(1.d0)), intent(out) :: res(:)

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, ix, temp(2), ijx, m, imax, i, k, idamax
  real(kind(1.d0))         :: locmax(2), amax
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_dmamaxs'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  locmax(:)=0.d0
  amax=0.d0

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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

  m = desc_a%matrix_data(psb_m_)
  k  = min(size(x,2),size(res,1))

  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  ! compute local max
  if ((desc_a%matrix_data(psb_n_row_).gt.0).and.(m.ne.0)) then
     do i=1,k
        imax=idamax(desc_a%matrix_data(psb_n_row_)-iix+1,x(iix,jjx),1)
        res(i)=abs(x(iix+imax-1,jjx+i-1))
     end do
  end if
  
  ! compute global max
  call dgamx2d(icontxt, 'A', ' ', ione, ione, amax, ione,&
           &temp ,temp,-ione ,-ione,-ione)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_dmamaxs
