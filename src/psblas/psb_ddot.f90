! File: psb_ddot.f90
!
! Function: psb_ddot
!    psb_ddot forms the dot product of two distributed vectors,
!
!    dot := sub( X )**T * sub( Y )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Parameters:
!    x      -  real,dimension(:,:).       The input vector containing the entries of sub( X ).
!    y      -  real,dimension(:,:).       The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset for sub( X ).
!    jy     -  integer(optional).         The column offset for sub( Y ).
!
function psb_ddot(x, y,desc_a, info, jx, jy)  
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:,:), y(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(in), optional    :: jx, jy
  integer, intent(out)             :: info
  real(kind(1.D0))                 :: f90_psddot

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, me, mycol,&
       & err_act, n, iix, jjx, temp(2)
  real(kind(1.d0)),pointer :: tmpx(:)
  real(kind(1.D0))         :: dot_local
  character(len=20)        :: name, ch_err

  name='psb_ddot'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (nprow == -ione) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= ione) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = ione
  if (present(jx)) then
     ijx = jx
  else
     ijx = ione
  endif

  iy = ione
  if (present(jy)) then
     ijy = jy
  else
     ijy = ione
  endif

  if(ijx.ne.ijy) then
     info=3050
     call psb_errpush(info,name)
     goto 9999
  end if

  m = desc_a%matrix_data(m_)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ijx,desc_data%matrix_data,info,iix,jjx)
  call psb_chkvect(m,ione,size(y,1),iy,ijy,desc_data%matrix_data,info,iiy,jjy)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((iix.ne.ione).or.(iiy.ne.ione)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  if(m.ne.0) then
     if(desc_a%matrix_data(psb_n_row_).gt.0) then
        dot = ddot(desc_a%matrix_data(psb_n_row_),&
             & x(iix,jjx),ione,y(iiy,jjy),ione)
        ! adjust dot because overlapped elements are computed more than once
        i=1
        do while (desc_a%ovrlap_elem(i).ne.-ione)
           dot = dot -&
                & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
                & x(iix+desc_a%ovrlap_elem(i)-1,jjx)*
                & y(iiy+desc_a%ovrlap_elem(i)-1,jjy)
           i = i+2
        end do
     else
        dot=0.d0
     end if
  else
     dot=0.d0
  end if

  ! compute global sum
  call dgsum2d(icontxt, 'A', ' ', ione, ione, dot,&
       & ione, mone ,mycol)
  
  psb_ddot = dot

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_ddot




! Function: psb_ddotv
!    psb_ddot forms the dot product of two distributed vectors,
!
!    dot := X**T * Y
!
! Parameters:
!    x      -  real,dimension(:).         The input vector containing the entries of X.
!    y      -  real,dimension(:).         The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
function psb_ddotv(x, y,desc_a, info)  
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:), y(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  real(kind(1.D0))                 :: psb_ddotv

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, me, mycol,&
       & err_act, n, iix, jjx, temp(2)
  real(kind(1.d0)),pointer :: tmpx(:)
  real(kind(1.D0))         :: dot_local
  character(len=20)        :: name, ch_err

  name='psb_ddot'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (nprow == -ione) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= ione) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = ione
  iy = ione
  m = desc_a%matrix_data(m_)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ijx,desc_data%matrix_data,info,iix,jjx)
  call psb_chkvect(m,ione,size(y,1),iy,ijy,desc_data%matrix_data,info,iiy,jjy)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((iix.ne.ione).or.(iiy.ne.ione)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  if(m.ne.0) then
     if(desc_a%matrix_data(psb_n_row_).gt.0) then
        dot = ddot(desc_a%matrix_data(psb_n_row_),&
             & x,ione,y,ione)
        ! adjust dot because overlapped elements are computed more than once
        i=1
        do while (desc_a%ovrlap_elem(i).ne.-ione)
           dot = dot -&
                & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
                & x(desc_a%ovrlap_elem(i))*
                & y(desc_a%ovrlap_elem(i))
           i = i+2
        end do
     else
        dot=0.d0
     end if
  else
     dot=0.d0
  end if

  ! compute global sum
  call dgsum2d(icontxt, 'A', ' ', ione, ione, dot,&
       & ione, mone ,mycol)
  
  psb_ddotv = dot

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_ddotv



! Subroutine: psb_ddotvs
!    psb_ddot forms the dot product of two distributed vectors,
!
!    dot := X**T * Y
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:).         The input vector containing the entries of X.
!    y      -  real,dimension(:).         The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
subroutine psb_ddotvs(res, x, y,desc_a, info)  
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:), y(:)
  real(kind(1.d0)), intent(out)    :: res
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, me, mycol,&
       & err_act, n, iix, jjx, temp(2)
  real(kind(1.d0)),pointer :: tmpx(:)
  real(kind(1.D0))         :: dot_local
  character(len=20)        :: name, ch_err

  name='psb_ddot'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (nprow == -ione) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= ione) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = ione
  iy = ione
  m = desc_a%matrix_data(m_)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ijx,desc_data%matrix_data,info,iix,jjx)
  call psb_chkvect(m,ione,size(y,1),iy,ijy,desc_data%matrix_data,info,iiy,jjy)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((iix.ne.ione).or.(iiy.ne.ione)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  if(m.ne.0) then
     if(desc_a%matrix_data(psb_n_row_).gt.0) then
        dot = ddot(desc_a%matrix_data(psb_n_row_),&
             & x,ione,y,ione)
        ! adjust dot because overlapped elements are computed more than once
        i=1
        do while (desc_a%ovrlap_elem(i).ne.-ione)
           dot = dot -&
                & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
                & x(desc_a%ovrlap_elem(i))*
                & y(desc_a%ovrlap_elem(i))
           i = i+2
        end do
     else
        dot=0.d0
     end if
  else
     dot=0.d0
  end if

  ! compute global sum
  call dgsum2d(icontxt, 'A', ' ', ione, ione, dot,&
       & ione, mone ,mycol)
  
  res = dot

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_ddotvs




! Subroutine: psb_dmdots
!    psb_ddot forms the dot product of two distributed vectors,
!
!    dot := sub( X )**T * sub( Y )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:,:).       The input vector containing the entries of sub( X ).
!    y      -  real,dimension(:,:).       The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
subroutine psb_dmdots(res, x, y, desc_a, info)  
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)     :: x(:,:), y(:,:)
  real(kind(1.d0)), intent(out)    :: res(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, me, mycol,&
       & err_act, n, iix, jjx, temp(2)
  real(kind(1.d0)),pointer :: dot(:)
  real(kind(1.D0))         :: dot_local
  character(len=20)        :: name, ch_err

  name='psb_dmdots'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (nprow == -ione) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= ione) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  ix = ione
  iy = ione

  m = desc_a%matrix_data(m_)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ijx,desc_data%matrix_data,info,iix,jjx)
  call psb_chkvect(m,ione,size(y,1),iy,ijy,desc_data%matrix_data,info,iiy,jjy)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((iix.ne.ione).or.(iiy.ne.ione)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  k = min(size(x,2),size(y,2))
  allocate(dot(k))

  if(m.ne.0) then
     if(desc_a%matrix_data(psb_n_row_).gt.0) then
        do j=1,k
           dot(j) = ddot(desc_a%matrix_data(psb_n_row_),&
                & x(iix,jjx+j-1),ione,y(iiy,jjy+j-1),ione)
           ! adjust dot because overlapped elements are computed more than once
           i=1
           do while (desc_a%ovrlap_elem(i).ne.-ione)
              dot(j) = dot(j) -&
                   & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
                   & x(iix+desc_a%ovrlap_elem(i)-1,jjx+j-1)*
                   & y(iiy+desc_a%ovrlap_elem(i)-1,jjy+j-1)
              i = i+2
           end do
        end do
     else
        dot(:)=0.d0
     end if
  else
     dot(:)=0.d0
  end if

  ! compute global sum
  call dgsum2d(icontxt, 'A', ' ', ione, ione, dot,&
       & ione, mone ,mycol)
  
  res(1:k) = dot(1:k)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_dmdots
