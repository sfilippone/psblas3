! File: psb_daxpby.f90
!
! Subroutine: psb_daxpby
!    Adds one distributed matrix to another,
!
!    sub( Y ) := beta * sub( Y ) + alpha * sub( X )
!
!    where sub( X ) denotes X(:,JX)
!
!    sub( Y ) denotes Y(:,JY).
!
! Parameters:
!    alpha  -  real.                      The scalar used to multiply each component of sub( X ).
!    x      -  real,dimension(:,:).       The input vector containing the entries of sub( X ).
!    beta   -  real.                      The scalar used to multiply each component of sub( Y ).
!    y      -  real,dimension(:,:).       The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset for sub( X ).
!    jy     -  integer(optional).         The column offset for sub( Y ).
!
subroutine  psb_daxpby(alpha, x, beta,y,desc_a,info, n, jx, jy)
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none                    

  integer, intent(in), optional   :: n, jx, jy
  integer, intent(out)            :: info
  type(psb_desc_type), intent(in) :: desc_a
  real(kind(1.D0)), intent(in)    :: alpha, beta
  real(kind(1.D0)), intent(in)    :: x(:,:)
  real(kind(1.D0)), intent(inout) :: y(:,:)

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, iix, jjx, temp(2), ix, iy, ijx, ijy, m, iiy, in, jjy
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_daxpby'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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

  if (present(n)) then
     if(((ijx+n).le.size(x,2)).and.&
          & ((ijy+n).le.size(y,2))) then 
        in = n
     else
        in = min(size(x,2),size(y,2))
     end if
  else
     in = min(size(x,2),size(y,2))
  endif

  if(ijx.ne.ijy) then
     info=3050
     call psb_errpush(info,name)
     goto 9999
  end if

  m = desc_a%matrix_data(psb_m_)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
  call psb_chkvect(m,ione,size(y,1),iy,ijy,desc_a%matrix_data,info,iiy,jjy)
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
  
  if ((in.ne.0)) then
     if(desc_a%matrix_data(psb_n_row_).gt.0) then
        call daxpby(desc_a%matrix_data(psb_n_col_),in,&
             & alpha,x(iix,jjx),size(x,1),beta,&
             & y(iiy,jjy),size(y,1),info)
     end if
  end if
        
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_daxpby





!
! Subroutine: psb_daxpbyv
!    Adds one distributed matrix to another,
!
!    Y := beta * Y + alpha * X
!
! Parameters:
!    alpha  -  real.                      The scalar used to multiply each component of X.
!    x      -  real,dimension(:).         The input vector containing the entries of X.
!    beta   -  real.                      The scalar used to multiply each component of Y.
!    y      -  real,dimension(:).         The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
subroutine  psb_daxpbyv(alpha, x, beta,y,desc_a,info)
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none                    

  integer, intent(out)            :: info
  type(psb_desc_type), intent(in) :: desc_a
  real(kind(1.D0)), intent(in)    :: alpha, beta
  real(kind(1.D0)), intent(in)    :: x(:)
  real(kind(1.D0)), intent(inout) :: y(:)

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, temp(2), ix, iy, ijx, m, iiy, in, jjy
  character(len=20)        :: name, ch_err

  name='psb_daxpby'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
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

  m = desc_a%matrix_data(psb_m_)

  ! check vector correctness
  call psb_chkvect(m,ione,size(x),ix,ione,desc_a%matrix_data,info,iix,jjx)
  call psb_chkvect(m,ione,size(y),iy,ione,desc_a%matrix_data,info,iiy,jjy)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
  end if

  if ((iix.ne.ione).or.(iiy.ne.ione)) then
     info=3040
     call psb_errpush(info,name)
  end if
  
  write(0,'(i2," before daxpby",2(i6,2x),2(f10.2,2x))')myrow,desc_a%matrix_data(psb_n_row_),&
       & desc_a%matrix_data(psb_n_col_),alpha,beta
  if(desc_a%matrix_data(psb_n_row_).gt.0) then
     call daxpby(desc_a%matrix_data(psb_n_col_),ione,&
          & alpha,x,size(x),beta,&
          & y,size(y),info)
  end if
  
  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_daxpbyv
