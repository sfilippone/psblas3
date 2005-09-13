! File: psb_dasum.f90
!
! Function: psb_dasum 
!    Computes norm1 of X
!
!    norm1 := sum(sub( X )(i))
!
!    where sub( X ) denotes X(1:N,JX:).
!
! Parameters:
!    x      -  real,dimension(:,:).       The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset.
!
function psb_dasum (x,desc_a, info, jx)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  integer, optional, intent(in)     :: jx
  real(kind(1.d0))                  :: psb_dasum

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, temp(2), ix, ijx, m, i
  real(kind(1.d0))         :: asum, dasum
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_dasum'
  info=0
  call psb_erractionsave(err_act)

  asum=0.d0

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

  ! check vector correctness
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
  if ((m.ne.0)) then
     if(desc_a%matrix_data(psb_n_row_).gt.0) then
        asum=dasum(desc_a%matrix_data(psb_n_row_)-iix+1,x(iix,jjx),ione)

        ! adjust asum because overlapped elements are computed more than once
        i=1
        do while (desc_a%ovrlap_elem(i).ne.-ione)
           asum = asum -&
                & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
                & tmpx(desc_a%ovrlap_elem(i))
           i = i+2
        end do

        ! compute global sum
        call dgsum2d(icontxt, 'A', ' ', ione, ione, asum,&
             & ione, mone ,mycol)

     else
        asum=0.d0
        ! compute global sum
        call dgsum2d(icontxt, 'A', ' ', ione, ione, asum,&
             & ione, mone ,mycol)
     end if
  else
     asum=0.d0
  end if
  

  psb_dasum=asum

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_dasum


! Function: psb_dasumv 
!    Computes norm1 of X
!
!    norm1 := sum(X(i))
!
! Parameters:
!    x      -  real,dimension(:).       The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!
function psb_dasumv (x,desc_a, info)

  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info
  real(kind(1.d0))                  :: psb_dasumv

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, temp(2), jx, ix, ijx, m, i
  real(kind(1.d0))         :: asum, dasum
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_dasumv'
  info=0
  call psb_erractionsave(err_act)

  asum=0.d0

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
  jx=1

  m = desc_a%matrix_data(psb_m_)

  ! check vector correctness
  call psb_chkvect(m,1,size(x),ix,jx,desc_a%matrix_data,info,iix,jjx)
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
  if ((m.ne.0)) then
     if(desc_a%matrix_data(psb_n_row_).gt.0) then
        asum=dasum(desc_a%matrix_data(psb_n_row_),x,ione)

        ! adjust asum because overlapped elements are computed more than once
        i=1
        do while (desc_a%ovrlap_elem(i).ne.-ione)
           asum = asum -&
                & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
                & tmpx(desc_a%ovrlap_elem(i))
           i = i+2
        end do

        ! compute global sum
        call dgsum2d(icontxt, 'A', ' ', ione, ione, asum,&
             & ione, mone ,mycol)

     else
        asum=0.d0
        ! compute global sum
        call dgsum2d(icontxt, 'A', ' ', ione, ione, asum,&
             & ione, mone ,mycol)
     end if
  else
     asum=0.d0
  end if
  
  psb_dasumv=asum

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end function psb_dasumv


! Subroutine: psb_dasum vs
!    Computes norm1 of X
!
!    norm1 := sum(X(i))
!
! Parameters:
!    res    -  real.                      The result.
!    x      -  real,dimension(:).         The input vector.
!    desc_a -  type(<psb_desc_type>).     The communication descriptor.
!    info   -  integer.                   Eventually returns an error code.
!    jx     -  integer(optional).         The column offset.
!
subroutine psb_dasumvs (res,x,desc_a, info)
  use psb_serial_mod
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)      :: x(:)
  real(kind(1.d0)), intent(out)     :: res
  type(psb_desc_type), intent(in)   :: desc_a
  integer, intent(out)              :: info

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, temp(2), ix, jx, ijx, m, i
  real(kind(1.d0))         :: asum, dasum
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_dasumvs'
  info=0
  call psb_erractionsave(err_act)

  asum=0.d0

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

  ! check vector correctness
  call psb_chkvect(m,1,size(x),ix,jx,desc_a%matrix_data,info,iix,jjx)
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
  if ((m.ne.0)) then
     if(desc_a%matrix_data(psb_n_row_).gt.0) then
        asum=dasum(desc_a%matrix_data(psb_n_row_),x,ione)

        ! adjust asum because overlapped elements are computed more than once
        i=1
        do while (desc_a%ovrlap_elem(i).ne.-ione)
           asum = asum -&
                & (desc_a%ovrlap_elem(i+1)-1)/desc_a%ovrlap_elem(i+1)*&
                & tmpx(desc_a%ovrlap_elem(i))
           i = i+2
        end do

        ! compute global sum
        call dgsum2d(icontxt, 'A', ' ', ione, ione, asum,&
             & ione, mone ,mycol)

     else
        asum=0.d0
        ! compute global sum
        call dgsum2d(icontxt, 'A', ' ', ione, ione, asum,&
             & ione, mone ,mycol)
     end if
  else
     asum=0.d0
  end if
  

  res = asum

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_dasumvs



