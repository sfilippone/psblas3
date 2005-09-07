! File: psb_dallc.f90
!
! Function: psb_dalloc
!    Allocates dense matrix for PSBLAS routines
! 
! Parameters: 
!    m      - number of rows.
!    n      - number of columns.
!    x      - the matrix to be allocated.
!    desc_a - the communication descriptor.
!    info   - eventually returns an error code
!    js     - (optional) the starting column
subroutine psb_dalloc(m, n, x, desc_a, info, js)
  !....allocate dense  matrix for psblas routines.....
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod

  implicit none
  
  !....parameters...
  integer, intent(in)                   :: m,n
  real(kind(1.d0)), pointer             :: x(:,:)
  type(psb_desc_type), intent(inout)    :: desc_a
  integer                               :: info
  integer, optional, intent(in)         :: js

  !locals
  integer             :: nprow,npcol,myrow,mycol,err,n_col,n_row,i,j,jj,err_act
  integer             :: icontxt,dectype
  integer             :: int_err(5),temp(1),exch(3)
  real(kind(1.d0))    :: real_err(5)
  integer, allocatable:: prc_v(:)
  character(len=20)   :: name, ch_err

  info=0
  err=0
  int_err(1)=0
  call psb_erractionsave(err_act)
  
  icontxt=desc_a%matrix_data(psb_ctxt_)
  
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  !     ....verify blacs grid correctness..
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,int_err)
     goto 9999
  endif
  
  dectype=desc_a%matrix_data(psb_dec_type_)
  !... check m and n parameters....
  if (m.lt.0) then
     info = 10
     int_err(1) = 1
     int_err(2) = m
     call psb_errpush(info,name,int_err)
     goto 9999
  else if (n.lt.0) then
     info = 10
     int_err(1) = 2
     int_err(2) = n
     call psb_errpush(info,name,int_err)
  else if (.not.psb_is_ok_dec(dectype)) then 
     info = 3110
     call psb_errpush(info,name)
     goto 9999
  else if (m.ne.desc_a%matrix_data(psb_n_)) then
     info = 300
     int_err(1) = 1
     int_err(2) = m
     int_err(3) = 4
     int_err(4) = psb_n_
     int_err(5) = desc_a%matrix_data(psb_n_)
     call psb_errpush(info,name,int_err)
     goto 9999
  endif
  
  if (present(js)) then 
    j=js
  else
    j=1
  endif
  !global check on m and n parameters
  if (myrow.eq.root) then
     exch(1)=m
     exch(2)=n
     exch(3)=j
     call igebs2d(icontxt,psb_all_,psb_topdef_, ithree,ione, exch, ithree)
  else
     call igebr2d(icontxt,psb_all_,psb_topdef_, ithree,ione, exch, ithree, root, 0)
     if (exch(1).ne.m) then
	info=550
	int_err(1)=1
        call psb_errpush(info,name,int_err)
        goto 9999
     else if (exch(2).ne.n) then
	info=550
	int_err(1)=2
        call psb_errpush(info,name,int_err)
        goto 9999
     else if (exch(3).ne.j) then
	info=550
	int_err(1)=3
        call psb_errpush(info,name,int_err)
        goto 9999
     endif
  endif

  !....allocate x .....
  if (psb_is_asb_dec(dectype).or.psb_is_upd_dec(dectype)) then
     n_col = max(1,desc_a%matrix_data(psb_n_col_))
     allocate(x(n_col,j:j+n-1),stat=info)
!     call sprealloc(n_col,j:j+n-1,x,info)
    if (info.ne.0) then
        info=4010
        ch_err='psb_sprealloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     endif
     do jj=j,j+n-1
       do i=1,n_col
         x(i,j) = 0.0d0
       end do
     end do
  else if (psb_is_bld_dec(dectype)) then
     n_row = max(1,desc_a%matrix_data(psb_n_row_))
     allocate(x(n_row,j:j+n-1),stat=info)
!     call sprealloc(n_row,j:j+n-1,x,info)
     if (info.ne.0) then
        info=4010
        ch_err='psb_sprealloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     endif
     do jj=j,j+n-1
       do i=1,n_row
         x(i,j) = 0.0d0
       end do
     end do
  endif

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dalloc



! Function: psb_dallocv
!    Allocates dense matrix for PSBLAS routines
! 
! Parameters: 
!    m      - number of rows.
!    x      - the matrix to be allocated.
!    desc_a - the communication descriptor.
!    info   - eventually returns an error code
subroutine psb_dallocv(m, x, desc_a,info)
  !....allocate sparse matrix structure for psblas routines.....
  use psb_descriptor_type
  use psb_const_mod
  use psb_realloc_mod
  use psb_error_mod

  implicit none
  
  !....parameters...
  integer, intent(in)            :: m
  real(kind(1.d0)), pointer      :: x(:)
  type(psb_desc_type), intent(in):: desc_a
  integer                        :: info

  !locals
  integer             :: nprow,npcol,myrow,mycol,err,n_col,n_row,dectype,i,err_act
  integer             :: icontxt
  integer             :: int_err(5),temp(1),exch(2)
  real(kind(1.d0))    :: real_err(5)
  logical, parameter  :: debug=.false. 
  character(len=20)   :: name, ch_err

  info=0
  call psb_erractionsave(err_act)
   
  icontxt=desc_a%matrix_data(psb_ctxt_)
  
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  !     ....verify blacs grid correctness..
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,int_err)
     goto 9999
  endif
  
  dectype=desc_a%matrix_data(psb_dec_type_)
  if (debug) write(0,*) 'dall: dectype',dectype
  if (debug) write(0,*) 'dall: is_ok? dectype',psb_is_ok_dec(dectype)
  !... check m and n parameters....
  if (m.lt.0) then
     info = 10
     int_err(1) = 1
     int_err(2) = m
     call psb_errpush(info,name,int_err)
     goto 9999
  else if (.not.psb_is_ok_dec(dectype)) then 
     info = 3110
     call psb_errpush(info,name)
     goto 9999
  else if (m.ne.desc_a%matrix_data(psb_n_)) then
     info = 300
     int_err(1) = 1
     int_err(2) = m
     int_err(3) = 4
     int_err(4) = psb_n_
     int_err(5) = desc_a%matrix_data(psb_n_)
     call psb_errpush(info,name,int_err)
     goto 9999
  endif

  !global check on m and n parameters
  if (myrow.eq.root) then
     exch(1) = m
     call igebs2d(icontxt,psb_all_,psb_topdef_, ione,ione, exch, ione)
  else
     call igebr2d(icontxt,psb_all_,psb_topdef_, ione,ione, exch, ione, root, 0)
     if (exch(1) .ne. m) then
	info = 550
	int_err(1) = 1
        call psb_errpush(info,name,int_err)
        goto 9999
     endif
  endif

  !....allocate x .....
  if (psb_is_asb_dec(dectype).or.psb_is_upd_dec(dectype)) then
     n_col = max(1,desc_a%matrix_data(psb_n_col_))
     call psb_realloc(n_col,x,info)
     if (info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     endif
     do i=1,n_col
       x(i) = 0.0d0
     end do

  else if (psb_is_bld_dec(dectype)) then
     n_row = max(1,desc_a%matrix_data(psb_n_row_))
     call psb_realloc(n_row,x,info)
     if (info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     endif
     do i=1,n_row
       x(i) = 0.0d0
     end do
  endif

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dallocv

