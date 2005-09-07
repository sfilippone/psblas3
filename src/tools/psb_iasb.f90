! File: psb_iasb.f90
!
! Subroutine: psb_iasb
!    Assembles a dense matrix for PSBLAS routines
! 
! Parameters: 
!    x       - integer,pointer,dimension(:,:).    The matrix to be assembled.
!    desc_a  - type(<psb_desc_type>).             The communication descriptor.
!    info    - integer.                           Eventually returns an error code
subroutine psb_iasb(x, desc_a, info)
  !....assembly dense matrix x .....
  use psb_descriptor_type
  use psb_const_mod
  use psb_psblas_mod
  use psb_error_mod
  implicit none

  type(psb_desc_type), intent(in) ::  desc_a
  integer, pointer                ::  x(:,:)
  integer, intent(out)            ::  info

  ! local variables
  integer :: icontxt,nprow,npcol,me,mypcol,temp,lwork,nrow,ncol,err_act
  integer, pointer ::  itemp(:,:)
  integer :: int_err(5), i1sz, i2sz, dectype, i
  real(kind(1.d0)) :: real_err(5)
  integer, parameter  :: ione=1
  real(kind(1.d0)),parameter    :: one=1
  logical, parameter :: debug=.false.
  character(len=20)   :: name, char_err

  info=0
  name='psb_iasb'
  call psb_erractionsave(err_act)

  if ((.not.associated(desc_a%matrix_data))) then
     info=3110
     call psb_errpush(info,name)
     return
  endif
  
  icontxt=desc_a%matrix_data(psb_ctxt_)
  dectype=desc_a%matrix_data(psb_dec_type_)
  
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
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

  ! check size
  icontxt=desc_a%matrix_data(psb_ctxt_)
  nrow=desc_a%matrix_data(psb_n_row_)
  ncol=desc_a%matrix_data(psb_n_col_)
  i1sz = size(x,dim=1)
  i2sz = size(x,dim=2)
  if (debug) write(*,*) 'asb: ',i1sz,i2sz,nrow,ncol
  if (i1sz.lt.ncol) then
     allocate(itemp(ncol,i2sz),stat=info)
     if (info.ne.0) then
        info=2025
        int_err(1)=ncol
        call psb_errpush(info,name,int_err)
        goto 9999
     endif
     itemp(nrow+1:,:) = 0
     itemp(1:nrow,:) = x(1:nrow,:)
     deallocate(x)
     x => itemp
  endif
  
  ! ..update halo elements..
  call psb_halo(x,desc_a,info,alpha=one)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
  
end subroutine psb_iasb



! Subroutine: psb_iasbv
!    Assembles a dense matrix for PSBLAS routines
! 
! Parameters: 
!    x       - integer,pointer,dimension(:).      The matrix to be assembled.
!    desc_a  - type(<psb_desc_type>).             The communication descriptor.
!    info    - integer.                           Eventually returns an error code
subroutine psb_iasbv(x, desc_a, info)
  !....assembly dense matrix x .....
  use psb_descriptor_type
  use psb_const_mod
  use psb_psblas_mod
  use psb_error_mod
  implicit none

  type(psb_desc_type), intent(in) ::  desc_a
  integer, pointer                ::  x(:)
  integer, intent(out)            ::  info

  ! local variables
  integer :: icontxt,nprow,npcol,me,mypcol,temp,lwork, err_act
  integer :: int_err(5), i1sz,nrow,ncol, dectype, i
  integer, pointer ::  itemp(:)
  real(kind(1.d0)) :: real_err(5)
  integer, parameter  :: ione=1  
  real(kind(1.d0)),parameter    :: one=1
  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  info=0
  call psb_erractionsave(err_act)
  name = 'psb_iasbv'
  
  
  icontxt=desc_a%matrix_data(psb_ctxt_)
  dectype=desc_a%matrix_data(psb_dec_type_)

  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
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

  nrow=desc_a%matrix_data(psb_n_row_)
  ncol=desc_a%matrix_data(psb_n_col_)
  if (debug) write(*,*) name,' sizes: ',nrow,ncol
  i1sz = size(x)
  if (debug) write(*,*) 'dasb: sizes ',i1sz,ncol
  if (i1sz.lt.ncol) then
    allocate(itemp(ncol),stat=info)  
    if (info.ne.0) then           
      info=2025
      int_err(1)=ncol
      call psb_errpush(info,name,int_err)
      goto 9999
    endif
    itemp(nrow+1:) = 0
    itemp(1:nrow) = x(1:nrow)
    deallocate(x)
    x => itemp
  endif  
  
  ! ..update halo elements..
  call psb_halo(x,desc_a,info,alpha=one)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
  
end subroutine psb_iasbv

