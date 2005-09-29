! File:  psb_dovrl.f90
!
! Subroutine: psb_dovrlm
!   This subroutine performs the exchange of the overlap elements in a distributed dense matrix between all the processes.
!
! Parameters:
!   x           -  real,dimension(:,:).          The local part of the dense matrix.
!   desc_a      -  type(<psb_desc_type>).        The communication descriptor.
!   info        -  integer.                      Eventually returns an error code.
!   jx          -  integer(optional).            The starting column of the global matrix. 
!   ik          -  integer(optional).            The number of columns to gather. 
!   work        -  real(optional).               A working area.
!   choice      -  logical(optional).            ???.
!   update_type -  integer(optional).            ???.
!
subroutine  psb_dovrlm(x,desc_a,info,jx,ik,work,choice,update_type)
  use psb_descriptor_type
  use psb_const_mod
  use psi_mod
  use psb_realloc_mod
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(inout), target   :: x(:,:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), optional, target        :: work(:)
  logical, intent(in), optional             :: choice
  integer, intent(in), optional             :: update_type,jx,ik

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, m, n, iix, jjx, temp(2), ix, ijx, nrow, ncol, k, maxk, iupdate,&
       & imode, err, liwork, i
  real(kind(1.d0)),pointer :: iwork(:), xp(:,:)
  logical                  :: ichoice
  character(len=20)        :: name, ch_err

  name='psb_dovrlm'
  info=0
  call psb_erractionsave(err_act)

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
  n = desc_a%matrix_data(psb_n_)
  nrow = desc_a%matrix_data(psb_n_row_)
  ncol = desc_a%matrix_data(psb_n_col_)

  maxk=size(x,2)-ijx+1

  if(present(ik)) then
     if(ik.gt.maxk) then
        k=maxk
     else
        k=ik
     end if
  else
     k = maxk
  end if

  if (present(choice)) then     
     ichoice = choice
  else
     ichoice = .true.
  endif
  if (present(update_type)) then 
     iupdate = update_type
  else
     iupdate = psb_none_
  endif

  imode = IOR(psb_swap_send_,psb_swap_recv_)

  ! check vector correctness
  call psb_chkvect(m,1,size(x,1),ix,ijx,desc_a%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
  end if

  err=info
  call psb_errcomm(icontxt,err)
  if(err.ne.0) goto 9999

  ! check for presence/size of a work area
  liwork=ncol
  if (present(work)) then
     if(size(work).ge.liwork) then
        iwork => work
     else
        call psb_realloc(liwork,iwork,info)
        if(info.ne.0) then
           info=4010
           ch_err='psb_realloc'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if
     end if
  else
     call psb_realloc(liwork,iwork,info)
     if(info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
  end if

  ! exchange overlap elements
  if(ichoice) then
     xp => x(iix:size(x,1),jjx:jjx+k-1)
     call psi_swapdata(imode,k,1.d0,xp,&
          & desc_a,iwork,info,data=psb_comm_ovr_)
  end if

  if(info.ne.0) then
     call psb_errpush(4010,name,a_err='psi_swapdata')
     goto 9999
  end if

  i=0
  ! switch on update type
  select case (iupdate)
  case(psb_square_root_)
     do while(desc_a%ovrlap_elem(i).ne.-ione)
        x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:) =&
             & x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:)/&
             & sqrt(real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_)))
        i = i+2
     end do
  case(psb_avg_)
     do while(desc_a%ovrlap_elem(i).ne.-ione)
        x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:) =&
             & x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:)/&
             & real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_))
        i = i+2
     end do
  case(psb_sum_)
     ! do nothing
  case default 
     ! wrong value for choice argument
     info = 70
     int_err=(/10,iupdate,0,0,0/)
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  end select



  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_dovrlm






! Subroutine: psb_dovrlv
!   This subroutine performs the exchange of the overlap elements in a distributed dense vector between all the processes.
!
! Parameters:
!   x           -  real,dimension(:).          The local part of the dense vector.
!   desc_a      -  type(<psb_desc_type>).        The communication descriptor.
!   info        -  integer.                      Eventually returns an error code.
!   work        -  real(optional).               A working area.
!   choice      -  logical(optional).            ???.
!   update_type -  integer(optional).            ???.
!
subroutine  psb_dovrlv(x,desc_a,info,work,choice,update_type)
  use psb_descriptor_type
  use psi_mod
  use psb_const_mod
  use psb_realloc_mod
  use psb_check_mod
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(inout), target   :: x(:)
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), optional, target        :: work(:)
  logical, intent(in), optional             :: choice
  integer, intent(in), optional             :: update_type

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, m, n, iix, jjx, temp(2), ix, ijx, nrow, ncol, k, maxk, iupdate,&
       & imode, err, liwork, i
  real(kind(1.d0)),pointer :: iwork(:)
  logical                  :: ichoice
  character(len=20)        :: name, ch_err

  name='psb_dovrlv'
  info=0
  call psb_erractionsave(err_act)

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
  ijx = 1

  m = desc_a%matrix_data(psb_m_)
  n = desc_a%matrix_data(psb_n_)
  nrow = desc_a%matrix_data(psb_n_row_)
  ncol = desc_a%matrix_data(psb_n_col_)

  k = 1

  if (present(choice)) then     
     ichoice = choice
  else
     ichoice = .true.
  endif
  if (present(update_type)) then 
     iupdate = update_type
  else
     iupdate = psb_none_
  endif

  imode = IOR(psb_swap_send_,psb_swap_recv_)

  ! check vector correctness
  call psb_chkvect(m,1,size(x),ix,ijx,desc_a%matrix_data,info,iix,jjx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chkvect'
     call psb_errpush(info,name,a_err=ch_err)
  end if

  if (iix.ne.1) then
     info=3040
     call psb_errpush(info,name)
  end if

  err=info
  call psb_errcomm(icontxt,err)
  if(err.ne.0) goto 9999

  ! check for presence/size of a work area
  liwork=ncol
  if (present(work)) then
     if(size(work).ge.liwork) then
        iwork => work
     else
        call psb_realloc(liwork,iwork,info)
        if(info.ne.0) then
           info=4010
           ch_err='psb_realloc'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        end if
     end if
  else
     call psb_realloc(liwork,iwork,info)
     if(info.ne.0) then
        info=4010
        ch_err='psb_realloc'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
  end if

  ! exchange overlap elements
  if(ichoice) then
     call psi_swapdata(imode,1.d0,x(iix:size(x)),&
          & desc_a,iwork,info,data=psb_comm_ovr_)
  end if

  if(info.ne.0) then
     call psb_errpush(4010,name,a_err='PSI_SwapData')
     goto 9999
  end if

  i=0
  ! switch on update type
  select case (iupdate)
  case(psb_square_root_)
     do while(desc_a%ovrlap_elem(i).ne.-ione)
        x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_)) =&
             & x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_))/&
             & sqrt(real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_)))
        i = i+2
     end do
  case(psb_avg_)
     do while(desc_a%ovrlap_elem(i).ne.-ione)
        x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_)) =&
             & x(desc_a%ovrlap_elem(i+psb_ovrlp_elem_))/&
             & real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_))
        i = i+2
     end do
  case(psb_sum_)
     ! do nothing
  case default 
     ! wrong value for choice argument
     info = 70
     int_err=(/10,iupdate,0,0,0/)
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  end select



  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_dovrlv
