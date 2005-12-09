
! File:  psb_ihalo.f90
!
! Subroutine: psb_ihalom
!   This subroutine performs the exchange of the halo elements in a distributed dense matrix between all the processes.
!
! Parameters:
!   x         -  integer,dimension(:,:).       The local part of the dense matrix.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   alpha     -  real(optional).               ???.
!   jx        -  integer(optional).            The starting column of the global matrix. 
!   ik        -  integer(optional).            The number of columns to gather. 
!   work      -  integer(optional).            A working area.
!   tran      -  character(optional).          ???.
!   mode      -  integer(optional).
!
subroutine  psb_ihalom(x,desc_a,info,alpha,jx,ik,work,tran,mode)
  use psb_descriptor_type
  use psb_const_mod
  use psi_mod
  use psb_realloc_mod
  use psb_check_mod
  use psb_error_mod
  implicit none

  integer, intent(inout), target           :: x(:,:)
  type(psb_desc_type), intent(in)          :: desc_a
  integer, intent(out)                     :: info
  real(kind(1.d0)), intent(in), optional   :: alpha
  integer, intent(inout), optional, target :: work(:)
  integer, intent(in), optional            :: mode,jx,ik
  character, intent(in), optional          :: tran

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, m, n, iix, jjx, temp(2), ix, ijx, nrow, ncol, k, maxk, liwork,&
       & imode, err
  integer, pointer         :: xp(:,:), iwork(:)
  character                :: ltran
  character(len=20)        :: name, ch_err

  name='psb_ihalom'
  if(psb_get_errstatus().ne.0) return 
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

  if (present(tran)) then     
     ltran = tran
  else
     ltran = 'N'
  endif
  if (present(mode)) then 
     imode = mode
  else
     imode = IOR(psb_swap_send_,psb_swap_recv_)
  endif

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


  ! we should write an "iscal"
!!$  if(present(alpha)) then
!!$     if(alpha.ne.1.d0) then
!!$        do i=0, k-1
!!$           call iscal(nrow,alpha,x(1,jjx+i),1)
!!$        end do
!!$     end if
!!$  end if

  liwork=nrow
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

  xp => x(iix:size(x,1),jjx:jjx+k-1)
  ! exchange halo elements
  if(ltran.eq.'N') then
     call psi_swapdata(imode,k,0,xp,&
          & desc_a,iwork,info)
  else if((ltran.eq.'T').or.(ltran.eq.'H')) then
     call psi_swaptran(imode,k,1,xp,&
          & desc_a,iwork,info)
  end if

  if(info.ne.0) then
     call psb_errpush(4010,name,a_err='PSI_iSwap...')
     goto 9999
  end if

  if(.not.present(work)) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_ihalom





! Subroutine: psb_ihalov
!   This subroutine performs the exchange of the halo elements in a distributed dense matrix between all the processes.
!
! Parameters:
!   x         -  integer,dimension(:).         The local part of the dense matrix.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   alpha     -  real(optional).               ???.
!   work      -  integer(optional).            A working area.
!   tran      -  character(optional).          ???.
!   mode      -  integer(optional).
!
subroutine  psb_ihalov(x,desc_a,info,alpha,work,tran,mode)
  use psb_descriptor_type
  use psb_const_mod
  use psi_mod
  use psb_realloc_mod
  use psb_check_mod
  use psb_error_mod
  implicit none

  integer, intent(inout)                   :: x(:)
  type(psb_desc_type), intent(in)          :: desc_a
  integer, intent(out)                     :: info
  real(kind(1.d0)), intent(in), optional   :: alpha
  integer, intent(inout), optional, target :: work(:)
  integer, intent(in), optional            :: mode
  character, intent(in), optional          :: tran

  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, m, n, iix, jjx, temp(2), ix, ijx, nrow, ncol, k, maxk, imode,&
       & err, liwork
  integer,pointer          :: iwork(:)
  character                :: ltran
  character(len=20)        :: name, ch_err

  name='psb_ihalov'
  if(psb_get_errstatus().ne.0) return 
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
!  ncol = desc_a%matrix_data(psb_n_col_)
  

  if (present(tran)) then     
     ltran = tran
  else
     ltran = 'N'
  endif
  if (present(mode)) then 
     imode = mode
  else
     imode = IOR(psb_swap_send_,psb_swap_recv_)
  endif

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

!!$  if(present(alpha)) then
!!$     if(alpha.ne.1.d0) then
!!$        call dscal(nrow,alpha,x,1)
!!$     end if
!!$  end if

  liwork=nrow
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

  ! exchange halo elements
  if(ltran.eq.'N') then
     call psi_swapdata(imode,0,x(iix:size(x)),&
          & desc_a,iwork,info)
  else if((ltran.eq.'T').or.(ltran.eq.'H')) then
     call psi_swaptran(imode,1,x(iix:size(x)),&
          & desc_a,iwork,info)
  end if

  if(info.ne.0) then
     call psb_errpush(4010,name,a_err='PSI_iSwap...')
     goto 9999
  end if

  if(.not.present(work)) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psb_ihalov



