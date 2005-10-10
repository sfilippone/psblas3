subroutine psb_dprecfree(p,info)
  !...free sparse matrix structure...
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psb_prec_type
  use psb_tools_mod
  use psb_error_mod
  implicit none
  !....parameters...

  type(psb_dprec_type), intent(inout) :: p
  integer, intent(out)                :: info

  !...locals....
  integer             :: int_err(5)
  integer             :: temp(1), me
  real(kind(1.d0))    :: real_err(5)
  integer             :: icontxt,err_act,i
  integer,parameter   :: ione=1
  character(len=20)   :: name, ch_err

  info=0
  name = 'psdprecfree'
  call psb_erractionsave(err_act)

  me=-1

!!$  if (associated(p%baseprec)) then 
!!$    call base_precfree(p%baseprec,info)
!!$        if (info /= 0) then
!!$           info=4010
!!$           ch_err='base_precfree'
!!$           call psb_errpush(info,name,a_err=ch_err)
!!$           goto 9999
!!$        end if
!!$    deallocate(p%baseprec,stat=info)
!!$    nullify(p%baseprec) 
!!$  endif
!!$
!!$  if (associated(p%mlprec)) then 
!!$    ! Check this !!!!! 
!!$    call base_precfree(p%mlprec,info)
!!$    if (info /= 0) then
!!$      write(0,*) 'From Base_precfree',info
!!$    end if
!!$    deallocate(p%mlprec,stat=info)
!!$    nullify(p%mlprec) 
!!$  endif

  if (associated(p%baseprecv)) then 
    do i=1,size(p%baseprecv) 
      call psb_base_precfree(p%baseprecv(i),info)
    end do
    deallocate(p%baseprecv)
    p%baseprecv => null()
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

end subroutine psb_dprecfree
