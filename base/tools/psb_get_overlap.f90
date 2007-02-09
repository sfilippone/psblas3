subroutine psb_get_ovrlap(ovrel,desc,info)
  use psb_descriptor_type
  use psb_realloc_mod
  use psb_error_mod
  implicit none 
  integer, allocatable            :: ovrel(:)
  type(psb_desc_type), intent(in) :: desc
  integer, intent(out)            :: info

  integer  :: i,j, err_act
  character(len=20)    :: name

  info = 0
  name='psi_get_overlap'
  call psb_erractionsave(err_act)

  if (.not.psb_is_asb_desc(desc)) then
    info = 1122
    call psb_errpush(info,name)
    goto 9999
  end if

  i=0
  j=1
  do while(desc%ovrlap_elem(j) /= -1) 
    i  = i +1 
    j  = j + 2
  enddo

  if (i > 0) then 

    allocate(ovrel(i),stat=info)
    if (info /= 0 ) then 
      info = 4000
      call psb_errpush(info,name)
      goto 9999
    end if
    
    i=0
    j=1
    do while(desc%ovrlap_elem(j) /= -1) 
      i  = i +1 
      ovrel(i) = desc%ovrlap_elem(j) 
      j  = j + 2
    enddo

  else

    if (allocated(ovrel)) then 
      deallocate(ovrel,stat=info)
      if (info /= 0) then 
        call psb_errpush(4010,name,a_err='Deallocate')
        goto 9999      
      end if
    end if

  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_get_ovrlap
