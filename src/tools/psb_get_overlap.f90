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
  name='psb_get_overlap'
  call psb_erractionsave(err_act)

  if (psb_is_ovl_asb(desc)) then 
    i=0
    j=1
    do while(desc%ovrlap_elem(j) /= -1) 
      i  = i +1 
      j  = j + 2
    enddo

    call psb_realloc(i,ovrel,info)
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
    info = 1122
    call psb_errpush(info,name)
    goto 9999
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

end subroutine psb_get_ovrlap
