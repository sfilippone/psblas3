! File:  psb_dipcsr2coo.f90 
! Subroutine: 
! Parameters:

Subroutine psb_dipcsr2coo(a,info)
  use psb_spmat_type
  use psb_error_mod
  implicit none

  !....Parameters...
  Type(psb_dspmat_type), intent(inout) :: A
  Integer, intent(out)                 :: info

  integer, pointer    :: iaux(:), itemp(:)
  !locals
  Integer             :: nza, nr
  integer             :: i,j,err_act
  logical, parameter  :: debug=.false.
  character(len=20)                 :: name, ch_err

  name='psb_dipcsr2coo'
  info  = 0
  call psb_erractionsave(err_act)

  if (a%fida /= 'CSR') then 
    info = 5
    call psb_errpush(info,name)
    goto 9999
  end if

  nr  = a%m 
  nza = a%ia2(nr+1) - 1
  allocate(iaux(nza),stat=info)
  if (info /=0) then 
    write(0,*) 'Failed allocation ',info, nza
    return
  end if
!!$  write(0,*) 'ipcsr2coo ',a%m      
  itemp => a%ia2
  a%ia2 => a%ia1
  a%ia1 => iaux
  
  do i=1, nr
    do j=itemp(i),itemp(i+1)-1
      a%ia1(j) = i
    end do
  end do
  
  a%fida='COO'
  a%infoa(nnz_) = nza
  a%infoa(srtd_) = isrtdcoo
  deallocate(itemp)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return

end Subroutine psb_dipcsr2coo
