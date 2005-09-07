! File:  psb_dspscal.f90 
! Subroutine: 
! Parameters:

!*****************************************************************************
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine psb_dspscal(a,d,info)
  ! the input format 
  use psb_spmat_type
  use psb_error_mod
  use psb_const_mod
  implicit none

  type(psb_dspmat_type), intent(inout) :: a
  integer, intent(out)                 :: info
  real(kind(1.d0)), intent(in)         :: d(:) 

  integer :: i,j,k,nr, nz,err_act
  character(len=20)                 :: name, ch_err

  name='psb_dspscal'
  info  = 0
  call psb_erractionsave(err_act)


  if (a%fida == 'CSR') then 

     do i=1, a%m
        do j=a%ia2(i),a%ia2(i+1)-1
           a%aspk(j) = a%aspk(j) * d(i)
        end do
     end do

  else if (a%fida == 'COO') then 

     do i=1,a%infoa(nnz_)
        j=a%ia1(i)
        a%aspk(i) = a%aspk(i) * d(j)
     enddo

  else if (a%fida == 'JAD') then 
     info=135
     ch_err=a%fida(1:3)
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  else
     info=136
     ch_err=a%fida(1:3)
     call psb_errpush(info,name,a_err=ch_err)
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

end subroutine psb_dspscal

