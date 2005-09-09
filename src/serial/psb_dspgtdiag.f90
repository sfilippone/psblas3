! File:  psb_dspgtdiag.f90 
! Subroutine: 
! Parameters:

!*****************************************************************************
!*                                                                           *
!* Takes a specified row from matrix A and copies into matrix B (possibly    *
!*  appending to B). Output is always COO. Input might be anything, once     *
!* we get to actually write the code.....                                    *
!*                                                                           *
!*****************************************************************************
subroutine psb_dspgtdiag(a,d,info)
  ! Output is always in  COO format  into B, irrespective of 
  ! the input format 
  use psb_spmat_type
  use psb_error_mod
  use psb_const_mod
  implicit none

  type(psb_dspmat_type), intent(in)     :: a
  real(kind(1.d0)), intent(inout)       :: d(:) 
  integer, intent(out)                  :: info

  integer :: i,j,k,nr, nz, err_act
  character(len=20)                 :: name, ch_err

  name='psb_dspgtdiag'
  info  = 0
  call psb_erractionsave(err_act)

  if (size(d) < min(a%k,a%m)) then 
    write(0,*) 'Insufficient space in DSPGTDIAG ', size(d),min(a%m,a%k)
  end if
  d(:) = 0.d0
  if (a%fida == 'CSR') then 
    
    do i=1, min(a%m,a%k)
      do j=a%ia2(i),a%ia2(i+1)-1
        if (a%ia1(j) == i) then 
          d(i) = a%aspk(j)
        end if
      end do
    end do

  else if (a%fida == 'COO') then 

    do i=1,a%infoa(psb_nnz_)
      j=a%ia1(i)
      if ((j==a%ia2(i)).and.(j <= min(a%k,a%m)) .and.(j>0)) then 
        d(j) = a%aspk(i)
      endif
    enddo
    
 else if (a%fida == 'JAD') then 
    info=135
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
 
end subroutine psb_dspgtdiag

