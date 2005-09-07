! File:  psb_drwextd.f90 
! Subroutine: 
! Parameters:

subroutine psb_drwextd(nr,a,info,b)
  use psb_spmat_type
  use psb_error_mod
  implicit none

  ! Extend matrix A up to NR rows with empty ones (i.e.: all zeroes)
  integer, intent(in)                            :: nr
  type(psb_dspmat_type), intent(inout)           :: a
  integer,intent(out)                            :: info
  type(psb_dspmat_type), intent(in), optional    :: b
  integer :: i,j,ja,jb,err_act
  character(len=20)                 :: name, ch_err

  name='psb_drwextd'
  info  = 0
  call psb_erractionsave(err_act)

  if (nr > a%m) then 

    if (a%fida == 'CSR') then 
      call psb_realloc(nr+1,a%ia2,info)
      if (present(b)) then 
        jb = b%ia2(b%m+1)-1
        call psb_realloc(size(a%ia1)+jb,a%ia1,info)
        call psb_realloc(size(a%aspk)+jb,a%aspk,info)
        do i=1, min(nr-a%m,b%m)
          ! Should use spgtrow. 
          ! Don't care for the time being.
          a%ia2(a%m+i+1) =  a%ia2(a%m+i) + b%ia2(i+1) - b%ia2(i)
          ja = a%ia2(a%m+i)
          jb = b%ia2(i)
          do 
            if (jb >=  b%ia2(i+1)) exit
            a%aspk(ja) = b%aspk(jb)
            a%ia1(ja) = b%ia1(jb)
            ja = ja + 1
            jb = jb + 1
          end do
        end do
        do j=i,nr-a%m
          a%ia2(a%m+i+1) = a%ia2(a%m+i)
        end do

      else
        do i=a%m+2,nr+1
          a%ia2(i) = a%ia2(i-1)
        end do
      end if
      a%m = nr
    else if (a%fida == 'COO') then 
      if (present(b)) then 
      else
      endif
      a%m = nr
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

end subroutine psb_drwextd
