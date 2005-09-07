! File:  psb_dcsrws.f90 
! Subroutine: 
! Parameters:

subroutine psb_dcsrws(rw,a,info,trans)
  use psb_spmat_type
  use psb_error_mod
  implicit none 

  type(psb_dspmat_type)      :: a
  real(kind(1.d0)), pointer  :: rw(:) 
  integer                    :: info
  character, optional        :: trans

  Interface dcsrws
    subroutine  dcsrws(trans,m,n,fida,descra,a,ia1,ia2,&
         &                infoa,rowsum,ierror)
      integer, intent(in)        :: m,n
      integer, intent(out)       :: ierror
      double precision, intent(in) :: a(*)
      double precision, intent(out) :: rowsum(*)
      integer, intent(in)          :: ia1(*), ia2(*), infoa(*)
      character, intent(in)        :: descra*11,fida*5,trans*1
    end subroutine dcsrws
  end interface

  character :: trans_
  integer   :: iwsz,m,n,k,lb,lc,err_act
  character(len=20)                 :: name, ch_err

  name='psb_dcsrws'
  info  = 0
  call psb_erractionsave(err_act)

  if (present(trans)) then 
    trans_ = trans
  else
    trans_ = 'N'
  endif

  if (trans_=='N') then 
    m = a%m
    k = a%k
  else
    k = a%m
    m = a%k
  end if

  if (size(rw) < m) then 
    call psb_realloc(m,rw,info)
    if (info /= 0) then
       info = 4000
       call psb_errpush(info,name)
       goto 9999
    end if
  end if

  call  dcsrws(trans,m,k,a%fida,a%descra,&
       & a%aspk,a%ia1,a%ia2,a%infoa,rw,info)


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return

end subroutine psb_dcsrws
