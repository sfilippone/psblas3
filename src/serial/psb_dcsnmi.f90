! File:  psb_dcsnmi.f90 
! Subroutine: 
! Parameters:

real(kind(1.d0)) function psb_dcsnmi(a,info,trans)

  use psb_spmat_type
  use psb_error_mod
  implicit none

  type(psb_dspmat_type), intent(in)  :: a
  integer, intent(out)               :: info
  character, optional                :: trans

  interface
     real(kind(1.d0)) function dcsnmi(trans,m,n,fida,descra,a,ia1,ia2,&
          &                 infoa,ierror)
       integer          ::  m,n, ierror
       character        ::  trans
       integer          ::  ia1(*),ia2(*),infoa(*)
       character        ::  descra*11, fida*5
       real(kind(1.d0)) ::  a(*)
     end function dcsnmi
  end interface

  integer             :: err_act
  character           :: itrans
  character(len=20)   :: name, ch_err

  name='psb_dcsnmi'
  call psb_erractionsave(err_act)

  if(present(trans)) then
     itrans=trans
  else
     itrans='N'
  end if

  dcsnmi90 = dcsnmi(itrans,a%m,a%k,a%fida,a%descra,a%aspk,a%ia1,a%ia2,a%infoa,info)
  if(info/=0) then
     dcsnmi90 = -1
     info=4010
     ch_err='dcsnmi'
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

end function psb_dcsnmi
