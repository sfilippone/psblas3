! File:  psb_dcsmv.f90 
! Subroutine: 
! Parameters:

subroutine psb_dcsmv(alpha,a,b,beta,c,info,trans)
  use psb_spmat_type
  implicit none 

  type(psb_dspmat_type) :: a
  real(kind(1.d0))      :: alpha, beta, b(:), c(:)
  integer               :: info
  character, optional   :: trans
  
  real(kind(1.d0)), allocatable :: work(:)
  character :: trans_
  integer   :: iwsz,m,n,k,lb,lc, err_act
  character(len=20)                 :: name, ch_err

  name='psb_dcsmv'
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
  n  = 1
  lb = size(b,1)
  lc = size(c,1)
  iwsz = 2*m*n
  allocate(work(iwsz))
  
  call dcsmm(trans_,m,n,k,alpha,&
       & a%pl,a%fida,a%descra,a%aspk,a%ia1,a%ia2,a%infoa,a%pr,&
       & b,lb,beta,c,lc,work,iwsz,info)
  
  deallocate(work)
  call psb_erractionrestore(err_act)

  if(info.ne.0) then
     if (err_act.eq.act_abort) then
        call psb_error()
        return
     end if
  end if
     
  return

end subroutine psb_dcsmv
