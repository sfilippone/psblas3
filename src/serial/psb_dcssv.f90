! File:  psb_dcssv.f90 
! Subroutine: 
! Parameters:

subroutine psb_dcssv(alpha,t,b,beta,c,info,trans,unitd,d)
  use psb_spmat_type
  use psb_error_mod
  implicit none

  type(psb_dspmat_type) :: t
  real(kind(1.d0))      :: alpha, beta, b(:), c(:)
  integer               :: info
  character, optional   :: trans, unitd
  real(kind(1.d0)), optional, target :: d(:)
  
  real(kind(1.d0)), allocatable :: work(:)
  real(kind(1.d0)), pointer :: ddl(:)
  character :: lt, lu
  integer   :: iwsz,m,n,lb,lc,err_act
  character(len=20)                 :: name, ch_err

  name='psb_dcssv'
  info  = 0
  call psb_erractionsave(err_act)
  
  
  if (present(trans)) then 
    lt = trans
  else
    lt = 'N'
  endif
  if (present(unitd)) then 
    lu = unitd
  else
    lu = 'U'
  endif
  if (present(d)) then 
    ddl => d
  else
    allocate(ddl(1))
  endif

  m = t%m
  n = 1
  lb = size(b,1)
  lc = size(c,1)
  iwsz = 2*m*n
  allocate(work(iwsz))
  
  call dcssm(lt,m,n,alpha,lu,ddl,&
       & t%pl,t%fida,t%descra,t%aspk,t%ia1,t%ia2,t%infoa,t%pr,&
       & b,lb,beta,c,lc,work,iwsz,info)
  
  if (.not.present(d)) then 
    deallocate(ddl)
  endif
  deallocate(work)
  call psb_erractionrestore(err_act)

  if(info.ne.0) then
     if (err_act.eq.act_abort) then
        call psb_error()
        return
     end if
  end if
     
  return

end subroutine psb_dcssv
