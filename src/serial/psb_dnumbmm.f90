! File:  psb_dnumbmm.f90 
! Subroutine: 
! Parameters:

subroutine psb_dnumbmm(a,b,c)
  use psb_spmat_type
  implicit none

  type(psb_dspmat_type)         :: a,b,c
  real(kind(1.d0)), allocatable :: temp(:)
  integer                       :: info

  allocate(temp(max(a%m,a%k,b%m,b%k)),stat=info)

  call psb_realloc(size(c%ia1),c%aspk,info)
  call numbmm(a%m,a%k,b%k,a%ia2,a%ia1,0,a%aspk,&
       & b%ia2,b%ia1,0,b%aspk,&
       & c%ia2,c%ia1,0,c%aspk,temp)
  deallocate(temp) 
  return
end subroutine psb_dnumbmm
