! File:  psb_dsymbmm.f90 
! Subroutine: 
! Parameters:

subroutine psb_dsymbmm(a,b,c)
  use psb_spmat_type
  implicit none 

  type(psb_dspmat_type) :: a,b,c
  integer, allocatable  :: itemp(:)
  integer               :: nze,info

  interface 
    subroutine symbmm (n, m, l, ia, ja, diaga, &
         & ib, jb, diagb, ic, jc, diagc, index)
      integer  n,m,l,  ia(*), ja(*), diaga, ib(*), jb(*), diagb,&
           & diagc,  index(*)
      integer, pointer :: ic(:),jc(:)
    end subroutine symbmm
  end interface

  if (b%m /= a%k) then 
    write(0,*) 'Mismatch in SYMBMM: ',a%m,a%k,b%m,b%k
  endif
  allocate(itemp(max(a%m,a%k,b%m,b%k)),stat=info)    
  nze = max(a%m+1,2*a%m)
  call psb_spreall(c,nze,info)
!!$  write(0,*) 'SYMBMM90 ',size(c%pl),size(c%pr)
  call symbmm(a%m,a%k,b%k,a%ia2,a%ia1,0,&
       & b%ia2,b%ia1,0,&
       & c%ia2,c%ia1,0,itemp)
  c%pl(1) = 0
  c%pr(1) = 0
  c%m=a%m
  c%k=b%k
  c%fida='CSR'
  deallocate(itemp) 
  return
end subroutine psb_dsymbmm
