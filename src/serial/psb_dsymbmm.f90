!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
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
