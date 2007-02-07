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
!
!
! Note: This subroutine performs the symbolic product of two sparse matrices.
!       It is modeled after the SMMP package by R. Bank and C. Douglas, but is 
!       rewritten in Fortran 95 making use of our sparse matrix facilities.
!

subroutine psb_dsymbmm(a,b,c,info)
  use psb_spmat_type
  use psb_string_mod
  use psb_serial_mod, only : psb_msort
  implicit none 

  type(psb_dspmat_type) :: a,b,c
  integer, allocatable  :: itemp(:)
  integer               :: nze,info

  interface 
    subroutine symbmm (n, m, l, ia, ja, diaga, &
         & ib, jb, diagb, ic, jc, diagc, index)
      integer  n,m,l,  ia(*), ja(*), diaga, ib(*), jb(*), diagb,&
           & diagc,  index(*)
      integer, allocatable :: ic(:),jc(:)
    end subroutine symbmm 
  end interface
  interface psb_sp_getrow
    subroutine psb_dspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw)
      use psb_spmat_type
      type(psb_dspmat_type), intent(in) :: a
      integer, intent(in)       :: irw
      integer, intent(out)      :: nz
      integer, intent(inout)    :: ia(:), ja(:)
      real(kind(1.d0)),  intent(inout)    :: val(:)
      integer, intent(in), target, optional :: iren(:)
      integer, intent(in), optional :: lrw
      integer, intent(out)  :: info
    end subroutine psb_dspgetrow
  end interface

  character(len=20)         :: name, ch_err
  integer                   :: err_act
  name='psb_symbmm'
  call psb_erractionsave(err_act)

  select case(toupper(a%fida(1:3)))
  case  ('CSR')
  case default
    info=135
    ch_err=a%fida(1:3)
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end select
  select case(toupper(b%fida(1:3)))
  case  ('CSR')
  case default
    info=136
    ch_err=b%fida(1:3)
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end select

  if (b%m /= a%k) then 
    write(0,*) 'Mismatch in SYMBMM: ',a%m,a%k,b%m,b%k
  endif
  allocate(itemp(max(a%m,a%k,b%m,b%k)),stat=info)    
  if (info /= 0) then 
    return
  endif
  nze = max(a%m+1,2*a%m)
  call psb_sp_reall(c,nze,info)
  !
  ! Note: we need to test whether there is a performance impact 
  !       in not using the original Douglas & Bank code. 
  !
  if (.true.) then 
    call symbmm(a%m,a%k,b%k,a%ia2,a%ia1,0,&
         & b%ia2,b%ia1,0,&
         & c%ia2,c%ia1,0,itemp)
  else 
    call inner_symbmm(a,b,c,itemp,info)
  endif
  call psb_realloc(size(c%ia1),c%aspk,info)

  c%pl(1) = 0
  c%pr(1) = 0
  c%m=a%m
  c%k=b%k
  c%fida='CSR'
  c%descra='GUN'
  deallocate(itemp) 
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
     call psb_error()
     return
  end if
  return

contains

  subroutine inner_symbmm(a,b,c,index,info)
    type(psb_dspmat_type) :: a,b,c
    integer               :: index(:),info
    integer, allocatable  :: iarw(:), iacl(:),ibrw(:),ibcl(:)
    real(kind(1.d0)), allocatable :: aval(:),bval(:)
    integer  :: maxlmn,i,j,m,n,k,l,istart,length,nazr,nbzr,jj,ii,minlm,minmn


    n = a%m
    m = a%k 
    l = b%k 
    maxlmn = max(l,m,n)

    allocate(iarw(maxlmn),iacl(maxlmn),ibrw(maxlmn),ibcl(maxlmn),&
         & aval(maxlmn),bval(maxlmn), stat=info)
    if (info /= 0) then 
      return
    endif


    if (size(c%ia2) < n+1) then 

      call psb_realloc(n+1,c%ia2,info)
    endif
    do i=1,maxlmn
      index(i)=0
    end do

      c%ia2(1)=1
      minlm = min(l,m)
      minmn = min(m,n)

      main: do  i=1,n
        istart=-1
        length=0
        call psb_sp_getrow(i,a,nazr,iarw,iacl,aval,info)
        do jj=1, nazr
          
          j=iacl(jj)
          
          if ((j<1).or.(j>m)) then 
            write(0,*) ' SymbMM: Problem with A ',i,jj,j,m
            info = 1
            return
          endif
          call psb_sp_getrow(j,b,nbzr,ibrw,ibcl,bval,info)
          do k=1,nbzr
            if ((ibcl(k)<1).or.(ibcl(k)>maxlmn)) then 
                write(0,*) 'Problem in SYMBMM 1:',j,k,ibcl(k),maxlmn
                info=2
                return
            else
              if(index(ibcl(k)).eq.0) then
                index(ibcl(k))=istart
                istart=ibcl(k)
                length=length+1
              endif
            endif
          end do
        end do

        c%ia2(i+1)=c%ia2(i)+length
        
        if (c%ia2(i+1) > size(c%ia1)) then 
          if (n > (2*i)) then 
            nze = max(c%ia2(i+1), c%ia2(i)*((n+i-1)/i))
          else
            nze = max(c%ia2(i+1), nint((dble(c%ia2(i))*(dble(n)/i)))   )
          endif 
          call psb_realloc(nze,c%ia1,info)
        end if 
        do j= c%ia2(i),c%ia2(i+1)-1
          c%ia1(j)=istart
          istart=index(istart)
          index(c%ia1(j))=0
        end do
        call psb_msort(c%ia1(c%ia2(i):c%ia2(i)+length-1))
        index(i) = 0
      end do main

  end subroutine inner_symbmm

end subroutine psb_dsymbmm
