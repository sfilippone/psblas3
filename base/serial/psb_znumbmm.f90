!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File:  psb_znumbmm.f90 
! Subroutine: 
! Arguments:
!
!
! Note: This subroutine performs the numerical product of two sparse matrices.
!       It is modeled after the SMMP package by R. Bank and C. Douglas, but is 
!       rewritten in Fortran 95 making use of our sparse matrix facilities.
!
!

subroutine psb_znumbmm(a,b,c)
  use psb_sparse_mod, psb_protect_name => psb_znumbmm
  implicit none 

  type(psb_z_sparse_mat), intent(in) :: a,b
  type(psb_z_sparse_mat), intent(inout)  :: c
  integer               :: info
  integer               :: err_act
  character(len=*), parameter ::  name='psb_numbmm'

  call psb_erractionsave(err_act)
  info = psb_success_

  if ((a%is_null()) .or.(b%is_null()).or.(c%is_null())) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  select type(aa=>c%a)
  type is (psb_z_csr_sparse_mat)
    call psb_numbmm(a%a,b%a,aa)
  class default
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end select

  call c%set_asb()

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
     call psb_error()
     return
  end if
  return

end subroutine psb_znumbmm

subroutine psb_zbase_numbmm(a,b,c)
  use psb_mat_mod
  use psb_string_mod
  use psb_serial_mod, psb_protect_name => psb_zbase_numbmm
  implicit none 

  class(psb_z_base_sparse_mat), intent(in) :: a,b
  type(psb_z_csr_sparse_mat), intent(inout)  :: c
  integer, allocatable  :: itemp(:)
  integer               :: nze, ma,na,mb,nb  
  character(len=20)     :: name
  complex(psb_dpk_), allocatable :: temp(:)
  integer               :: info
  integer               :: err_act
  name='psb_numbmm'
  call psb_erractionsave(err_act)
  info = psb_success_


  ma = a%get_nrows()
  na = a%get_ncols()
  mb = b%get_nrows()
  nb = b%get_ncols()


  if ( mb /= na ) then 
    write(psb_err_unit,*) 'Mismatch in SYMBMM: ',ma,na,mb,nb
  endif
  allocate(temp(max(ma,na,mb,nb)),stat=info)    
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_Errpush(info,name)
    goto 9999
  endif

  !
  ! Note: we still have to test about possible performance hits. 
  !
  !
  call psb_ensure_size(size(c%ja),c%val,info)
  select type(a)
  type is (psb_z_csr_sparse_mat) 
    select type(b)
    type is (psb_z_csr_sparse_mat) 
      call csr_numbmm(a,b,c,temp,info)
    class default
      call gen_numbmm(a,b,c,temp,info)
    end select
  class default
    call gen_numbmm(a,b,c,temp,info)
  end select
  
  if (info /= psb_success_) then 
    call psb_errpush(info,name)
    goto 9999
  end if

  call c%set_asb()
  deallocate(temp) 

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
     call psb_error()
     return
  end if
  return

contains 

  subroutine csr_numbmm(a,b,c,temp,info)
    type(psb_z_csr_sparse_mat), intent(in)  :: a,b
    type(psb_z_csr_sparse_mat), intent(inout) :: c
    complex(psb_dpk_)                          :: temp(:)
    integer, intent(out)                    :: info
    integer               :: nze, ma,na,mb,nb

    info = psb_success_
    ma = a%get_nrows()
    na = a%get_ncols()
    mb = b%get_nrows()
    nb = b%get_ncols()
    
    call cnumbmm(ma,na,nb,a%irp,a%ja,0,a%val,&
         & b%irp,b%ja,0,b%val,&
         & c%irp,c%ja,0,c%val,temp)

    
  end subroutine csr_numbmm

  subroutine gen_numbmm(a,b,c,temp,info)
    class(psb_z_base_sparse_mat), intent(in)  :: a,b
    type(psb_z_csr_sparse_mat), intent(inout) :: c
    integer               :: info
    complex(psb_dpk_)      :: temp(:)
    integer, allocatable  :: iarw(:), iacl(:),ibrw(:),ibcl(:)
    complex(psb_dpk_), allocatable :: aval(:),bval(:)
    integer  :: maxlmn,i,j,m,n,k,l,nazr,nbzr,jj,minlm,minmn,minln
    complex(psb_dpk_)      :: ajj

    n = a%get_nrows()
    m = a%get_ncols() 
    l = b%get_ncols()
    maxlmn = max(l,m,n)
    allocate(iarw(maxlmn),iacl(maxlmn),ibrw(maxlmn),ibcl(maxlmn),&
         & aval(maxlmn),bval(maxlmn), stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      return
    endif

    do i = 1,maxlmn
      temp(i) = dzero
    end do
    minlm = min(l,m)
    minln = min(l,n)
    minmn = min(m,n)
    do  i = 1,n

      call a%csget(i,i,nazr,iarw,iacl,aval,info)
      do jj=1, nazr
        j=iacl(jj)
        ajj = aval(jj)
        if ((j<1).or.(j>m)) then 
          write(psb_err_unit,*) ' NUMBMM: Problem with A ',i,jj,j,m
            info = 1
            return
          
        endif
        call b%csget(j,j,nbzr,ibrw,ibcl,bval,info)
        do k=1,nbzr
          if ((ibcl(k)<1).or.(ibcl(k)>maxlmn)) then 
            write(psb_err_unit,*) 'Problem in NUMBM 1:',j,k,ibcl(k),maxlmn
            info = psb_err_pivot_too_small_
            return
          else
            temp(ibcl(k)) = temp(ibcl(k)) + ajj * bval(k)
          endif
        enddo
      end do
      do  j = c%irp(i),c%irp(i+1)-1
        if((c%ja(j)<1).or. (c%ja(j) > maxlmn))  then 
          write(psb_err_unit,*) ' NUMBMM: output problem',i,j,c%ja(j),maxlmn
            info = psb_err_invalid_ovr_num_
            return
        else
          c%val(j) = temp(c%ja(j))
          temp(c%ja(j)) = dzero
        endif
      end do
    end do

    
  end subroutine gen_numbmm

end subroutine psb_zbase_numbmm

