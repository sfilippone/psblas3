!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
! File:  psb_zsymbmm.f90 
! Subroutine: 
! Arguments:
!
!
! Note: This subroutine performs the symbolic product of two sparse matrices.
!       It is modeled after the SMMP package by R. Bank and C. Douglas, but is 
!       rewritten in Fortran 95/2003 making use of our sparse matrix facilities.
!

subroutine psb_zsymbmm(a,b,c,info)
  use psb_base_mod, psb_protect_name => psb_zsymbmm
  implicit none 

  type(psb_zspmat_type), intent(in)    :: a,b
  type(psb_zspmat_type), intent(out)   :: c
  integer(psb_ipk_), intent(out)                  :: info
  type(psb_z_csr_sparse_mat), allocatable :: ccsr
  integer(psb_ipk_) :: err_act
  character(len=*), parameter ::  name='psb_symbmm'
  call psb_erractionsave(err_act)
  info = psb_success_

  if ((a%is_null()) .or.(b%is_null())) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(ccsr,stat=info)    

  if (info == psb_success_) then 
      call psb_symbmm(a%a,b%a,ccsr,info)
  else
    info = psb_err_alloc_dealloc_
  end if
  
  if (info /= psb_success_) then 
    call psb_errpush(info,name) 
    goto 9999
  end if
  call move_alloc(ccsr,c%a)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_zsymbmm

subroutine psb_zbase_symbmm(a,b,c,info)
  use psb_mat_mod
  use psb_serial_mod, psb_protect_name => psb_zbase_symbmm
  implicit none 

  class(psb_z_base_sparse_mat), intent(in) :: a,b
  type(psb_z_csr_sparse_mat), intent(out)  :: c
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_), allocatable  :: itemp(:)
  integer(psb_ipk_) :: nze, ma,na,mb,nb
  character(len=20)     :: name
  integer(psb_ipk_) :: err_act
  name='psb_symbmm'
  call psb_erractionsave(err_act)
  info = psb_success_

  ma = a%get_nrows()
  na = a%get_ncols()
  mb = b%get_nrows()
  nb = b%get_ncols()


  if ( mb /= na ) then 
    write(psb_err_unit,*) 'Mismatch in SYMBMM: ',ma,na,mb,nb
    info = psb_err_invalid_matrix_sizes_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(itemp(max(ma,na,mb,nb)),stat=info)    
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_Errpush(info,name)
    goto 9999
  endif
  !
  ! Note: we need to test whether there is a performance impact 
  !       in not using the original Douglas & Bank code. 
  !
  select type(a)
  type is (psb_z_csr_sparse_mat) 
    select type(b)
    type is (psb_z_csr_sparse_mat) 
      call csr_symbmm(a,b,c,itemp,info)
    class default
      call gen_symbmm(a,b,c,itemp,info)
    end select
  class default
    call gen_symbmm(a,b,c,itemp,info)
  end select

  if (info /= psb_success_) then 
    call psb_errpush(info,name) 
    goto 9999
  end if

  call psb_realloc(size(c%ja),c%val,info)
  deallocate(itemp) 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

contains

  subroutine csr_symbmm(a,b,c,itemp,info)
    type(psb_z_csr_sparse_mat), intent(in)  :: a,b
    type(psb_z_csr_sparse_mat), intent(out) :: c
    integer(psb_ipk_) :: itemp(:) 
    integer(psb_ipk_), intent(out)                    :: info
    interface 
      subroutine symbmm (n, m, l, ia, ja, diaga, &
           & ib, jb, diagb, ic, jc, diagc, index)
        import :: psb_ipk_
        integer(psb_ipk_) :: n,m,l,  ia(*), ja(*), diaga, ib(*), jb(*), diagb,&
             & diagc,  index(*)
        integer(psb_ipk_), allocatable :: ic(:),jc(:)
      end subroutine symbmm
    end interface
    integer(psb_ipk_) :: nze, ma,na,mb,nb

    info = psb_success_
    ma = a%get_nrows()
    na = a%get_ncols()
    mb = b%get_nrows()
    nb = b%get_ncols()
    
    nze = max(ma+1,2*ma)
    call c%allocate(ma,nb,nze)
    call symbmm(ma,na,nb,a%irp,a%ja,izero,&
         & b%irp,b%ja,izero,&
         & c%irp,c%ja,izero,itemp)
    
  end subroutine csr_symbmm
  subroutine gen_symbmm(a,b,c,index,info)
    class(psb_z_base_sparse_mat), intent(in)  :: a,b
    type(psb_z_csr_sparse_mat), intent(out) :: c
    integer(psb_ipk_) :: index(:),info
    integer(psb_ipk_), allocatable  :: iarw(:), iacl(:),ibrw(:),ibcl(:)
    integer(psb_ipk_) :: maxlmn,i,j,m,n,k,l,istart,length,nazr,nbzr,jj,minlm,minmn
    integer(psb_ipk_) :: nze, ma,na,mb,nb

    ma = a%get_nrows()
    na = a%get_ncols()
    mb = b%get_nrows()
    nb = b%get_ncols()

    nze = max(ma+1,2*ma)
    call c%allocate(ma,nb,nze)

    n = ma
    m = na 
    l = nb 
    maxlmn = max(l,m,n)

    allocate(iarw(maxlmn),iacl(maxlmn),ibrw(maxlmn),ibcl(maxlmn),&
         & stat=info)
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      return
    endif

    do i=1,maxlmn
      index(i)=0
    end do

    c%irp(1)=1
    minlm = min(l,m)
    minmn = min(m,n)

    main: do  i=1,n
      istart=-1
      length=0
      call a%csget(i,i,nazr,iarw,iacl,info)
      do jj=1, nazr

        j=iacl(jj)

        if ((j<1).or.(j>m)) then 
          write(psb_err_unit,*) ' SymbMM: Problem with A ',i,jj,j,m
          info = 1
          return
        endif
        call b%csget(j,j,nbzr,ibrw,ibcl,info)
        do k=1,nbzr
          if ((ibcl(k)<1).or.(ibcl(k)>maxlmn)) then 
            write(psb_err_unit,*) 'Problem in SYMBMM 1:',j,k,ibcl(k),maxlmn
            info=psb_err_pivot_too_small_
            return
          else
            if(index(ibcl(k)) == 0) then
              index(ibcl(k))=istart
              istart=ibcl(k)
              length=length+1
            endif
          endif
        end do
      end do

      c%irp(i+1)=c%irp(i)+length

      if (c%irp(i+1) > size(c%ja)) then 
        if (n > (2*i)) then 
          nze = max(c%irp(i+1), c%irp(i)*((n+i-1)/i))
        else
          nze = max(c%irp(i+1), nint((dble(c%irp(i))*(dble(n)/i)))   )
        endif
        call psb_realloc(nze,c%ja,info)
      end if
      do j= c%irp(i),c%irp(i+1)-1
        c%ja(j)=istart
        istart=index(istart)
        index(c%ja(j))=0
      end do
      call psb_msort(c%ja(c%irp(i):c%irp(i)+length-1))
      index(i) = 0
    end do main

  end subroutine gen_symbmm

end subroutine psb_zbase_symbmm
