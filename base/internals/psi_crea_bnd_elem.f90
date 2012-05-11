!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
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
!
! File: psi_crea_bnd_elem.f90
!
! Subroutine: psi_crea_bnd_elem
!    Extracts a list of boundary indices. If no boundary is present in 
!    the distribution the output vector is put in the unallocated state,
!    otherwise its size is equal to the number of boundary indices on the 
!    current (calling) process. 
! 
! Arguments: 
!    bndel(:) - integer(psb_ipk_), allocatable      Array containing the output list              
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
! 
subroutine psi_crea_bnd_elem(bndel,desc_a,info)
  use psi_mod, psb_protect_name => psi_crea_bnd_elem
  use psb_realloc_mod
  use psb_descriptor_type
  use psb_error_mod
  use psb_serial_mod
  implicit none
  
  integer(psb_ipk_), allocatable :: bndel(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out) :: info

  integer(psb_ipk_), allocatable :: work(:)
  integer(psb_ipk_) :: i, j, nr, ns, k, err_act
  character(len=20)    :: name

  info = psb_success_
  name='psi_crea_bnd_elem'
  call psb_erractionsave(err_act)

  allocate(work(size(desc_a%halo_index)),stat=info)
  if (info /= psb_success_ ) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  i=0
  j=1
  do while(desc_a%halo_index(j) /= -1) 

    nr = desc_a%halo_index(j+1)
    ns = desc_a%halo_index(j+1+nr+1)
    do k=1, ns
      i = i + 1
      work(i) = desc_a%halo_index(j+1+nr+1+k)
    enddo
    j  = j + 1 + ns + 1 + nr + 1
  enddo

  call psb_msort_unique(work(1:i),j)

  if (.true.) then 
    if (j>=0) then 
      call psb_realloc(j,bndel,info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999      
      end if
      bndel(1:j) = work(1:j)
    else
      if (allocated(bndel)) then 
        deallocate(bndel)
      end if
    end if
  else
    call psb_realloc(j+1,bndel,info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if
    bndel(1:j) = work(1:j)
    bndel(j+1) = -1
  endif

  deallocate(work)
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
end subroutine psi_crea_bnd_elem
