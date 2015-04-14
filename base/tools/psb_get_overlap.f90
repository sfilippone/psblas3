!!$ 
!!$              Parallel Sparse BLAS  version 3.4
!!$    (C) Copyright 2006, 2010, 2015
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
! File: psb_get_overlap.f90
!
! Subroutine: psb_get_overlap
!    Extracts a list of overlap indices. If no overlap is present in 
!    the distribution the output vector is put in the unallocated state,
!    otherwise its size is equal to the number of overlap indices on the 
!    current (calling) process. 
! 
! Arguments: 
!    ovrel(:) - integer(psb_ipk_), allocatable      Array containing the output list              
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                  return code.
!
subroutine psb_get_ovrlap(ovrel,desc,info)
  use psb_base_mod, psb_protect_name => psb_get_ovrlap
  implicit none 
  integer(psb_ipk_), allocatable, intent(out) :: ovrel(:)
  type(psb_desc_type), intent(in) :: desc
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: i,j, err_act
  character(len=20)    :: name

  info = psb_success_
  name='psi_get_overlap'
  call psb_erractionsave(err_act)

  if (.not.psb_is_asb_desc(desc)) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (allocated(desc%ovrlap_elem) ) then 

    i=size(desc%ovrlap_elem,1) 

    allocate(ovrel(i),stat=info)
    if (info /= psb_success_ ) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if

    do j=1,i
      ovrel(j) = desc%ovrlap_elem(j,1)       
    end do

  else

    if (allocated(ovrel)) then 
      deallocate(ovrel,stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Deallocate')
        goto 9999      
      end if
    end if

  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_get_ovrlap
