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
! File: psb_zspfree.f90
!
! Subroutine: psb_zspfree
!    Frees a sparse matrix structure.
! 
! Parameters: 
!    a        - type(<psb_zspmat_type>).          The sparse matrix to be freed.      
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code.
!
subroutine psb_zspfree(a, desc_a,info)
  !...free sparse matrix structure...
  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(in)      :: desc_a
  type(psb_zspmat_type), intent(inout) :: a
  integer, intent(out)        :: info
  !...locals....
  integer             :: ictxt, err_act
  character(len=20)   :: name

  if(psb_get_errstatus().ne.0) return 
  info=0
  name = 'psb_zspfree'
  call psb_erractionsave(err_act)

  if (.not.allocated(desc_a%matrix_data)) then 
    info = 295
    call psb_errpush(info,name)
    return
  else
    ictxt = psb_cd_get_context(desc_a)
  end if

  !...deallocate a....
  call psb_sp_free(a,info)


  if(info /= 0) then
    info=2045
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_zspfree



