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
! File: psb_dfree.f90
!
! Subroutine: psb_dfree
!    frees a dense matrix structure
! 
! Parameters: 
!    x        - real, pointer, dimension(:,:).    The dense matrix to be freed.
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code
subroutine psb_dfree(x, desc_a, info)
  !...free dense matrix structure...
  use psb_const_mod
  use psb_descriptor_type
  use psb_error_mod
  use psb_penv_mod
  implicit none

  !....parameters...
  real(kind(1.d0)),allocatable, intent(inout)  :: x(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer                     :: info

  !...locals....
  integer             :: ictxt,np,me,err_act
  character(len=20)   :: name


  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)
  name='psb_dfree'
  if (.not.allocated(desc_a%matrix_data)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)
  !     ....verify blacs grid correctness..
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.allocated(x)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

  !deallocate x
  deallocate(x,stat=info)
  if (info /= psb_no_err_) then
    info=4000
    call psb_errpush(info,name)
    goto 9999
  endif


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return

end subroutine psb_dfree



! Subroutine: psb_dfreev
!    frees a dense matrix structure
! 
! Parameters: 
!    x        - real, pointer, dimension(:).    The dense matrix to be freed.
!    desc_a   - type(<psb_desc_type>).          The communication descriptor.
!    info     - integer.                        Eventually returns an error code
subroutine psb_dfreev(x, desc_a, info)
  !...free dense matrix structure...
  use psb_const_mod
  use psb_descriptor_type
  use psb_error_mod
  use psb_penv_mod

  implicit none
  !....parameters...
  real(kind(1.d0)),allocatable, intent(inout)    :: x(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info 
  !...locals....
  integer             :: ictxt,np,me,err_act
  character(len=20)   :: name


  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)
  name='psb_dfreev'

  if (.not.allocated(desc_a%matrix_data)) then
     info=295
     call psb_errpush(info,name)
     return
  end if
  ictxt = psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.allocated(x)) then
     info=295
     call psb_errpush(info,name)
     goto 9999
  end if

  !deallocate x
  deallocate(x,stat=info)
  if (info /= psb_no_err_) then
    info=4000
    call psb_errpush(info,name)
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return

end subroutine psb_dfreev
