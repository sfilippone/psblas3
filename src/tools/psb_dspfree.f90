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
! File: psb_dspfree.f90
!
! Subroutine: psb_dspfree
!    Frees a sparse matrix structure.
! 
! Parameters: 
!    a        - type(<psb_dspmat_type>).          The sparse matrix to be freed.      
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code.
!
subroutine psb_dspfree(a, desc_a,info)
  !...free sparse matrix structure...
  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(in) :: desc_a
  type(psb_dspmat_type), intent(inout)       ::a
  integer, intent(out)        :: info
  !...locals....
  integer             :: int_err(5)
  integer             :: temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol,err, err_act
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name = 'psb_dspfree'
  call psb_erractionsave(err_act)

  if (.not.associated(desc_a%matrix_data)) then 
     info=295
     call psb_errpush(info,name)
     return
  else
     icontxt=desc_a%matrix_data(psb_ctxt_)
  end if

  !...deallocate a....

  if ((info.eq.0).and.(.not.associated(a%pr))) info=2951
  if (info.eq.0) then
     !deallocate pr field
     deallocate(a%pr,stat=info)
     if (info.ne.0) info=2045
  end if
  if ((info.eq.0).and.(.not.associated(a%pl))) info=2952
  !deallocate pl  field
  if (info.eq.0) then 
     deallocate(a%pl,stat=info)
     if (info.ne.0) info=2046
  end if
  if ((info.eq.0).and.(.not.associated(a%ia2))) info=2953
  if (info.eq.0) then
     !deallocate ia2 field
     deallocate(a%ia2,stat=info)
     if (info.ne.0) info=2047
  end if
  if ((info.eq.0).and.(.not.associated(a%ia1))) info=2954
  if (info.eq.0) then
     !deallocate ia1  field
     deallocate(a%ia1,stat=info)
     if (info.ne.0) info=2048
  endif
  if ((info.eq.0).and.(.not.associated(a%aspk))) info=2955
  if (info.eq.0) then
     !deallocate aspk field
     deallocate(a%aspk,stat=info)
     if (info.ne.0) info=2049
  endif
  if (info.eq.0) call psb_nullify_sp(a)

  if(info.ne.0) then
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dspfree



subroutine psb_dspfrees(a, info)
  !...free sparse matrix structure...
  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_dspmat_type), intent(inout)       ::a
  integer, intent(out)        :: info
  !...locals....
  integer             :: int_err(5)
  integer             :: temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer             :: icontxt,nprow,npcol,me,mypcol,err, err_act
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name = 'psb_dspfrees'
  call psb_erractionsave(err_act)

  !...deallocate a....

!  if ((info.eq.0).and.(.not.associated(a%pr))) info=2951
  if ((info.eq.0).and.(associated(a%pr))) then
     !deallocate pr field
     deallocate(a%pr,stat=info)
     if (info.ne.0) info=2045
  end if
!  if ((info.eq.0).and.(.not.associated(a%pl))) info=2952
  !deallocate pl  field
  if ((info.eq.0).and.(associated(a%pl))) then 
     deallocate(a%pl,stat=info)
     if (info.ne.0) info=2046
  end if
!  if ((info.eq.0).and.(.not.associated(a%ia2))) info=2953
  if ((info.eq.0).and.(associated(a%ia2))) then
     !deallocate ia2 field
     deallocate(a%ia2,stat=info)
     if (info.ne.0) info=2047
  end if
!  if ((info.eq.0).and.(.not.associated(a%ia1))) info=2954
  if ((info.eq.0).and.(associated(a%ia1))) then
     !deallocate ia1  field
     deallocate(a%ia1,stat=info)
     if (info.ne.0) info=2048
  endif
!  if ((info.eq.0).and.(.not.associated(a%aspk))) info=2955
  if ((info.eq.0).and.(associated(a%aspk))) then
     !deallocate aspk field
     deallocate(a%aspk,stat=info)
     if (info.ne.0) info=2049
  endif
  if (info.eq.0) call psb_nullify_sp(a)

  if(info.ne.0) then
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return

end subroutine psb_dspfrees
