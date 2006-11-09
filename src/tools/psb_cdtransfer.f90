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
! File: psb_cdtransfer.f90
!
! Subroutine: psb_cdtransfer
!   Transfers data and allocation from in to out (just like MOVE_ALLOC). 
! 
! Parameters: 
!    desc_in  - type(<psb_desc_type>).         The communication descriptor to be transferred.
!    desc_out - type(<psb_desc_type>).         The output communication descriptor.
!    info     - integer.                       Eventually returns an error code.
subroutine psb_cdtransfer(desc_in, desc_out, info)

  use psb_descriptor_type
  use psb_serial_mod
  use psb_realloc_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod

  implicit none
  !....parameters...

  type(psb_desc_type), intent(inout)  :: desc_in
  type(psb_desc_type), intent(out)  :: desc_out
  integer, intent(out)                :: info

  !locals
  integer             :: np,me,ictxt, err_act
  logical, parameter  :: debug=.false.,debugprt=.false.
  character(len=20)   :: name

  if (debug) write(0,*) me,'Entered CDTRANSFER'
  if (psb_get_errstatus().ne.0) return 
  info = 0
  call psb_erractionsave(err_act)
  name = 'psb_cdtransfer'
  
  ictxt=desc_in%matrix_data(psb_ctxt_)

  call psb_info(ictxt, me, np)
  if (debug) write(0,*) me,'Entered CDTRANSFER'
  if (np == -1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  endif

!!$  call psb_nullify_desc(desc_out)
!!$  

  call psb_transfer( desc_in%matrix_data ,    desc_out%matrix_data  , info)
  call psb_transfer( desc_in%halo_index  ,    desc_out%halo_index   , info)
  call psb_transfer( desc_in%bnd_elem    ,    desc_out%bnd_elem     , info)
  call psb_transfer( desc_in%ovrlap_elem ,    desc_out%ovrlap_elem  , info)
  call psb_transfer( desc_in%ovrlap_index,    desc_out%ovrlap_index , info)
  call psb_transfer( desc_in%loc_to_glob ,    desc_out%loc_to_glob  , info)
  call psb_transfer( desc_in%glob_to_loc ,    desc_out%glob_to_loc  , info)
  call psb_transfer( desc_in%lprm        ,    desc_out%lprm         , info)
  call psb_transfer( desc_in%idx_space   ,    desc_out%idx_space    , info)


  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act == act_ret) then
     return
  else
     call psb_error(ictxt)
  end if
  return

end subroutine psb_cdtransfer
