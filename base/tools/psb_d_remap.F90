!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
!
! Subroutine: psb_cdcpy
!   Produces a clone of a descriptor.
! 
! Arguments: 
!    desc_in  - type(psb_desc_type).       The communication descriptor to be cloned.
!    desc_out - type(psb_desc_type).       The output communication descriptor.
!    info     - integer.                       Return code.
subroutine psb_d_remap(np_remap, desc_in, a_in, desc_out, a_out, info)

  use psb_base_mod, psb_protect_name => psb_d_remap

  implicit none
  !....parameters...
  integer(psb_ipk_), intent(in)        :: np_remap
  type(psb_desc_type), intent(inout)   :: desc_in
  type(psb_dspmat_type), intent(inout) :: a_in
  type(psb_dspmat_type), intent(out)   :: a_out
  type(psb_desc_type), intent(out)     :: desc_out
  integer(psb_ipk_), intent(out)       :: info


  !locals
  integer(psb_ipk_) :: np, me, ictxt, err_act
  integer(psb_ipk_) :: newctxt, rnp, rme
  integer(psb_ipk_) :: ipdest, id1, id2, imd, i
  integer(psb_ipk_), allocatable :: newnl(:)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cdcpy'

  ictxt = desc_in%get_context()

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': Entered'
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  write(0,*) ' Remapping from ',np,' onto ', np_remap
 
  if (desc_in%get_fmt() == 'BLOCK') then
    ! OK
    call psb_init(newctxt,np=np_remap,basectxt=ictxt)
    call psb_info(newctxt,rme,rnp)
    write(0,*) 'Old context: ',me,np,' New context: ',rme,rnp
    call psb_bcast(ictxt,rnp)
    allocate(newnl(rnp),stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif    
    if (rnp >= np) then 
      write(0,*) ' No remapping on larger proc count now'
      info = psb_err_internal_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    id2 = np/rnp
    id1 = id2+1
    imd = mod(np,rnp)
    if (me < (imd*id1)) then
      ipdest = (me/id1)
    else
      ipdest = ( ((me-imd*id1)/id2) +  imd)
    end if

    write(0,*) ' Sending my data from ',me,' to ', &
         & ipdest, 'out of ',rnp,rnp-1
    newnl = 0
    newnl(ipdest+1) = desc_in%get_local_rows()
    call psb_sum(ictxt,newnl)
    if (rme>=0) then
      call psb_cdall(newctxt,desc_out,info,nl=newnl(rme+1))
      call psb_cdasb(desc_out,info)
      write(0,*) me,rme,'In ',desc_in%get_local_rows(),desc_in%get_global_rows(),&
           & ' out ',desc_out%get_local_rows(),desc_out%get_global_rows()
    else 
      write(0,*) me,rme,'In ',desc_in%get_local_rows(),desc_in%get_global_rows(),&
           & ' out ',0,0
    end if
    call psb_exit(newctxt,close=.false.)
  else
    write(0,*) 'Right now only BLOCK on input '
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  endif
    
  !call psb_cdall()
  
  ! For the time being cleanup
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return  
  
end subroutine psb_d_remap
