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
subroutine psb_zprecinit(ictxt,p,ptype,info)

  use psb_base_mod
  use psb_z_prec_type, psb_protect_name => psb_zprecinit
  use psb_z_nullprec, only : psb_z_null_prec_type
  use psb_z_diagprec, only : psb_z_diag_prec_type
  use psb_z_bjacprec, only : psb_z_bjac_prec_type
  implicit none
  integer(psb_ipk_), intent(in)          :: ictxt
  class(psb_zprec_type), intent(inout) :: p
  character(len=*), intent(in)           :: ptype
  integer(psb_ipk_), intent(out)         :: info

  info = psb_success_

  if (allocated(p%prec) ) then
    call p%prec%free(info)
    if (info == psb_success_) deallocate(p%prec,stat=info) 
    if (info /= psb_success_) return
  end if
  
  p%ictxt = ictxt
  
  select case(psb_toupper(ptype(1:len_trim(ptype))))
  case ('NONE','NOPREC') 
    
    allocate(psb_z_null_prec_type :: p%prec, stat=info)       
    
  case ('DIAG')
    allocate(psb_z_diag_prec_type :: p%prec, stat=info)       
    
  case ('BJAC') 
    allocate(psb_z_bjac_prec_type :: p%prec, stat=info)       
    
  case default
    write(psb_err_unit,*) 'Unknown preconditioner type request "',ptype,'"'
    info = psb_err_invalid_preca_
    
  end select
  if (info == psb_success_)  call p%prec%precinit(info)
  
end subroutine psb_zprecinit
