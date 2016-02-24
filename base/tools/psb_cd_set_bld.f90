!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
subroutine psb_cd_set_ovl_bld(desc,info)
  use psb_base_mod, psb_protect_name => psb_cd_set_ovl_bld
  implicit none
  type(psb_desc_type), intent(inout) :: desc
  integer(psb_ipk_) :: info

  call psb_cd_set_bld(desc,info) 
  if (info == psb_success_) then 
    if (desc%indxmap%row_extendable()) then 
      call desc%indxmap%set_state(psb_desc_ovl_bld_)
!!$      desc%matrix_data(psb_dec_type_) = psb_desc_ovl_bld_ 
    else
      info = psb_err_invalid_cd_state_
    end if
  end if
    
end subroutine psb_cd_set_ovl_bld

subroutine psb_cd_set_bld(desc,info)
  use psb_base_mod, psb_protect_name => psb_cd_set_bld
  use psi_mod
  implicit none
  type(psb_desc_type), intent(inout) :: desc
  integer(psb_ipk_) :: info
  !locals
  integer(psb_ipk_) :: np,me,ictxt, err_act,idx,gidx,lidx,nc
  logical, parameter  :: debug=.false.,debugprt=.false.
  character(len=20)   :: name
  if (debug) write(psb_err_unit,*) me,'Entered CDCPY'
  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cd_set_bld'

  ictxt = desc%get_context()

  if (debug) write(psb_err_unit,*)'Entered CDSETBLD',ictxt
  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (debug) write(psb_err_unit,*) me,'Entered CDSETBLD'


  if (desc%is_asb())  call psb_cd_reinit(desc,info) 
  
  call desc%indxmap%set_state(psb_desc_bld_)

  if (debug) write(psb_err_unit,*) me,'SET_BLD: done'
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_cd_set_bld
