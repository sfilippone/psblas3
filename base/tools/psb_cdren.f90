!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006, 2010, 2015, 2017
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: psb_cdren.f90
!  
!  WARNING: this routine is almost certainly obsolete. Must be reviewed.
! 
! Subroutine: psb_cdren
!    Updates a communication descriptor according to a renumbering scheme.
! 
! Arguments: 
!    trans    - character.                     Whether iperm or its transpose 
!                                              should be applied.
!    iperm    - integer(psb_ipk_),dimension(:).          The renumbering scheme.
!    desc_a   - type(psb_desc_type).         The communication descriptor
!                                              to be updated.
!    info     - integer.                       Return code
!
subroutine psb_cdren(trans,iperm,desc_a,info)
  use psb_base_mod, psb_protect_name => psb_cdren
  use psi_mod
  implicit none


  !...parameters....
  type(psb_desc_type), intent(inout)  :: desc_a
  integer(psb_ipk_), intent(inout)              :: iperm(:)
  character, intent(in)               :: trans
  integer(psb_ipk_), intent(out)                :: info
  !....locals....
  integer(psb_ipk_) :: i,j,np,me, n_col, kh, nh
  integer(psb_ipk_) :: dectype
  integer(psb_ipk_) :: ictxt,n_row, int_err(5), err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cdren'
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt   = desc_a%get_context()
  dectype = desc_a%get_dectype()
  n_row   = desc_a%get_local_rows()
  n_col   = desc_a%get_local_cols()

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_asb_desc(desc_a)) then 
    info = psb_err_invalid_cd_state_
    int_err(1) = dectype
    call psb_errpush(info,name,int_err)
    goto 9999
  endif

  if (iperm(1) /= 0) then 
    if (.not.psb_isaperm(n_row,iperm)) then
      info = 610
      int_err(1) = iperm(1)
      call psb_errpush(info,name,int_err)
      goto 9999
    endif
  endif


  !check on errors encountered in psdspins

  if ((iperm(1) /= 0))   then 

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': here we go with ',iperm(1) 
    call psb_ensure_size(n_col,desc_a%lprm,info)
    if (psb_toupper(trans) == 'N') then 
      do i=1, n_row
        desc_a%lprm(iperm(i)) = i
      enddo
      do i=n_row+1,n_col
        desc_a%lprm(i) = i
      enddo
    else if (psb_toupper(trans) == 'T') then 
      do i=1, n_row
        desc_a%lprm(i) = iperm(i)
      enddo
      do i=n_row+1,n_col
        desc_a%lprm(i) = i
      enddo
    endif
    ! crossed fingers.....
    ! fix glob_to_loc/loc_to_glob  mappings, then indices lists
    ! hmm, maybe we should just move all of this onto a different level,
    ! have a specialized subroutine, and do it in the solver context???? 
    if (allocated(desc_a%halo_index)) &
         & call psi_renum_index(desc_a%lprm,desc_a%halo_index,info)
    if (allocated(desc_a%ovrlap_index)) &
         & call psi_renum_index(desc_a%lprm,desc_a%ovrlap_index,info)
    if (allocated(desc_a%ovr_mst_idx)) &
         & call psi_renum_index(desc_a%lprm,desc_a%ovr_mst_idx,info)
    if (allocated(desc_a%ext_index)) &
         & call psi_renum_index(desc_a%lprm,desc_a%ext_index,info)
          
    do i=1, size(desc_a%ovrlap_elem,1)
      desc_a%ovrlap_elem(i,1) = desc_a%lprm(desc_a%ovrlap_elem(i,1))
    end do
    do i=1, size(desc_a%bnd_elem)
      desc_a%bnd_elem(i) = desc_a%lprm(desc_a%bnd_elem(i))
    end do

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': done renumbering'
  else 
    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': nothing to be done'
  endif
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_cdren
