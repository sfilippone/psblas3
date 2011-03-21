!!$  
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
! File: psb_icdasb.f90
!
! Subroutine: psb_icdasb
!   Assemble the psblas communications descriptor: inner part.
!   The user callable routine is defined in the psb_tools_mod module.
! 
! Arguments: 
!    desc_a  - type(psb_desc_type).    The communication descriptor.
!    info    - integer.                return code.
!    ext_hv  - logical                 Essentially this distinguishes a call 
!                                      coming from the build of an extended
!                                      halo descriptor with respect to a normal call. 
!
subroutine psb_icdasb(desc_a,info,ext_hv)
  use psb_base_mod, psb_protect_name => psb_icdasb
  use psi_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  !...Parameters....
  type(psb_desc_type), intent(inout) :: desc_a
  integer, intent(out)               :: info
  logical, intent(in), optional      :: ext_hv

  !....Locals....
  integer          ::  int_err(5)
  integer,allocatable ::  ovrlap_index(:),halo_index(:), ext_index(:)

  integer             ::  i,np,me, n_col, dectype, err_act, icomm
  integer             :: ictxt,n_row
  logical             :: ext_hv_
  integer             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  int_err(1) = 0
  name = 'psb_cdasb'

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt   = psb_cd_get_context(desc_a)
  dectype = psb_cd_get_dectype(desc_a)
  n_row   = psb_cd_get_local_rows(desc_a)
  n_col   = psb_cd_get_local_cols(desc_a)
  call psb_get_mpicomm(ictxt,icomm )

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_ok_desc(desc_a)) then 
    info = psb_err_spmat_invalid_state_
    int_err(1) = dectype
    call psb_errpush(info,name)
    goto 9999
  endif

  info = psb_get_errstatus()
  if (info /= psb_success_) then 
    ! Something went wrong in cdins/spins
    ! signal and exit
    info = psb_err_wrong_ins_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(ext_hv)) then 
    ext_hv_ = ext_hv
  else
    ext_hv_ = .false.
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit, *) me,' ',trim(name),': start'

  if (allocated(desc_a%indxmap)) then 
    call psi_ldsc_pre_halo(desc_a,ext_hv_,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='ldsc_pre_halo')
      goto 9999
    end if

    ! Take out the lists for ovrlap, halo and ext...
    call psb_move_alloc(desc_a%ovrlap_index,ovrlap_index,info)
    call psb_move_alloc(desc_a%halo_index,halo_index,info)
    call psb_move_alloc(desc_a%ext_index,ext_index,info)

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': Final conversion'
    ! Then convert and put them back where they belong.    
    call psi_cnv_dsc(halo_index,ovrlap_index,ext_index,desc_a,info) 

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_cnv_dsc')
      goto 9999
    end if

    deallocate(ovrlap_index, halo_index, ext_index, stat=info)
    if (info /= psb_success_) then
      info =psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if

    call desc_a%indxmap%asb(info)
    if (info /= psb_success_) then 
      write(0,*) 'Error from internal indxmap asb ',info
      info = psb_success_
    end if

!!$    desc_a%matrix_data(psb_n_row_) = desc_a%indxmap%get_lr()    
!!$    desc_a%matrix_data(psb_n_col_) = desc_a%indxmap%get_lc()
!!$    ! Ok, register into MATRIX_DATA 
!!$    desc_a%matrix_data(psb_dec_type_) = psb_desc_asb_

  else
    info = psb_err_spmat_invalid_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_icdasb
