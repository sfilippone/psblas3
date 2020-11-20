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
! File: psb_icdasb.f90
!
! Subroutine: psb_icdasb
!   Assemble the psblas communications descriptor: inner part.
!   The user callable routine is defined in the psb_tools_mod module.
! 
! Arguments: 
!    desc  - type(psb_desc_type).    The communication descriptor.
!    info    - integer.                return code.
!    ext_hv  - logical                 Essentially this distinguishes a call 
!                                      coming from the build of an extended
!                                      halo descriptor with respect to a normal call. 
!
subroutine psb_icdasb(desc,info,ext_hv,mold)
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
  type(psb_desc_type), intent(inout) :: desc
  integer(psb_ipk_), intent(out)               :: info
  logical, intent(in), optional      :: ext_hv
  class(psb_i_base_vect_type), optional, intent(in) :: mold

  !....Locals....
  integer(psb_ipk_) ::  int_err(5)
  integer(psb_ipk_),allocatable ::  ovrlap_index(:),halo_index(:), ext_index(:)

  integer(psb_ipk_)  ::  i, n_col, dectype, err_act, n_row
  type(psb_ctxt_type) :: ctxt
  integer(psb_mpk_) ::  icomm
  integer(psb_ipk_) ::  np,me
  logical             :: ext_hv_
  logical, parameter  :: do_timings=.true.
  integer(psb_ipk_), save  :: idx_phase1=-1, idx_phase2=-1, idx_phase3=-1
  integer(psb_ipk_), save  :: idx_phase11=-1, idx_phase12=-1, idx_phase13=-1
  integer(psb_ipk_), save  :: idx_total=-1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
  int_err(1) = 0
  name = 'psb_cdasb'

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt   = desc%get_context()
  dectype = desc%get_dectype()
  n_row   = desc%get_local_rows()
  n_col   = desc%get_local_cols()
  icomm   = desc%get_mpic()
  if ((do_timings).and.(idx_total==-1))       &
       & idx_total = psb_get_timer_idx("ICDASB: total ")
  if ((do_timings).and.(idx_phase1==-1))       &
       & idx_phase1 = psb_get_timer_idx("ICDASB: phase1 ")
  if ((do_timings).and.(idx_phase2==-1))       &
       & idx_phase2 = psb_get_timer_idx("ICDASB: phase2")
  if ((do_timings).and.(idx_phase3==-1))       &
       & idx_phase3 = psb_get_timer_idx("ICDASB: phase3")
!!$  if ((do_timings).and.(idx_phase11==-1))       &
!!$       & idx_phase11 = psb_get_timer_idx("ICDASB: phase11 ")
!!$  if ((do_timings).and.(idx_phase12==-1))       &
!!$       & idx_phase12 = psb_get_timer_idx("ICDASB: phase12")
!!$  if ((do_timings).and.(idx_phase13==-1))       &
!!$       & idx_phase13 = psb_get_timer_idx("ICDASB: phase13")

  call psb_tic(idx_total)
  ! check on blacs grid 
  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.desc%is_ok()) then 
    info = psb_err_invalid_cd_state_
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

  if (allocated(desc%indxmap)) then 
    if (do_timings) call psb_tic(idx_phase1)    
    if (.not.ext_hv_) then 
      call psi_bld_tmphalo(desc,info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='bld_tmphalo')
        goto 9999
      end if
    end if
    if (do_timings) call psb_toc(idx_phase1)
    if (do_timings) call psb_tic(idx_phase2)    
    ! Take out the lists for ovrlap, halo and ext...
    call psb_move_alloc(desc%ovrlap_index,ovrlap_index,info)
    call psb_move_alloc(desc%halo_index,halo_index,info)
    call psb_move_alloc(desc%ext_index,ext_index,info)

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': Final conversion'
    ! Then convert and put them back where they belong.    
    call psi_cnv_dsc(halo_index,ovrlap_index,ext_index,desc,info,mold=mold) 

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
    if (do_timings) call psb_toc(idx_phase2)
    if (do_timings) call psb_tic(idx_phase3)    

    call desc%indxmap%asb(info)
    if (info == psb_success_) then 
      if (allocated(desc%indxmap%tempvg)) &
           & deallocate(desc%indxmap%tempvg,stat=info)
    end if
    if (info /= psb_success_) then 
      write(0,*) 'Error from internal indxmap asb ',info
      info = psb_success_
    end if
    if (do_timings) call psb_toc(idx_phase3)    
  else
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_toc(idx_total)  
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_icdasb
