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
! Subroutine: psb_cd_renum_block
!   Produces a clone of a descriptor.
! 
! Arguments: 
!    desc_in  - type(psb_desc_type).         The communication descriptor to be cloned.
!    desc_out - type(psb_desc_type).         The output communication descriptor.
!    info     - integer.                       Return code.
subroutine psb_cd_renum_block(desc_in, desc_out, info)

  use psb_base_mod, psb_protect_name => psb_cd_renum_block

  implicit none
  !....parameters...

  type(psb_desc_type), intent(inout) :: desc_in
  type(psb_desc_type), intent(out)   :: desc_out
  integer(psb_ipk_), intent(out)     :: info

  !locals
  type(psb_gen_block_map), allocatable :: blck_map
  integer(psb_ipk_), allocatable :: lidx(:),gidx(:),reflidx(:),vnl(:)
  integer(psb_ipk_) :: i,n_row, n_col, n_glob_row, n_glob_col
  integer(psb_ipk_) :: np,me,ictxt, err_act
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name
  
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cd_renum_block'

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

  call desc_in%clone(desc_out,info)
  if (desc_in%is_ovl()) then
    write(0,*) 'Warning: descriptor with overlap, not going to clone into BLOCK'
  else
    !
    ! Ok, convert into a GEN BLOCK. 
    !
    allocate(psb_gen_block_map :: blck_map, stat=info) 
    if (info == 0) allocate(vnl(0:np),stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif
    n_row = desc_in%get_local_rows()
    n_col = desc_in%get_local_cols()
    n_glob_row = desc_in%get_global_rows()
    n_glob_col = desc_in%get_global_cols()
    vnl = 0
    vnl(me) = n_row
    call psb_sum(ictxt,vnl)
    vnl(1:np) = vnl(0:np-1)
    vnl(0) = 0
    do i=1,np
      vnl(i) = vnl(i-1)+vnl(i)
    end do
    allocate(lidx(n_col),reflidx(n_col),gidx(n_col),stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif
    !
    ! A GEN_BLOCK distribution would just assign indices 1:N_R0
    ! to process 0, N_R0+1:N_R0+N_R1 to process 1 and so on;
    ! once these are set you can get the others by invoking
    ! HALO on the old descriptor, then put the map into the new
    ! one. Halo lists are stored in local indices, so they will
    ! continue to work since local indices stay the same, it's
    ! only the global indices that were reshuffled.
    !
    reflidx(1:n_col) = [(i,i=1,n_col)]
    gidx(1:n_row) = reflidx(1:n_row) + vnl(me)
    call psb_halo(gidx,desc_in,info)
    if (info == 0) call blck_map%gen_block_map_init(ictxt,vnl(me),info)
    if (info == 0) call blck_map%g2l_ins(gidx,lidx,info,lidx=reflidx)
    if (info == 0) call blck_map%asb(info)
  
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    endif
    call move_alloc(blck_map,desc_out%indxmap)
  end if
  
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_cd_renum_block
