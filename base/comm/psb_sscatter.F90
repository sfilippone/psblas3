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
! File:  psb_sscatter.f90
!
! Subroutine: psb_sscatter_vect
!   This subroutine scatters a global vector locally owned by one process
!   into pieces that are local to all the processes.
!
! Arguments:
!   globx     -  real,dimension(:)          The global matrix to scatter.
!   locx      -  type(psb_s_vect_type)      The local piece of the distributed matrix.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Error code.
!   iroot     -  integer(optional).            The process that owns the global matrix. 
!                                              If -1 all the processes have a copy. 
!                                              Default -1
subroutine  psb_sscatter_vect(globx, locx, desc_a, info, root, mold)
  use psb_base_mod, psb_protect_name => psb_sscatter_vect
  implicit none
  type(psb_s_vect_type), intent(inout) :: locx
  real(psb_spk_), intent(in)     :: globx(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_), intent(in), optional    :: root
  class(psb_s_base_vect_type), intent(in), optional :: mold
  
  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_mpk_) :: np, me, icomm, myrank, rootrank
  integer(psb_ipk_) :: ierr(5), err_act, m, n, i, j, idx, nrow, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, k, pos, ilx, jlx
  real(psb_spk_), allocatable  :: vlocx(:)
  character(len=20)        :: name, ch_err
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_scatter_vect'
  info=psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  ctxt=desc_a%get_context()
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()


  ! check on blacs grid 
  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  
  if (info == psb_success_) call psb_scatter(globx, vlocx, desc_a, info, root=root)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_scatterv')
    goto 9999
  endif
  
  call locx%bld(vlocx,mold=mold)
  
  call psb_erractionrestore(err_act)
  return  
  
9999 call psb_error_handler(ctxt,err_act)
  
  return
  
end subroutine psb_sscatter_vect
