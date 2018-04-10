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
! File: psb_dspasb.f90
!
! Subroutine: psb_dspasb
!    Assemble sparse matrix
!
! Arguments: 
!    a        - type(psb_dspmat_type).     The sparse matrix to be allocated.      
!    desc_a   - type(psb_desc_type).       The communication descriptor.
!    info     - integer.                     return code.
!    afmt     - character(optional)          The desired output storage format.
!    upd      - character(optional).         How will the matrix be updated? 
!                                            psb_upd_srch_    Simple strategy  
!                                            psb_upd_perm_    Permutation(more memory)
!    dupl     - integer(optional).           Duplicate coefficient handling:
!                                            psb_dupl_ovwrt_     overwrite
!                                            psb_dupl_add_       add 
!                                            psb_dupl_err_       raise an error. 
! 
!
subroutine psb_dspasb(a,desc_a, info, afmt, upd, dupl, mold)
  use psb_base_mod, psb_protect_name => psb_dspasb
  use psi_mod
  implicit none


  !...Parameters....
  type(psb_dspmat_type), intent (inout)  :: a
  type(psb_desc_type), intent(in)         :: desc_a
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_ipk_),optional, intent(in)            :: dupl, upd
  character(len=*), optional, intent(in)         :: afmt
  class(psb_d_base_sparse_mat), intent(in), optional :: mold
  !....Locals....
  integer(psb_ipk_) :: int_err(5)
  integer(psb_ipk_) :: np,me,n_col, err_act
  integer(psb_ipk_) :: spstate
  integer(psb_ipk_) :: ictxt,n_row
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)     :: name, ch_err

  info = psb_success_
  int_err(1)=0
  name = 'psb_spasb'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt    = desc_a%get_context()
  n_row    = desc_a%get_local_rows()
  n_col    = desc_a%get_local_cols()

  ! check on BLACS grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.desc_a%is_asb()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if


  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit, *) me,' ',trim(name),&
       & '   Begin matrix assembly...'

  !check on errors encountered in psdspins


  if (a%is_bld()) then 
    !
    ! First case: we come from a fresh build. 
    ! 

    n_row = desc_a%get_local_rows()
    n_col = desc_a%get_local_cols()
    call a%set_nrows(n_row)
    call a%set_ncols(n_col)
  end if

  if (a%is_bld()) then 
    call a%cscnv(info,type=afmt,dupl=dupl, mold=mold)
  else if (a%is_upd()) then 
    call a%asb(mold=mold)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
    
  end if

  
  IF (debug_level >= psb_debug_ext_) then 
    ch_err=a%get_fmt()
    write(debug_unit, *) me,' ',trim(name),':  From SPCNV',&
         & info,' ',ch_err
  end IF
  
  if (psb_errstatus_fatal()) then    
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='cscnv')
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_dspasb
