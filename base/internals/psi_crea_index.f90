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
!
! File: psi_crea_index.f90
!
! Subroutine: psb_crea_index
!    Converts a list of data exchanges from build format to assembled format. 
!    See psi_desc_index for a description of the formats. 
!    Works by first finding a suitable ordering for the data exchanges, 
!    then doing the actual conversion. 
!
! Arguments:
! desc_a       - type(psb_desc_type)   The descriptor; in this context only the index 
!                                       mapping parts are used.
! index_in(:)  - integer               The index list, build format  
! index_out(:) - integer(psb_ipk_), allocatable  The index list, assembled format
! nxch         - integer               The number of data exchanges on the calling process
! nsnd         - integer               Total send buffer size       on the calling process
! nrcv         - integer               Total receive buffer size    on the calling process
!
!  
subroutine psi_crea_index(desc_a,index_in,index_out,nxch,nsnd,nrcv,info)
  use psb_realloc_mod
  use psb_desc_mod
  use psb_error_mod
  use psb_penv_mod
  use psi_mod, psb_protect_name => psi_crea_index
  implicit none

  type(psb_desc_type), intent(in)     :: desc_a
  integer(psb_ipk_), intent(out)                :: info,nxch,nsnd,nrcv
  integer(psb_ipk_), intent(in)                 :: index_in(:)
  integer(psb_ipk_), allocatable, intent(inout) :: index_out(:)

  !         ....local scalars...      
  integer(psb_ipk_) :: ictxt, me, np, mode, err_act, dl_lda
  !         ...parameters...
  integer(psb_ipk_), allocatable :: dep_list(:,:), length_dl(:)
  integer(psb_ipk_),parameter    :: root=psb_root_,no_comm=-1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name

  info = psb_success_
  name='psi_crea_index'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_ctxt()

  call psb_info(ictxt,me,np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! ...extract dependence list (ordered list of identifer process
  !    which every process must communcate with...
  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': calling extract_dep_list'
  mode = 1

  call psi_extract_dep_list(ictxt,&
       & desc_a%is_bld(), desc_a%is_upd(),&
       & index_in, dep_list,length_dl,dl_lda,mode,info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='extrct_dl')
    goto 9999
  end if

  ! ...now process root contains dependence list of all processes...
  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': root sorting dep list'

  call psi_dl_check(dep_list,dl_lda,np,length_dl)

  ! ....now i can sort dependency lists.
  call psi_sort_dl(dep_list,length_dl,np,info)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_sort_dl')
    goto 9999
  end if
  
  if(debug_level >= psb_debug_inner_)&
       & write(debug_unit,*) me,' ',trim(name),': calling psi_desc_index'
  ! Do the actual format conversion. 
  call psi_desc_index(desc_a,index_in,dep_list(1:,me),&
       & length_dl(me),nsnd,nrcv, index_out,info)
  if(debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': out of  psi_desc_index',&
       & size(index_out)
  nxch = length_dl(me)
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_desc_index')
    goto 9999
  end if

  deallocate(dep_list,length_dl)
  if(debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
end subroutine psi_crea_index
