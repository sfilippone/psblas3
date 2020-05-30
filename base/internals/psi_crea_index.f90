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
! File: psi_i_crea_index.f90
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
subroutine psi_i_crea_index(desc_a,index_in,index_out,nxch,nsnd,nrcv,info)
  use psb_realloc_mod
  use psb_desc_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_timers_mod
  use psi_mod, psb_protect_name => psi_i_crea_index
  implicit none

  type(psb_desc_type), intent(in)     :: desc_a
  integer(psb_ipk_), intent(out)                :: info,nxch,nsnd,nrcv
  integer(psb_ipk_), intent(in)                 :: index_in(:)
  integer(psb_ipk_), allocatable, intent(inout) :: index_out(:)

  !         ....local scalars...      
  integer(psb_ipk_) :: ictxt, me, np, mode, err_act, dl_lda, ldl
  !         ...parameters...
  integer(psb_ipk_), allocatable :: dep_list(:,:), length_dl(:), loc_dl(:), c_dep_list(:), dl_ptr(:)
  integer(psb_ipk_) :: dlmax, dlavg
  integer(psb_ipk_),parameter    :: root=psb_root_,no_comm=-1
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name
  logical, parameter  :: do_timings=.false.
  integer(psb_ipk_), save  :: idx_phase1=-1, idx_phase2=-1, idx_phase3=-1
  integer(psb_ipk_), save  :: idx_phase11=-1, idx_phase12=-1, idx_phase13=-1

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
  if ((do_timings).and.(idx_phase1==-1))       &
       & idx_phase1 = psb_get_timer_idx("PSI_CREA_INDEX: phase1 ")
  if ((do_timings).and.(idx_phase2==-1))       &
       & idx_phase2 = psb_get_timer_idx("PSI_CREA_INDEX: phase2")
  if ((do_timings).and.(idx_phase3==-1))       &
       & idx_phase3 = psb_get_timer_idx("PSI_CREA_INDEX: phase3")
!!$  if ((do_timings).and.(idx_phase11==-1))       &
!!$       & idx_phase11 = psb_get_timer_idx("PSI_CREA_INDEX: phase11 ")
!!$  if ((do_timings).and.(idx_phase12==-1))       &
!!$       & idx_phase12 = psb_get_timer_idx("PSI_CREA_INDEX: phase12")
!!$  if ((do_timings).and.(idx_phase13==-1))       &
!!$       & idx_phase13 = psb_get_timer_idx("PSI_CREA_INDEX: phase13")


  ! ...extract dependence list (ordered list of identifer process
  !    which every process must communcate with...
  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': calling extract_dep_list'
  mode = 1
  if (.false.) then 
    if (do_timings) call psb_tic(idx_phase1)

    call psi_extract_dep_list(ictxt,&
         & desc_a%is_bld(), desc_a%is_upd(),&
         & index_in, dep_list,length_dl,dl_lda,mode,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='extrct_dl')
      goto 9999
    end if

    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),': from extract_dep_list',&
         &     me,length_dl(0),index_in(1), ':',dep_list(:length_dl(me),me)
    ! ...now process root contains dependence list of all processes...
    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),': root sorting dep list'
    if (do_timings) call psb_toc(idx_phase1)
    if (do_timings) call psb_tic(idx_phase2)

    call psi_dl_check(dep_list,dl_lda,np,length_dl)

    ! ....now i can sort dependency lists.
    call psi_sort_dl(dep_list,length_dl,np,info)
    if(info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_sort_dl')
      goto 9999
    end if
    if (do_timings) call psb_toc(idx_phase2)
    ldl = length_dl(me)
    loc_dl = dep_list(1:ldl,me)

  else

    if (do_timings) call psb_tic(idx_phase1)

    call psi_extract_loc_dl(ictxt,&
         & desc_a%is_bld(), desc_a%is_upd(),&
         & index_in, loc_dl,length_dl,info)

    dlmax = maxval(length_dl(:))
    dlavg = (sum(length_dl(:))+np-1)/np
!!$    if ((dlmax>0).and.(me==0)) write(0,*) 'Dependency list : max:',dlmax,&
!!$         & '  avg:',dlavg, choose_sorting(dlmax,dlavg,np)

    if (choose_sorting(dlmax,dlavg,np)) then 
      if (.true.) then 
        call psi_bld_glb_dep_list(ictxt,&
             & loc_dl,length_dl,dep_list,dl_lda,info)

        if (info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_,name,a_err='extrct_dl')
          goto 9999
        end if

        if (debug_level >= psb_debug_inner_) &
             & write(debug_unit,*) me,' ',trim(name),': from extract_dep_list',&
             &     me,length_dl(0),index_in(1), ':',dep_list(:length_dl(me),me)
        ! ...now process root contains dependence list of all processes...
        if (debug_level >= psb_debug_inner_) &
             & write(debug_unit,*) me,' ',trim(name),': root sorting dep list'
        if (do_timings) call psb_toc(idx_phase1)
        if (do_timings) call psb_tic(idx_phase2)

        !
        ! The dependency list has been symmetrized inside xtract_loc_dl
        !       
!!$       call psi_dl_check(dep_list,dl_lda,np,length_dl)
        
        ! ....now i can sort dependency lists.
        call psi_sort_dl(dep_list,length_dl,np,info)
        if(info /= psb_success_) then
          call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_sort_dl')
          goto 9999
        end if
        if (do_timings) call psb_toc(idx_phase2)
        ldl = length_dl(me)
        loc_dl = dep_list(1:ldl,me)
      else
        if (do_timings) call psb_toc(idx_phase1)
        if (do_timings) call psb_tic(idx_phase2)
        call psi_bld_glb_dep_list(ictxt,&
             & loc_dl,length_dl,c_dep_list,dl_ptr,info)

!!$        call psi_dl_check(dep_list,dl_lda,np,length_dl)
!!$
!!$        ! ....now i can sort dependency lists.
        call psi_sort_dl(dl_ptr,c_dep_list,length_dl,np,info)
!!$        if(info /= psb_success_) then
!!$          call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_sort_dl')
!!$          goto 9999
!!$        end if
        if (do_timings) call psb_toc(idx_phase2)
        
        
      end if
    else
      ! Do nothing
      ldl    = length_dl(me)
      loc_dl = loc_dl(1:ldl)      
    end if

  end if

  if (do_timings) call psb_tic(idx_phase3)
  if(debug_level >= psb_debug_inner_)&
       & write(debug_unit,*) me,' ',trim(name),': calling psi_desc_index',ldl,':',loc_dl(1:ldl)
  ! Do the actual format conversion. 
  call psi_desc_index(desc_a,index_in,loc_dl,ldl,nsnd,nrcv,index_out,info)
  if(debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': out of  psi_desc_index',&
       & size(index_out)
  nxch = ldl
  if(info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_desc_index')
    goto 9999
  end if
  if (do_timings) call psb_toc(idx_phase3)

  if (allocated(dep_list)) deallocate(dep_list,stat=info)
  if ((info==0).and.allocated(length_dl)) deallocate(length_dl,stat=info)
  if (info /= 0) then
    info = psb_err_alloc_dealloc_
    goto 9999
  end if
  if(debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': done'

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return
contains
  function choose_sorting(dlmax,dlavg,np) result(val)
    implicit none
    integer(psb_ipk_), intent(in) :: dlmax,dlavg,np
    logical                       :: val

    val = .not.(((dlmax>(26*4)).or.((dlavg>=(26*2)).and.(np>=128))))
  end function choose_sorting

end subroutine psi_i_crea_index
