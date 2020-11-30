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
! Purpose: 
!  Provide a set of timers 
! 
!
module psb_timers_mod
  use psb_const_mod
  use psb_realloc_mod
  use psb_error_mod
  use psb_penv_mod
  
  public psb_init_timers, psb_get_timer_idx, psb_reset_timers,&
       & psb_tic, psb_toc, psb_print_timers, psb_get_timer, psb_free_timers
  private

  ! Reallocation
  integer(psb_ipk_), parameter :: ntchunk = 50
  integer(psb_ipk_), save      :: active_timers = 0
  ! Indices   
  integer(psb_ipk_), parameter :: timer_tic_ = 1, timer_toc_ = 2
  integer(psb_ipk_), parameter :: timer_x_ = 3, timer_sum_ = 4
  integer(psb_ipk_), parameter :: timer_max_ = 5
  integer(psb_ipk_), parameter :: timer_avg_ = 6, timer_min_ = 7
  integer(psb_ipk_), parameter :: timer_entries_ = 7

  ! The data itself 
  type psb_string_item
    character(len=40) :: data
  end type psb_string_item
  integer(psb_ipk_), allocatable  :: nsamples(:) 
  real(psb_dpk_), allocatable     :: timers(:,:)
  type(psb_string_item), allocatable :: timers_descr(:)
  logical                         :: wanted(timer_entries_)
  type(psb_string_item)           :: entries_descr(timer_entries_)
  save  :: nsamples, timers, timers_descr, wanted, entries_descr
  
  interface psb_realloc
    module procedure psb_string_item_realloc
  end interface psb_realloc
  
contains

  subroutine  print_timer(me, timer, timer_descr, iout)
    implicit none
    integer(psb_ipk_), intent(in)     :: me
    real(psb_dpk_), intent(in)        :: timer(timer_entries_) 
    type(psb_string_item), intent(in) :: timer_descr
    integer(psb_ipk_), optional       :: iout
    character(len=36) :: tmpname
    integer(psb_ipk_) :: iout_, i
    if (present(iout)) then
      iout_ = iout
    else
      iout_ = psb_out_unit
    end if

    write(tmpname,'(a)') trim(timer_descr%data)//":"
    
    write(iout_,'(a36,4(1x,a,f10.2))') tmpname,&
         & "Sum: ",timer(timer_sum_), &
         & "Avg: ",timer(timer_avg_), &
         & "Max: ",timer(timer_max_), &
         & "Min: ",-timer(timer_min_)
    
  end subroutine print_timer
  
  subroutine psb_print_timers(ctxt, idx, proc, global, iout)
    implicit none
    type(psb_ctxt_type), intent(in)         :: ctxt
    integer(psb_ipk_), intent(in), optional :: idx, proc, iout
    logical, optional :: global
    !
    !
    !
    integer(psb_ipk_) :: me,np,info,i,j, idxmin_, idxmax_, proc_
    real(psb_dpk_) :: gtimers(timer_entries_)
    real(psb_dpk_), allocatable :: ptimers(:,:) 
    logical :: global_
    
    call psb_info(ctxt,me,np)
    if (present(global)) then
      global_ = global
    else
      global_ = .true.
    end if
    if (present(proc)) then
      proc_ = proc
    else
      proc_ = -1
    end if
    if (present(idx)) then
      idxmin_ = idx
      idxmax_ = idx
    else
      idxmin_ = 1
      idxmax_ = active_timers
    end if
      
    if (global_) then
      if (allocated(timers)) then 
        allocate(ptimers(timer_entries_,size(timers,2)),stat=info)
        if (info /= 0) then
          write(0,*) 'Error while trying to allocate temporary ',info
          call psb_abort(ctxt)
        end if
        ptimers = timers
        call psb_max(ctxt,ptimers)
        if (me == psb_root_) then
          do i=idxmin_, idxmax_
            call print_timer(me, ptimers(:,i), timers_descr(i), iout)
          end do
        end if
      end if
    else
      if ((proc_ == -1).or.(me==proc_)) then
        do i=idxmin_, idxmax_
          call print_timer(me, ptimers(:,i), timers_descr(i), iout)
        end do
      end if
    end if
      
  end subroutine psb_print_timers
        
  subroutine psb_reset_timers()
    active_timers = 0
    if (allocated(nsamples)) nsamples = 0
    if (allocated(timers))  then
      timers   = dzero
      timers(timer_min_,:) = -huge(dzero)
    end if
    wanted        = .true.
    wanted(timer_tic_) = .false.
    wanted(timer_toc_) = .false.
    entries_descr(timer_tic_)%data = "tic"
    entries_descr(timer_toc_)%data = "toc"
    entries_descr(timer_x_)%data   = "Time"
    entries_descr(timer_sum_)%data = "Total time"
    entries_descr(timer_max_)%data = "Max time "
    entries_descr(timer_avg_)%data = "Average time"
    entries_descr(timer_min_)%data = "Min time"
    
  end subroutine psb_reset_timers
  
  subroutine psb_init_timers(ntimers)
    implicit none 
    integer(psb_ipk_), optional :: ntimers
    integer(psb_ipk_)  :: ntimers_, info
    
    if (present(ntimers)) then
      ntimers_ = ntimers
    else
      ntimers_ = ntchunk
    end if
    call reallocate_timers(ntimers_,info)
    if (info == 0) call psb_reset_timers()
    
  end subroutine psb_init_timers

  subroutine psb_free_timers()
    implicit none 

    integer(psb_ipk_)  :: info

    if (allocated(nsamples)) deallocate(nsamples,stat=info)
    if (allocated(timers)) deallocate(timers,stat=info)
    if (allocated(timers_descr)) deallocate(timers_descr,stat=info)
    
  end subroutine psb_free_timers

  subroutine reallocate_timers(tsz,info)
    implicit none 
    integer(psb_ipk_), intent(in)  :: tsz
    integer(psb_ipk_), intent(out) :: info

    call psb_realloc(timer_entries_,tsz,timers,info)
    if (info == 0) call psb_realloc(tsz,nsamples,info)
    if (info == 0) call psb_realloc(tsz,timers_descr,info)
    
  end subroutine reallocate_timers
  
  function psb_get_timer_idx(string) result(val)
    implicit none     
    integer(psb_ipk_) :: val
    character(len=*), intent(in), optional :: string
    integer(psb_ipk_) :: info
    
    val = -1
    if (.not.allocated(timers)) call psb_init_timers()
    if (active_timers >= size(timers,2)) then 
      call reallocate_timers((active_timers+ntchunk),info)
      if (info /= 0) return
      nsamples(active_timers+1:) = 0
      timers(:,active_timers+1:) = dzero
      timers(timer_min_,active_timers+1:) = -huge(dzero)
    end if
    active_timers = active_timers + 1
    if (present(string)) then
      timers_descr(active_timers)%data = string
    end if
    val          = active_timers
  end function psb_get_timer_idx

  function psb_get_timer(idx) result(val)
    implicit none     
    real(psb_dpk_) :: val
    integer(psb_ipk_), intent(in) :: idx
    !
    integer(psb_ipk_) :: info

    val = dzero
    if ((1<=idx).and.(idx <= active_timers)) then
      val = timers(timer_x_,idx)
    end if
  end function psb_get_timer

  subroutine psb_tic(idx)
    implicit none     
    integer(psb_ipk_), intent(in) :: idx
    
    if ((1<=idx).and.(idx <= active_timers)) &
         & timers(timer_tic_,idx) = psb_wtime()

  end subroutine psb_tic
  
  subroutine psb_toc(idx)
    implicit none     
    integer(psb_ipk_), intent(in) :: idx
    
    if ((1<=idx).and.(idx <= active_timers)) then 
      timers(timer_toc_,idx) = psb_wtime()
      timers(timer_x_,idx) = &
           & timers(timer_toc_,idx) - timers(timer_tic_,idx)
      nsamples(idx) = nsamples(idx) + 1
      timers(timer_sum_,idx) = &
           & timers(timer_sum_,idx) + timers(timer_x_,idx)
      timers(timer_avg_,idx) = timers(timer_sum_,idx) / nsamples(idx)
      timers(timer_max_,idx) = &
           & max(timers(timer_max_,idx), timers(timer_x_,idx))
      ! Trick: keep the MAX of negative times, so that
      ! a MAX over all processes for an entire section
      ! will give the MIN for the MIN entry. 
      timers(timer_min_,idx) = &
           & max(timers(timer_min_,idx), (-timers(timer_x_,idx)))
    end if
  end subroutine psb_toc
 
  Subroutine psb_string_item_realloc(len,rrax,info,lb)
    use psb_const_mod
    use psb_error_mod
    implicit none 
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                     :: len
    type(psb_string_item),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: lb

    ! ...Local Variables
    type(psb_string_item),allocatable  :: tmp(:)
    integer(psb_ipk_) :: dim,err_act,err, lb_, lbi,ub_
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_string_item_realloc'
    call psb_erractionsave(err_act)
    info=psb_success_ 

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
      goto 9999
    end if
    ub_ = lb_ + len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1)
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call move_alloc(tmp,rrax)
      End If
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
        goto 9999
      end if
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_string_item_realloc

  
end module psb_timers_mod

