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
! File: psb_cdprt.f90
!
! Subroutine: psb_cdprt
!    Prints the descriptor to an output file
! 
! Arguments: 
!    iout          - integer.                The output unit to print to.
!    desc_p        - type(psb_desc_type).  The communication descriptor to be printed.
!    glob          - logical(otpional).      Wheter to print out global or local data.
!    short         - logical(optional).      Used to choose a verbose output.
!    verbosity     - integer(optional)       Choose a verbosity level 
!
subroutine psb_cdprt(iout,desc_p,glob,short, verbosity)
  use psb_base_mod, psb_protect_name => psb_cdprt
  implicit none 
  type(psb_desc_type), intent(in)  :: desc_p
  integer(psb_ipk_), intent(in)                :: iout
  integer(psb_ipk_), intent(in), optional     :: verbosity
  logical, intent(in), optional      :: glob,short
  logical :: short_, glob_

  integer(psb_ipk_) :: m, n_row, n_col,counter,idx,&
       & n_elem_recv,n_elem_send,proc,i, verb_
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: me, np
  integer(psb_ipk_)   :: total_snd, total_rcv, total_xhcg, global_halo, global_points
  integer(psb_ipk_)   :: local_snd, local_rcv, local_xhcg, local_halo, local_points
  integer(psb_ipk_)   :: max_xchg
  real(psb_dpk_)      :: av2s, v2s, anxchg
  
  if (present(glob)) then 
    glob_ = glob
  else
    glob_ = .false.
  endif
  if (present(short)) then 
    short_ = short
  else
    short_ = .true.
  endif
  if (present(verbosity)) then 
    verb_ = verbosity
  else
    verb_ = 1
  endif
  
  ctxt = desc_p%get_ctxt()
  call psb_info(ctxt, me,np)
  call psb_min(ctxt,verb_)
  
  !
  ! Level 1: Print global info
  !
  global_points = desc_p%get_global_rows()
  local_points  = desc_p%get_local_rows()
  local_halo    = desc_p%get_local_cols() - desc_p%get_local_rows()
  if (local_halo>0) then
    v2s = real(local_points,psb_dpk_)/real(local_halo,psb_dpk_)
  else
    v2s = 0.0_psb_dpk_
  end if

  global_halo   = local_halo
  av2s          = v2s
  call psb_sum(ctxt, global_halo)
  call psb_sum(ctxt, av2s)
  av2s = av2s / np
  max_xchg = get_nxchg(desc_p)
  anxchg   = max_xchg
  call psb_max(ctxt, max_xchg)
  call psb_sum(ctxt, anxchg)
  anxchg = anxchg/np
  
  if (me == psb_root_) then
    write(iout,*) ' Communication descriptor details '
    write(iout,*) '               Descriptor format:       ',desc_p%get_fmt()    
    write(iout,*) ' Global descriptor data: points:',global_points,' halo:',global_halo
    write(iout,*) ' Maximum number of exchanges     :',max_xchg
    write(iout,*) ' Average number of exchanges     :',anxchg
    write(iout,*)
    write(iout,*) ' Average volume to surface ratio :',av2s
    write(iout,*)
  end if
  flush(iout)
  call psb_barrier(ctxt)
  
  if (verb_ <= 1) return
  
  !
  ! Level 2: Statistics at process level
  !
    do i=0, np-1
    if (me == i)  then
      write(iout,*) me,': Local descriptor data: points:',local_points,&
           & ' halo:',local_halo
      if (local_halo>0) then 
        write(iout,*) me,': Volume to surface ratio:',real(local_points,psb_dpk_)/real(local_halo,psb_dpk_)
      else
        write(iout,*) me,': Volume to surface ratio:',0.0_psb_dpk_
      end if
    end if
    call psb_barrier(ctxt)
  end do

  if (verb_ <= 2) return

  if (me==psb_root_) write(iout,*) 'Communication data for : comm_halo'
  do i=0, np-1
    if (me == i)  &
       &  call  print_my_xchg(iout,desc_p,verbosity=verb_,data=psb_comm_halo_,glob=glob_)
    call psb_barrier(ctxt)
  end do

  if (verb_ <= 3) return
  
  if (me==psb_root_) write(iout,*) 'Communication data for : comm_ext'
  do i=0, np-1
    if (me == i)  &
       &  call  print_my_xchg(iout,desc_p,verbosity=verb_,data=psb_comm_ext_,glob=glob_)
    call psb_barrier(ctxt)
  end do
  
  return

contains

  function get_nxchg(desc,data) result(res)
    implicit none
    type(psb_desc_type), intent(in), target :: desc
    integer(psb_ipk_), intent(in), optional :: data
    integer(psb_ipk_) :: res
    
    class(psb_i_base_vect_type), pointer :: vpnt
    integer(psb_ipk_) :: ip, nerv, nesd, data_, totxch, idxr, idxs, info
    
    if (present(data)) then
      data_ = data
    else
      data_ = psb_comm_halo_
    end if
    call psb_cd_v_get_list(data_,desc_p,vpnt,totxch,idxr,idxs,info)
    res = totxch
  end function get_nxchg
  
    
  subroutine print_my_xchg(iout,desc_p,data,glob,short, verbosity)
    implicit none 
    type(psb_desc_type), intent(in), target   :: desc_p
    integer(psb_ipk_), intent(in)             :: iout
    integer(psb_ipk_), intent(in), optional   :: verbosity, data
    logical, intent(in), optional      :: glob,short
    logical :: short_, glob_
    
    integer(psb_ipk_)   :: ip, nerv, nesd, totxch,idxr,idxs
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)   :: me, np, data_, info, verb_
    integer(psb_lpk_), allocatable :: gidx(:)
    class(psb_i_base_vect_type), pointer :: vpnt
    
    ctxt = desc_p%get_ctxt()
    call psb_info(ctxt, me,np)
    if (present(data)) then
      data_ = data
    else
      data_ = psb_comm_halo_
    end if
    if (present(verbosity)) then
      verb_ = verbosity
    else
      verb_ = 1
    end if
    
    call psb_cd_v_get_list(data_,desc_p,vpnt,totxch,idxr,idxs,info)
    if (glob) &
         &   call psb_realloc(max(idxr,idxs,1),gidx,info)

    select case(verb_)
    case (1) 
      write(iout,*) me,': Total exchanges  :',totxch
      write(iout,*) me,':      Total sends    :',idxs
      write(iout,*) me,':      Total receives :', idxr
    case (2)
      write(iout,*) me,': Total exchanges  :',totxch
      write(iout,*) me,':      Total sends    :',idxs
      write(iout,*) me,':      Total receives :', idxr
      if (totxch == 0) return
      if (.not.associated(vpnt)) return
      if (.not.allocated(vpnt%v)) return
      associate(idx => vpnt%v)
        ip     = 1        
        do 
          if (ip > size(idx)) then 
            write(psb_err_unit,*) ': Warning: out of size of input vector '
            exit
          end if
          if (idx(ip) == -1) exit
          totxch = totxch+1
          nerv   = idx(ip+psb_n_elem_recv_)
          nesd   = idx(ip+nerv+psb_n_elem_send_)
          write(iout,*) '   ',me,': Exchanging with:',idx(ip),' Sends:',nesd,' Receives:', nerv
          idxs = idxs + nesd
          idxr = idxr + nerv
          ip   = ip+nerv+nesd+3
        end do
      end associate
    case (3)
      write(iout,*) me,': Total exchanges  :',totxch
      write(iout,*) me,':      Total sends    :',idxs
      write(iout,*) me,':      Total receives :', idxr
      if (totxch == 0) return
      if (.not.associated(vpnt)) return
      if (.not.allocated(vpnt%v)) return
      associate(idx => vpnt%v)
        ip     = 1        
        do 
          if (ip > size(idx)) then 
            write(psb_err_unit,*) ': Warning: out of size of input vector '
            exit
          end if
          if (idx(ip) == -1) exit
          totxch = totxch+1
          nerv   = idx(ip+psb_n_elem_recv_)
          nesd   = idx(ip+nerv+psb_n_elem_send_)
          write(iout,*) '   ',me,': Exchanging with:',idx(ip),' Sends:',nesd,' Receives:', nerv
          if (glob) then
            call desc_p%l2g(idx(ip+nerv+psb_n_elem_send_+1:ip+nerv+psb_n_elem_send_+nesd),gidx,info)
            write(iout,*) '   ',me,': sending to:',idx(ip),' :',gidx(1:nesd)
            call desc_p%l2g(idx(ip+psb_n_elem_recv_+1:ip+psb_n_elem_recv_+nerv),gidx,info)
            write(iout,*) '   ',me,': rcvng from:',idx(ip),' :',gidx(1:nerv)
          else
            write(iout,*) '   ',me,': sending to:',idx(ip),' :',&
                 & idx(ip+nerv+psb_n_elem_send_+1:ip+nerv+psb_n_elem_send_+nesd)
            write(iout,*) '   ',me,': rcvng from:',idx(ip),' :',&
                 & idx(ip+psb_n_elem_recv_+1:ip+psb_n_elem_recv_+nerv)
          end if
          idxs = idxs + nesd
          idxr = idxr + nerv
          ip   = ip+nerv+nesd+3
        end do
      end associate
    case default
      ! Do nothing
    end select
    flush(iout)
  end subroutine print_my_xchg

end subroutine psb_cdprt
