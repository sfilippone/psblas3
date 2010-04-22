!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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

function psb_c_linmap(map_kind,desc_X, desc_Y, map_X2Y, map_Y2X,iaggr,naggr) result(this)

  use psb_sparse_mod, psb_protect_name => psb_c_linmap

  implicit none 
  type(psb_clinmap_type)         :: this
  type(psb_desc_type), target       :: desc_X, desc_Y
  type(psb_c_sparse_mat), intent(in) :: map_X2Y, map_Y2X
  integer, intent(in)               :: map_kind
  integer, intent(in), optional     :: iaggr(:), naggr(:)
  !
  integer                      :: info
  character(len=20), parameter :: name='psb_linmap'

  info = psb_success_
  select case(map_kind) 
  case (psb_map_aggr_)
    ! OK    
    if (psb_is_ok_desc(desc_X)) then 
      this%p_desc_X=>desc_X
    else
      info = psb_err_pivot_too_small_
    endif
    if (psb_is_ok_desc(desc_Y)) then 
      this%p_desc_Y=>desc_Y
    else
      info = psb_err_invalid_ovr_num_
    endif
    if (present(iaggr)) then 
      if (.not.present(naggr)) then 
        info = 7
      else
        allocate(this%iaggr(size(iaggr)),&
             & this%naggr(size(naggr)), stat=info) 
        if (info == psb_success_) then 
          this%iaggr = iaggr
          this%naggr = naggr
        end if
      end if
    else
      allocate(this%iaggr(0), this%naggr(0), stat=info) 
    end if

  case(psb_map_gen_linear_)    
    
    if (psb_is_ok_desc(desc_X)) then 
      call psb_cdcpy(desc_X, this%desc_X,info) 
    else
      info = psb_err_pivot_too_small_
    endif
    if (psb_is_ok_desc(desc_Y)) then 
      call psb_cdcpy(desc_Y, this%desc_Y,info) 
    else
      info = psb_err_invalid_ovr_num_
    endif
    ! For a general linear map ignore iaggr,naggr
    allocate(this%iaggr(0), this%naggr(0), stat=info) 

  case default
    write(0,*) 'Bad map kind into psb_linmap ',map_kind
    info = 1
  end select

  if (info == psb_success_) call psb_clone(map_X2Y,this%map_X2Y,info)
  if (info == psb_success_) call psb_clone(map_Y2X,this%map_Y2X,info)
  if (info == psb_success_) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == psb_success_) then
    call psb_set_map_kind(map_kind, this)
  end if
  if (info /= psb_success_) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

end function psb_c_linmap

function psb_d_linmap(map_kind,desc_X, desc_Y, map_X2Y, map_Y2X,iaggr,naggr) result(this)

  use psb_sparse_mod, psb_protect_name => psb_d_linmap

  implicit none 
  type(psb_dlinmap_type)         :: this
  type(psb_desc_type), target       :: desc_X, desc_Y
  type(psb_d_sparse_mat), intent(in) :: map_X2Y, map_Y2X
  integer, intent(in)               :: map_kind
  integer, intent(in), optional     :: iaggr(:), naggr(:)
  !
  integer                      :: info
  character(len=20), parameter :: name='psb_linmap'
  logical, parameter           :: debug=.false.

  info = psb_success_
  select case(map_kind) 
  case (psb_map_aggr_)
    ! OK
    
    if (psb_is_ok_desc(desc_X)) then 
      this%p_desc_X=>desc_X
    else
      info = psb_err_pivot_too_small_
    endif
    if (psb_is_ok_desc(desc_Y)) then 
      this%p_desc_Y=>desc_Y
    else
      info = psb_err_invalid_ovr_num_
    endif
    if (present(iaggr)) then 
      if (.not.present(naggr)) then 
        info = 7
      else
        allocate(this%iaggr(size(iaggr)),&
             & this%naggr(size(naggr)), stat=info) 
        if (info == psb_success_) then 
          this%iaggr = iaggr
          this%naggr = naggr
        end if
      end if
    else
      allocate(this%iaggr(0), this%naggr(0), stat=info) 
    end if

  case(psb_map_gen_linear_)    
    
    if (psb_is_ok_desc(desc_X)) then 
      call psb_cdcpy(desc_X, this%desc_X,info) 
    else
      info = psb_err_pivot_too_small_
    endif
    if (psb_is_ok_desc(desc_Y)) then 
      call psb_cdcpy(desc_Y, this%desc_Y,info) 
    else
      info = psb_err_invalid_ovr_num_
    endif
    ! For a general linear map ignore iaggr,naggr
    allocate(this%iaggr(0), this%naggr(0), stat=info) 

  case default
    write(0,*) 'Bad map kind into psb_linmap ',map_kind
    info = 1
  end select

  if (info == psb_success_) call psb_clone(map_X2Y,this%map_X2Y,info)
  if (info == psb_success_) call psb_clone(map_Y2X,this%map_Y2X,info)
  if (info == psb_success_) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == psb_success_) then
    call psb_set_map_kind(map_kind, this)
  end if
  if (info /= psb_success_) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if
  if (debug) then 
!!$    write(0,*) trim(name),'  forward map:',allocated(this%map_X2Y%aspk)
!!$    write(0,*) trim(name),' backward map:',allocated(this%map_Y2X%aspk)
  end if

end function psb_d_linmap

function psb_s_linmap(map_kind,desc_X, desc_Y, map_X2Y, map_Y2X,iaggr,naggr) result(this)

  use psb_sparse_mod, psb_protect_name => psb_s_linmap

  implicit none 
  type(psb_slinmap_type)             :: this
  type(psb_desc_type), target        :: desc_X, desc_Y
  type(psb_s_sparse_mat), intent(in) :: map_X2Y, map_Y2X
  integer, intent(in)                :: map_kind
  integer, intent(in), optional      :: iaggr(:), naggr(:)
  !
  integer                      :: info
  character(len=20), parameter :: name='psb_linmap'

  info = psb_success_

  select case(map_kind) 
  case (psb_map_aggr_)
    ! OK
    
    if (psb_is_ok_desc(desc_X)) then 
      this%p_desc_X=>desc_X
    else
      info = psb_err_pivot_too_small_
    endif
    if (psb_is_ok_desc(desc_Y)) then 
      this%p_desc_Y=>desc_Y
    else
      info = psb_err_invalid_ovr_num_
    endif
    if (present(iaggr)) then 
      if (.not.present(naggr)) then 
        info = 7
      else
        allocate(this%iaggr(size(iaggr)),&
             & this%naggr(size(naggr)), stat=info) 
        if (info == psb_success_) then 
          this%iaggr = iaggr
          this%naggr = naggr
        end if
      end if
    else
      allocate(this%iaggr(0), this%naggr(0), stat=info) 
    end if

  case(psb_map_gen_linear_)    
    
    if (psb_is_ok_desc(desc_X)) then 
      call psb_cdcpy(desc_X, this%desc_X,info) 
    else
      info = psb_err_pivot_too_small_
    endif
    if (psb_is_ok_desc(desc_Y)) then 
      call psb_cdcpy(desc_Y, this%desc_Y,info) 
    else
      info = psb_err_invalid_ovr_num_
    endif
    ! For a general linear map ignore iaggr,naggr
    allocate(this%iaggr(0), this%naggr(0), stat=info) 

  case default
    write(0,*) 'Bad map kind into psb_linmap ',map_kind
    info = 1
  end select
  

  if (info == psb_success_) call psb_clone(map_X2Y,this%map_X2Y,info)
  if (info == psb_success_) call psb_clone(map_Y2X,this%map_Y2X,info)
  if (info == psb_success_) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == psb_success_) then
    call psb_set_map_kind(map_kind, this)
  end if
  if (info /= psb_success_) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

end function psb_s_linmap

function psb_z_linmap(map_kind,desc_X, desc_Y, map_X2Y, map_Y2X,iaggr,naggr) result(this)

  use psb_sparse_mod, psb_protect_name => psb_z_linmap

  implicit none 
  type(psb_zlinmap_type)         :: this
  type(psb_desc_type), target       :: desc_X, desc_Y
  type(psb_z_sparse_mat), intent(in) :: map_X2Y, map_Y2X
  integer, intent(in)               :: map_kind
  integer, intent(in), optional     :: iaggr(:), naggr(:)
  !
  integer                      :: info
  character(len=20), parameter :: name='psb_linmap'

  info = psb_success_
  select case(map_kind) 
  case (psb_map_aggr_)
    ! OK
    
    if (psb_is_ok_desc(desc_X)) then 
      this%p_desc_X=>desc_X
    else
      info = psb_err_pivot_too_small_
    endif
    if (psb_is_ok_desc(desc_Y)) then 
      this%p_desc_Y=>desc_Y
    else
      info = psb_err_invalid_ovr_num_
    endif
    if (present(iaggr)) then 
      if (.not.present(naggr)) then 
        info = 7
      else
        allocate(this%iaggr(size(iaggr)),&
             & this%naggr(size(naggr)), stat=info) 
        if (info == psb_success_) then 
          this%iaggr = iaggr
          this%naggr = naggr
        end if
      end if
    else
      allocate(this%iaggr(0), this%naggr(0), stat=info) 
    end if

  case(psb_map_gen_linear_)    
    
    if (psb_is_ok_desc(desc_X)) then 
      call psb_cdcpy(desc_X, this%desc_X,info) 
    else
      info = psb_err_pivot_too_small_
    endif
    if (psb_is_ok_desc(desc_Y)) then 
      call psb_cdcpy(desc_Y, this%desc_Y,info) 
    else
      info = psb_err_invalid_ovr_num_
    endif
    ! For a general linear map ignore iaggr,naggr
    allocate(this%iaggr(0), this%naggr(0), stat=info) 

  case default
    write(0,*) 'Bad map kind into psb_linmap ',map_kind
    info = 1
  end select

  if (info == psb_success_) call psb_clone(map_X2Y,this%map_X2Y,info)
  if (info == psb_success_) call psb_clone(map_Y2X,this%map_Y2X,info)
  if (info == psb_success_) call psb_realloc(psb_itd_data_size_,this%itd_data,info) 
  if (info == psb_success_) then
    call psb_set_map_kind(map_kind, this)
  end if
  if (info /= psb_success_) then
    write(0,*) trim(name),' Invalid descriptor input'
    return
  end if

end function psb_z_linmap
