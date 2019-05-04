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
!
!
! Takes a vector x from space map%p_desc_X and maps it onto
! map%p_desc_Y under map%map_X2Y possibly with communication
! due to exch_fw_idx
!
subroutine psb_c_map_X2Y(alpha,x,beta,y,map,info,work)
  use psb_base_mod, psb_protect_name => psb_c_map_X2Y

  implicit none 
  type(psb_clinmap_type), intent(in) :: map
  complex(psb_spk_), intent(in)       :: alpha,beta
  complex(psb_spk_), intent(inout)    :: x(:)
  complex(psb_spk_), intent(out)      :: y(:)
  integer(psb_ipk_), intent(out)                  :: info 
  complex(psb_spk_), optional         :: work(:)

  !
  complex(psb_spk_), allocatable :: xt(:), yt(:)
  integer(psb_ipk_) :: i, j, nr1, nc1,nr2, nc2,&
       & map_kind, nr, ictxt
  character(len=20), parameter  :: name='psb_map_X2Y'

  info = psb_success_
  if (.not.map%is_asb()) then 
    write(psb_err_unit,*) trim(name),' Invalid map input: unassembled'
    info = 1
    return 
  end if

  map_kind = map%get_kind()

  select case(map_kind)
  case(psb_map_aggr_)

    ictxt = map%p_desc_Y%get_context()
    nr2   = map%p_desc_Y%get_global_rows()
    nc2   = map%p_desc_Y%get_local_cols() 
    allocate(yt(nc2),stat=info) 
    if (info == psb_success_) call psb_halo(x,map%p_desc_X,info,work=work)
    if (info == psb_success_) call psb_csmm(cone,map%map_X2Y,x,czero,yt,info)
    if ((info == psb_success_) .and. psb_is_repl_desc(map%p_desc_Y)) then
      call psb_sum(ictxt,yt(1:nr2))
    end if
    if (info == psb_success_) call psb_geaxpby(alpha,yt,beta,y,map%p_desc_Y,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case(psb_map_gen_linear_)

    ictxt = map%desc_Y%get_context()
    nr1   = map%desc_X%get_local_rows() 
    nc1   = map%desc_X%get_local_cols() 
    nr2   = map%desc_Y%get_global_rows()
    nc2   = map%desc_Y%get_local_cols() 
    allocate(xt(nc1),yt(nc2),stat=info) 
    xt(1:nr1) = x(1:nr1) 
    if (info == psb_success_) call psb_halo(xt,map%desc_X,info,work=work)
    if (info == psb_success_) call psb_csmm(cone,map%map_X2Y,xt,czero,yt,info)
    if ((info == psb_success_) .and. psb_is_repl_desc(map%desc_Y)) then
      call psb_sum(ictxt,yt(1:nr2))
    end if
    if (info == psb_success_) call psb_geaxpby(alpha,yt,beta,y,map%desc_Y,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) trim(name),' Error from inner routines',info
      info = -1
    end if
   

  case default
    write(psb_err_unit,*) trim(name),' Invalid descriptor input', &
         & map_kind, psb_map_aggr_, psb_map_gen_linear_
    info = 1
    return 
  end select

end subroutine psb_c_map_X2Y

subroutine psb_c_map_X2Y_vect(alpha,x,beta,y,map,info,work,vtx,vty)
  use psb_base_mod, psb_protect_name => psb_c_map_X2Y_vect
  implicit none 
  type(psb_clinmap_type), intent(in)   :: map
  complex(psb_spk_), intent(in)        :: alpha,beta
  type(psb_c_vect_type), intent(inout) :: x,y
  integer(psb_ipk_), intent(out)                 :: info 
  complex(psb_spk_), optional          :: work(:)
  type(psb_c_vect_type), optional, target, intent(inout)  :: vtx,vty
  ! Local
  type(psb_c_vect_type), target  :: xt, yt
  type(psb_c_vect_type),pointer  :: ptx, pty
  complex(psb_spk_), allocatable :: xta(:), yta(:)
  integer(psb_ipk_) :: i, j, nr1, nc1,nr2, nc2 ,&
       &  map_kind, nr, ictxt, iam, np
  character(len=20), parameter   :: name='psb_map_X2Yv'

  info = psb_success_
  if (.not.map%is_asb()) then 
    write(psb_err_unit,*) trim(name),' Invalid map input: unassembled'
    info = 1
    return 
  end if

  map_kind = map%get_kind()

  select case(map_kind)
  case(psb_map_aggr_)

    ictxt = map%p_desc_Y%get_context()
    call psb_info(ictxt,iam,np)
    nr2   = map%p_desc_Y%get_global_rows()
    nc2   = map%p_desc_Y%get_local_cols()
    if (present(vty)) then
      pty => vty
    else
      call psb_geasb(yt,map%p_desc_Y,info,scratch=.true.,mold=x%v)
      pty => yt
    end if
    if (info == psb_success_) call psb_halo(x,map%p_desc_X,info,work=work)
    if (info == psb_success_) call psb_csmm(cone,map%map_X2Y,x,czero,pty,info)
    if ((info == psb_success_) .and. map%p_desc_Y%is_repl().and.(np>1)) then
      yta = pty%get_vect()
      call psb_sum(ictxt,yta(1:nr2))
      call pty%set(yta)
    end if
    if (info == psb_success_) call psb_geaxpby(alpha,pty,beta,y,map%p_desc_Y,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) trim(name),' Error from inner routines',info
      info = -1
    else 
      if (.not.present(vty)) call yt%free(info)
    end if

  case(psb_map_gen_linear_)

    ictxt = map%desc_Y%get_context()
    call psb_info(ictxt,iam,np)
    nr1   = map%desc_X%get_local_rows() 
    nc1   = map%desc_X%get_local_cols() 
    nr2   = map%desc_Y%get_global_rows()
    nc2   = map%desc_Y%get_local_cols() 

    if (present(vtx).and.present(vty)) then
      ptx => vtx
      pty => vty
    else
      call psb_geasb(xt,map%desc_X,info,scratch=.true.,mold=x%v)
      call psb_geasb(yt,map%desc_Y,info,scratch=.true.,mold=x%v)
      ptx => xt
      pty => yt
    end if

    call psb_geaxpby(cone,x,czero,ptx,map%desc_X,info)
    if (info == psb_success_) call psb_halo(ptx,map%desc_X,info,work=work)
    if (info == psb_success_) call psb_csmm(cone,map%map_X2Y,ptx,czero,pty,info)
    if ((info == psb_success_) .and. map%desc_Y%is_repl().and.(np>1)) then
      yta = pty%get_vect()
      call psb_sum(ictxt,yta(1:nr2))
      call pty%set(yta)
    end if
    if (info == psb_success_) call psb_geaxpby(alpha,pty,beta,y,map%desc_Y,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) trim(name),' Error from inner routines',info
      info = -1
    else
      if (.not.(present(vtx).and.present(vty) )) then 
        call xt%free(info)
        call yt%free(info)
      end if
    end if
   
  case default
    write(psb_err_unit,*) trim(name),' Invalid descriptor input', &
         & map_kind, psb_map_aggr_, psb_map_gen_linear_
    info = 1
    return 
  end select

  return
end subroutine psb_c_map_X2Y_vect


!
! Takes a vector x from space map%p_desc_Y and maps it onto
! map%p_desc_X under map%map_Y2X possibly with communication
! due to exch_bk_idx
!
subroutine psb_c_map_Y2X(alpha,x,beta,y,map,info,work)
  use psb_base_mod, psb_protect_name => psb_c_map_Y2X

  implicit none 
  type(psb_clinmap_type), intent(in) :: map
  complex(psb_spk_), intent(in)       :: alpha,beta
  complex(psb_spk_), intent(inout)    :: x(:)
  complex(psb_spk_), intent(out)      :: y(:)
  integer(psb_ipk_), intent(out)                :: info 
  complex(psb_spk_), optional         :: work(:)

  !
  complex(psb_spk_), allocatable :: xt(:), yt(:)
  integer(psb_ipk_) :: i, j, nr1, nc1,nr2, nc2,&
       & map_kind, nr, ictxt
  character(len=20), parameter  :: name='psb_map_Y2X'

  info = psb_success_
  if (.not.map%is_asb()) then 
    write(psb_err_unit,*) trim(name),' Invalid map input: unassembled'
    info = 1
    return 
  end if

  map_kind = map%get_kind()

  select case(map_kind)
  case(psb_map_aggr_)

    ictxt = map%p_desc_X%get_context()
    nr2   = map%p_desc_X%get_global_rows()
    nc2   = map%p_desc_X%get_local_cols() 
    allocate(yt(nc2),stat=info) 
    if (info == psb_success_) call psb_halo(x,map%p_desc_Y,info,work=work)
    if (info == psb_success_) call psb_csmm(cone,map%map_Y2X,x,czero,yt,info)
    if ((info == psb_success_) .and. psb_is_repl_desc(map%p_desc_X)) then
      call psb_sum(ictxt,yt(1:nr2))
    end if
    if (info == psb_success_) call psb_geaxpby(alpha,yt,beta,y,map%p_desc_X,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) trim(name),' Error from inner routines',info
      info = -1
    end if

  case(psb_map_gen_linear_)

    ictxt = map%desc_X%get_context()
    nr1   = map%desc_Y%get_local_rows() 
    nc1   = map%desc_Y%get_local_cols() 
    nr2   = map%desc_X%get_global_rows()
    nc2   = map%desc_X%get_local_cols() 
    allocate(xt(nc1),yt(nc2),stat=info) 
    xt(1:nr1) = x(1:nr1) 
    if (info == psb_success_) call psb_halo(xt,map%desc_Y,info,work=work)
    if (info == psb_success_) call psb_csmm(cone,map%map_Y2X,xt,czero,yt,info)
    if ((info == psb_success_) .and. psb_is_repl_desc(map%desc_X)) then
      call psb_sum(ictxt,yt(1:nr2))
    end if
    if (info == psb_success_) call psb_geaxpby(alpha,yt,beta,y,map%desc_X,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) trim(name),' Error from inner routines',info
      info = -1
    end if
   

  case default
    write(psb_err_unit,*) trim(name),' Invalid descriptor input'
    info = 1
    return 
  end select

end subroutine psb_c_map_Y2X

subroutine psb_c_map_Y2X_vect(alpha,x,beta,y,map,info,work,vtx,vty)
  use psb_base_mod, psb_protect_name => psb_c_map_Y2X_vect
  implicit none 
  type(psb_clinmap_type), intent(in)   :: map
  complex(psb_spk_), intent(in)        :: alpha,beta
  type(psb_c_vect_type), intent(inout) :: x,y
  integer(psb_ipk_), intent(out)                 :: info 
  complex(psb_spk_), optional          :: work(:)
  type(psb_c_vect_type), optional, target, intent(inout)  :: vtx,vty
  ! Local
  type(psb_c_vect_type), target  :: xt, yt
  type(psb_c_vect_type),pointer  :: ptx, pty
  complex(psb_spk_), allocatable :: xta(:), yta(:)
  integer(psb_ipk_) :: i, j, nr1, nc1,nr2, nc2,&
       & map_kind, nr, ictxt, iam, np
  character(len=20), parameter   :: name='psb_map_Y2Xv'

  info = psb_success_
  if (.not.map%is_asb()) then 
    write(psb_err_unit,*) trim(name),' Invalid map input: unassembled'
    info = 1
    return 
  end if

  map_kind = map%get_kind()

  select case(map_kind)
  case(psb_map_aggr_)

    ictxt = map%p_desc_X%get_context()
    call psb_info(ictxt,iam,np)
    nr2   = map%p_desc_X%get_global_rows()
    nc2   = map%p_desc_X%get_local_cols() 
    if (present(vty)) then
      pty => vty
    else
      call psb_geasb(yt,map%p_desc_X,info,scratch=.true.,mold=x%v)
      pty => yt
    end if
    if (info == psb_success_) call psb_halo(x,map%p_desc_Y,info,work=work)
    if (info == psb_success_) call psb_csmm(cone,map%map_Y2X,x,czero,pty,info)
    if ((info == psb_success_) .and. map%p_desc_X%is_repl().and.(np>1)) then
      yta = pty%get_vect()
      call psb_sum(ictxt,yta(1:nr2))
      call pty%set(yta)
    end if
    if (info == psb_success_) call psb_geaxpby(alpha,pty,beta,y,map%p_desc_X,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) trim(name),' Error from inner routines',info
      info = -1
    else
      if (.not.present(vty)) call yt%free(info)
    end if

  case(psb_map_gen_linear_)

    ictxt = map%desc_X%get_context()
    call psb_info(ictxt,iam,np)
    nr1   = map%desc_Y%get_local_rows() 
    nc1   = map%desc_Y%get_local_cols() 
    nr2   = map%desc_X%get_global_rows()
    nc2   = map%desc_X%get_local_cols() 
    if (present(vtx).and.present(vty)) then
      ptx => vtx
      pty => vty
    else
      call psb_geasb(xt,map%desc_Y,info,scratch=.true.,mold=x%v)
      call psb_geasb(yt,map%desc_X,info,scratch=.true.,mold=x%v)
      ptx => xt
      pty => yt
    end if

    call psb_geaxpby(cone,x,czero,ptx,map%desc_Y,info)

    if (info == psb_success_) call psb_halo(ptx,map%desc_Y,info,work=work)
    if (info == psb_success_) call psb_csmm(cone,map%map_Y2X,ptx,czero,pty,info)
    if ((info == psb_success_) .and. map%desc_X%is_repl().and.(np>1)) then
      yta = pty%get_vect()
      call psb_sum(ictxt,yta(1:nr2))
      call pty%set(yta)
    end if
    if (info == psb_success_) call psb_geaxpby(alpha,pty,beta,y,map%desc_X,info)
    if (info /= psb_success_) then 
      write(psb_err_unit,*) trim(name),' Error from inner routines',info
      info = -1
    else
      if (.not.(present(vtx).and.present(vty) )) then 
        call xt%free(info)
        call yt%free(info)
      end if
    end if

  case default
    write(psb_err_unit,*) trim(name),' Invalid descriptor input'
    info = 1
    return 
  end select

end subroutine psb_c_map_Y2X_vect

function psb_c_linmap(map_kind,desc_X, desc_Y, map_X2Y, map_Y2X,iaggr,naggr) &
     & result(this)

  use psb_base_mod, psb_protect_name => psb_c_linmap
  implicit none 
  type(psb_clinmap_type)         :: this
  type(psb_desc_type), target       :: desc_X, desc_Y
  type(psb_cspmat_type), intent(inout) :: map_X2Y, map_Y2X
  integer(psb_ipk_), intent(in)               :: map_kind
  integer(psb_ipk_), intent(in), optional     :: iaggr(:), naggr(:)
  !
  integer(psb_ipk_) :: info
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
          this%iaggr(:) = iaggr(:)
          this%naggr(:) = naggr(:)
        end if
      end if
    else
      allocate(this%iaggr(0), this%naggr(0), stat=info) 
    end if

  case(psb_map_gen_linear_)    
    
    if (desc_X%is_ok()) then 
      call desc_X%clone(this%desc_X,info) 
    else
      info = psb_err_pivot_too_small_
    endif
    if (desc_Y%is_ok()) then 
      call desc_Y%clone(this%desc_Y,info) 
    else
      info = psb_err_invalid_ovr_num_
    endif
    ! If iaggr/naggr are present, copy them anyway.
    if (present(iaggr)) then 
      if (.not.present(naggr)) then 
        info = 7
      else
        allocate(this%iaggr(size(iaggr)),&
             & this%naggr(size(naggr)), stat=info) 
        if (info == psb_success_) then 
          this%iaggr(:) = iaggr(:)
          this%naggr(:) = naggr(:)
        end if
      end if
    else
      allocate(this%iaggr(0), this%naggr(0), stat=info) 
    end if

  case default
    write(psb_err_unit,*) 'Bad map kind into psb_linmap ',map_kind
    info = 1
  end select

  if (info == psb_success_) call map_X2Y%clone(this%map_X2Y,info)
  if (info == psb_success_) call map_Y2X%clone(this%map_Y2X,info)
  if (info == psb_success_) then
    call this%set_kind(map_kind)
  end if
  if (info /= psb_success_) then
    write(psb_err_unit,*) trim(name),' Invalid descriptor input'
    return
  end if

end function psb_c_linmap
