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
! package: psb_glist_map_mod
!    Defines the GLIST_MAP type.
!
! This is almost identical to the LIST_MAP type, but it has an additional
! vector of size GLOB_ROWS giving, for each index, the owning process.
! This implies that:
! 1. We have room for such an additional vector;
! 2. There are no overlap (only one process owns a given index). 
!
!
module psb_glist_map_mod
  use psb_const_mod
  use psb_desc_const_mod
  use psb_list_map_mod
  
  type, extends(psb_list_map) :: psb_glist_map
    integer(psb_ipk_), allocatable :: vgp(:)
  contains
    procedure, pass(idxmap)  :: glist_map_init   => glist_initvg
    procedure, pass(idxmap)  :: sizeof  => glist_sizeof
    procedure, pass(idxmap)  :: free    => glist_free
    procedure, pass(idxmap)  :: clone   => glist_clone
    procedure, nopass        :: get_fmt => glist_get_fmt
    procedure, pass(idxmap)  :: fnd_owner => glist_fnd_owner

  end type psb_glist_map

  private :: glist_initvg, glist_sizeof, glist_free, glist_get_fmt


contains

    
  function glist_sizeof(idxmap) result(val)
    implicit none 
    class(psb_glist_map), intent(in) :: idxmap
    integer(psb_epk_) :: val
    
    val = idxmap%psb_list_map%sizeof()

    if (allocated(idxmap%vgp)) &
         & val = val + size(idxmap%vgp)*psb_sizeof_ip

  end function glist_sizeof


  subroutine glist_free(idxmap)
    implicit none 
    class(psb_glist_map), intent(inout) :: idxmap
    
    if (allocated(idxmap%vgp)) &
         & deallocate(idxmap%vgp)
    
    call idxmap%psb_list_map%free()
    
  end subroutine glist_free




  subroutine glist_initvg(idxmap,ictxt,vg,info)
    use psb_penv_mod
    use psb_error_mod
    implicit none 
    class(psb_glist_map), intent(inout) :: idxmap
    integer(psb_mpk_), intent(in)  :: ictxt    
    integer(psb_ipk_), intent(in)  :: vg(:)
    integer(psb_ipk_), intent(out) :: info
    !  To be implemented
    integer(psb_mpk_) :: iam, np
    integer(psb_ipk_) :: i, n, nl
    

    info = 0
    call psb_info(ictxt,iam,np) 
    if (np < 0) then 
      write(psb_err_unit,*) 'Invalid ictxt:',ictxt
      info = -1
      return
    end if
    n = size(vg) 
    
    idxmap%global_rows  = n
    idxmap%global_cols  = n

    allocate(idxmap%loc_to_glob(n),idxmap%glob_to_loc(n),&
         & idxmap%vgp(n),stat=info) 
    if (info /= 0)  then
      info = -2
      return
    end if

    idxmap%ictxt        = ictxt
    idxmap%state        = psb_desc_bld_
    call psb_get_mpicomm(ictxt,idxmap%mpic)

    nl = 0 
    do i=1, n 
      if ((vg(i)  > np-1).or.(vg(i) < 0)) then
        info=psb_err_partfunc_wrong_pid_
        exit
      end if
      idxmap%vgp(i) = vg(i)
      if (vg(i) == iam) then
        ! this point belongs to me
        nl = nl + 1
        idxmap%glob_to_loc(i)  = nl
        idxmap%loc_to_glob(nl) = i
      else
        idxmap%glob_to_loc(i) = -(np+vg(i)+1)
      end if
    end do
    
    call idxmap%set_lr(nl)
    call idxmap%set_lc(nl)
   
  end subroutine glist_initvg

  subroutine glist_fnd_owner(idx,iprc,idxmap,info)
    use psb_penv_mod
    use psb_sort_mod
    implicit none 
    integer(psb_ipk_), intent(in) :: idx(:)
    integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
    class(psb_glist_map), intent(in) :: idxmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_mpk_) :: ictxt, iam, np
    integer(psb_ipk_) :: nv, i, ngp
    
    ictxt = idxmap%get_ctxt()
    call psb_info(ictxt,iam,np)
    nv = size(idx)
    allocate(iprc(nv),stat=info) 
    if (info /= 0) then 
      write(0,*) 'Memory allocation failure in repl_map_fnd-owner'
      return
    end if

    ngp = size(idxmap%vgp)
    do i=1, nv 
      if ((1<=idx(i)).and.(idx(i)<=ngp)) then
        iprc(i) = idxmap%vgp(idx(i))
      else
        iprc(i) = -1
      end if
    end do

  end subroutine glist_fnd_owner

  function glist_get_fmt() result(res)
    implicit none 
    character(len=5) :: res
    res = 'GLIST'
  end function glist_get_fmt



  subroutine glist_clone(idxmap,outmap,info)
    use psb_penv_mod
    use psb_error_mod
    use psb_realloc_mod
    implicit none 
    class(psb_glist_map), intent(inout)    :: idxmap
    class(psb_indx_map), allocatable, intent(out) :: outmap
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='glist_clone'
    logical, parameter :: debug=.false.

    info = psb_success_
    call psb_get_erraction(err_act)
    if (allocated(outmap)) then 
      write(0,*) 'Error: should not be allocated on input'
      info = -87
      goto 9999
    end if
    
    allocate(psb_glist_map :: outmap, stat=info) 
    if (info /= psb_success_) then 
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if

    select type (outmap)
    type is (psb_glist_map) 

      if (info == psb_success_) then 
        outmap%psb_indx_map = idxmap%psb_indx_map
        outmap%pnt_h        = idxmap%pnt_h
      end if
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%loc_to_glob,outmap%loc_to_glob,info)
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%glob_to_loc,outmap%glob_to_loc,info)
      if (info == psb_success_)&
           &  call psb_safe_ab_cpy(idxmap%vgp,outmap%vgp,info)
    class default
      ! This should be impossible 
      info = -1
    end select
      
    if (info /= psb_success_) then 
      info = psb_err_from_subroutine_
      call psb_errpush(info,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return


9999 call psb_error_handler(err_act)

    return
  end subroutine glist_clone

end module psb_glist_map_mod
