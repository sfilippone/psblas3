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
!  d_mat_impl:
!   implementation of the outer matrix methods.
!   Most of the methods rely on the STATE design pattern:
!   the inner class(psb_d_base_sparse_mat) is responsbile
!   for actually executing the method.
!
!


! == ===================================
!
!
!
! Setters 
!
!
!
!
!
!
! == ===================================


subroutine  psb_d_set_nrows(m,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_nrows
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in) :: m
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='set_nrows'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_nrows(m)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_nrows


subroutine  psb_d_set_ncols(n,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_ncols
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in) :: n
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call a%a%set_ncols(n)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_ncols



!
!  Valid values for DUPL: 
!  psb_dupl_ovwrt_ 
!  psb_dupl_add_   
!  psb_dupl_err_   
!

subroutine  psb_d_set_dupl(n,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_dupl
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in) :: n
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_dupl(n)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_dupl


!
! Set the STATE of the internal matrix object
!

subroutine  psb_d_set_null(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_null
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_null()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_null


subroutine  psb_d_set_bld(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_bld
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_bld()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_bld


subroutine  psb_d_set_upd(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_upd
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_upd()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return


end subroutine psb_d_set_upd


subroutine  psb_d_set_asb(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_asb
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_asb()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_asb


subroutine psb_d_set_sorted(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_sorted
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_sorted(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_sorted


subroutine psb_d_set_triangle(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_triangle
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_triangle(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_triangle


subroutine psb_d_set_unit(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_unit
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_unit(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_unit


subroutine psb_d_set_lower(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_lower
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_lower(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_lower


subroutine psb_d_set_upper(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_d_set_upper
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_upper(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_set_upper



! == ===================================
!
!
!
! Data management
!
!
!
!
!
! == ===================================  


subroutine psb_d_sparse_print(iout,a,iv,head,ivr,ivc)
  use psb_d_mat_mod, psb_protect_name => psb_d_sparse_print
  use psb_error_mod
  implicit none 

  integer(psb_ipk_), intent(in)               :: iout
  class(psb_dspmat_type), intent(in) :: a   
  integer(psb_ipk_), intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%print(iout,iv,head,ivr,ivc)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_sparse_print


subroutine psb_d_n_sparse_print(fname,a,iv,head,ivr,ivc)
  use psb_d_mat_mod, psb_protect_name => psb_d_n_sparse_print
  use psb_error_mod
  implicit none 

  character(len=*), intent(in)  :: fname   
  class(psb_dspmat_type), intent(in) :: a   
  integer(psb_ipk_), intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act, info, iout
  logical :: isopen
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  iout = max(psb_inp_unit,psb_err_unit,psb_out_unit) + 1
  do 
    inquire(unit=iout, opened=isopen)
    if (.not.isopen) exit
    iout = iout + 1
    if (iout > 99) exit
  end do
  if (iout > 99) then 
    write(psb_err_unit,*) 'Error: could not find a free unit for I/O'
    return
  end if
  open(iout,file=fname,iostat=info)
  if (info == psb_success_) then 
    call a%a%print(iout,iv,head,ivr,ivc)
    close(iout)
  else
    write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
  end if

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_n_sparse_print


subroutine psb_d_get_neigh(a,idx,neigh,n,info,lev)
  use psb_d_mat_mod, psb_protect_name => psb_d_get_neigh
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a   
  integer(psb_ipk_), intent(in)                :: idx 
  integer(psb_ipk_), intent(out)               :: n   
  integer(psb_ipk_), allocatable, intent(out)  :: neigh(:)
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_), optional, intent(in)      :: lev 

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='get_neigh'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%get_neigh(idx,neigh,n,info,lev)

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_get_neigh



subroutine psb_d_csall(nr,nc,a,info,nz) 
  use psb_d_mat_mod, psb_protect_name => psb_d_csall
  use psb_d_base_mat_mod
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in)             :: nr,nc
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: nz

  integer(psb_ipk_) :: err_act 
  character(len=20)  :: name='csall'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)

  call a%free()

  info = psb_success_
  allocate(psb_d_coo_sparse_mat :: a%a, stat=info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  call a%a%allocate(nr,nc,nz)
  call a%set_bld() 

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csall


subroutine  psb_d_reallocate_nz(nz,a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_reallocate_nz
  use psb_error_mod
  implicit none 
  integer(psb_ipk_), intent(in) :: nz
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%reallocate(nz)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_reallocate_nz


subroutine  psb_d_free(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_free
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a

  if (allocated(a%a)) then 
    call a%a%free()
    deallocate(a%a) 
  endif

end subroutine psb_d_free


subroutine  psb_d_trim(a) 
  use psb_d_mat_mod, psb_protect_name => psb_d_trim
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%trim()

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_trim



subroutine psb_d_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_d_mat_mod, psb_protect_name => psb_d_csput_a
  use psb_d_base_mat_mod
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: gtl(:)

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csput_a'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.(a%is_bld().or.a%is_upd())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info,gtl) 
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csput_a

subroutine psb_d_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_d_mat_mod, psb_protect_name => psb_d_csput_v
  use psb_d_base_mat_mod
  use psb_d_vect_mod, only : psb_d_vect_type
  use psb_i_vect_mod, only : psb_i_vect_type
  use psb_error_mod
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  type(psb_d_vect_type), intent(inout)  :: val
  type(psb_i_vect_type), intent(inout)  :: ia, ja
  integer(psb_ipk_), intent(in)             :: nz, imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: gtl(:)

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csput_v'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.(a%is_bld().or.a%is_upd())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  
  if (allocated(val%v).and.allocated(ia%v).and.allocated(ja%v)) then
    call a%a%csput(nz,ia%v,ja%v,val%v,imin,imax,jmin,jmax,info,gtl) 
  else
    info = psb_err_invalid_mat_state_
  endif

  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csput_v


subroutine psb_d_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csgetptn
  implicit none

  class(psb_dspmat_type), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csget(imin,imax,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csgetptn


subroutine psb_d_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csgetrow
  implicit none

  class(psb_dspmat_type), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csget(imin,imax,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csgetrow




subroutine psb_d_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csgetblk
  implicit none

  class(psb_dspmat_type), intent(in)    :: a
  class(psb_dspmat_type), intent(inout) :: b
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.
  logical            :: append_
  type(psb_d_coo_sparse_mat), allocatable  :: acoo


  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (present(append))  then 
    append_ = append
  else
    append_ = .false.
  end if

  allocate(acoo,stat=info)    
  if (append_.and.(info==psb_success_)) then 
    if (allocated(b%a)) &
         & call b%a%mv_to_coo(acoo,info)
  end if

  if (info == psb_success_) then 
    call a%a%csget(imin,imax,acoo,info,&
         & jmin,jmax,iren,append,rscale,cscale)
  else
    info = psb_err_alloc_dealloc_
  end if
  if (info == psb_success_) call move_alloc(acoo,b%a)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csgetblk


subroutine psb_d_tril(a,l,info,diag,imin,imax,&
     & jmin,jmax,rscale,cscale,u)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_tril
  implicit none 
  class(psb_dspmat_type), intent(in)      :: a
  class(psb_dspmat_type), intent(inout)   :: l
  integer(psb_ipk_),intent(out)           :: info
  integer(psb_ipk_), intent(in), optional :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional           :: rscale,cscale
  class(psb_dspmat_type), optional, intent(inout)   :: u

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='tril'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: lcoo, ucoo

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(lcoo,stat=info)    
  call l%free()
  if (present(u)) then
    if (info == psb_success_) allocate(ucoo,stat=info)    
    call u%free()
    if (info == psb_success_) call a%a%tril(lcoo,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,ucoo)
    if (info == psb_success_) call move_alloc(ucoo,u%a)
    if (info == psb_success_) call u%cscnv(info,mold=a%a)
  else
    if (info == psb_success_) then 
      call a%a%tril(lcoo,info,diag,imin,imax,&
           & jmin,jmax,rscale,cscale)
    else
      info = psb_err_alloc_dealloc_
    end if
  end if
  if (info == psb_success_) call move_alloc(lcoo,l%a)
  if (info == psb_success_) call l%cscnv(info,mold=a%a)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return


end subroutine psb_d_tril

subroutine psb_d_triu(a,u,info,diag,imin,imax,&
     & jmin,jmax,rscale,cscale,l)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_triu
  implicit none 
  class(psb_dspmat_type), intent(in)      :: a
  class(psb_dspmat_type), intent(inout)   :: u
  integer(psb_ipk_),intent(out)           :: info
  integer(psb_ipk_), intent(in), optional :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional           :: rscale,cscale
  class(psb_dspmat_type), optional, intent(inout)   :: l

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='triu'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: lcoo, ucoo


  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(ucoo,stat=info)    
  call u%free()

  if (present(l)) then
    if (info == psb_success_) allocate(lcoo,stat=info)    
    call l%free()
    if (info == psb_success_) call a%a%triu(ucoo,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,lcoo)
    if (info == psb_success_) call move_alloc(lcoo,l%a)
    if (info == psb_success_) call l%cscnv(info,mold=a%a)
  else
    if (info == psb_success_) then 
      call a%a%triu(ucoo,info,diag,imin,imax,&
           & jmin,jmax,rscale,cscale)
    else
      info = psb_err_alloc_dealloc_
    end if
  end if
  if (info == psb_success_) call move_alloc(ucoo,u%a)
  if (info == psb_success_) call u%cscnv(info,mold=a%a)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return


end subroutine psb_d_triu


subroutine psb_d_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csclip
  implicit none

  class(psb_dspmat_type), intent(in) :: a
  class(psb_dspmat_type), intent(inout) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: acoo

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  call b%free()
  if (info == psb_success_) then 
    call a%a%csclip(acoo,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else
    info = psb_err_alloc_dealloc_
  end if
 
  if (info == psb_success_) call move_alloc(acoo,b%a)
  if (info /= psb_success_) goto 9999 
  
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csclip


subroutine psb_d_b_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_b_csclip
  implicit none

  class(psb_dspmat_type), intent(in) :: a
  type(psb_d_coo_sparse_mat), intent(out) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%csclip(b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_b_csclip




subroutine psb_d_cscnv(a,b,info,type,mold,upd,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cscnv
  implicit none 
  class(psb_dspmat_type), intent(in)      :: a
  class(psb_dspmat_type), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psb_d_base_sparse_mat), intent(in), optional :: mold


  class(psb_d_base_sparse_mat), allocatable  :: altmp
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call b%free()
  if (count( (/present(mold),present(type) /)) > 1) then
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='TYPE, MOLD')
    goto 9999
  end if

  if (present(mold)) then 

    allocate(altmp, mold=mold,stat=info) 

  else if (present(type)) then 

    select case (psb_toupper(type))
    case ('CSR')
      allocate(psb_d_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psb_d_coo_sparse_mat :: altmp, stat=info) 
    case ('CSC')
      allocate(psb_d_csc_sparse_mat :: altmp, stat=info) 
    case default
      info = psb_err_format_unknown_ 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(altmp, mold=psb_get_mat_default(a),stat=info) 
  end if

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  
  if (present(dupl)) then 
    call altmp%set_dupl(dupl)
  else if (a%is_bld()) then 
    ! Does this make sense at all?? Who knows..
    call altmp%set_dupl(psb_dupl_def_)
  end if

  if (debug) write(psb_err_unit,*) 'Converting from ',&
       & a%get_fmt(),' to ',altmp%get_fmt()

  call altmp%cp_from_fmt(a%a, info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call move_alloc(altmp,b%a)
  call b%trim()
  call b%asb() 
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_cscnv



subroutine psb_d_cscnv_ip(a,info,type,mold,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cscnv_ip
  implicit none 

  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl
  character(len=*), optional, intent(in) :: type
  class(psb_d_base_sparse_mat), intent(in), optional :: mold


  class(psb_d_base_sparse_mat), allocatable  :: altmp
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cscnv_ip'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(dupl)) then 
    call a%set_dupl(dupl)
  else if (a%is_bld()) then 
    call a%set_dupl(psb_dupl_def_)
  end if

  if (count( (/present(mold),present(type) /)) > 1) then
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='TYPE, MOLD')
    goto 9999
  end if

  if (present(mold)) then 

    allocate(altmp, mold=mold,stat=info) 

  else if (present(type)) then 

    select case (psb_toupper(type))
    case ('CSR')
      allocate(psb_d_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psb_d_coo_sparse_mat :: altmp, stat=info) 
    case ('CSC')
      allocate(psb_d_csc_sparse_mat :: altmp, stat=info) 
    case default
      info = psb_err_format_unknown_ 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(altmp, mold=psb_get_mat_default(a),stat=info) 
  end if

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug) write(psb_err_unit,*) 'Converting in-place from ',&
       & a%get_fmt(),' to ',altmp%get_fmt()

  call altmp%mv_from_fmt(a%a, info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call move_alloc(altmp,a%a)
  call a%set_asb() 
  call a%trim()
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_cscnv_ip



subroutine psb_d_cscnv_base(a,b,info,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cscnv_base
  implicit none 
  class(psb_dspmat_type), intent(in)       :: a
  class(psb_d_base_sparse_mat), intent(out) :: b
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl


  type(psb_d_coo_sparse_mat)  :: altmp
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%cp_to_coo(altmp,info )
  if ((info == psb_success_).and.present(dupl)) then 
    call altmp%set_dupl(dupl)
  end if
  call altmp%fix(info)
  if (info == psb_success_) call altmp%trim()
  if (info == psb_success_) call altmp%set_asb() 
  if (info == psb_success_) call b%mv_from_coo(altmp,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_cscnv_base



subroutine psb_d_clip_d(a,b,info)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_clip_d
  implicit none

  class(psb_dspmat_type), intent(in)    :: a
  class(psb_dspmat_type), intent(inout) :: b
  integer(psb_ipk_),intent(out)                  :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: acoo
  integer(psb_ipk_) :: i, j, nz

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  if (info == psb_success_) call a%a%cp_to_coo(acoo,info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  nz = acoo%get_nzeros()
  j = 0
  do i=1, nz
    if (acoo%ia(i) /= acoo%ja(i)) then 
      j = j + 1 
      acoo%ia(j)  = acoo%ia(i)
      acoo%ja(j)  = acoo%ja(i)
      acoo%val(j) = acoo%val(i)
    end if
  end do
  call acoo%set_nzeros(j)
  call acoo%trim()
  call b%mv_from(acoo)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_clip_d



subroutine psb_d_clip_d_ip(a,info)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_clip_d_ip
  implicit none

  class(psb_dspmat_type), intent(inout) :: a
  integer(psb_ipk_),intent(out)                  :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_d_coo_sparse_mat), allocatable  :: acoo
  integer(psb_ipk_) :: i, j, nz

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  if (info == psb_success_) call a%a%mv_to_coo(acoo,info)
  if (info /= psb_success_) then
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  endif

  nz = acoo%get_nzeros()
  j = 0
  do i=1, nz
    if (acoo%ia(i) /= acoo%ja(i)) then 
      j = j + 1 
      acoo%ia(j)  = acoo%ia(i)
      acoo%ja(j)  = acoo%ja(i)
      acoo%val(j) = acoo%val(i)
    end if
  end do
  call acoo%set_nzeros(j)
  call acoo%trim()
  call a%mv_from(acoo)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_clip_d_ip


subroutine psb_d_mv_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_mv_from
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  call a%free()
  allocate(a%a,mold=b, stat=info)
  call a%a%mv_from_fmt(b,info)
  call b%free()

  return
end subroutine psb_d_mv_from


subroutine psb_d_cp_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cp_from
  implicit none 
  class(psb_dspmat_type), intent(out)      :: a
  class(psb_d_base_sparse_mat), intent(in) :: b
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='cp_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  
  call a%free()
  !
  ! Note: it is tempting to use SOURCE allocation below;
  ! however this would run the risk of messing up with data
  ! allocated externally (e.g. GPU-side data).
  !
  allocate(a%a,mold=b,stat=info)
  if (info /= psb_success_) info = psb_err_alloc_dealloc_
  if (info == psb_success_) call a%a%cp_from_fmt(b, info)    
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return
end subroutine psb_d_cp_from


subroutine psb_d_mv_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_mv_to
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  call b%mv_from_fmt(a%a,info)

  return
end subroutine psb_d_mv_to


subroutine psb_d_cp_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cp_to
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  call b%cp_from_fmt(a%a,info)

  return
end subroutine psb_d_cp_to

subroutine psb_d_mold(a,b)
  use psb_d_mat_mod, psb_protect_name => psb_d_mold
  class(psb_dspmat_type), intent(inout)     :: a
  class(psb_d_base_sparse_mat), allocatable, intent(out) :: b
  integer(psb_ipk_) :: info

  allocate(b,mold=a%a, stat=info)
  
end subroutine psb_d_mold

subroutine psb_dspmat_type_move(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_dspmat_type_move
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  class(psb_dspmat_type), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                   :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='move_alloc'
  logical, parameter :: debug=.false.

  info = psb_success_
  call b%free()
  call move_alloc(a%a,b%a)

  return
end subroutine psb_dspmat_type_move


subroutine psb_dspmat_clone(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_dspmat_clone
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  class(psb_dspmat_type), intent(inout) :: b
  integer(psb_ipk_), intent(out)        :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='clone'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  call b%free()
  if (allocated(a%a)) then 
    call a%a%clone(b%a,info)
  end if
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_dspmat_clone


subroutine psb_d_transp_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_transp_1mat
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%transp()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_transp_1mat



subroutine psb_d_transp_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_transp_2mat
  implicit none 
  class(psb_dspmat_type), intent(in)  :: a
  class(psb_dspmat_type), intent(inout) :: b

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call b%free()
  allocate(b%a,mold=a%a,stat=info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    goto 9999
  end if
  call a%a%transp(b%a)    

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_transp_2mat


subroutine psb_d_transc_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_transc_1mat
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%transc()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_transc_1mat



subroutine psb_d_transc_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_transc_2mat
  implicit none 
  class(psb_dspmat_type), intent(in)    :: a
  class(psb_dspmat_type), intent(inout) :: b

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call b%free()
  allocate(b%a,mold=a%a,stat=info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    goto 9999
  end if
  call a%a%transc(b%a)    

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_transc_2mat


subroutine psb_d_asb(a,mold)
  use psb_d_mat_mod, psb_protect_name => psb_d_asb
  use psb_error_mod
  implicit none 

  class(psb_dspmat_type), intent(inout) :: a   
  class(psb_d_base_sparse_mat), optional, intent(in) :: mold
  class(psb_d_base_sparse_mat), allocatable :: tmp
  class(psb_d_base_sparse_mat), pointer :: mld
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='d_asb'

  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%asb()
  if (present(mold)) then 
    if (.not.same_type_as(a%a,mold)) then 
      allocate(tmp,mold=mold)
      call tmp%mv_from_fmt(a%a,info)
      call a%a%free()
      call move_alloc(tmp,a%a)
    end if
  else
    mld => psb_d_get_base_mat_default()
    if (.not.same_type_as(a%a,mld)) &
         & call a%cscnv(info)
  end if
  

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_asb

subroutine psb_d_reinit(a,clear)
  use psb_d_mat_mod, psb_protect_name => psb_d_reinit
  use psb_error_mod
  implicit none 

  class(psb_dspmat_type), intent(inout) :: a   
  logical, intent(in), optional :: clear
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='reinit'

  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (a%a%has_update()) then 
    call a%a%reinit(clear)
  else
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_reinit




! == ===================================
!
!
!
! Computational routines
!
!
!
!
!
!
! == ===================================


subroutine psb_d_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csmm
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='psb_csmm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%spmm(alpha,x,beta,y,info,trans) 
  if (info /= psb_success_) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csmm


subroutine psb_d_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csmv
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='psb_csmv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%spmm(alpha,x,beta,y,info,trans) 
  if (info /= psb_success_) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csmv

subroutine psb_d_csmv_vect(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_d_vect_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_csmv_vect
  implicit none 
  class(psb_dspmat_type), intent(in)   :: a
  real(psb_dpk_), intent(in)        :: alpha, beta
  type(psb_d_vect_type), intent(inout) :: x
  type(psb_d_vect_type), intent(inout) :: y
  integer(psb_ipk_), intent(out)                 :: info
  character, optional, intent(in)      :: trans
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='psb_csmv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(y%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%spmm(alpha,x%v,beta,y%v,info,trans) 
  if (info /= psb_success_) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_csmv_vect



subroutine psb_d_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_error_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cssm
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  real(psb_dpk_), intent(in), optional :: d(:)
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='psb_cssm'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%spsm(alpha,x,beta,y,info,trans,scale,d) 
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_cssm


subroutine psb_d_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_error_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cssv
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  real(psb_dpk_), intent(in), optional :: d(:)
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='psb_cssv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%spsm(alpha,x,beta,y,info,trans,scale,d) 

  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_cssv


subroutine psb_d_cssv_vect(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_error_mod
  use psb_d_vect_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cssv_vect
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), intent(in)         :: alpha, beta
  type(psb_d_vect_type), intent(inout)   :: x
  type(psb_d_vect_type), intent(inout)   :: y
  integer(psb_ipk_), intent(out)               :: info
  character, optional, intent(in)    :: trans, scale
  type(psb_d_vect_type), optional, intent(inout)   :: d
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='psb_cssv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(x%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(y%v)) then 
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (present(d)) then 
    if (.not.allocated(d%v)) then 
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
    endif
    call a%a%spsm(alpha,x%v,beta,y%v,info,trans,scale,d%v) 
  else
    call a%a%spsm(alpha,x%v,beta,y%v,info,trans,scale) 
  end if

  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_cssv_vect

function psb_d_maxval(a) result(res)
  use psb_d_mat_mod, psb_protect_name => psb_d_maxval
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_)         :: res


  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='maxval'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  res = a%a%maxval()
  return


9999 call psb_error_handler(err_act)

  return

end function psb_d_maxval

function psb_d_csnmi(a) result(res)
  use psb_d_mat_mod, psb_protect_name => psb_d_csnmi
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='csnmi'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  res = a%a%spnmi()
  return


9999 call psb_error_handler(err_act)

  return

end function psb_d_csnmi


function psb_d_csnm1(a) result(res)
  use psb_d_mat_mod, psb_protect_name => psb_d_csnm1
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='csnm1'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  res = a%a%spnm1()
  return


9999 call psb_error_handler(err_act)

  return

end function psb_d_csnm1


function psb_d_rowsum(a,info) result(d)
  use psb_d_mat_mod, psb_protect_name => psb_d_rowsum
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), allocatable     :: d(:)
  integer(psb_ipk_), intent(out)               :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(d(max(1,a%a%get_nrows())), stat=info)
  if (info /= psb_success_) goto 9999
  call a%a%rowsum(d)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_d_rowsum

function psb_d_arwsum(a,info) result(d)
  use psb_d_mat_mod, psb_protect_name => psb_d_arwsum
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), allocatable           :: d(:)
  integer(psb_ipk_), intent(out)       :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='arwsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(d(max(1,a%a%get_nrows())), stat=info)
  if (info /= psb_success_) goto 9999

  call a%a%arwsum(d)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_d_arwsum

function psb_d_colsum(a,info) result(d)
  use psb_d_mat_mod, psb_protect_name => psb_d_colsum
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), allocatable         :: d(:)
  integer(psb_ipk_), intent(out)       :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(d(max(1,a%a%get_ncols())), stat=info)
  if (info /= psb_success_) goto 9999

  call a%a%colsum(d)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_d_colsum

function psb_d_aclsum(a,info) result(d)
  use psb_d_mat_mod, psb_protect_name => psb_d_aclsum
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), allocatable           :: d(:)
  integer(psb_ipk_), intent(out)       :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='aclsum'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(d(max(1,a%a%get_ncols())), stat=info)
  if (info /= psb_success_) goto 9999

  call a%a%aclsum(d)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_d_aclsum


function psb_d_get_diag(a,info) result(d)
  use psb_d_mat_mod, psb_protect_name => psb_d_get_diag
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  real(psb_dpk_), allocatable         :: d(:)
  integer(psb_ipk_), intent(out)       :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(d(max(1,min(a%a%get_nrows(),a%a%get_ncols()))), stat=info)
  if (info /= 0) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info,name)
    goto 9999
  end if
  call a%a%get_diag(d,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_d_get_diag


subroutine psb_d_scal(d,a,info,side)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_scal
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  real(psb_dpk_), intent(in)              :: d(:)
  integer(psb_ipk_), intent(out)                    :: info
  character, intent(in), optional :: side

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scal(d,info,side=side)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_scal


subroutine psb_d_scals(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_scals
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  real(psb_dpk_), intent(in)              :: d
  integer(psb_ipk_), intent(out)                    :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scal(d,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_d_scals

subroutine psb_d_mv_from_lb(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_mv_from_lb
  implicit none 
  
  class(psb_dspmat_type), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  if (.not.allocated(a%a)) allocate(psb_d_csr_sparse_mat :: a%a, stat=info)
  if (info == psb_success_) call a%a%mv_from_lfmt(b,info)
  
end subroutine psb_d_mv_from_lb

  
subroutine psb_d_cp_from_lb(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cp_from_lb
  implicit none 
  
  class(psb_dspmat_type), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (.not.allocated(a%a)) allocate(psb_d_csr_sparse_mat :: a%a, stat=info)
  if (info == psb_success_) call a%a%cp_from_lfmt(b,info)

end subroutine psb_d_cp_from_lb

subroutine psb_d_mv_to_lb(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_mv_to_lb
  implicit none 
  
  class(psb_dspmat_type), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  if (.not.allocated(a%a)) then
    call b%free()
  else
    call a%a%mv_to_lfmt(b,info)
    call a%free()
  end if
  
end subroutine psb_d_mv_to_lb

subroutine psb_d_cp_to_lb(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cp_to_lb
  implicit none 
  class(psb_dspmat_type), intent(in) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  if (.not.allocated(a%a)) then
    call b%free()
  else
    call a%a%cp_to_lfmt(b,info)
  end if

end subroutine psb_d_cp_to_lb

subroutine psb_d_mv_from_l(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_mv_from_l
  implicit none 
  class(psb_dspmat_type), intent(inout) :: a
  class(psb_ldspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (allocated(b%a)) then
    if (.not.allocated(a%a)) allocate(psb_d_csr_sparse_mat :: a%a, stat=info)
    call a%a%mv_from_lfmt(b%a,info)
  else
    call a%free()
  end if
  call b%free()
  
end subroutine psb_d_mv_from_l

  
subroutine psb_d_cp_from_l(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cp_from_l
  implicit none 

  class(psb_dspmat_type), intent(out) :: a
  class(psb_ldspmat_type), intent(in) :: b
  integer(psb_ipk_) :: info
  
  if (allocated(b%a)) then
    if (.not.allocated(a%a)) allocate(psb_d_csr_sparse_mat :: a%a, stat=info)
    call a%a%cp_from_lfmt(b%a,info)
  else
    call a%free()
  end if
end subroutine psb_d_cp_from_l

subroutine psb_d_mv_to_l(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_mv_to_l
  implicit none 

  class(psb_dspmat_type), intent(inout) :: a
  class(psb_ldspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  if (allocated(a%a)) then
    if (.not.allocated(b%a)) allocate(psb_ld_csr_sparse_mat :: b%a, stat=info)
    call a%a%mv_to_lfmt(b%a,info)
  else
    call b%free()
  end if
  call a%free()
  
end subroutine psb_d_mv_to_l

subroutine psb_d_cp_to_l(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_d_cp_to_l
  implicit none 
  
  class(psb_dspmat_type), intent(in) :: a
  class(psb_ldspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  if (allocated(a%a)) then
    if (.not.allocated(b%a)) allocate(psb_ld_csr_sparse_mat :: b%a, stat=info)
    call a%a%cp_to_lfmt(b%a,info)
  else
    call b%free()
  end if
  
end subroutine psb_d_cp_to_l


!
!
! ld versions
!


subroutine  psb_ld_set_nrows(m,a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_nrows
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_lpk_), intent(in) :: m
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='set_nrows'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_nrows(m)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_nrows


subroutine  psb_ld_set_ncols(n,a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_ncols
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_lpk_), intent(in) :: n
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call a%a%set_ncols(n)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_ncols



!
!  Valid values for DUPL: 
!  psb_dupl_ovwrt_ 
!  psb_dupl_add_   
!  psb_dupl_err_   
!

subroutine  psb_ld_set_dupl(n,a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_dupl
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in) :: n
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_dupl(n)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_dupl


!
! Set the STATE of the internal matrix object
!

subroutine  psb_ld_set_null(a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_null
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_null()

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_null


subroutine  psb_ld_set_bld(a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_bld
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_bld()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_bld


subroutine  psb_ld_set_upd(a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_upd
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_upd()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return


end subroutine psb_ld_set_upd


subroutine  psb_ld_set_asb(a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_asb
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_asb()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_asb


subroutine psb_ld_set_sorted(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_sorted
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_sorted(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_sorted


subroutine psb_ld_set_triangle(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_triangle
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_triangle(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_triangle


subroutine psb_ld_set_unit(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_unit
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_unit(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_unit


subroutine psb_ld_set_lower(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_lower
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_lower(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_lower


subroutine psb_ld_set_upper(a,val) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_set_upper
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%set_upper(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_set_upper



! == ===================================
!
!
!
! Data management
!
!
!
!
!
! == ===================================  


subroutine psb_ld_sparse_print(iout,a,iv,head,ivr,ivc)
  use psb_d_mat_mod, psb_protect_name => psb_ld_sparse_print
  use psb_error_mod
  implicit none 

  integer(psb_ipk_), intent(in)               :: iout
  class(psb_ldspmat_type), intent(in) :: a   
  integer(psb_lpk_), intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%print(iout,iv,head,ivr,ivc)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_sparse_print


subroutine psb_ld_n_sparse_print(fname,a,iv,head,ivr,ivc)
  use psb_d_mat_mod, psb_protect_name => psb_ld_n_sparse_print
  use psb_error_mod
  implicit none 

  character(len=*), intent(in)  :: fname   
  class(psb_ldspmat_type), intent(in) :: a   
  integer(psb_lpk_), intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act, info, iout
  logical :: isopen
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  iout = max(psb_inp_unit,psb_err_unit,psb_out_unit) + 1
  do 
    inquire(unit=iout, opened=isopen)
    if (.not.isopen) exit
    iout = iout + 1
    if (iout > 99) exit
  end do
  if (iout > 99) then 
    write(psb_err_unit,*) 'Error: could not find a free unit for I/O'
    return
  end if
  open(iout,file=fname,iostat=info)
  if (info == psb_success_) then 
    call a%a%print(iout,iv,head,ivr,ivc)
    close(iout)
  else
    write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
  end if

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_n_sparse_print


subroutine psb_ld_get_neigh(a,idx,neigh,n,info,lev)
  use psb_d_mat_mod, psb_protect_name => psb_ld_get_neigh
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(in) :: a   
  integer(psb_lpk_), intent(in)                :: idx 
  integer(psb_lpk_), intent(out)               :: n   
  integer(psb_lpk_), allocatable, intent(out)  :: neigh(:)
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_lpk_), optional, intent(in)      :: lev 

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='get_neigh'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%get_neigh(idx,neigh,n,info,lev)

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_get_neigh



subroutine psb_ld_csall(nr,nc,a,info,nz) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_csall
  use psb_d_base_mat_mod
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_lpk_), intent(in)             :: nr,nc
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_lpk_), intent(in), optional   :: nz

  integer(psb_ipk_) :: err_act 
  character(len=20)  :: name='csall'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)

  call a%free()

  info = psb_success_
  allocate(psb_ld_coo_sparse_mat :: a%a, stat=info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  call a%a%allocate(nr,nc,nz)
  call a%set_bld() 

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_csall


subroutine  psb_ld_reallocate_nz(nz,a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_reallocate_nz
  use psb_error_mod
  implicit none 
  integer(psb_lpk_), intent(in) :: nz
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%reallocate(nz)

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_reallocate_nz


subroutine  psb_ld_free(a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_free
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a

  if (allocated(a%a)) then 
    call a%a%free()
    deallocate(a%a) 
  endif

end subroutine psb_ld_free


subroutine  psb_ld_trim(a) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_trim
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%trim()

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_trim



subroutine psb_ld_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_csput_a
  use psb_d_base_mat_mod
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer(psb_lpk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_lpk_), intent(in), optional   :: gtl(:)

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csput_a'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.(a%is_bld().or.a%is_upd())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info,gtl) 
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_csput_a

subroutine psb_ld_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_d_mat_mod, psb_protect_name => psb_ld_csput_v
  use psb_d_base_mat_mod
  use psb_d_vect_mod, only : psb_d_vect_type
  use psb_l_vect_mod, only : psb_l_vect_type
  use psb_error_mod
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  type(psb_d_vect_type), intent(inout)  :: val
  type(psb_l_vect_type), intent(inout)  :: ia, ja
  integer(psb_lpk_), intent(in)             :: nz, imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_lpk_), intent(in), optional   :: gtl(:)

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csput_v'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.(a%is_bld().or.a%is_upd())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  
  if (allocated(val%v).and.allocated(ia%v).and.allocated(ja%v)) then
    call a%a%csput(nz,ia%v,ja%v,val%v,imin,imax,jmin,jmax,info,gtl) 
  else
    info = psb_err_invalid_mat_state_
  endif

  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_csput_v


subroutine psb_ld_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_csgetptn
  implicit none

  class(psb_ldspmat_type), intent(in) :: a
  integer(psb_lpk_), intent(in)                  :: imin,imax
  integer(psb_lpk_), intent(out)                 :: nz
  integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_lpk_), intent(in), optional        :: iren(:)
  integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csget(imin,imax,nz,ia,ja,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_csgetptn


subroutine psb_ld_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_csgetrow
  implicit none

  class(psb_ldspmat_type), intent(in) :: a
  integer(psb_lpk_), intent(in)                  :: imin,imax
  integer(psb_lpk_), intent(out)                 :: nz
  integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_lpk_), intent(in), optional        :: iren(:)
  integer(psb_lpk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  call a%a%csget(imin,imax,nz,ia,ja,val,info,&
       & jmin,jmax,iren,append,nzin,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_csgetrow




subroutine psb_ld_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_csgetblk
  implicit none

  class(psb_ldspmat_type), intent(in)    :: a
  class(psb_ldspmat_type), intent(inout) :: b
  integer(psb_lpk_), intent(in)                  :: imin,imax
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_lpk_), intent(in), optional        :: iren(:)
  integer(psb_lpk_), intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.
  logical            :: append_
  type(psb_ld_coo_sparse_mat), allocatable  :: acoo


  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (present(append))  then 
    append_ = append
  else
    append_ = .false.
  end if

  allocate(acoo,stat=info)    
  if (append_.and.(info==psb_success_)) then 
    if (allocated(b%a)) &
         & call b%a%mv_to_coo(acoo,info)
  end if

  if (info == psb_success_) then 
    call a%a%csget(imin,imax,acoo,info,&
         & jmin,jmax,iren,append,rscale,cscale)
  else
    info = psb_err_alloc_dealloc_
  end if
  if (info == psb_success_) call move_alloc(acoo,b%a)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_csgetblk


subroutine psb_ld_tril(a,l,info,diag,imin,imax,&
     & jmin,jmax,rscale,cscale,u)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_tril
  implicit none 
  class(psb_ldspmat_type), intent(in)      :: a
  class(psb_ldspmat_type), intent(inout)   :: l
  integer(psb_ipk_),intent(out)           :: info
  integer(psb_lpk_), intent(in), optional :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional           :: rscale,cscale
  class(psb_ldspmat_type), optional, intent(inout)   :: u

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='tril'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat), allocatable  :: lcoo, ucoo

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(lcoo,stat=info)    
  call l%free()
  if (present(u)) then
    if (info == psb_success_) allocate(ucoo,stat=info)    
    call u%free()
    if (info == psb_success_) call a%a%tril(lcoo,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,ucoo)
    if (info == psb_success_) call move_alloc(ucoo,u%a)
    if (info == psb_success_) call u%cscnv(info,mold=a%a)
  else
    if (info == psb_success_) then 
      call a%a%tril(lcoo,info,diag,imin,imax,&
           & jmin,jmax,rscale,cscale)
    else
      info = psb_err_alloc_dealloc_
    end if
  end if
  if (info == psb_success_) call move_alloc(lcoo,l%a)
  if (info == psb_success_) call l%cscnv(info,mold=a%a)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return


end subroutine psb_ld_tril

subroutine psb_ld_triu(a,u,info,diag,imin,imax,&
     & jmin,jmax,rscale,cscale,l)
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_triu
  implicit none 
  class(psb_ldspmat_type), intent(in)      :: a
  class(psb_ldspmat_type), intent(inout)   :: u
  integer(psb_ipk_),intent(out)           :: info
  integer(psb_lpk_), intent(in), optional :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional           :: rscale,cscale
  class(psb_ldspmat_type), optional, intent(inout)   :: l

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='triu'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat), allocatable  :: lcoo, ucoo


  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(ucoo,stat=info)    
  call u%free()

  if (present(l)) then
    if (info == psb_success_) allocate(lcoo,stat=info)    
    call l%free()
    if (info == psb_success_) call a%a%triu(ucoo,info,diag,imin,imax,&
         & jmin,jmax,rscale,cscale,lcoo)
    if (info == psb_success_) call move_alloc(lcoo,l%a)
    if (info == psb_success_) call l%cscnv(info,mold=a%a)
  else
    if (info == psb_success_) then 
      call a%a%triu(ucoo,info,diag,imin,imax,&
           & jmin,jmax,rscale,cscale)
    else
      info = psb_err_alloc_dealloc_
    end if
  end if
  if (info == psb_success_) call move_alloc(ucoo,u%a)
  if (info == psb_success_) call u%cscnv(info,mold=a%a)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return


end subroutine psb_ld_triu


subroutine psb_ld_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_csclip
  implicit none

  class(psb_ldspmat_type), intent(in) :: a
  class(psb_ldspmat_type), intent(inout) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_lpk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_ld_coo_sparse_mat), allocatable  :: acoo

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)    
  call b%free()
  if (info == psb_success_) then 
    call a%a%csclip(acoo,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else
    info = psb_err_alloc_dealloc_
  end if
 
  if (info == psb_success_) call move_alloc(acoo,b%a)
  if (info /= psb_success_) goto 9999 
  
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_csclip


subroutine psb_ld_b_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_d_base_mat_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_b_csclip
  implicit none

  class(psb_ldspmat_type), intent(in) :: a
  type(psb_ld_coo_sparse_mat), intent(out) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_lpk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%csclip(b,info,&
       & imin,imax,jmin,jmax,rscale,cscale)
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_b_csclip




subroutine psb_ld_cscnv(a,b,info,type,mold,upd,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cscnv
  implicit none 
  class(psb_ldspmat_type), intent(in)      :: a
  class(psb_ldspmat_type), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psb_ld_base_sparse_mat), intent(in), optional :: mold


  class(psb_ld_base_sparse_mat), allocatable  :: altmp
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call b%free()
  if (count( (/present(mold),present(type) /)) > 1) then
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='TYPE, MOLD')
    goto 9999
  end if

  if (present(mold)) then 

    allocate(altmp, mold=mold,stat=info) 

  else if (present(type)) then 

    select case (psb_toupper(type))
    case ('CSR')
      allocate(psb_ld_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psb_ld_coo_sparse_mat :: altmp, stat=info) 
    case ('CSC')
      allocate(psb_ld_csc_sparse_mat :: altmp, stat=info) 
    case default
      info = psb_err_format_unknown_ 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(altmp, mold=psb_get_mat_default(a),stat=info) 
  end if

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  
  if (present(dupl)) then 
    call altmp%set_dupl(dupl)
  else if (a%is_bld()) then 
    ! Does this make sense at all?? Who knows..
    call altmp%set_dupl(psb_dupl_def_)
  end if

  if (debug) write(psb_err_unit,*) 'Converting from ',&
       & a%get_fmt(),' to ',altmp%get_fmt()

  call altmp%cp_from_fmt(a%a, info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call move_alloc(altmp,b%a)
  call b%trim()
  call b%asb() 
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_cscnv



subroutine psb_ld_cscnv_ip(a,info,type,mold,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cscnv_ip
  implicit none 

  class(psb_ldspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl
  character(len=*), optional, intent(in) :: type
  class(psb_ld_base_sparse_mat), intent(in), optional :: mold


  class(psb_ld_base_sparse_mat), allocatable  :: altmp
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cscnv_ip'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(dupl)) then 
    call a%set_dupl(dupl)
  else if (a%is_bld()) then 
    call a%set_dupl(psb_dupl_def_)
  end if

  if (count( (/present(mold),present(type) /)) > 1) then
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='TYPE, MOLD')
    goto 9999
  end if

  if (present(mold)) then 

    allocate(altmp, mold=mold,stat=info) 

  else if (present(type)) then 

    select case (psb_toupper(type))
    case ('CSR')
      allocate(psb_ld_csr_sparse_mat :: altmp, stat=info) 
    case ('COO')
      allocate(psb_ld_coo_sparse_mat :: altmp, stat=info) 
    case ('CSC')
      allocate(psb_ld_csc_sparse_mat :: altmp, stat=info) 
    case default
      info = psb_err_format_unknown_ 
      call psb_errpush(info,name,a_err=type)
      goto 9999
    end select
  else
    allocate(altmp, mold=psb_get_mat_default(a),stat=info) 
  end if

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug) write(psb_err_unit,*) 'Converting in-place from ',&
       & a%get_fmt(),' to ',altmp%get_fmt()

  call altmp%mv_from_fmt(a%a, info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call move_alloc(altmp,a%a)
  call a%set_asb() 
  call a%trim()
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_cscnv_ip



subroutine psb_ld_cscnv_base(a,b,info,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cscnv_base
  implicit none 
  class(psb_ldspmat_type), intent(in)       :: a
  class(psb_ld_base_sparse_mat), intent(out) :: b
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl


  type(psb_ld_coo_sparse_mat)  :: altmp
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%cp_to_coo(altmp,info )
  if ((info == psb_success_).and.present(dupl)) then 
    call altmp%set_dupl(dupl)
  end if
  call altmp%fix(info)
  if (info == psb_success_) call altmp%trim()
  if (info == psb_success_) call altmp%set_asb() 
  if (info == psb_success_) call b%mv_from_coo(altmp,info)

  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_cscnv_base



!!$subroutine psb_ld_clip_d(a,b,info)
!!$  ! Output is always in  COO format 
!!$  use psb_error_mod
!!$  use psb_const_mod
!!$  use psb_d_base_mat_mod
!!$  use psb_d_mat_mod, psb_protect_name => psb_ld_clip_d
!!$  implicit none
!!$
!!$  class(psb_ldspmat_type), intent(in)    :: a
!!$  class(psb_ldspmat_type), intent(inout) :: b
!!$  integer(psb_ipk_),intent(out)                  :: info
!!$
!!$  integer(psb_ipk_) :: err_act
!!$  character(len=20)  :: name='clip_diag'
!!$  logical, parameter :: debug=.false.
!!$  type(psb_ld_coo_sparse_mat), allocatable  :: acoo
!!$  integer(psb_lpk_) :: i, j, nz
!!$
!!$  info = psb_success_
!!$  call psb_erractionsave(err_act)
!!$  if (a%is_null()) then 
!!$    info = psb_err_invalid_mat_state_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  endif
!!$
!!$  allocate(acoo,stat=info)    
!!$  if (info == psb_success_) call a%a%cp_to_coo(acoo,info)
!!$  if (info /= psb_success_) then 
!!$    info = psb_err_alloc_dealloc_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  endif
!!$
!!$  nz = acoo%get_nzeros()
!!$  j = 0
!!$  do i=1, nz
!!$    if (acoo%ia(i) /= acoo%ja(i)) then 
!!$      j = j + 1 
!!$      acoo%ia(j)  = acoo%ia(i)
!!$      acoo%ja(j)  = acoo%ja(i)
!!$      acoo%val(j) = acoo%val(i)
!!$    end if
!!$  end do
!!$  call acoo%set_nzeros(j)
!!$  call acoo%trim()
!!$  call b%mv_from(acoo)
!!$
!!$  call psb_erractionrestore(err_act)
!!$  return
!!$
!!$
!!$9999 call psb_error_handler(err_act)
!!$
!!$  return
!!$
!!$end subroutine psb_ld_clip_d
!!$
!!$
!!$
!!$subroutine psb_ld_clip_d_ip(a,info)
!!$  ! Output is always in  COO format 
!!$  use psb_error_mod
!!$  use psb_const_mod
!!$  use psb_d_base_mat_mod
!!$  use psb_d_mat_mod, psb_protect_name => psb_ld_clip_d_ip
!!$  implicit none
!!$
!!$  class(psb_ldspmat_type), intent(inout) :: a
!!$  integer(psb_ipk_),intent(out)                  :: info
!!$
!!$  integer(psb_ipk_) :: err_act
!!$  character(len=20)  :: name='clip_diag'
!!$  logical, parameter :: debug=.false.
!!$  type(psb_ld_coo_sparse_mat), allocatable  :: acoo
!!$  integer(psb_lpk_) :: i, j, nz
!!$
!!$  info = psb_success_
!!$  call psb_erractionsave(err_act)
!!$  if (a%is_null()) then 
!!$    info = psb_err_invalid_mat_state_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  endif
!!$
!!$  allocate(acoo,stat=info)    
!!$  if (info == psb_success_) call a%a%mv_to_coo(acoo,info)
!!$  if (info /= psb_success_) then
!!$    info = psb_err_alloc_dealloc_
!!$    call psb_errpush(info,name)
!!$    goto 9999
!!$  endif
!!$
!!$  nz = acoo%get_nzeros()
!!$  j = 0
!!$  do i=1, nz
!!$    if (acoo%ia(i) /= acoo%ja(i)) then 
!!$      j = j + 1 
!!$      acoo%ia(j)  = acoo%ia(i)
!!$      acoo%ja(j)  = acoo%ja(i)
!!$      acoo%val(j) = acoo%val(i)
!!$    end if
!!$  end do
!!$  call acoo%set_nzeros(j)
!!$  call acoo%trim()
!!$  call a%mv_from(acoo)
!!$
!!$  call psb_erractionrestore(err_act)
!!$  return
!!$
!!$
!!$9999 call psb_error_handler(err_act)
!!$
!!$  return
!!$
!!$end subroutine psb_ld_clip_d_ip
!!$

subroutine psb_ld_mv_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_mv_from
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  call a%free()
  allocate(a%a,mold=b, stat=info)
  call a%a%mv_from_fmt(b,info)
  call b%free()

  return
end subroutine psb_ld_mv_from


subroutine psb_ld_cp_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cp_from
  implicit none 
  class(psb_ldspmat_type), intent(out)      :: a
  class(psb_ld_base_sparse_mat), intent(in) :: b
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='cp_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  
  call a%free()
  !
  ! Note: it is tempting to use SOURCE allocation below;
  ! however this would run the risk of messing up with data
  ! allocated externally (e.g. GPU-side data).
  !
  allocate(a%a,mold=b,stat=info)
  if (info /= psb_success_) info = psb_err_alloc_dealloc_
  if (info == psb_success_) call a%a%cp_from_fmt(b, info)    
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return
end subroutine psb_ld_cp_from


subroutine psb_ld_mv_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_mv_to
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  call b%mv_from_fmt(a%a,info)

  return
end subroutine psb_ld_mv_to


subroutine psb_ld_cp_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cp_to
  implicit none 
  class(psb_ldspmat_type), intent(in) :: a
  class(psb_ld_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  call b%cp_from_fmt(a%a,info)

  return
end subroutine psb_ld_cp_to

subroutine psb_ld_mold(a,b)
  use psb_d_mat_mod, psb_protect_name => psb_ld_mold
  class(psb_ldspmat_type), intent(inout)     :: a
  class(psb_ld_base_sparse_mat), allocatable, intent(out) :: b
  integer(psb_ipk_) :: info

  allocate(b,mold=a%a, stat=info)
  
end subroutine psb_ld_mold

subroutine psb_ldspmat_type_move(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ldspmat_type_move
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_ldspmat_type), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                   :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='move_alloc'
  logical, parameter :: debug=.false.

  info = psb_success_
  call b%free()
  call move_alloc(a%a,b%a)

  return
end subroutine psb_ldspmat_type_move


subroutine psb_ldspmat_clone(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ldspmat_clone
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_ldspmat_type), intent(inout) :: b
  integer(psb_ipk_), intent(out)        :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='clone'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  call b%free()
  if (allocated(a%a)) then 
    call a%a%clone(b%a,info)
  end if
  if (info /= psb_success_) goto 9999 

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ldspmat_clone


subroutine psb_ld_transp_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_transp_1mat
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%transp()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_transp_1mat



subroutine psb_ld_transp_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_transp_2mat
  implicit none 
  class(psb_ldspmat_type), intent(in)  :: a
  class(psb_ldspmat_type), intent(inout) :: b

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call b%free()
  allocate(b%a,mold=a%a,stat=info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    goto 9999
  end if
  call a%a%transp(b%a)    

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_transp_2mat


subroutine psb_ld_transc_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_transc_1mat
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%transc()

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_transc_1mat



subroutine psb_ld_transc_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_transc_2mat
  implicit none 
  class(psb_ldspmat_type), intent(in)    :: a
  class(psb_ldspmat_type), intent(inout) :: b

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.


  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call b%free()
  allocate(b%a,mold=a%a,stat=info)
  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_
    goto 9999
  end if
  call a%a%transc(b%a)    

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_transc_2mat


subroutine psb_ld_asb(a,mold)
  use psb_d_mat_mod, psb_protect_name => psb_ld_asb
  use psb_error_mod
  implicit none 

  class(psb_ldspmat_type), intent(inout) :: a   
  class(psb_ld_base_sparse_mat), optional, intent(in) :: mold
  class(psb_ld_base_sparse_mat), allocatable :: tmp
  class(psb_ld_base_sparse_mat), pointer :: mld
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='ld_asb'

  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%asb()
  if (present(mold)) then 
    if (.not.same_type_as(a%a,mold)) then 
      allocate(tmp,mold=mold)
      call tmp%mv_from_fmt(a%a,info)
      call a%a%free()
      call move_alloc(tmp,a%a)
    end if
  else
    mld => psb_ld_get_base_mat_default()
    if (.not.same_type_as(a%a,mld)) &
         & call a%cscnv(info)
  end if
  

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_asb

subroutine psb_ld_reinit(a,clear)
  use psb_d_mat_mod, psb_protect_name => psb_ld_reinit
  use psb_error_mod
  implicit none 

  class(psb_ldspmat_type), intent(inout) :: a   
  logical, intent(in), optional :: clear
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='reinit'

  call psb_erractionsave(err_act)
  if (a%is_null()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (a%a%has_update()) then 
    call a%a%reinit(clear)
  else
    info = psb_err_missing_override_method_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_reinit




function psb_ld_get_diag(a,info) result(d)
  use psb_d_mat_mod, psb_protect_name => psb_ld_get_diag
  use psb_error_mod
  use psb_const_mod
  implicit none 
  class(psb_ldspmat_type), intent(in) :: a
  real(psb_dpk_), allocatable         :: d(:)
  integer(psb_ipk_), intent(out)       :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(d(max(1,min(a%a%get_nrows(),a%a%get_ncols()))), stat=info)
  if (info /= 0) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info,name)
    goto 9999
  end if
  call a%a%get_diag(d,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_ld_get_diag


subroutine psb_ld_scal(d,a,info,side)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_scal
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  real(psb_dpk_), intent(in)              :: d(:)
  integer(psb_ipk_), intent(out)                    :: info
  character, intent(in), optional :: side

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scal(d,info,side=side)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_scal


subroutine psb_ld_scals(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_scals
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  real(psb_dpk_), intent(in)              :: d
  integer(psb_ipk_), intent(out)                    :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scal(d,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_ld_scals

subroutine psb_ld_mv_from_ib(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_mv_from_ib
  implicit none 
  
  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  if (.not.allocated(a%a)) allocate(psb_ld_csr_sparse_mat :: a%a, stat=info)
  if (info == psb_success_) call a%a%mv_from_ifmt(b,info)
  
end subroutine psb_ld_mv_from_ib
  
subroutine psb_ld_cp_from_ib(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cp_from_ib
  implicit none 
  
  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (.not.allocated(a%a)) allocate(psb_ld_csr_sparse_mat :: a%a, stat=info)
  if (info == psb_success_) call a%a%cp_from_ifmt(b,info)

end subroutine psb_ld_cp_from_ib

subroutine psb_ld_mv_to_ib(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_mv_to_ib
  implicit none 
  
  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  if (.not.allocated(a%a)) then
    call b%free()
  else
    call a%a%mv_to_ifmt(b,info)
    call a%free()
  end if
  
end subroutine psb_ld_mv_to_ib

subroutine psb_ld_cp_to_ib(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cp_to_ib
  implicit none 
  class(psb_ldspmat_type), intent(in) :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  if (.not.allocated(a%a)) then
    call b%free()
  else
    call a%a%cp_to_ifmt(b,info)
  end if

end subroutine psb_ld_cp_to_ib

subroutine psb_ld_mv_from_i(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_mv_from_i
  implicit none 
  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_dspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (allocated(b%a)) then
    if (.not.allocated(a%a)) allocate(psb_ld_csr_sparse_mat :: a%a, stat=info)
    call a%a%mv_from_ifmt(b%a,info)
  else
    call a%free()
  end if
  call b%free()
  
end subroutine psb_ld_mv_from_i

  
subroutine psb_ld_cp_from_i(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cp_from_i
  implicit none 

  class(psb_ldspmat_type), intent(out) :: a
  class(psb_dspmat_type), intent(in) :: b
  integer(psb_ipk_) :: info
  
  if (allocated(b%a)) then
    if (.not.allocated(a%a)) allocate(psb_ld_csr_sparse_mat :: a%a, stat=info)
    call a%a%cp_from_ifmt(b%a,info)
  else
    call a%free()
  end if
end subroutine psb_ld_cp_from_i

subroutine psb_ld_mv_to_i(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_mv_to_i
  implicit none 

  class(psb_ldspmat_type), intent(inout) :: a
  class(psb_dspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  if (allocated(a%a)) then
    if (.not.allocated(b%a)) allocate(psb_d_csr_sparse_mat :: b%a, stat=info)
    call a%a%mv_to_ifmt(b%a,info)
  else
    call b%free()
  end if
  call a%free()
  
end subroutine psb_ld_mv_to_i

subroutine psb_ld_cp_to_i(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_d_mat_mod, psb_protect_name => psb_ld_cp_to_i
  implicit none 
  
  class(psb_ldspmat_type), intent(in) :: a
  class(psb_dspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info
  
  if (allocated(a%a)) then
    if (.not.allocated(b%a)) allocate(psb_d_csr_sparse_mat :: b%a, stat=info)
    call a%a%cp_to_ifmt(b%a,info)
  else
    call b%free()
  end if
  
end subroutine psb_ld_cp_to_i



  
