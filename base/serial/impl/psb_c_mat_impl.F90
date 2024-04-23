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
!  c_mat_impl:
!   implementation of the outer matrix methods.
!   Most of the methods rely on the STATE design pattern:
!   the inner class(psb_c_base_sparse_mat) is responsbile
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


subroutine  psb_c_set_nrows(m,a)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_nrows
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in) :: m
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='set_nrows'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_nrows(m)
  else if (allocated(a%ad)) then
    call a%ad%set_nrows(m)
    call a%and%set_nrows(m)
  else 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_nrows


subroutine  psb_c_set_ncols(n,a)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_ncols
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in) :: n
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.
  integer(psb_ipk_) :: nr

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_ncols(n)
  else if (allocated(a%ad)) then
    nr = a%get_nrows()
    call a%ad%set_ncols(nr)
    call a%and%set_ncols(max(0,n-nr))
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_ncols



!
!  Valid values for DUPL:
!  psb_dupl_ovwrt_
!  psb_dupl_add_
!  psb_dupl_err_
!

subroutine  psb_c_set_dupl(n,a)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_dupl
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in) :: n
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_dupl(n)
  else if (allocated(a%ad)) then
    call a%ad%set_dupl(n)
    call a%and%set_dupl(n)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_dupl


!
! Set the STATE of the internal matrix object
!

subroutine  psb_c_set_null(a)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_null
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_null()
  else if (allocated(a%ad)) then
    call a%ad%set_null()
    call a%and%set_null()
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_null


subroutine  psb_c_set_bld(a)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_bld
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_bld
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_bld


subroutine  psb_c_set_upd(a)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_upd
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_upd()
  else if (allocated(a%ad)) then
    call a%ad%set_upd()
    call a%and%set_upd()
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
end subroutine psb_c_set_upd

subroutine  psb_c_set_asb(a)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_asb
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_asb()
  else if (allocated(a%ad)) then
    call a%ad%set_asb()
    call a%and%set_asb()
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_asb


subroutine psb_c_set_sorted(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_sorted
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_sorted(val)
  else if (allocated(a%ad)) then
    call a%ad%set_sorted(val)
    call a%and%set_sorted(val)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_sorted


subroutine psb_c_set_triangle(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_triangle
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_triangle(val)
  else if (allocated(a%ad)) then
    call a%ad%set_triangle(val)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_triangle

subroutine psb_c_set_symmetric(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_symmetric
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_symmetric(val)
  else if (allocated(a%ad)) then
    call a%ad%set_symmetric(val)    
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_symmetric

subroutine psb_c_set_unit(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_unit
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_unit(val)
  else if (allocated(a%ad)) then
    call a%ad%set_unit(val)    
  else      

    call psb_errpush(info,name)
    goto 9999
  endif
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_unit

subroutine psb_c_set_lower(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_lower
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_lower(val)
  else if (allocated(a%ad)) then
    call a%ad%set_lower(val)    
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_lower


subroutine psb_c_set_upper(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_c_set_upper
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: val
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='get_nzeros'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    call a%a%set_lower(val)
  else if (allocated(a%ad)) then
    call a%ad%set_lower(val)    
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_set_upper



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
subroutine psb_c_sparse_print(iout,a,iv,head,ivr,ivc)
  use psb_c_mat_mod, psb_protect_name => psb_c_sparse_print
  use psb_error_mod
  implicit none

  integer(psb_ipk_), intent(in)               :: iout
  class(psb_cspmat_type), intent(in) :: a
  integer(psb_lpk_), intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat) :: acoo, ac1,ac2
  
  info = psb_success_
  call psb_get_erraction(err_act)
  if (allocated(a%a)) then
    call a%a%print(iout,iv,head,ivr,ivc)    
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,acoo,a%get_nrows(),a%get_ncols(),info)
    call acoo%print(iout,iv,head,ivr,ivc)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_sparse_print


subroutine psb_c_n_sparse_print(fname,a,iv,head,ivr,ivc)
  use psb_c_mat_mod, psb_protect_name => psb_c_n_sparse_print
  use psb_error_mod
  implicit none

  character(len=*), intent(in)  :: fname
  class(psb_cspmat_type), intent(in) :: a
  integer(psb_lpk_), intent(in), optional     :: iv(:)
  character(len=*), optional        :: head
  integer(psb_lpk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act, info, iout
  logical :: isopen
  character(len=20)  :: name='sparse_print'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_get_erraction(err_act)
  
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
    call a%print(iout,iv,head,ivr,ivc)
    close(iout)
  else
    write(psb_err_unit,*) 'Error: could not open ',fname,' for output'
  end if

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_n_sparse_print

subroutine psb_c_get_neigh(a,idx,neigh,n,info,lev)
  use psb_c_mat_mod, psb_protect_name => psb_c_get_neigh
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  integer(psb_ipk_), intent(in)                :: idx
  integer(psb_ipk_), intent(out)               :: n
  integer(psb_ipk_), allocatable, intent(out)  :: neigh(:)
  integer(psb_ipk_), intent(out)               :: info
  integer(psb_ipk_), optional, intent(in)      :: lev

  integer(psb_ipk_) :: err_act, n1
  character(len=20)  :: name='get_neigh'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)

  if (allocated(a%a)) then
    call a%a%get_neigh(idx,neigh,n,info,lev)
  else if (allocated(a%ad)) then
    call a%ad%get_neigh(idx,neigh,n1,info,lev)
    call a%ad%get_neigh(idx,neigh,n,info,lev,nin=n1)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_get_neigh



subroutine psb_c_csall(nr,nc,a,info,nz,type,mold)
  use psb_c_mat_mod, psb_protect_name => psb_c_csall
  use psb_c_base_mat_mod
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in)             :: nr,nc
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: nz
  character(len=*), intent(in), optional    :: type
  class(psb_c_base_sparse_mat), optional, intent(in) :: mold

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csall'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)

  call a%free()

  info = psb_success_
  if (present(mold)) then
    allocate(a%a, stat=info, mold=mold)
  else if (present(type)) then
    select case (type)
    case('CSR')
      allocate(psb_c_csr_sparse_mat :: a%a, stat=info)
    case('COO')
      allocate(psb_c_coo_sparse_mat :: a%a, stat=info)
    case('CSC')
      allocate(psb_c_csc_sparse_mat :: a%a, stat=info)
    case default
      allocate(psb_c_coo_sparse_mat :: a%a, stat=info)
    end select
  else
    allocate(psb_c_coo_sparse_mat :: a%a, stat=info)
  end if

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

end subroutine psb_c_csall


subroutine  psb_c_reallocate_nz(nz,a)
  use psb_c_mat_mod, psb_protect_name => psb_c_reallocate_nz
  use psb_error_mod
  implicit none
  integer(psb_ipk_), intent(in) :: nz
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (allocated(a%a)) then
    call a%a%reallocate(nz)
  else if (allocated(a%ad)) then
    call a%ad%reallocate(nz)
    call a%and%reallocate(nz)    
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_reallocate_nz


subroutine  psb_c_free(a)
  use psb_c_mat_mod, psb_protect_name => psb_c_free
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a

  if (allocated(a%a)) then
    call a%a%free()
    deallocate(a%a)
  endif
  if (allocated(a%ad)) then
    call a%ad%free()
    deallocate(a%ad)
  endif
  if (allocated(a%and)) then
    call a%and%free()
    deallocate(a%and)
  endif
  if (allocated(a%rmta)) then
    call a%rmta%free()
    deallocate(a%rmta)
  end if
  a%remote_build = psb_matbld_noremote_
  
end subroutine psb_c_free


subroutine  psb_c_trim(a)
  use psb_c_mat_mod, psb_protect_name => psb_c_trim
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  if (allocated(a%a)) then
    call a%a%trim()
  else if (allocated(a%ad)) then
    call a%ad%trim()
    call a%and%trim()    
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_trim

subroutine psb_c_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
  use psb_c_mat_mod, psb_protect_name => psb_c_csput_a
  use psb_c_base_mat_mod
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)      :: val(:)
  integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info

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
  
  if (allocated(a%a)) then 
    call a%a%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info)
  else if (allocated(a%ad)) then
    call a%ad%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info)
    call a%and%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
    
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csput_a

subroutine psb_c_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
  use psb_c_mat_mod, psb_protect_name => psb_c_csput_v
  use psb_c_base_mat_mod
  use psb_c_vect_mod, only : psb_c_vect_type
  use psb_i_vect_mod, only : psb_i_vect_type
  use psb_error_mod
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  type(psb_c_vect_type), intent(inout)  :: val
  type(psb_i_vect_type), intent(inout)  :: ia, ja
  integer(psb_ipk_), intent(in)             :: nz, imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info

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
    call a%csput(nz,ia%v%v,ja%v%v,val%v%v,imin,imax,jmin,jmax,info)
  else
    info = psb_err_invalid_mat_state_
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csput_v

subroutine psb_c_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csgetptn
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act, nz1
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (allocated(a%a)) then
    call a%a%csget(imin,imax,nz,ia,ja,info,&
         & jmin=jmin,jmax=jmax,iren=iren,append=append,nzin=nzin,&
         & rscale=rscale,cscale=cscale)

  else if (allocated(a%ad)) then
    call a%ad%csget(imin,imax,nz1,ia,ja,info,&
         & jmin=jmin,jmax=jmax,iren=iren,append=append,nzin=nzin,&
         & rscale=rscale,cscale=cscale)
    call a%and%csget(imin,imax,nz,ia,ja,info,&
         & jmin=jmin,jmax=jmax,iren=iren,append=.true.,nzin=nz1,&
         & rscale=rscale,cscale=cscale)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csgetptn


subroutine psb_c_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csgetrow
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale,chksz

  integer(psb_ipk_) :: err_act, nz1
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    call a%a%csget(imin,imax,nz,ia,ja,val,info,&
         & jmin=jmin,jmax=jmax,iren=iren,append=append,nzin=nzin,&
         & rscale=rscale,cscale=cscale,chksz=chksz)
  else if (allocated(a%ad)) then
    call a%ad%csget(imin,imax,nz1,ia,ja,val,info,&
         & jmin=jmin,jmax=jmax,iren=iren,append=append,nzin=nzin,&
         & rscale=rscale,cscale=cscale,chksz=chksz)
    call a%and%csget(imin,imax,nz,ia,ja,val,info,&
         & jmin=jmin,jmax=jmax,iren=iren,append=.true.,nzin=nz1,&
         & rscale=rscale,cscale=cscale,chksz=chksz)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csgetrow


subroutine psb_c_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csgetblk
  implicit none

  class(psb_cspmat_type), intent(in)    :: a
  class(psb_cspmat_type), intent(inout) :: b
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
  type(psb_c_coo_sparse_mat), allocatable  :: acoo


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
    if (allocated(a%a)) then
      call a%a%csget(imin,imax,acoo,info,&
           & jmin=jmin,jmax=jmax,iren=iren,append=append,&
           & rscale=rscale,cscale=cscale)
    else if (allocated(a%ad)) then
      call a%ad%csget(imin,imax,acoo,info,&
           & jmin=jmin,jmax=jmax,iren=iren,append=append,&
           & rscale=rscale,cscale=cscale)
      call a%and%csget(imin,imax,acoo,info,&
           & jmin=jmin,jmax=jmax,iren=iren,append=.true.,&
           & rscale=rscale,cscale=cscale)
    else
      info = psb_err_invalid_mat_state_
      call psb_errpush(info,name)
      goto 9999
    endif      
  else
    info = psb_err_alloc_dealloc_
  end if
  if (info == psb_success_) call move_alloc(acoo,b%a)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csgetblk


subroutine psb_c_tril(a,l,info,diag,imin,imax,&
     & jmin,jmax,rscale,cscale,u)
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_tril
  implicit none
  class(psb_cspmat_type), intent(in)      :: a
  class(psb_cspmat_type), intent(inout)   :: l
  integer(psb_ipk_),intent(out)           :: info
  integer(psb_ipk_), intent(in), optional :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional           :: rscale,cscale
  class(psb_cspmat_type), optional, intent(inout)   :: u

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='tril'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: lcoo, ucoo
  type(psb_c_coo_sparse_mat) :: acoo, ac1,ac2

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  allocate(lcoo,stat=info)
  call l%free()
  if (allocated(a%a)) then

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

  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,acoo,a%get_nrows(),a%get_ncols(),info)

    if (present(u)) then
      if (info == psb_success_) allocate(ucoo,stat=info)
      call u%free()
      if (info == psb_success_) call acoo%tril(lcoo,info,diag,imin,imax,&
           & jmin,jmax,rscale,cscale,ucoo)
      if (info == psb_success_) call move_alloc(ucoo,u%a)
      if (info == psb_success_) call u%cscnv(info,mold=a%a)
    else
      if (info == psb_success_) then
        call acoo%tril(lcoo,info,diag,imin,imax,&
             & jmin,jmax,rscale,cscale)
      else
        info = psb_err_alloc_dealloc_
      end if
    end if
    if (info == psb_success_) call move_alloc(lcoo,l%a)
    if (info == psb_success_) call l%cscnv(info,mold=a%a)
    if (info /= psb_success_) goto 9999

  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return


end subroutine psb_c_tril

subroutine psb_c_triu(a,u,info,diag,imin,imax,&
     & jmin,jmax,rscale,cscale,l)
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_triu
  implicit none
  class(psb_cspmat_type), intent(in)      :: a
  class(psb_cspmat_type), intent(inout)   :: u
  integer(psb_ipk_),intent(out)           :: info
  integer(psb_ipk_), intent(in), optional :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional           :: rscale,cscale
  class(psb_cspmat_type), optional, intent(inout)   :: l

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='triu'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: lcoo, ucoo
  type(psb_c_coo_sparse_mat) :: acoo, ac1,ac2


  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(ucoo,stat=info)
  call u%free()

  if (allocated(a%a)) then

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

  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,acoo,a%get_nrows(),a%get_ncols(),info)

    if (present(l)) then
      if (info == psb_success_) allocate(lcoo,stat=info)
      call l%free()
      if (info == psb_success_) call acoo%triu(ucoo,info,diag,imin,imax,&
           & jmin,jmax,rscale,cscale,lcoo)
      if (info == psb_success_) call move_alloc(lcoo,l%a)
      if (info == psb_success_) call l%cscnv(info,mold=a%a)
    else
      if (info == psb_success_) then
        call acoo%triu(ucoo,info,diag,imin,imax,&
             & jmin,jmax,rscale,cscale)
      else
        info = psb_err_alloc_dealloc_
      end if
    end if
    if (info == psb_success_) call move_alloc(ucoo,u%a)
    if (info == psb_success_) call u%cscnv(info,mold=a%a)
    if (info /= psb_success_) goto 9999
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_triu


subroutine psb_c_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csclip
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  class(psb_cspmat_type), intent(inout) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: acoo
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)
  if (info /= 0) then
    info = psb_err_alloc_dealloc_
     goto 9999
  end if
  call b%free()
  if (allocated(a%a)) then
    call a%a%csclip(acoo,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call aa%csclip(acoo,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (info == psb_success_) call move_alloc(acoo,b%a)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csclip

subroutine psb_c_csclip_ip(a,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csclip_ip
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: acoo
  type(psb_c_coo_sparse_mat)  :: ac1,ac2,aa

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    call a%a%csclip(acoo,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call aa%csclip(acoo,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (info == psb_success_) call a%free()
  if (info == psb_success_) call move_alloc(acoo,a%a)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_csclip_ip

subroutine psb_c_b_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_b_csclip
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  type(psb_c_coo_sparse_mat), intent(out) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    call a%a%csclip(b,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call aa%csclip(b,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_b_csclip

subroutine psb_c_split_nd(a,n_rows,n_cols,info)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_split_nd
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in)           :: n_rows, n_cols
  integer(psb_ipk_), intent(out)          :: info
!!$      integer(psb_ipk_),optional, intent(in)           :: dupl
!!$      character(len=*), optional, intent(in) :: type
!!$      class(psb_c_base_sparse_mat), intent(in), optional :: mold
  type(psb_c_coo_sparse_mat) :: acoo
  type(psb_c_csr_sparse_mat), allocatable :: aclip
  type(psb_c_ecsr_sparse_mat), allocatable :: andclip
  logical, parameter :: use_ecsr=.true.
  character(len=20)     :: name, ch_err
  integer(psb_ipk_) :: err_act

  info = psb_success_
  name = 'psb_split'
  call psb_erractionsave(err_act)
  if (allocated(a%a)) then
    allocate(aclip)
    call a%a%csclip(acoo,info,jmax=n_rows,rscale=.false.,cscale=.false.)
    allocate(a%ad,mold=a%a)
    call a%ad%mv_from_coo(acoo,info)
    call a%a%csclip(acoo,info,jmin=n_rows+1,jmax=n_cols,rscale=.false.,cscale=.false.)
    if (use_ecsr) then
      allocate(andclip)
      call andclip%mv_from_coo(acoo,info)
      call move_alloc(andclip,a%and)
    else
      allocate(a%and,mold=a%a)
      call a%and%mv_from_coo(acoo,info)
    end if
    call a%a%free()
    deallocate(a%a)
  end if
  if (psb_errstatus_fatal()) then    
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='cscnv')
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_split_nd

subroutine psb_c_merge_nd(a,n_rows,n_cols,info,acoo)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_merge_nd
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(in)           :: n_rows, n_cols
  integer(psb_ipk_), intent(out)          :: info
      class(psb_c_coo_sparse_mat), intent(out), optional :: acoo
!!$      integer(psb_ipk_),optional, intent(in)           :: dupl
!!$      character(len=*), optional, intent(in) :: type
!!$      class(psb_c_base_sparse_mat), intent(in), optional :: mold
  type(psb_c_coo_sparse_mat) :: acoo1
  integer(psb_ipk_) :: nz
  logical, parameter :: use_ecsr=.true.
  character(len=20)     :: name, ch_err
  integer(psb_ipk_) :: err_act

  info = psb_success_
  name = 'psb_split'
  call psb_erractionsave(err_act)

  call a%ad%csmerge(a%and,acoo1,n_rows,n_cols,info)

  if (present(acoo)) then
    call acoo%mv_from_coo(acoo1,info)
  else
    call a%ad%free()
    call a%and%free()
    if (allocated(a%a)) then
      call a%a%free()
      deallocate(a%a)
    end if
    allocate(a%a,mold=a%ad)
    call a%a%mv_from_coo(acoo1,info)
    deallocate(a%ad,a%and)
  end if

  if (psb_errstatus_fatal()) then    
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='cscnv')
    goto 9999
  endif
  
  call psb_erractionrestore(err_act)
  return
  
9999 call psb_error_handler(err_act)
  
  return
  
end subroutine psb_c_merge_nd

subroutine psb_c_cscnv(a,b,info,type,mold,upd,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cscnv
  implicit none
  class(psb_cspmat_type), intent(in)      :: a
  class(psb_cspmat_type), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psb_c_base_sparse_mat), intent(in), optional :: mold


  class(psb_c_base_sparse_mat), allocatable  :: altmp
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

  if (allocated(a%a))   call inner_cp_fmt(a%a,b%a,info,type,mold,dupl)
  if (allocated(a%ad))  call inner_cp_fmt(a%ad,b%ad,info,type,mold,dupl)
  if (allocated(a%and)) call inner_cp_fmt(a%and,b%and,info,type,mold,dupl)

  call b%trim()
  call b%set_asb()
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return
contains
  subroutine inner_cp_fmt(a,b,info,type,mold,dupl)
    class(psb_c_base_sparse_mat), intent(in) :: a
    class(psb_c_base_sparse_mat), intent(inout), allocatable  :: b
    integer(psb_ipk_), intent(out)         :: info
    integer(psb_ipk_),optional, intent(in) :: dupl
    character(len=*), optional, intent(in) :: type
    class(psb_c_base_sparse_mat), intent(in), optional :: mold

    class(psb_c_base_sparse_mat), allocatable  :: altmp
    integer(psb_ipk_) :: err_act

    info = psb_success_
    call psb_erractionsave(err_act)

    if (present(mold)) then

      allocate(altmp, mold=mold,stat=info)

    else if (present(type)) then

      select case (psb_toupper(type))
      case ('CSR')
        allocate(psb_c_csr_sparse_mat :: altmp, stat=info)
      case ('COO')
        allocate(psb_c_coo_sparse_mat :: altmp, stat=info)
      case ('CSC')
        allocate(psb_c_csc_sparse_mat :: altmp, stat=info)
      case default
        info = psb_err_format_unknown_
        call psb_errpush(info,name,a_err=type)
        goto 9999
      end select
    else
      allocate(psb_c_csr_sparse_mat :: altmp, stat=info)
      !allocate(altmp, mold=psb_get_mat_default(a),stat=info)
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

    call altmp%cp_from_fmt(a, info)

    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err="mv_from")
      goto 9999
    end if

    call move_alloc(altmp,b)

    call psb_erractionrestore(err_act)
    return


9999 call psb_error_handler(err_act)

    return
  end subroutine inner_cp_fmt
end subroutine psb_c_cscnv

subroutine psb_c_cscnv_ip(a,info,type,mold,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cscnv_ip
  implicit none

  class(psb_cspmat_type), intent(inout)   :: a
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_),optional, intent(in)  :: dupl
  character(len=*), optional, intent(in)  :: type
  class(psb_c_base_sparse_mat), intent(in), optional :: mold

  class(psb_c_base_sparse_mat), allocatable  :: altmp
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

  if (allocated(a%a))   call inner_mv_fmt(a%a,info,type,mold,dupl)    
  if (allocated(a%ad))  call inner_mv_fmt(a%ad,info,type,mold,dupl)
  if (allocated(a%and)) call inner_mv_fmt(a%and,info,type,mold,dupl)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err="mv_from")
    goto 9999
  end if

  call a%trim()
  call a%set_asb()
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return
contains
  subroutine inner_mv_fmt(a,info,type,mold,dupl)
    class(psb_c_base_sparse_mat), intent(inout), allocatable  :: a
    integer(psb_ipk_), intent(out)                   :: info
    integer(psb_ipk_),optional, intent(in)           :: dupl
    character(len=*), optional, intent(in) :: type
    class(psb_c_base_sparse_mat), intent(in), optional :: mold
    class(psb_c_base_sparse_mat), allocatable  :: altmp
    integer(psb_ipk_) :: err_act

    info = psb_success_
    call psb_erractionsave(err_act)

    if (present(mold)) then

      allocate(altmp, mold=mold,stat=info)

    else if (present(type)) then

      select case (psb_toupper(type))
      case ('CSR')
        allocate(psb_c_csr_sparse_mat :: altmp, stat=info)
      case ('COO')
        allocate(psb_c_coo_sparse_mat :: altmp, stat=info)
      case ('CSC')
        allocate(psb_c_csc_sparse_mat :: altmp, stat=info)
      case default
        info = psb_err_format_unknown_
        call psb_errpush(info,name,a_err=type)
        goto 9999
      end select
    else
      allocate(psb_c_csr_sparse_mat :: altmp, stat=info)
      !allocate(altmp, mold=psb_get_mat_default(a),stat=info)
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

    call altmp%mv_from_fmt(a, info)

    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err="mv_from")
      goto 9999
    end if

    call move_alloc(altmp,a)

    call psb_erractionrestore(err_act)
    return


9999 call psb_error_handler(err_act)

    return
  end subroutine inner_mv_fmt

end subroutine psb_c_cscnv_ip


subroutine psb_c_cscnv_base(a,b,info,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cscnv_base
  implicit none
  class(psb_cspmat_type), intent(in)       :: a
  class(psb_c_base_sparse_mat), intent(out) :: b
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl


  type(psb_c_coo_sparse_mat)  :: altmp
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cscnv'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2

  info = psb_success_
  call psb_erractionsave(err_act)

  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    call a%a%cp_to_coo(altmp,info )
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call aa%cp_to_coo(altmp,info )    
  else      
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
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

end subroutine psb_c_cscnv_base


subroutine psb_c_clip_d(a,b,info)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_clip_d
  implicit none

  class(psb_cspmat_type), intent(in)    :: a
  class(psb_cspmat_type), intent(inout) :: b
  integer(psb_ipk_),intent(out)                  :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: acoo
  integer(psb_ipk_) :: i, j, nz

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)
  if (info == psb_success_) call a%cp_to(acoo)
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

end subroutine psb_c_clip_d



subroutine psb_c_clip_d_ip(a,info)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_clip_d_ip
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  integer(psb_ipk_),intent(out)                  :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat), allocatable  :: acoo
  integer(psb_ipk_) :: i, j, nz

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)
  if (info == psb_success_) call a%mv_to(acoo)
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

end subroutine psb_c_clip_d_ip


subroutine psb_c_mv_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_mv_from
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  logical :: do_split

  do_split = allocated(a%ad)
  
  call a%free()
  allocate(a%a,mold=b, stat=info)
  call a%a%mv_from_fmt(b,info)
  if (do_split) call a%split_nd(a%get_nrows(),a%get_ncols(),info) 
  call b%free()

  return
end subroutine psb_c_mv_from


subroutine psb_c_cp_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cp_from
  implicit none
  class(psb_cspmat_type), intent(out)      :: a
  class(psb_c_base_sparse_mat), intent(in) :: b
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='cp_from'
  logical, parameter :: debug=.false.
  logical :: do_split

  do_split = allocated(a%ad)
  

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
  if (do_split) call a%split_nd(a%get_nrows(),a%get_ncols(),info) 
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return
end subroutine psb_c_cp_from


subroutine psb_c_mv_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_mv_to
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2

  if (allocated(a%a)) then
    call b%mv_from_fmt(a%a,info)
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call b%mv_from_coo(aa,info)
  end if

  return
end subroutine psb_c_mv_to


subroutine psb_c_cp_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cp_to
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2

  if (allocated(a%a)) then
    call b%cp_from_fmt(a%a,info)
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call b%cp_from_coo(aa,info)
  end if
  return
end subroutine psb_c_cp_to

subroutine psb_c_mold(a,b)
  use psb_c_mat_mod, psb_protect_name => psb_c_mold
  class(psb_cspmat_type), intent(inout)     :: a
  class(psb_c_base_sparse_mat), allocatable, intent(out) :: b
  integer(psb_ipk_) :: info

  if (allocated(a%a)) then
    allocate(b,mold=a%a, stat=info)
  else if (allocated(a%ad)) then
    allocate(b,mold=a%ad, stat=info)
  end if

end subroutine psb_c_mold

subroutine psb_cspmat_type_move(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_cspmat_type_move
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_cspmat_type), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                   :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='move_alloc'
  logical, parameter :: debug=.false.

  info = psb_success_
  call b%free()
  call move_alloc(a%a,b%a)
  call move_alloc(a%ad,b%ad)
  call move_alloc(a%and,b%and)

  return
end subroutine psb_cspmat_type_move

subroutine psb_cspmat_clone(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_cspmat_clone
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_cspmat_type), intent(inout) :: b
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
  if (allocated(a%ad)) then
    call a%ad%clone(b%ad,info)
  end if
  if (allocated(a%and)) then
    call a%and%clone(b%and,info)
  end if
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_cspmat_clone


subroutine psb_c_transp_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_transp_1mat
  implicit none
  class(psb_cspmat_type), intent(inout) :: a

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2


  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    call a%a%transp()
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call aa%transp()
    call a%mv_from(aa)
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_transp_1mat

subroutine psb_c_transp_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_transp_2mat
  implicit none
  class(psb_cspmat_type), intent(in)  :: a
  class(psb_cspmat_type), intent(inout) :: b

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transp'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2


  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call b%free()
  if (allocated(a%a)) then
    allocate(b%a,mold=a%a,stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      goto 9999
    end if
    call a%a%transp(b%a)
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call aa%transp()
    call b%mv_from(aa)
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_transp_2mat

subroutine psb_c_transc_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_transc_1mat
  implicit none
  class(psb_cspmat_type), intent(inout) :: a

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2


  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    call a%a%transc()
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call aa%transc()
    call a%mv_from(aa)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_transc_1mat

subroutine psb_c_transc_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_transc_2mat
  implicit none
  class(psb_cspmat_type), intent(in)    :: a
  class(psb_cspmat_type), intent(inout) :: b

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='transc'
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat) :: aa, ac1,ac2

  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call b%free()
  if (allocated(a%a)) then
    allocate(b%a,mold=a%a,stat=info)
    if (info /= psb_success_) then
      info = psb_err_alloc_dealloc_
      goto 9999
    end if
    call a%a%transc(b%a)
  else if (allocated(a%ad)) then
    call a%ad%cp_to_coo(ac1,info)
    call a%and%cp_to_coo(ac2,info)
    call ac1%csmerge(ac2,aa,a%get_nrows(),a%get_ncols(),info)
    call aa%transc()
    call b%mv_from(aa)
  end if
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_transc_2mat

subroutine psb_c_asb(a,mold)
  use psb_c_mat_mod, psb_protect_name => psb_c_asb
  use psb_error_mod
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  class(psb_c_base_sparse_mat), optional, intent(in) :: mold
  class(psb_c_base_sparse_mat), allocatable :: tmp
  class(psb_c_base_sparse_mat), pointer :: mld
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='c_asb'

  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then 
    call a%a%asb()
    if (present(mold)) then
      if (.not.same_type_as(a%a,mold)) then
        allocate(tmp,mold=mold)
        call tmp%mv_from_fmt(a%a,info)
        call a%a%free()
        call move_alloc(tmp,a%a)
      end if
    else
      mld => psb_c_get_base_mat_default()
      if (.not.same_type_as(a%a,mld)) &
           & call a%cscnv(info)
    end if
    call a%split_nd(a%get_nrows(),a%get_ncols(),info) 
    
  else if (allocated(a%ad)) then 
    call a%ad%asb()
    call a%and%asb()
    if (present(mold)) then
      if (.not.same_type_as(a%ad,mold)) then
        allocate(tmp,mold=mold)
        call tmp%mv_from_fmt(a%ad,info)
        call a%ad%free()
        call move_alloc(tmp,a%ad)
        allocate(tmp,mold=mold)
        call tmp%mv_from_fmt(a%and,info)
        call a%and%free()
        call move_alloc(tmp,a%and)
      end if
    else
      mld => psb_c_get_base_mat_default()
      if (.not.same_type_as(a%a,mld)) &
           & call a%cscnv(info)
    end if
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_asb

subroutine psb_c_reinit(a,clear)
  use psb_c_mat_mod, psb_protect_name => psb_c_reinit
  use psb_error_mod
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  logical, intent(in), optional :: clear
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='reinit'

  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (allocated(a%a)) then
    call inner_reinit(a%a,name,info)
  else if (allocated(a%ad)) then
    call inner_reinit(a%ad,name,info)
    call inner_reinit(a%and,name,info)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif
  if (info /= 0) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return
contains
  subroutine inner_reinit(aa,name,info)
    class(psb_c_base_sparse_mat) :: aa
    character(len=*) :: name
    integer(psb_ipk_) :: info
    info = 0
    if (aa%has_update()) then
      call aa%reinit(clear)
    else
      info = psb_err_missing_override_method_
      call psb_errpush(info,name)
    endif
  end subroutine inner_reinit
end subroutine psb_c_reinit




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
!
!
! What do we do here??????
!

subroutine psb_c_csmm(alpha,a,x,beta,y,info,trans)
  use psb_error_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csmm
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
  complex(psb_spk_), intent(inout) :: y(:,:)
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

end subroutine psb_c_csmm


subroutine psb_c_csmv(alpha,a,x,beta,y,info,trans)
  use psb_error_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csmv
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)    :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout) :: y(:)
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

end subroutine psb_c_csmv

subroutine psb_c_csmv_vect(alpha,a,x,beta,y,info,trans)
  use psb_error_mod
  use psb_c_vect_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_csmv_vect
  implicit none
  class(psb_cspmat_type), intent(in)   :: a
  complex(psb_spk_), intent(in)        :: alpha, beta
  type(psb_c_vect_type), intent(inout) :: x
  type(psb_c_vect_type), intent(inout) :: y
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

end subroutine psb_c_csmv_vect



subroutine psb_c_cssm(alpha,a,x,beta,y,info,trans,scale,d)
  use psb_error_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cssm
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
  complex(psb_spk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  complex(psb_spk_), intent(in), optional :: d(:)
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

  if (allocated(a%a)) then
    call a%a%spsm(alpha,x,beta,y,info,trans,scale,d)
  else if (allocated(a%ad)) then 
    call a%ad%spsm(alpha,x,beta,y,info,trans,scale,d)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_cssm


subroutine psb_c_cssv(alpha,a,x,beta,y,info,trans,scale,d)
  use psb_error_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cssv
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)    :: alpha, beta, x(:)
  complex(psb_spk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  complex(psb_spk_), intent(in), optional :: d(:)
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

  if (allocated(a%a)) then
    call a%a%spsm(alpha,x,beta,y,info,trans,scale,d)
  else if (allocated(a%ad)) then 
    call a%ad%spsm(alpha,x,beta,y,info,trans,scale,d)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif
    
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_cssv


subroutine psb_c_cssv_vect(alpha,a,x,beta,y,info,trans,scale,d)
  use psb_error_mod
  use psb_c_vect_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cssv_vect
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), intent(in)         :: alpha, beta
  type(psb_c_vect_type), intent(inout)   :: x
  type(psb_c_vect_type), intent(inout)   :: y
  integer(psb_ipk_), intent(out)               :: info
  character, optional, intent(in)    :: trans, scale
  type(psb_c_vect_type), optional, intent(inout)   :: d
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
  if (allocated(a%a)) then
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
  else if (allocated(a%ad)) then 
    if (present(d)) then
      if (.not.allocated(d%v)) then
        info = psb_err_invalid_vect_state_
        call psb_errpush(info,name)
        goto 9999
      endif
      call a%ad%spsm(alpha,x%v,beta,y%v,info,trans,scale,d%v)
    else
      call a%ad%spsm(alpha,x%v,beta,y%v,info,trans,scale)
    end if
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_cssv_vect

function psb_c_maxval(a) result(res)
  use psb_c_mat_mod, psb_protect_name => psb_c_maxval
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  real(psb_spk_)         :: res


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

  if (allocated(a%a)) then
    res = a%a%maxval()
  else if (allocated(a%ad)) then 
    res = max(a%ad%maxval(),a%and%maxval())
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  return


9999 call psb_error_handler(err_act)

  return

end function psb_c_maxval

function psb_c_csnmi(a) result(res)
  use psb_c_mat_mod, psb_protect_name => psb_c_csnmi
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  real(psb_spk_)         :: res

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

  res = maxval(a%arwsum(info))

  return


9999 call psb_error_handler(err_act)

  return

end function psb_c_csnmi


function psb_c_csnm1(a) result(res)
  use psb_c_mat_mod, psb_protect_name => psb_c_csnm1
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  real(psb_spk_)         :: res

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

  res = maxval(a%aclsum(info))

  return

9999 call psb_error_handler(err_act)

  return

end function psb_c_csnm1


function psb_c_rowsum(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_c_rowsum
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), allocatable     :: d(:),d1(:)
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
  if (allocated(a%a)) then
    call a%a%rowsum(d)
  else if (allocated(a%ad)) then 
    call a%ad%rowsum(d)
    call a%and%rowsum(d1)
    d=d+d1
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_c_rowsum

function psb_c_arwsum(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_c_arwsum
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  real(psb_spk_), allocatable           :: d(:),d1(:)
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

  if (allocated(a%a)) then
    call a%a%arwsum(d)
  else if (allocated(a%ad)) then 
    call a%ad%arwsum(d)
    call a%and%arwsum(d1)
    d=d+d1
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_c_arwsum

function psb_c_colsum(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_c_colsum
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), allocatable         :: d(:), d1(:)
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

  if (allocated(a%a)) then
    call a%a%colsum(d)
  else if (allocated(a%ad)) then 
    call a%ad%colsum(d)
    call a%and%colsum(d1)
    d = [d,d1]
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
    
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_c_colsum

function psb_c_aclsum(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_c_aclsum
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  real(psb_spk_), allocatable           :: d(:),d1(:)
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

  if (allocated(a%a)) then
    call a%a%aclsum(d)
  else if (allocated(a%ad)) then 
    call a%ad%aclsum(d)
    call a%and%aclsum(d1)
    d = [d,d1]
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_c_aclsum


function psb_c_get_diag(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_c_get_diag
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_cspmat_type), intent(in) :: a
  complex(psb_spk_), allocatable         :: d(:)
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
  if (allocated(a%a)) then
    call a%a%get_diag(d,info)
  else if (allocated(a%ad)) then 
    call a%ad%get_diag(d,info)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_c_get_diag


subroutine psb_c_scal(d,a,info,side)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_scal
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)              :: d(:)
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

  if (allocated(a%a)) then
    call a%a%scal(d,info,side=side)
  else if (allocated(a%ad)) then
    call a%ad%scal(d,info,side=side)
    !
    ! FIXME
    !
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_scal


subroutine psb_c_scals(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_scals
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)              :: d
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

  if (allocated(a%a)) then
    call a%a%scal(d,info)
  else if (allocated(a%ad)) then 
    call a%ad%scal(d,info)
    call a%and%scal(d,info)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_scals

subroutine psb_c_scalplusidentity(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_scalplusidentity
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)              :: d
  integer(psb_ipk_), intent(out)                    :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='scalplusidentity'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    call a%a%scalpid(d,info)
  else if (allocated(a%ad)) then 
    call a%ad%scalpid(d,info)
    call a%and%scal(d,info)    
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_scalplusidentity

subroutine psb_c_spaxpby(alpha,a,beta,b,info)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_spaxpby
  implicit none
  complex(psb_spk_), intent(in)             :: alpha
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)             :: beta
  class(psb_cspmat_type), intent(inout) :: b
  integer(psb_ipk_), intent(out)          :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='spaxby'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    call a%a%spaxpby(alpha,beta,b%a,info)
  else if (allocated(a%ad)) then 
    call a%ad%spaxpby(alpha,beta,b%a,info)
    call a%and%spaxpby(alpha,cone,b%a,info)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif
  
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_spaxpby


function psb_c_cmpval(a,val,tol,info) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cmpval
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)             :: val
  real(psb_spk_), intent(in)            :: tol
  logical                                 :: res
  integer(psb_ipk_), intent(out)          :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cmpval'
  logical, parameter :: debug=.false.

  res = .false.
  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    res = a%a%spcmp(val,tol,info)
  else if (allocated(a%ad)) then 
    res = a%ad%spcmp(val,tol,info) .and. a%and%spcmp(val,tol,info)
1  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif
    
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_c_cmpval

function psb_c_cmpmat(a,b,tol,info) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cmpmat
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_cspmat_type), intent(inout) :: b
  real(psb_spk_), intent(in)            :: tol
  logical                                 :: res
  integer(psb_ipk_), intent(out)          :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cmpmat'
  logical, parameter :: debug=.false.

  res = .false.
  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (allocated(a%a)) then
    res = a%a%spcmp(b%a,tol,info)
  else if (allocated(a%ad)) then 
    res = a%ad%spcmp(b%ad,tol,info) .and. a%and%spcmp(b%and,tol,info)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
  endif
    
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_c_cmpmat

subroutine psb_c_mv_from_lb(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_mv_from_lb
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  class(psb_lc_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  info = psb_success_
  call a%free()
  if (.not.allocated(a%a)) allocate(psb_c_csr_sparse_mat :: a%a, stat=info)
  if (info == psb_success_) call a%a%mv_from_lfmt(b,info)

end subroutine psb_c_mv_from_lb


subroutine psb_c_cp_from_lb(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cp_from_lb
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  class(psb_lc_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  info = psb_success_
  call a%free()
  if (.not.allocated(a%a)) allocate(psb_c_csr_sparse_mat :: a%a, stat=info)
  if (info == psb_success_) call a%a%cp_from_lfmt(b,info)

end subroutine psb_c_cp_from_lb

subroutine psb_c_mv_to_lb(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_mv_to_lb
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  class(psb_lc_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  type(psb_c_coo_sparse_mat) :: acoo

  if (.not.allocated(a%a)) then
    if (allocated(a%ad)) then 
      call a%merge_nd(a%get_nrows(),a%get_ncols(),info,acoo=acoo)
      call acoo%mv_to_lfmt(b,info)
    else
      call b%free()
    end if
  else
    call a%a%mv_to_lfmt(b,info)
  end if
  call a%free()

end subroutine psb_c_mv_to_lb

subroutine psb_c_cp_to_lb(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cp_to_lb
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_lc_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info
  type(psb_c_coo_sparse_mat) :: acoo
  
  if (.not.allocated(a%a)) then
    if (allocated(a%ad)) then 
      call a%merge_nd(a%get_nrows(),a%get_ncols(),info,acoo=acoo)
      call acoo%mv_to_lfmt(b,info)
    else
      call b%free()
    end if
  else
    call a%a%cp_to_lfmt(b,info)
  end if

end subroutine psb_c_cp_to_lb

subroutine psb_c_mv_from_l(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_mv_from_l
  implicit none
  class(psb_cspmat_type), intent(inout) :: a
  class(psb_lcspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info

  info = psb_success_
  if (allocated(b%a)) then
    if (.not.allocated(a%a)) allocate(psb_c_csr_sparse_mat :: a%a, stat=info)
    call a%a%mv_from_lfmt(b%a,info)
  else
    call a%free()
  end if
  call b%free()

end subroutine psb_c_mv_from_l


subroutine psb_c_cp_from_l(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cp_from_l
  implicit none

  class(psb_cspmat_type), intent(out) :: a
  class(psb_lcspmat_type), intent(in) :: b
  integer(psb_ipk_) :: info

  info = psb_success_
  if (allocated(b%a)) then
    if (.not.allocated(a%a)) allocate(psb_c_csr_sparse_mat :: a%a, stat=info)
    call a%a%cp_from_lfmt(b%a,info)
  else
    call a%free()
  end if
end subroutine psb_c_cp_from_l

subroutine psb_c_mv_to_l(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_mv_to_l
  implicit none

  class(psb_cspmat_type), intent(inout) :: a
  class(psb_lcspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (allocated(a%a)) then
    if (.not.allocated(b%a)) allocate(psb_lc_csr_sparse_mat :: b%a, stat=info)
    call a%a%mv_to_lfmt(b%a,info)
  else
    call b%free()
  end if
  call a%free()

end subroutine psb_c_mv_to_l

subroutine psb_c_cp_to_l(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_c_cp_to_l
  implicit none

  class(psb_cspmat_type), intent(in) :: a
  class(psb_lcspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (allocated(a%a)) then
    if (.not.allocated(b%a)) allocate(psb_lc_csr_sparse_mat :: b%a, stat=info)
    call a%a%cp_to_lfmt(b%a,info)
  else
    call b%free()
  end if

end subroutine psb_c_cp_to_l


!
!
! lc versions
!


subroutine  psb_lc_set_lnrows(m,a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_lnrows
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_lnrows

#if defined(IPK4) && defined(LPK8)
subroutine  psb_lc_set_inrows(m,a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_inrows
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_inrows
#endif

subroutine  psb_lc_set_lncols(n,a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_lncols
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_lncols

#if defined(IPK4) && defined(LPK8)
subroutine  psb_lc_set_incols(n,a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_incols
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_incols
#endif

!
!  Valid values for DUPL:
!  psb_dupl_ovwrt_
!  psb_dupl_add_
!  psb_dupl_err_
!

subroutine  psb_lc_set_dupl(n,a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_dupl
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_dupl


!
! Set the STATE of the internal matrix object
!

subroutine  psb_lc_set_null(a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_null
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_null


subroutine  psb_lc_set_bld(a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_bld
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_bld


subroutine  psb_lc_set_upd(a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_upd
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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


end subroutine psb_lc_set_upd


subroutine  psb_lc_set_asb(a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_asb
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_asb


subroutine psb_lc_set_sorted(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_sorted
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_sorted


subroutine psb_lc_set_triangle(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_triangle
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_triangle

subroutine psb_lc_set_symmetric(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_symmetric
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

  call a%a%set_symmetric(val)

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_lc_set_symmetric

subroutine psb_lc_set_unit(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_unit
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_unit


subroutine psb_lc_set_lower(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_lower
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_lower


subroutine psb_lc_set_upper(a,val)
  use psb_c_mat_mod, psb_protect_name => psb_lc_set_upper
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_set_upper



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


subroutine psb_lc_sparse_print(iout,a,iv,head,ivr,ivc)
  use psb_c_mat_mod, psb_protect_name => psb_lc_sparse_print
  use psb_error_mod
  implicit none

  integer(psb_ipk_), intent(in)               :: iout
  class(psb_lcspmat_type), intent(in) :: a
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

end subroutine psb_lc_sparse_print


subroutine psb_lc_n_sparse_print(fname,a,iv,head,ivr,ivc)
  use psb_c_mat_mod, psb_protect_name => psb_lc_n_sparse_print
  use psb_error_mod
  implicit none

  character(len=*), intent(in)  :: fname
  class(psb_lcspmat_type), intent(in) :: a
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

end subroutine psb_lc_n_sparse_print


subroutine psb_lc_get_neigh(a,idx,neigh,n,info,lev)
  use psb_c_mat_mod, psb_protect_name => psb_lc_get_neigh
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
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

end subroutine psb_lc_get_neigh



subroutine psb_lc_csall(nr,nc,a,info,nz,type,mold)
  use psb_c_mat_mod, psb_protect_name => psb_lc_csall
  use psb_c_base_mat_mod
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  integer(psb_lpk_), intent(in)             :: nr,nc
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_lpk_), intent(in), optional   :: nz
  character(len=*), intent(in), optional    :: type
  class(psb_lc_base_sparse_mat), optional, intent(in) :: mold

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csall'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)

  call a%free()

  info = psb_success_
  if (present(mold)) then
    allocate(a%a, stat=info, mold=mold)
  else if (present(type)) then
    select case (type)
    case('CSR')
      allocate(psb_lc_csr_sparse_mat :: a%a, stat=info)
    case('COO')
      allocate(psb_lc_coo_sparse_mat :: a%a, stat=info)
    case('CSC')
      allocate(psb_lc_csc_sparse_mat :: a%a, stat=info)
    case default
      allocate(psb_lc_coo_sparse_mat :: a%a, stat=info)
    end select
  else
    allocate(psb_lc_coo_sparse_mat :: a%a, stat=info)
  end if
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

end subroutine psb_lc_csall


subroutine  psb_lc_reallocate_nz(nz,a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_reallocate_nz
  use psb_error_mod
  implicit none
  integer(psb_lpk_), intent(in) :: nz
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_reallocate_nz


subroutine  psb_lc_free(a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_free
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a

  if (allocated(a%a)) then
    call a%a%free()
    deallocate(a%a)
  endif

end subroutine psb_lc_free


subroutine  psb_lc_trim(a)
  use psb_c_mat_mod, psb_protect_name => psb_lc_trim
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_trim



subroutine psb_lc_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
  use psb_c_mat_mod, psb_protect_name => psb_lc_csput_a
  use psb_c_base_mat_mod
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)      :: val(:)
  integer(psb_lpk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info

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


  call a%a%csput(nz,ia,ja,val,imin,imax,jmin,jmax,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_lc_csput_a

subroutine psb_lc_csput_v(nz,ia,ja,val,a,imin,imax,jmin,jmax,info)
  use psb_c_mat_mod, psb_protect_name => psb_lc_csput_v
  use psb_c_base_mat_mod
  use psb_c_vect_mod, only : psb_c_vect_type
  use psb_l_vect_mod, only : psb_l_vect_type
  use psb_error_mod
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  type(psb_c_vect_type), intent(inout)  :: val
  type(psb_l_vect_type), intent(inout)  :: ia, ja
  integer(psb_lpk_), intent(in)             :: nz, imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info

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
    call a%a%csput(nz,ia%v,ja%v,val%v,imin,imax,jmin,jmax,info)
  else
    info = psb_err_invalid_mat_state_
  endif

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_lc_csput_v


subroutine psb_lc_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_csgetptn
  implicit none

  class(psb_lcspmat_type), intent(in) :: a
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

end subroutine psb_lc_csgetptn


subroutine psb_lc_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_csgetrow
  implicit none

  class(psb_lcspmat_type), intent(in) :: a
  integer(psb_lpk_), intent(in)                  :: imin,imax
  integer(psb_lpk_), intent(out)                 :: nz
  integer(psb_lpk_), allocatable, intent(inout)  :: ia(:), ja(:)
  complex(psb_spk_), allocatable,  intent(inout)    :: val(:)
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

end subroutine psb_lc_csgetrow




subroutine psb_lc_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_csgetblk
  implicit none

  class(psb_lcspmat_type), intent(in)    :: a
  class(psb_lcspmat_type), intent(inout) :: b
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
  type(psb_lc_coo_sparse_mat), allocatable  :: acoo


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

end subroutine psb_lc_csgetblk


subroutine psb_lc_tril(a,l,info,diag,imin,imax,&
     & jmin,jmax,rscale,cscale,u)
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_tril
  implicit none
  class(psb_lcspmat_type), intent(in)      :: a
  class(psb_lcspmat_type), intent(inout)   :: l
  integer(psb_ipk_),intent(out)           :: info
  integer(psb_lpk_), intent(in), optional :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional           :: rscale,cscale
  class(psb_lcspmat_type), optional, intent(inout)   :: u

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='tril'
  logical, parameter :: debug=.false.
  type(psb_lc_coo_sparse_mat), allocatable  :: lcoo, ucoo

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


end subroutine psb_lc_tril

subroutine psb_lc_triu(a,u,info,diag,imin,imax,&
     & jmin,jmax,rscale,cscale,l)
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_triu
  implicit none
  class(psb_lcspmat_type), intent(in)      :: a
  class(psb_lcspmat_type), intent(inout)   :: u
  integer(psb_ipk_),intent(out)           :: info
  integer(psb_lpk_), intent(in), optional :: diag,imin,imax,jmin,jmax
  logical, intent(in), optional           :: rscale,cscale
  class(psb_lcspmat_type), optional, intent(inout)   :: l

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='triu'
  logical, parameter :: debug=.false.
  type(psb_lc_coo_sparse_mat), allocatable  :: lcoo, ucoo


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


end subroutine psb_lc_triu


subroutine psb_lc_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_csclip
  implicit none

  class(psb_lcspmat_type), intent(in) :: a
  class(psb_lcspmat_type), intent(inout) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_lpk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_lc_coo_sparse_mat), allocatable  :: acoo

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

end subroutine psb_lc_csclip

subroutine psb_lc_csclip_ip(a,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_csclip_ip
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_lpk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='csclip'
  logical, parameter :: debug=.false.
  type(psb_lc_coo_sparse_mat), allocatable  :: acoo

  info = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_null()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  allocate(acoo,stat=info)
  if (info == psb_success_) then
    call a%a%csclip(acoo,info,&
         & imin,imax,jmin,jmax,rscale,cscale)
  else
    info = psb_err_alloc_dealloc_
  end if
  if (info == psb_success_) call a%free()
  if (info == psb_success_) call move_alloc(acoo,a%a)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_lc_csclip_ip

subroutine psb_lc_b_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_b_csclip
  implicit none

  class(psb_lcspmat_type), intent(in) :: a
  type(psb_lc_coo_sparse_mat), intent(out) :: b
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

end subroutine psb_lc_b_csclip




subroutine psb_lc_cscnv(a,b,info,type,mold,upd,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cscnv
  implicit none
  class(psb_lcspmat_type), intent(in)      :: a
  class(psb_lcspmat_type), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl, upd
  character(len=*), optional, intent(in) :: type
  class(psb_lc_base_sparse_mat), intent(in), optional :: mold


  class(psb_lc_base_sparse_mat), allocatable  :: altmp
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
      allocate(psb_lc_csr_sparse_mat :: altmp, stat=info)
    case ('COO')
      allocate(psb_lc_coo_sparse_mat :: altmp, stat=info)
    case ('CSC')
      allocate(psb_lc_csc_sparse_mat :: altmp, stat=info)
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

end subroutine psb_lc_cscnv



subroutine psb_lc_cscnv_ip(a,info,type,mold,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cscnv_ip
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl
  character(len=*), optional, intent(in) :: type
  class(psb_lc_base_sparse_mat), intent(in), optional :: mold


  class(psb_lc_base_sparse_mat), allocatable  :: altmp
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
      allocate(psb_lc_csr_sparse_mat :: altmp, stat=info)
    case ('COO')
      allocate(psb_lc_coo_sparse_mat :: altmp, stat=info)
    case ('CSC')
      allocate(psb_lc_csc_sparse_mat :: altmp, stat=info)
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

end subroutine psb_lc_cscnv_ip



subroutine psb_lc_cscnv_base(a,b,info,dupl)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cscnv_base
  implicit none
  class(psb_lcspmat_type), intent(in)       :: a
  class(psb_lc_base_sparse_mat), intent(out) :: b
  integer(psb_ipk_), intent(out)                   :: info
  integer(psb_ipk_),optional, intent(in)           :: dupl


  type(psb_lc_coo_sparse_mat)  :: altmp
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

end subroutine psb_lc_cscnv_base



subroutine psb_lc_clip_d(a,b,info)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_clip_d
  implicit none

  class(psb_lcspmat_type), intent(in)    :: a
  class(psb_lcspmat_type), intent(inout) :: b
  integer(psb_ipk_),intent(out)                  :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_lc_coo_sparse_mat), allocatable  :: acoo
  integer(psb_lpk_) :: i, j, nz

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

end subroutine psb_lc_clip_d



subroutine psb_lc_clip_d_ip(a,info)
  ! Output is always in  COO format
  use psb_error_mod
  use psb_const_mod
  use psb_c_base_mat_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_clip_d_ip
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
  integer(psb_ipk_),intent(out)                  :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='clip_diag'
  logical, parameter :: debug=.false.
  type(psb_lc_coo_sparse_mat), allocatable  :: acoo
  integer(psb_lpk_) :: i, j, nz

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

end subroutine psb_lc_clip_d_ip


subroutine psb_lc_mv_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_mv_from
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_lc_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  call a%free()
  allocate(a%a,mold=b, stat=info)
  call a%a%mv_from_fmt(b,info)
  call b%free()

  return
end subroutine psb_lc_mv_from


subroutine psb_lc_cp_from(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cp_from
  implicit none
  class(psb_lcspmat_type), intent(out)      :: a
  class(psb_lc_base_sparse_mat), intent(in) :: b
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
end subroutine psb_lc_cp_from


subroutine psb_lc_mv_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_mv_to
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_lc_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  call b%mv_from_fmt(a%a,info)

  return
end subroutine psb_lc_mv_to


subroutine psb_lc_cp_to(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cp_to
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  class(psb_lc_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  call b%cp_from_fmt(a%a,info)

  return
end subroutine psb_lc_cp_to

subroutine psb_lc_mold(a,b)
  use psb_c_mat_mod, psb_protect_name => psb_lc_mold
  class(psb_lcspmat_type), intent(inout)     :: a
  class(psb_lc_base_sparse_mat), allocatable, intent(out) :: b
  integer(psb_ipk_) :: info

  allocate(b,mold=a%a, stat=info)

end subroutine psb_lc_mold

subroutine psb_lcspmat_type_move(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lcspmat_type_move
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_lcspmat_type), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                   :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='move_alloc'
  logical, parameter :: debug=.false.

  info = psb_success_
  call b%free()
  call move_alloc(a%a,b%a)

  return
end subroutine psb_lcspmat_type_move


subroutine psb_lcspmat_clone(a,b,info)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lcspmat_clone
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_lcspmat_type), intent(inout) :: b
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

end subroutine psb_lcspmat_clone


subroutine psb_lc_transp_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_transp_1mat
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a

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

end subroutine psb_lc_transp_1mat



subroutine psb_lc_transp_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_transp_2mat
  implicit none
  class(psb_lcspmat_type), intent(in)  :: a
  class(psb_lcspmat_type), intent(inout) :: b

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

end subroutine psb_lc_transp_2mat


subroutine psb_lc_transc_1mat(a)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_transc_1mat
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a

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

end subroutine psb_lc_transc_1mat



subroutine psb_lc_transc_2mat(a,b)
  use psb_error_mod
  use psb_string_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_transc_2mat
  implicit none
  class(psb_lcspmat_type), intent(in)    :: a
  class(psb_lcspmat_type), intent(inout) :: b

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

end subroutine psb_lc_transc_2mat


subroutine psb_lc_asb(a,mold)
  use psb_c_mat_mod, psb_protect_name => psb_lc_asb
  use psb_error_mod
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_lc_base_sparse_mat), optional, intent(in) :: mold
  class(psb_lc_base_sparse_mat), allocatable :: tmp
  class(psb_lc_base_sparse_mat), pointer :: mld
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='lc_asb'

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
    mld => psb_lc_get_base_mat_default()
    if (.not.same_type_as(a%a,mld)) &
         & call a%cscnv(info)
  end if


  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)

  return

end subroutine psb_lc_asb

subroutine psb_lc_reinit(a,clear)
  use psb_c_mat_mod, psb_protect_name => psb_lc_reinit
  use psb_error_mod
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
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

end subroutine psb_lc_reinit




function psb_lc_get_diag(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_lc_get_diag
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  complex(psb_spk_), allocatable         :: d(:)
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

end function psb_lc_get_diag


subroutine psb_lc_scal(d,a,info,side)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_scal
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)              :: d(:)
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

end subroutine psb_lc_scal


subroutine psb_lc_scals(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_scals
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)              :: d
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

end subroutine psb_lc_scals

subroutine psb_lc_scalplusidentity(d,a,info)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_scalplusidentity
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)              :: d
  integer(psb_ipk_), intent(out)                    :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='scalplusidentity'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%scalpid(d,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_lc_scalplusidentity

subroutine psb_lc_spaxpby(alpha,a,beta,b,info)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_spaxpby
  implicit none
  complex(psb_spk_), intent(in)             :: alpha
  class(psb_lcspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)             :: beta
  class(psb_lcspmat_type), intent(inout) :: b
  integer(psb_ipk_), intent(out)          :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='spaxby'
  logical, parameter :: debug=.false.

  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  call a%a%spaxpby(alpha,beta,b%a,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_lc_spaxpby

function psb_lc_cmpval(a,val,tol,info) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cmpval
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  complex(psb_spk_), intent(in)             :: val
  real(psb_spk_), intent(in)            :: tol
  logical                                 :: res
  integer(psb_ipk_), intent(out)          :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cmpval'
  logical, parameter :: debug=.false.

  res = .false.
  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  res = a%a%spcmp(val,tol,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_lc_cmpval

function psb_lc_cmpmat(a,b,tol,info) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cmpmat
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_lcspmat_type), intent(inout) :: b
  real(psb_spk_), intent(in)            :: tol
  logical                                 :: res
  integer(psb_ipk_), intent(out)          :: info

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='cmpmat'
  logical, parameter :: debug=.false.

  res = .false.
  info = psb_success_
  call psb_erractionsave(err_act)
  if (.not.allocated(a%a)) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  res = a%a%spcmp(b%a,tol,info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end function psb_lc_cmpmat

function psb_lc_maxval(a) result(res)
  use psb_c_mat_mod, psb_protect_name => psb_lc_maxval
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  real(psb_spk_)         :: res


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

end function psb_lc_maxval

function psb_lc_csnmi(a) result(res)
  use psb_c_mat_mod, psb_protect_name => psb_lc_csnmi
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  real(psb_spk_)         :: res

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

end function psb_lc_csnmi

function psb_lc_csnm1(a) result(res)
  use psb_c_mat_mod, psb_protect_name => psb_lc_csnm1
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  real(psb_spk_)         :: res

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

end function psb_lc_csnm1


function psb_lc_rowsum(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_lc_rowsum
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  complex(psb_spk_), allocatable     :: d(:)
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

end function psb_lc_rowsum

function psb_lc_arwsum(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_lc_arwsum
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  real(psb_spk_), allocatable           :: d(:)
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

end function psb_lc_arwsum

function psb_lc_colsum(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_lc_colsum
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  complex(psb_spk_), allocatable         :: d(:)
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

end function psb_lc_colsum

function psb_lc_aclsum(a,info) result(d)
  use psb_c_mat_mod, psb_protect_name => psb_lc_aclsum
  use psb_error_mod
  use psb_const_mod
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  real(psb_spk_), allocatable           :: d(:)
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

end function psb_lc_aclsum

subroutine psb_lc_mv_from_ib(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_mv_from_ib
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  info = psb_success_
  if (.not.allocated(a%a)) allocate(psb_lc_csr_sparse_mat :: a%a, stat=info)
  if (info == psb_success_) call a%a%mv_from_ifmt(b,info)

end subroutine psb_lc_mv_from_ib

subroutine psb_lc_cp_from_ib(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cp_from_ib
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  info = psb_success_
  if (.not.allocated(a%a)) allocate(psb_lc_csr_sparse_mat :: a%a, stat=info)
  if (info == psb_success_) call a%a%cp_from_ifmt(b,info)

end subroutine psb_lc_cp_from_ib

subroutine psb_lc_mv_to_ib(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_mv_to_ib
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (.not.allocated(a%a)) then
    call b%free()
  else
    call a%a%mv_to_ifmt(b,info)
    call a%free()
  end if

end subroutine psb_lc_mv_to_ib

subroutine psb_lc_cp_to_ib(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cp_to_ib
  implicit none
  class(psb_lcspmat_type), intent(in) :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (.not.allocated(a%a)) then
    call b%free()
  else
    call a%a%cp_to_ifmt(b,info)
  end if

end subroutine psb_lc_cp_to_ib

subroutine psb_lc_mv_from_i(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_mv_from_i
  implicit none
  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_cspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (allocated(b%a)) then
    if (.not.allocated(a%a)) allocate(psb_lc_csr_sparse_mat :: a%a, stat=info)
    call a%a%mv_from_ifmt(b%a,info)
  else
    call a%free()
  end if
  call b%free()

end subroutine psb_lc_mv_from_i


subroutine psb_lc_cp_from_i(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cp_from_i
  implicit none

  class(psb_lcspmat_type), intent(out) :: a
  class(psb_cspmat_type), intent(in) :: b
  integer(psb_ipk_) :: info

  if (allocated(b%a)) then
    if (.not.allocated(a%a)) allocate(psb_lc_csr_sparse_mat :: a%a, stat=info)
    call a%a%cp_from_ifmt(b%a,info)
  else
    call a%free()
  end if
end subroutine psb_lc_cp_from_i

subroutine psb_lc_mv_to_i(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_mv_to_i
  implicit none

  class(psb_lcspmat_type), intent(inout) :: a
  class(psb_cspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (allocated(a%a)) then
    if (.not.allocated(b%a)) allocate(psb_c_csr_sparse_mat :: b%a, stat=info)
    call a%a%mv_to_ifmt(b%a,info)
  else
    call b%free()
  end if
  call a%free()

end subroutine psb_lc_mv_to_i

subroutine psb_lc_cp_to_i(a,b)
  use psb_error_mod
  use psb_const_mod
  use psb_c_mat_mod, psb_protect_name => psb_lc_cp_to_i
  implicit none

  class(psb_lcspmat_type), intent(in) :: a
  class(psb_cspmat_type), intent(inout) :: b
  integer(psb_ipk_) :: info

  if (allocated(a%a)) then
    if (.not.allocated(b%a)) allocate(psb_c_csr_sparse_mat :: b%a, stat=info)
    call a%a%cp_to_ifmt(b%a,info)
  else
    call b%free()
  end if

end subroutine psb_lc_cp_to_i
