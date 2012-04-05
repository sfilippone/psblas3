! == ==================================
!
!
!
! Data management
!
!
!
!
!
! == ==================================

subroutine psb_z_base_cp_to_coo(a,b,info)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_cp_to_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  class(psb_z_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_cp_to_coo

subroutine psb_z_base_cp_from_coo(a,b,info)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_cp_from_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  class(psb_z_coo_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_cp_from_coo


subroutine psb_z_base_cp_to_fmt(a,b,info)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_cp_to_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  class(psb_z_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='to_fmt'
  logical, parameter :: debug=.false.
  type(psb_z_coo_sparse_mat)  :: tmp
  
  !
  ! Default implementation
  ! 
  info = psb_success_
  call psb_erractionsave(err_act)

  call a%cp_to_coo(tmp,info)
  if (info == psb_success_) call b%mv_from_coo(tmp,info)

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if
    
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
  return

end subroutine psb_z_base_cp_to_fmt

subroutine psb_z_base_cp_from_fmt(a,b,info)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_cp_from_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  class(psb_z_base_sparse_mat), intent(in) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='from_fmt'
  logical, parameter :: debug=.false.
  type(psb_z_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  ! 
  info  = psb_success_
  call psb_erractionsave(err_act)

  call b%cp_to_coo(tmp,info)
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='to/from coo')
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_z_base_cp_from_fmt


subroutine psb_z_base_mv_to_coo(a,b,info)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_mv_to_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  class(psb_z_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='to_coo'
  logical, parameter :: debug=.false.


  info  = psb_success_
  call psb_erractionsave(err_act)
  
  call a%cp_to_coo(b,info) 

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='to coo')
    goto 9999
  end if

  call a%free()
  
  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_z_base_mv_to_coo

subroutine psb_z_base_mv_from_coo(a,b,info)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_mv_from_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  class(psb_z_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='from_coo'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)
  
  call a%cp_from_coo(b,info) 

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='from coo')
    goto 9999
  end if

  call b%free()
  
  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_z_base_mv_from_coo


subroutine psb_z_base_mv_to_fmt(a,b,info)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_mv_to_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  class(psb_z_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='to_fmt'
  logical, parameter :: debug=.false.
  type(psb_z_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  ! 
  info = psb_success_
  call a%mv_to_coo(tmp,info)
  if (info == psb_success_) call b%mv_from_coo(tmp,info)

  return

end subroutine psb_z_base_mv_to_fmt

subroutine psb_z_base_mv_from_fmt(a,b,info)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_mv_from_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  class(psb_z_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='from_fmt'
  logical, parameter :: debug=.false.
  type(psb_z_coo_sparse_mat)  :: tmp

  !
  ! Default implementation
  ! 
  info = psb_success_
  call b%mv_to_coo(tmp,info)
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  return

end subroutine psb_z_base_mv_from_fmt

subroutine psb_z_base_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_error_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_csput
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  complex(psb_dpk_), intent(in)      :: val(:)
  integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: gtl(:)

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csput'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_csput

subroutine psb_z_base_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_csgetrow
  implicit none

  class(psb_z_base_sparse_mat), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  complex(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_csgetrow



!
! Here we have the base implementation of getblk and clip:
! this is just based on the getrow.
! If performance is critical it can be overridden.
!
subroutine psb_z_base_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_csgetblk
  implicit none

  class(psb_z_base_sparse_mat), intent(in) :: a
  class(psb_z_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale
  integer(psb_ipk_) :: err_act, nzin, nzout
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csget'
  logical :: append_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(append)) then 
    append_ = append
  else
    append_ = .false.
  endif
  if (append_) then 
    nzin = a%get_nzeros()
  else
    nzin = 0
  endif

  call a%csget(imin,imax,nzout,b%ia,b%ja,b%val,info,&
       & jmin=jmin, jmax=jmax, iren=iren, append=append_, &
       & nzin=nzin, rscale=rscale, cscale=cscale)

  if (info /= psb_success_) goto 9999

  call b%set_nzeros(nzin+nzout)
  call b%fix(info)
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_z_base_csgetblk


subroutine psb_z_base_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_csclip
  implicit none

  class(psb_z_base_sparse_mat), intent(in) :: a
  class(psb_z_coo_sparse_mat), intent(out) :: b
  integer(psb_ipk_),intent(out)                  :: info
  integer(psb_ipk_), intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  integer(psb_ipk_) :: err_act, nzin, nzout, imin_, imax_, jmin_, jmax_, mb,nb
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csget'
  logical :: rscale_, cscale_
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  nzin = 0
  if (present(imin)) then 
    imin_ = imin
  else
    imin_ = 1
  end if
  if (present(imax)) then 
    imax_ = imax
  else
    imax_ = a%get_nrows()
  end if
  if (present(jmin)) then 
    jmin_ = jmin
  else
    jmin_ = 1
  end if
  if (present(jmax)) then 
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  end if
  if (present(rscale)) then 
    rscale_ = rscale
  else
    rscale_ = .true.
  end if
  if (present(cscale)) then 
    cscale_ = cscale
  else
    cscale_ = .true.
  end if

  if (rscale_) then 
    mb = imax_ - imin_ +1
  else 
    mb = a%get_nrows() ! Should this be imax_ ?? 
  endif
  if (cscale_) then 
    nb = jmax_ - jmin_ +1
  else 
    nb = a%get_ncols()  ! Should this be jmax_ ?? 
  endif
  call b%allocate(mb,nb)
  call a%csget(imin_,imax_,nzout,b%ia,b%ja,b%val,info,&
       & jmin=jmin_, jmax=jmax_, append=.false., &
       & nzin=nzin, rscale=rscale_, cscale=cscale_)
  if (info /= psb_success_) goto 9999

  call b%set_nzeros(nzin+nzout)
  call b%fix(info)

  if (info /= psb_success_) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_z_base_csclip

subroutine psb_z_base_mold(a,b,info) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_mold
  use psb_error_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in)  :: a
  class(psb_z_base_sparse_mat), intent(out), allocatable  :: b
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_mold

subroutine psb_z_base_transp_2mat(a,b)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_transp_2mat
  use psb_error_mod
  implicit none 

  class(psb_z_base_sparse_mat), intent(in) :: a
  class(psb_base_sparse_mat), intent(out)    :: b

  type(psb_z_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=*), parameter :: name='z_base_transp'

  call psb_erractionsave(err_act)

  info = psb_success_
  select type(b)
  class is (psb_z_base_sparse_mat)
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call tmp%transp()
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  class default
    info = psb_err_invalid_dynamic_type_
  end select
  if (info /= psb_success_) then 
    ierr(1)=ione;
    call psb_errpush(info,name,a_err=b%get_fmt(),i_err=ierr)
    goto 9999
  end if
  call psb_erractionrestore(err_act) 

  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_z_base_transp_2mat

subroutine psb_z_base_transc_2mat(a,b)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_transc_2mat
  implicit none 

  class(psb_z_base_sparse_mat), intent(in) :: a
  class(psb_base_sparse_mat), intent(out)    :: b

  type(psb_z_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=*), parameter :: name='z_base_transc'

  call psb_erractionsave(err_act)

  info = psb_success_
  select type(b)
  class is (psb_z_base_sparse_mat)
    call a%cp_to_coo(tmp,info)
    if (info == psb_success_) call tmp%transc()
    if (info == psb_success_) call b%mv_from_coo(tmp,info)
  class default
    info = psb_err_invalid_dynamic_type_
  end select
  if (info /= psb_success_) then 
    ierr(1) = ione;
    call psb_errpush(info,name,a_err=b%get_fmt(),i_err=ierr)
    goto 9999
  end if
  call psb_erractionrestore(err_act) 

  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return
end subroutine psb_z_base_transc_2mat

subroutine psb_z_base_transp_1mat(a)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_transp_1mat
  use psb_error_mod
  implicit none 

  class(psb_z_base_sparse_mat), intent(inout) :: a

  type(psb_z_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=*), parameter :: name='z_base_transp'

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%mv_to_coo(tmp,info)
  if (info == psb_success_) call tmp%transp()
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then 
    info = psb_err_missing_override_method_ 
    call psb_errpush(info,name,a_err=a%get_fmt())
    goto 9999
  end if
  call psb_erractionrestore(err_act) 

  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_z_base_transp_1mat

subroutine psb_z_base_transc_1mat(a)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_transc_1mat
  implicit none 

  class(psb_z_base_sparse_mat), intent(inout) :: a

  type(psb_z_coo_sparse_mat) :: tmp
  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=*), parameter :: name='z_base_transc'

  call psb_erractionsave(err_act)
  info = psb_success_
  call a%mv_to_coo(tmp,info)
  if (info == psb_success_) call tmp%transc()
  if (info == psb_success_) call a%mv_from_coo(tmp,info)

  if (info /= psb_success_) then 
    info = psb_err_missing_override_method_ 
    call psb_errpush(info,name,a_err=a%get_fmt())
    goto 9999
  end if
  call psb_erractionrestore(err_act) 

  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_z_base_transc_1mat


! == ==================================
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
! == ==================================

subroutine psb_z_base_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_csmm
  use psb_error_mod

  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  complex(psb_dpk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_base_csmm'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_csmm


subroutine psb_z_base_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_csmv
  use psb_error_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  complex(psb_dpk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_base_csmv'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return


end subroutine psb_z_base_csmv


subroutine psb_z_base_inner_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_inner_cssm
  use psb_error_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  complex(psb_dpk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_base_inner_cssm'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_inner_cssm


subroutine psb_z_base_inner_cssv(alpha,a,x,beta,y,info,trans) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_inner_cssv
  use psb_error_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  complex(psb_dpk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_base_inner_cssv'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_inner_cssv


subroutine psb_z_base_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_cssm
  use psb_error_mod
  use psb_string_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)    :: alpha, beta, x(:,:)
  complex(psb_dpk_), intent(inout) :: y(:,:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  complex(psb_dpk_), intent(in), optional :: d(:)

  complex(psb_dpk_), allocatable :: tmp(:,:)
  integer(psb_ipk_) :: err_act, nar,nac,nc, i
  character(len=1) :: scale_
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  nar = a%get_nrows()
  nac = a%get_ncols()
  nc = min(size(x,2), size(y,2))
  if (size(x,1) < nac) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = nac; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (size(y,1) < nar) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = nar; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (.not. (a%is_triangle())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(d)) then 
    if (present(scale)) then 
      scale_ = scale
    else
      scale_ = 'L'
    end if

    if (psb_toupper(scale_) == 'R') then 
      if (size(d,1) < nac) then
        info = psb_err_input_asize_small_i_
        ierr(1) = 9; ierr(2) = nac; 
        call psb_errpush(info,name,i_err=ierr)
        goto 9999
      end if

      allocate(tmp(nac,nc),stat=info) 
      if (info /= psb_success_) info = psb_err_alloc_dealloc_ 
      if (info == psb_success_) then 
        do i=1, nac
          tmp(i,1:nc) = d(i)*x(i,1:nc) 
        end do
      end if
      if (info == psb_success_)&
           & call a%inner_cssm(alpha,tmp,beta,y,info,trans)

      if (info == psb_success_) then 
        deallocate(tmp,stat=info) 
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
      end if

    else if (psb_toupper(scale_) == 'L') then 

      if (size(d,1) < nar) then
        info = psb_err_input_asize_small_i_
             ierr(1) = 9; ierr(2) = nar; 
        call psb_errpush(info,name,i_err=ierr)
        goto 9999
      end if

      allocate(tmp(nar,nc),stat=info) 
      if (info /= psb_success_) info = psb_err_alloc_dealloc_ 
      if (info == psb_success_)&
           & call a%inner_cssm(zone,x,zzero,tmp,info,trans)

      if (info == psb_success_)then 
        do i=1, nar
          tmp(i,1:nc) = d(i)*tmp(i,1:nc) 
        end do
      end if
      if (info == psb_success_)&
           & call psb_geaxpby(nar,nc,alpha,tmp,beta,y,info)

      if (info == psb_success_) then 
        deallocate(tmp,stat=info) 
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
      end if

    else
      info = 31
      ierr(1) = 8; ierr(2) = izero; 
      call psb_errpush(info,name,i_err=ierr,a_err=scale_)
      goto 9999
    end if
  else 
    ! Scale is ignored in this case 
    call a%inner_cssm(alpha,x,beta,y,info,trans)
  end if

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='inner_cssm')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_z_base_cssm


subroutine psb_z_base_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_cssv
  use psb_error_mod
  use psb_string_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)    :: alpha, beta, x(:)
  complex(psb_dpk_), intent(inout) :: y(:)
  integer(psb_ipk_), intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  complex(psb_dpk_), intent(in), optional :: d(:)

  complex(psb_dpk_), allocatable :: tmp(:)
  integer(psb_ipk_) :: err_act, nar,nac,nc, i
  character(len=1) :: scale_
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  nar = a%get_nrows()
  nac = a%get_ncols()
  nc = 1
  if (size(x,1) < nac) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = nac; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (size(y,1) < nar) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = nar; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (.not. (a%is_triangle())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(d)) then 
    if (present(scale)) then 
      scale_ = scale
    else
      scale_ = 'L'
    end if

    if (psb_toupper(scale_) == 'R') then 
      if (size(d,1) < nac) then
        info = psb_err_input_asize_small_i_
        ierr(1) = 9; ierr(2) = nac; 
        call psb_errpush(info,name,i_err=ierr)
        goto 9999
      end if

      allocate(tmp(nac),stat=info) 
      if (info /= psb_success_) info = psb_err_alloc_dealloc_ 
      if (info == psb_success_) call inner_vscal(nac,d,x,tmp) 
      if (info == psb_success_)&
           & call a%inner_cssm(alpha,tmp,beta,y,info,trans)

      if (info == psb_success_) then 
        deallocate(tmp,stat=info) 
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
      end if

    else if (psb_toupper(scale_) == 'L') then 
      if (size(d,1) < nar) then
        info = psb_err_input_asize_small_i_
             ierr(1) = 9; ierr(2) = nar; 
        call psb_errpush(info,name,i_err=ierr)
        goto 9999
      end if

      if (beta == zzero) then 
        call a%inner_cssm(alpha,x,zzero,y,info,trans)
        if (info == psb_success_)  call inner_vscal1(nar,d,y)
      else
        allocate(tmp(nar),stat=info) 
        if (info /= psb_success_) info = psb_err_alloc_dealloc_ 
        if (info == psb_success_)&
             & call a%inner_cssm(alpha,x,zzero,tmp,info,trans)

        if (info == psb_success_)  call inner_vscal1(nar,d,tmp)
        if (info == psb_success_)&
             & call psb_geaxpby(nar,zone,tmp,beta,y,info)
        if (info == psb_success_) then 
          deallocate(tmp,stat=info) 
          if (info /= psb_success_) info = psb_err_alloc_dealloc_
        end if
      end if

    else
      info = 31
      ierr(1) = 8; ierr(2) = izero; 
      call psb_errpush(info,name,i_err=ierr,a_err=scale_)
      goto 9999
    end if
  else 
    ! Scale is ignored in this case 
    call a%inner_cssm(alpha,x,beta,y,info,trans)
  end if

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='inner_cssm')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return
contains
  subroutine inner_vscal(n,d,x,y)
    implicit none 
    integer(psb_ipk_), intent(in)         :: n
    complex(psb_dpk_), intent(in)  :: d(*),x(*)
    complex(psb_dpk_), intent(out) :: y(*)
    integer(psb_ipk_) :: i

    do i=1,n
      y(i) = d(i)*x(i) 
    end do
  end subroutine inner_vscal


  subroutine inner_vscal1(n,d,x)
    implicit none 
    integer(psb_ipk_), intent(in)         :: n
    complex(psb_dpk_), intent(in)  :: d(*)
    complex(psb_dpk_), intent(inout) :: x(*)
    integer(psb_ipk_) :: i

    do i=1,n
      x(i) = d(i)*x(i) 
    end do
  end subroutine inner_vscal1

end subroutine psb_z_base_cssv


subroutine psb_z_base_scals(d,a,info) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_scals
  use psb_error_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  complex(psb_dpk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_scals'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_scals



subroutine psb_z_base_scal(d,a,info,side) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_scal
  use psb_error_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(inout) :: a
  complex(psb_dpk_), intent(in)      :: d(:)
  integer(psb_ipk_), intent(out)            :: info
  character, intent(in), optional :: side

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_scal'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  return

end subroutine psb_z_base_scal



function psb_z_base_maxval(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_maxval

  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='maxval'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  res = dzero

  return

end function psb_z_base_maxval

function psb_z_base_csnmi(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_csnmi

  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csnmi'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  res = dzero

  return

end function psb_z_base_csnmi

function psb_z_base_csnm1(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_csnm1

  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='csnm1'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if
  res = dzero

  return

end function psb_z_base_csnm1

subroutine psb_z_base_rowsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_rowsum
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_z_base_rowsum

subroutine psb_z_base_arwsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_arwsum
  class(psb_z_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='arwsum'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_z_base_arwsum

subroutine psb_z_base_colsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_colsum
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_z_base_colsum

subroutine psb_z_base_aclsum(d,a) 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_aclsum
  class(psb_z_base_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: err_act, info
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='aclsum'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_z_base_aclsum


subroutine psb_z_base_get_diag(a,d,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_get_diag

  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(out)     :: d(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  ! This is the base version. If we get here
  ! it means the derived class is incomplete,
  ! so we throw an error.
  info = psb_err_missing_override_method_
  call psb_errpush(info,name,a_err=a%get_fmt())

  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_z_base_get_diag


! == ==================================
!
!
!
! Computational routines for z_VECT
! variables. If the actual data type is 
! a "normal" one, these are sufficient. 
! 
!
!
!
! == ==================================



subroutine psb_z_base_vect_mv(alpha,a,x,beta,y,info,trans) 
  use psb_error_mod
  use psb_const_mod
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_vect_mv
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)       :: alpha, beta
  class(psb_z_base_vect_type), intent(inout) :: x
  class(psb_z_base_vect_type), intent(inout) :: y
  integer(psb_ipk_), intent(out)             :: info
  character, optional, intent(in)  :: trans

  ! For the time being we just throw everything back
  ! onto the normal routines. 
  call x%sync()
  call y%sync()
  call a%csmm(alpha,x%v,beta,y%v,info,trans)
  
end subroutine psb_z_base_vect_mv

subroutine psb_z_base_vect_cssv(alpha,a,x,beta,y,info,trans,scale,d)
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_vect_cssv
  use psb_z_base_vect_mod
  use psb_error_mod
  use psb_string_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)       :: alpha, beta
  class(psb_z_base_vect_type), intent(inout) :: x,y
  integer(psb_ipk_), intent(out)             :: info
  character, optional, intent(in)  :: trans, scale
  class(psb_z_base_vect_type), intent(inout),optional  :: d

  complex(psb_dpk_), allocatable :: tmp(:)
  class(psb_z_base_vect_type), allocatable :: tmpv
  integer(psb_ipk_) :: err_act, nar,nac,nc, i
  character(len=1) :: scale_
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  nar = a%get_nrows()
  nac = a%get_ncols()
  nc = 1
  if (x%get_nrows() < nac) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = nac; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (y%get_nrows() < nar) then
    info = psb_err_input_asize_small_i_
    ierr(1) = 3; ierr(2) = nar; 
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if

  if (.not. (a%is_triangle())) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  call x%sync() 
  call y%sync()
  if (present(d)) then 
    call d%sync()
    if (present(scale)) then 
      scale_ = scale
    else
      scale_ = 'L'
    end if

    if (psb_toupper(scale_) == 'R') then 
      if (d%get_nrows() < nac) then
        info = psb_err_input_asize_small_i_
        ierr(1) = 9; ierr(2) = nac; 
        call psb_errpush(info,name,i_err=ierr)
        goto 9999
      end if
      
      allocate(tmpv, mold=y,stat=info)
      !    allocate(tmp(nac),stat=info) 
      if (info /= psb_success_) info = psb_err_alloc_dealloc_ 
      if (info == psb_success_) call tmpv%mlt(zone,d%v(1:nac),x,zzero,info) 
      if (info == psb_success_)&
           & call a%inner_cssm(alpha,tmpv,beta,y,info,trans)

      if (info == psb_success_) then 
        call tmpv%free(info)
        if (info == psb_success_) deallocate(tmpv,stat=info) 
        if (info /= psb_success_) info = psb_err_alloc_dealloc_
      end if

    else if (psb_toupper(scale_) == 'L') then 
      if (d%get_nrows() < nar) then
        info = psb_err_input_asize_small_i_
        ierr(1) = 9; ierr(2) = nar; 
        call psb_errpush(info,name,i_err=ierr)
        goto 9999
      end if

      if (beta == zzero) then 
        call a%inner_cssm(alpha,x,zzero,y,info,trans)
        if (info == psb_success_)  call y%mlt(d%v(1:nar),info)
!!$        if (info == psb_success_)  call inner_vscal1(nar,d,y)
      else
        !        allocate(tmp(nar),stat=info) 
        allocate(tmpv, mold=y,stat=info)
        if (info /= psb_success_) info = psb_err_alloc_dealloc_ 
        if (info == psb_success_)&
             & call a%inner_cssm(alpha,x,zzero,tmpv,info,trans)

        if (info == psb_success_)  call tmpv%mlt(d%v(1:nar),info)
        if (info == psb_success_)&
             & call y%axpby(nar,zone,tmpv,beta,info)
        if (info == psb_success_) then
          call tmpv%free(info)
          if (info == psb_success_) deallocate(tmpv,stat=info) 
          if (info /= psb_success_) info = psb_err_alloc_dealloc_
        end if
      end if

    else
      info = 31
      ierr(1) = 8; ierr(2) = izero; 
      call psb_errpush(info,name,i_err=ierr,a_err=scale_)
      goto 9999
    end if
  else 
    ! Scale is ignored in this case 
    call a%inner_cssm(alpha,x,beta,y,info,trans)
  end if

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='inner_cssm')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_z_base_vect_cssv


subroutine psb_z_base_inner_vect_sv(alpha,a,x,beta,y,info,trans) 
  use psb_z_base_mat_mod, psb_protect_name => psb_z_base_inner_vect_sv
  use psb_error_mod
  use psb_string_mod
  use psb_z_base_vect_mod
  implicit none 
  class(psb_z_base_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)       :: alpha, beta
  class(psb_z_base_vect_type), intent(inout) :: x, y
  integer(psb_ipk_), intent(out)             :: info
  character, optional, intent(in)  :: trans

  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ierr(5)
  character(len=20)  :: name='z_base_inner_vect_sv'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)

  call a%inner_cssm(alpha,x%v,beta,y%v,info,trans)  

  if (info /= psb_success_) then 
    info = psb_err_from_subroutine_ 
    call psb_errpush(info,name, a_err='inner_cssm')
    goto 9999
  end if


  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return

end subroutine psb_z_base_inner_vect_sv
