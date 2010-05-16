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

subroutine psb_s_base_cp_to_coo(a,b,info)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_cp_to_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  class(psb_s_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)            :: info

  Integer :: err_act
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

end subroutine psb_s_base_cp_to_coo

subroutine psb_s_base_cp_from_coo(a,b,info)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_cp_from_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  class(psb_s_coo_sparse_mat), intent(in) :: b
  integer, intent(out)            :: info

  Integer :: err_act
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

end subroutine psb_s_base_cp_from_coo


subroutine psb_s_base_cp_to_fmt(a,b,info)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_cp_to_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  class(psb_s_base_sparse_mat), intent(inout) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='to_fmt'
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

end subroutine psb_s_base_cp_to_fmt

subroutine psb_s_base_cp_from_fmt(a,b,info)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_cp_from_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  class(psb_s_base_sparse_mat), intent(in) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='from_fmt'
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

end subroutine psb_s_base_cp_from_fmt


subroutine psb_s_base_mv_to_coo(a,b,info)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_mv_to_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  class(psb_s_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)            :: info

  Integer :: err_act
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

end subroutine psb_s_base_mv_to_coo

subroutine psb_s_base_mv_from_coo(a,b,info)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_mv_from_coo
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  class(psb_s_coo_sparse_mat), intent(inout) :: b
  integer, intent(out)            :: info

  Integer :: err_act
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

end subroutine psb_s_base_mv_from_coo


subroutine psb_s_base_mv_to_fmt(a,b,info)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_mv_to_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  class(psb_s_base_sparse_mat), intent(inout) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='to_fmt'
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

end subroutine psb_s_base_mv_to_fmt

subroutine psb_s_base_mv_from_fmt(a,b,info)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_mv_from_fmt
  use psb_error_mod
  use psb_realloc_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  class(psb_s_base_sparse_mat), intent(inout) :: b
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='from_fmt'
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

end subroutine psb_s_base_mv_from_fmt

subroutine psb_s_base_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_error_mod
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_csput
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  real(psb_spk_), intent(in)      :: val(:)
  integer, intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer, intent(out)            :: info
  integer, intent(in), optional   :: gtl(:)

  Integer :: err_act
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

end subroutine psb_s_base_csput

subroutine psb_s_base_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_csgetrow
  implicit none

  class(psb_s_base_sparse_mat), intent(in) :: a
  integer, intent(in)                  :: imin,imax
  integer, intent(out)                 :: nz
  integer, allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_spk_), allocatable,  intent(inout)    :: val(:)
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale
  Integer :: err_act
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

end subroutine psb_s_base_csgetrow



subroutine psb_s_base_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_csgetblk
  implicit none

  class(psb_s_base_sparse_mat), intent(in) :: a
  class(psb_s_coo_sparse_mat), intent(inout) :: b
  integer, intent(in)                  :: imin,imax
  integer,intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer, intent(in), optional        :: iren(:)
  integer, intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale
  Integer :: err_act, nzin, nzout
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

end subroutine psb_s_base_csgetblk


subroutine psb_s_base_csclip(a,b,info,&
     & imin,imax,jmin,jmax,rscale,cscale)
  ! Output is always in  COO format 
  use psb_error_mod
  use psb_const_mod
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_csclip
  implicit none

  class(psb_s_base_sparse_mat), intent(in) :: a
  class(psb_s_coo_sparse_mat), intent(out) :: b
  integer,intent(out)                  :: info
  integer, intent(in), optional        :: imin,imax,jmin,jmax
  logical, intent(in), optional        :: rscale,cscale

  Integer :: err_act, nzin, nzout, imin_, imax_, jmin_, jmax_, mb,nb
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

end subroutine psb_s_base_csclip


subroutine psb_s_base_transp_2mat(a,b)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_transp_2mat
  use psb_error_mod
  implicit none 

  class(psb_s_base_sparse_mat), intent(out) :: a
  class(psb_base_sparse_mat), intent(in)   :: b

  type(psb_s_coo_sparse_mat) :: tmp
  integer err_act, info
  character(len=*), parameter :: name='s_base_transp'

  call psb_erractionsave(err_act)

  info = psb_success_
  select type(b)
    class is (psb_s_base_sparse_mat)
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call tmp%transp()
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
    class default
    info = psb_err_missing_override_method_
  end select
  if (info /= psb_success_) then 
    call psb_errpush(info,name,a_err=b%get_fmt())
    goto 9999
  end if
  call psb_erractionrestore(err_act) 

  return
9999 continue
  if (err_act /= psb_act_ret_) then
    call psb_error()
  end if

  return

end subroutine psb_s_base_transp_2mat

subroutine psb_s_base_transc_2mat(a,b)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_transc_2mat
  implicit none 

  class(psb_s_base_sparse_mat), intent(out) :: a
  class(psb_base_sparse_mat), intent(in)   :: b

  call a%transp(b) 
end subroutine psb_s_base_transc_2mat

subroutine psb_s_base_transp_1mat(a)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_transp_1mat
  use psb_error_mod
  implicit none 

  class(psb_s_base_sparse_mat), intent(inout) :: a

  type(psb_s_coo_sparse_mat) :: tmp
  integer :: err_act, info
  character(len=*), parameter :: name='s_base_transp'

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

end subroutine psb_s_base_transp_1mat

subroutine psb_s_base_transc_1mat(a)
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_transc_1mat
  implicit none 

  class(psb_s_base_sparse_mat), intent(inout) :: a

  call a%transp() 
end subroutine psb_s_base_transc_1mat


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

subroutine psb_s_base_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_csmm
  use psb_error_mod

  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_spk_), intent(inout) :: y(:,:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans

  Integer :: err_act
  character(len=20)  :: name='s_base_csmm'
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

end subroutine psb_s_base_csmm


subroutine psb_s_base_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_csmv
  use psb_error_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)    :: alpha, beta, x(:)
  real(psb_spk_), intent(inout) :: y(:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans

  Integer :: err_act
  character(len=20)  :: name='s_base_csmv'
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


end subroutine psb_s_base_csmv


subroutine psb_s_base_inner_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_inner_cssm
  use psb_error_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_spk_), intent(inout) :: y(:,:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans

  Integer :: err_act
  character(len=20)  :: name='s_base_inner_cssm'
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

end subroutine psb_s_base_inner_cssm


subroutine psb_s_base_inner_cssv(alpha,a,x,beta,y,info,trans) 
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_inner_cssv
  use psb_error_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)    :: alpha, beta, x(:)
  real(psb_spk_), intent(inout) :: y(:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans

  Integer :: err_act
  character(len=20)  :: name='s_base_inner_cssv'
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

end subroutine psb_s_base_inner_cssv


subroutine psb_s_base_cssm(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_cssm
  use psb_error_mod
  use psb_string_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)    :: alpha, beta, x(:,:)
  real(psb_spk_), intent(inout) :: y(:,:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  real(psb_spk_), intent(in), optional :: d(:)

  real(psb_spk_), allocatable :: tmp(:,:)
  Integer :: err_act, nar,nac,nc, i
  character(len=1) :: scale_
  character(len=20)  :: name='s_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then 
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  nar = a%get_nrows()
  nac = a%get_ncols()
  nc = min(size(x,2), size(y,2))
  if (size(x,1) < nac) then
    info = 36
    call psb_errpush(info,name,i_err=(/3,nac,0,0,0/))
    goto 9999
  end if
  if (size(y,1) < nar) then
    info = 36
    call psb_errpush(info,name,i_err=(/3,nar,0,0,0/))
    goto 9999
  end if

  if (.not. (a%is_triangle())) then 
    info = 1121
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
        info = 36
        call psb_errpush(info,name,i_err=(/9,nac,0,0,0/))
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
        info = 36
        call psb_errpush(info,name,i_err=(/9,nar,0,0,0/))
        goto 9999
      end if

      allocate(tmp(nar,nc),stat=info) 
      if (info /= psb_success_) info = psb_err_alloc_dealloc_ 
      if (info == psb_success_)&
           & call a%inner_cssm(sone,x,szero,tmp,info,trans)

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
      call psb_errpush(info,name,i_err=(/8,0,0,0,0/),a_err=scale_)
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


  return
  call psb_erractionrestore(err_act)
  return


9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error()
    return
  end if
  return


end subroutine psb_s_base_cssm


subroutine psb_s_base_cssv(alpha,a,x,beta,y,info,trans,scale,d) 
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_cssv
  use psb_error_mod
  use psb_string_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)    :: alpha, beta, x(:)
  real(psb_spk_), intent(inout) :: y(:)
  integer, intent(out)            :: info
  character, optional, intent(in) :: trans, scale
  real(psb_spk_), intent(in), optional :: d(:)

  real(psb_spk_), allocatable :: tmp(:)
  Integer :: err_act, nar,nac,nc, i
  character(len=1) :: scale_
  character(len=20)  :: name='s_cssm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  if (.not.a%is_asb()) then
    info = 1121
    call psb_errpush(info,name)
    goto 9999
  endif

  nar = a%get_nrows()
  nac = a%get_ncols()
  nc = 1
  if (size(x,1) < nac) then
    info = 36
    call psb_errpush(info,name,i_err=(/3,nac,0,0,0/))
    goto 9999
  end if
  if (size(y,1) < nar) then
    info = 36
    call psb_errpush(info,name,i_err=(/3,nar,0,0,0/))
    goto 9999
  end if

  if (.not. (a%is_triangle())) then 
    info = 1121
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
        info = 36
        call psb_errpush(info,name,i_err=(/9,nac,0,0,0/))
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
        info = 36
        call psb_errpush(info,name,i_err=(/9,nar,0,0,0/))
        goto 9999
      end if

      if (beta == szero) then 
        call a%inner_cssm(alpha,x,szero,y,info,trans)
        if (info == psb_success_)  call inner_vscal1(nar,d,y)
      else
        allocate(tmp(nar),stat=info) 
        if (info /= psb_success_) info = psb_err_alloc_dealloc_ 
        if (info == psb_success_)&
             & call a%inner_cssm(alpha,x,szero,tmp,info,trans)

        if (info == psb_success_)  call inner_vscal1(nar,d,tmp)
        if (info == psb_success_)&
             & call psb_geaxpby(nar,sone,tmp,beta,y,info)
        if (info == psb_success_) then 
          deallocate(tmp,stat=info) 
          if (info /= psb_success_) info = psb_err_alloc_dealloc_
        end if
      end if

    else
      info = 31
      call psb_errpush(info,name,i_err=(/8,0,0,0,0/),a_err=scale_)
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


  return
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
    integer, intent(in)         :: n
    real(psb_spk_), intent(in)  :: d(*),x(*)
    real(psb_spk_), intent(out) :: y(*)
    integer :: i

    do i=1,n
      y(i) = d(i)*x(i) 
    end do
  end subroutine inner_vscal


  subroutine inner_vscal1(n,d,x)
    implicit none 
    integer, intent(in)         :: n
    real(psb_spk_), intent(in)  :: d(*)
    real(psb_spk_), intent(inout) :: x(*)
    integer :: i

    do i=1,n
      x(i) = d(i)*x(i) 
    end do
  end subroutine inner_vscal1

end subroutine psb_s_base_cssv


subroutine psb_s_base_scals(d,a,info) 
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_scals
  use psb_error_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  real(psb_spk_), intent(in)      :: d
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='s_scals'
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

end subroutine psb_s_base_scals



subroutine psb_s_base_scal(d,a,info) 
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_scal
  use psb_error_mod
  implicit none 
  class(psb_s_base_sparse_mat), intent(inout) :: a
  real(psb_spk_), intent(in)      :: d(:)
  integer, intent(out)            :: info

  Integer :: err_act
  character(len=20)  :: name='s_scal'
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

end subroutine psb_s_base_scal



function psb_s_base_csnmi(a) result(res)
  use psb_error_mod
  use psb_const_mod
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_csnmi

  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  real(psb_spk_)         :: res

  Integer :: err_act, info
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
  res = -sone

  return

end function psb_s_base_csnmi

subroutine psb_s_base_get_diag(a,d,info) 
  use psb_error_mod
  use psb_const_mod
  use psb_s_base_mat_mod, psb_protect_name => psb_s_base_get_diag

  implicit none 
  class(psb_s_base_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(out)     :: d(:)
  integer, intent(out)            :: info

  Integer :: err_act
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

end subroutine psb_s_base_get_diag




