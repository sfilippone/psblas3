  
!> Function  csmv:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Product by a dense rank 1 array.
!!
!!        Compute
!!           Y = alpha*op(A)*X + beta*Y
!!
!! \param alpha  Scaling factor for Ax
!! \param A      the input sparse matrix
!! \param x(:)   the input dense X
!! \param beta   Scaling factor for y
!! \param y(:)   the input/output dense Y
!! \param info   return code
!! \param trans  [N] Whether to use A (N), its transpose (T)
!!               or its conjugate transpose (C)
!!
!
subroutine psb_s_dns_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_csmv
  implicit none 
  class(psb_s_dns_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)          :: alpha, beta, x(:)
  real(psb_spk_), intent(inout)       :: y(:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans
  !
  character :: trans_
  integer(psb_ipk_)  :: err_act, m, n, lda
  character(len=20)  :: name='s_dns_csmv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(trans)) then
    trans_ = psb_toupper(trans)
  else
    trans_ = 'N'
  end if

  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  
  if (a%is_dev()) call a%sync()
  if (trans_ == 'N') then
    m=a%get_nrows()
    n=a%get_ncols()
  else
    n=a%get_nrows()
    m=a%get_ncols()
  end if
  lda = size(a%val,1)


  call sgemv(trans_,a%get_nrows(),a%get_ncols(),alpha,&
       & a%val,size(a%val,1),x,1,beta,y,1)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_s_dns_csmv

  
!> Function  csmm:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Product by a dense rank 2 array.
!!
!!        Compute
!!           Y = alpha*op(A)*X + beta*Y
!!
!! \param alpha  Scaling factor for Ax
!! \param A      the input sparse matrix
!! \param x(:,:)   the input dense X
!! \param beta   Scaling factor for y
!! \param y(:,:)   the input/output dense Y
!! \param info   return code
!! \param trans  [N] Whether to use A (N), its transpose (T)
!!               or its conjugate transpose (C)
!!
!
subroutine psb_s_dns_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_csmm
  implicit none 
  class(psb_s_dns_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_spk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans
  !
  character :: trans_
  integer(psb_ipk_)  :: err_act,m,n,k, lda, ldx, ldy
  character(len=20)  :: name='s_dns_csmm'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_

  if (present(trans)) then
    trans_ = trans
  else
    trans_ = 'N'
  end if
 
  if (.not.a%is_asb()) then 
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (a%is_dev()) call a%sync()
  if (psb_toupper(trans_)=='N') then
    m = a%get_nrows()
    k = a%get_ncols()
    n = min(size(y,2),size(x,2))
  else
    k = a%get_nrows()
    m = a%get_ncols()
    n = min(size(y,2),size(x,2))    
  end if
  lda = size(a%val,1)
  ldx = size(x,1)
  ldy = size(y,1)
  call sgemm(trans_,'N',m,n,k,alpha,a%val,lda,x,ldx,beta,y,ldy)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine psb_s_dns_csmm



!
!
!> Function  csnmi:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Operator infinity norm
!!     CSNMI = MAXVAL(SUM(ABS(A(:,:)),dim=2))
!! 
!
function psb_s_dns_csnmi(a) result(res)
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_csnmi
  implicit none 
  class(psb_s_dns_sparse_mat), intent(in) :: a
  real(psb_spk_)         :: res
  !
  integer(psb_ipk_) :: i
  real(psb_spk_) :: acc

  res = szero 
  if (a%is_dev()) call a%sync()
 
  do i = 1, a%get_nrows()
    acc = sum(abs(a%val(i,:)))
    res = max(res,acc)
  end do

end function psb_s_dns_csnmi


!
!> Function  get_diag:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Extract the diagonal of A. 
!!        
!!   D(i) = A(i:i), i=1:min(nrows,ncols)
!!
!! \param d(:)  The output diagonal
!! \param info  return code. 
! 
subroutine psb_s_dns_get_diag(a,d,info) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_get_diag
  implicit none 
  class(psb_s_dns_sparse_mat), intent(in) :: a
  real(psb_spk_), intent(out)     :: d(:)
  integer(psb_ipk_), intent(out)            :: info
  !
  integer(psb_ipk_) :: err_act, mnm, i
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.

  info  = psb_success_
  call psb_erractionsave(err_act)
  if (a%is_dev()) call a%sync()

  mnm = min(a%get_nrows(),a%get_ncols())
  if (size(d) < mnm) then 
    info=psb_err_input_asize_invalid_i_
    call psb_errpush(info,name,i_err=(/2_psb_ipk_,size(d,kind=psb_ipk_)/))
    goto 9999
  end if


  do i=1, mnm
    d(i) = a%val(i,i)
  end do
  do i=mnm+1,size(d) 
    d(i) = szero
  end do

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine psb_s_dns_get_diag


!         
!
!> Function  reallocate_nz
!! \memberof  psb_s_dns_sparse_mat
!! \brief One--parameters version of (re)allocate
!!
!!  \param nz  number of nonzeros to allocate for
!!             i.e. makes sure that the internal storage
!!             allows for NZ coefficients and their indices. 
!  
subroutine  psb_s_dns_reallocate_nz(nz,a) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_reallocate_nz
  implicit none 
  integer(psb_ipk_), intent(in) :: nz
  class(psb_s_dns_sparse_mat), intent(inout) :: a
  !
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='s_dns_reallocate_nz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)

  !
  ! This is a no-op, allocation is fixed.
  ! 
  if (a%is_dev()) call a%sync()


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine psb_s_dns_reallocate_nz

!
!> Function  mold:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Allocate a class(psb_s_dns_sparse_mat) with the
!!     same dynamic type as the input.
!!     This is equivalent to allocate(  mold=  ) and is provided
!!     for those compilers not yet supporting mold.
!!   \param b The output variable
!!   \param info return code
! 
subroutine psb_s_dns_mold(a,b,info) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_mold
  implicit none 
  class(psb_s_dns_sparse_mat), intent(in)  :: a
  class(psb_s_base_sparse_mat), intent(inout), allocatable  :: b
  integer(psb_ipk_), intent(out)                  :: info
  !
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='dns_mold'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  
  allocate(psb_s_dns_sparse_mat :: b, stat=info)
  
  if (info /= 0) then 
    info = psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return
9999 call psb_error_handler(err_act)
  return
  
end subroutine psb_s_dns_mold

!         
!
!> Function  allocate_mnnz
!! \memberof  psb_s_dns_sparse_mat
!! \brief Three-parameters version of allocate
!!
!!  \param m  number of rows
!!  \param n  number of cols
!!  \param nz [estimated internally] number of nonzeros to allocate for
!
subroutine  psb_s_dns_allocate_mnnz(m,n,a,nz) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_allocate_mnnz
  implicit none 
  integer(psb_ipk_), intent(in) :: m,n
  class(psb_s_dns_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(in), optional :: nz
  !
  integer(psb_ipk_) :: err_act, info, nz_
  character(len=20)  :: name='allocate_mnz'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (m < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/1_psb_ipk_/))
    goto 9999
  endif
  if (n < 0) then 
    info = psb_err_iarg_neg_
    call psb_errpush(info,name,i_err=(/2_psb_ipk_/))
    goto 9999
  endif
  

  ! Basic stuff common to all formats
  call a%set_nrows(m)
  call a%set_ncols(n)
  call a%set_triangle(.false.)
  call a%set_unit(.false.)
  call a%set_dupl(psb_dupl_def_)
  call a%set_bld()
  call a%set_host()
  
  ! We ignore NZ in this case.

  call psb_realloc(m,n,a%val,info)
  if (info == psb_success_) then 
    a%val = szero
    a%nnz = 0
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine psb_s_dns_allocate_mnnz


!
!
!
!> Function  csgetrow:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Get a (subset of) row(s)
!!        
!!        getrow is the basic method by which the other (getblk, clip) can
!!        be implemented.
!!        
!!        Returns the set
!!           NZ, IA(1:nz), JA(1:nz), VAL(1:NZ)
!!         each identifying the position of a nonzero in A
!!         i.e.
!!           VAL(1:NZ) = A(IA(1:NZ),JA(1:NZ))
!!         with IMIN<=IA(:)<=IMAX
!!         with JMIN<=JA(:)<=JMAX
!!         IA,JA are reallocated as necessary.
!!         
!!  \param imin  the minimum row index we are interested in 
!!  \param imax  the minimum row index we are interested in 
!!  \param nz the number of output coefficients
!!  \param ia(:)  the output row indices
!!  \param ja(:)  the output col indices
!!  \param val(:)  the output coefficients
!!  \param info  return code
!!  \param jmin [1] minimum col index 
!!  \param jmax [a\%get_ncols()] maximum col index 
!!  \param iren(:) [none] an array to return renumbered indices (iren(ia(:)),iren(ja(:))
!!  \param rscale [false] map [min(ia(:)):max(ia(:))] onto [1:max(ia(:))-min(ia(:))+1]
!!  \param cscale [false] map [min(ja(:)):max(ja(:))] onto [1:max(ja(:))-min(ja(:))+1]
!!          ( iren cannot be specified with rscale/cscale)
!!  \param append [false] append to ia,ja 
!!  \param nzin [none]  if append, then first new entry should go in entry nzin+1
!!           
!
subroutine psb_s_dns_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale,chksz)
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_csgetrow
  implicit none

  class(psb_s_dns_sparse_mat), intent(in)      :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_spk_), allocatable,  intent(inout)   :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale,chksz
  !
  logical :: append_, rscale_, cscale_, chksz_ 
  integer(psb_ipk_) :: nzin_, jmin_, jmax_, err_act, i,j,k
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  if (a%is_dev()) call a%sync()
  
  if (present(jmin)) then
    jmin_ = jmin
  else
    jmin_ = 1
  endif
  if (present(jmax)) then
    jmax_ = jmax
  else
    jmax_ = a%get_ncols()
  endif

  if ((imax<imin).or.(jmax_<jmin_)) then 
    nz = 0
    return
  end if

  if (present(append)) then
    append_=append
  else
    append_=.false.
  endif
  if ((append_).and.(present(nzin))) then 
    nzin_ = nzin
  else
    nzin_ = 0
  endif
  if (present(rscale)) then 
    rscale_ = rscale
  else
    rscale_ = .false.
  endif
  if (present(cscale)) then 
    cscale_ = cscale
  else
    cscale_ = .false.
  endif
  if (present(chksz)) then 
    chksz_ = chksz
  else
    chksz_ = .true.
  endif

  if ((rscale_.or.cscale_).and.(present(iren))) then 
    info = psb_err_many_optional_arg_
    call psb_errpush(info,name,a_err='iren (rscale.or.cscale)')
    goto 9999
  end if

  if (append) then 
    write(0,*) 'APPEND=TRUE NOT IMPLEMENTED'
    info = -1
    call psb_errpush(info,name,a_err='not impl')
    goto 9999
  end if
  nz = count(a%val(imin:imax,jmin_:jmax_) /= szero)

  if (chksz_) then 
    call psb_ensure_size(nzin_+nz,ia,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nz,ja,info)
    if (info == psb_success_) call psb_ensure_size(nzin_+nz,val,info)
    if (info /= psb_success_) goto 9999
  end if

  k = 0
  do i=imin,imax
    do j=jmin_,jmax_
      if (a%val(i,j) /= szero) then 
        k = k + 1
        ia(k) = i
        ja(k) = j
        val(k) = a%val(i,j)
      end if
    end do
  end do

  if (rscale_) then 
    do i=nzin_+1, nzin_+nz
      ia(i) = ia(i) - imin + 1
    end do
  end if
  if (cscale_) then 
    do i=nzin_+1, nzin_+nz
      ja(i) = ja(i) - jmin_ + 1
    end do
  end if

  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine psb_s_dns_csgetrow


!> Function  trim 
!! \memberof  psb_s_dns_sparse_mat
!! \brief Memory trim
!! Make sure the memory allocation of the sparse matrix is as tight as
!! possible given the actual number of nonzeros it contains. 
!
subroutine  psb_s_dns_trim(a)
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_dns_trim
  implicit none 
  class(psb_s_dns_sparse_mat), intent(inout) :: a
  !
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  ! Do nothing, we are already at minimum memory.

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine psb_s_dns_trim

!
!> Function  cp_from_coo:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Copy and convert from psb_s_coo_sparse_mat
!!        Invoked from the target object.
!!   \param b The input variable
!!   \param info return code
!  

subroutine psb_s_cp_dns_from_coo(a,b,info) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_cp_dns_from_coo
  implicit none 

  class(psb_s_dns_sparse_mat), intent(inout) :: a
  class(psb_s_coo_sparse_mat), intent(in)    :: b
  integer(psb_ipk_), intent(out)               :: info
  !
  type(psb_s_coo_sparse_mat)   :: tmp
  integer(psb_ipk_)             :: nza, nr, i,err_act, nc
  integer(psb_ipk_), parameter  :: maxtry=8
  integer(psb_ipk_)             :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

  if (.not.b%is_by_rows()) then 
    ! This is to have fix_coo called behind the scenes
    call b%cp_to_coo(tmp,info)
    call tmp%fix(info)
    if (info /= psb_success_) return
    
    nr  = tmp%get_nrows()
    nc  = tmp%get_ncols()
    nza = tmp%get_nzeros()
    ! If it is sorted then we can lessen memory impact 
    a%psb_s_base_sparse_mat = tmp%psb_s_base_sparse_mat
    
    call psb_realloc(nr,nc,a%val,info) 
    if (info /= 0) goto 9999
    a%val = szero
    do i=1, nza
      a%val(tmp%ia(i),tmp%ja(i)) = tmp%val(i)
    end do
    a%nnz = nza
    call tmp%free()
  else
    if (b%is_dev()) call b%sync()
    nr  = b%get_nrows()
    nc  = b%get_ncols()
    nza = b%get_nzeros()
    ! If it is sorted then we can lessen memory impact 
    a%psb_s_base_sparse_mat = b%psb_s_base_sparse_mat
    
    call psb_realloc(nr,nc,a%val,info) 
    if (info /= 0) goto 9999
    a%val = szero
    do i=1, nza
      a%val(b%ia(i),b%ja(i)) = b%val(i)
    end do
    a%nnz = nza
  end if
  call a%set_host()

  return

9999 call psb_error_handler(err_act)
  return
  
end subroutine psb_s_cp_dns_from_coo


  
!
!> Function  cp_to_coo:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Copy and convert to psb_s_coo_sparse_mat
!!        Invoked from the source object.
!!   \param b The output variable
!!   \param info return code
!  

subroutine psb_s_cp_dns_to_coo(a,b,info) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_cp_dns_to_coo
  implicit none 

  class(psb_s_dns_sparse_mat), intent(in)  :: a
  class(psb_s_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)               :: info

  !locals
  Integer(Psb_Ipk_)             :: nza, nr, nc,i,j,k,err_act

  info = psb_success_

  if (a%is_dev())   call a%sync()
  nr  = a%get_nrows()
  nc  = a%get_ncols()
  nza = a%get_nzeros()

  call b%allocate(nr,nc,nza)
  b%psb_s_base_sparse_mat = a%psb_s_base_sparse_mat

  k = 0
  do i=1,a%get_nrows()
    do j=1,a%get_ncols()
      if (a%val(i,j) /= szero) then 
        k = k + 1
        b%ia(k) = i
        b%ja(k) = j
        b%val(k) = a%val(i,j)
      end if
    end do
  end do

  call b%set_nzeros(nza)
  call b%set_sort_status(psb_row_major_)
  call b%set_asb()
  call b%set_host()

end subroutine psb_s_cp_dns_to_coo



!
!> Function  mv_to_coo:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Convert to psb_s_coo_sparse_mat, freeing the source.
!!        Invoked from the source object.
!!   \param b The output variable
!!   \param info return code
!  
subroutine psb_s_mv_dns_to_coo(a,b,info) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_mv_dns_to_coo
  implicit none 

  class(psb_s_dns_sparse_mat), intent(inout) :: a
  class(psb_s_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)              :: info

  info = psb_success_

  call a%cp_to_coo(b,info)
  call a%free()
  return

end subroutine psb_s_mv_dns_to_coo


!
!> Function  mv_from_coo:
!! \memberof  psb_s_dns_sparse_mat
!! \brief Convert from psb_s_coo_sparse_mat, freeing the source.
!!        Invoked from the target object.
!!   \param b The input variable
!!   \param info return code
!  
!  
subroutine psb_s_mv_dns_from_coo(a,b,info) 
  use psb_base_mod
  use psb_s_dns_mat_mod, psb_protect_name => psb_s_mv_dns_from_coo
  implicit none 

  class(psb_s_dns_sparse_mat), intent(inout) :: a
  class(psb_s_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)               :: info

  info = psb_success_

  call a%cp_from_coo(b,info)
  call b%free()

  return

end subroutine psb_s_mv_dns_from_coo

