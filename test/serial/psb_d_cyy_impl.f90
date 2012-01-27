!
! A minimal, non functional implementation of a matrix type module. 
!  
subroutine psb_d_cyy_csmv(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_csmv
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_cyy_csmv'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
  return
end subroutine psb_d_cyy_csmv

subroutine psb_d_cyy_csmm(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_csmm
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_), allocatable  :: acc(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_cyy_csmm'
  logical, parameter :: debug=.false.

  info = psb_success_
end subroutine psb_d_cyy_csmm

subroutine psb_d_cyy_cssv(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_cssv
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:)
  real(psb_dpk_), intent(inout)       :: y(:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_cyy_cssv'
  logical, parameter :: debug=.false.

  info = psb_success_
end subroutine psb_d_cyy_cssv

subroutine psb_d_cyy_cssm(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_cssm
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(in)          :: alpha, beta, x(:,:)
  real(psb_dpk_), intent(inout)       :: y(:,:)
  integer(psb_ipk_), intent(out)                :: info
  character, optional, intent(in)     :: trans

  character :: trans_
  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: tmp(:,:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_cyy_cssm'
  logical, parameter :: debug=.false.

  info = psb_success_

end subroutine psb_d_cyy_cssm

function psb_d_cyy_csnmi(a) result(res)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_csnmi
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: i,j,k,m,n, nr, ir, jc, nc
  real(psb_dpk_) :: acc
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_csnmi'
  logical, parameter :: debug=.false.
  res = dzero 
end function psb_d_cyy_csnmi

function psb_d_cyy_csnm1(a) result(res)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_csnm1

  implicit none 
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_)         :: res

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc, info
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_cyy_csnm1'
  logical, parameter :: debug=.false.
  res = -done 
  return
end function psb_d_cyy_csnm1

subroutine psb_d_cyy_rowsum(d,a) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_rowsum
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)             :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.
  return
end subroutine psb_d_cyy_rowsum

subroutine psb_d_cyy_arwsum(d,a) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_arwsum
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
  character(len=20)  :: name='rowsum'
  logical, parameter :: debug=.false.
end subroutine psb_d_cyy_arwsum

subroutine psb_d_cyy_colsum(d,a) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_colsum
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
  character(len=20)  :: name='colsum'
  logical, parameter :: debug=.false.
end subroutine psb_d_cyy_colsum

subroutine psb_d_cyy_aclsum(d,a) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_aclsum
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)              :: d(:)

  integer(psb_ipk_) :: i,j,k,m,n, nnz, ir, jc, nc
  real(psb_dpk_) :: acc
  real(psb_dpk_), allocatable :: vt(:)
  logical   :: tra
  integer(psb_ipk_) :: err_act, info, int_err(5)
  character(len=20)  :: name='aclsum'
  logical, parameter :: debug=.false.
end subroutine psb_d_cyy_aclsum

subroutine psb_d_cyy_get_diag(a,d,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_get_diag
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(in) :: a
  real(psb_dpk_), intent(out)     :: d(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act, mnm, i, j, k
  character(len=20)  :: name='get_diag'
  logical, parameter :: debug=.false.
  info  = psb_success_
end subroutine psb_d_cyy_get_diag

subroutine psb_d_cyy_scal(d,a,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_scal
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d(:)
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act,mnm, i, j, m
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.
  info  = psb_success_
end subroutine psb_d_cyy_scal

subroutine psb_d_cyy_scals(d,a,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_scals
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: d
  integer(psb_ipk_), intent(out)            :: info

  integer(psb_ipk_) :: err_act,mnm, i, j, m
  character(len=20)  :: name='scal'
  logical, parameter :: debug=.false.
  info  = psb_success_
end subroutine psb_d_cyy_scals

subroutine  psb_d_cyy_reallocate_nz(nz,a) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_reallocate_nz
  implicit none 
  integer(psb_ipk_), intent(in) :: nz
  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='d_cyy_reallocate_nz'
  logical, parameter :: debug=.false.

end subroutine psb_d_cyy_reallocate_nz

subroutine psb_d_cyy_mold(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_mold
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(in)  :: a
  class(psb_d_base_sparse_mat), intent(out), allocatable  :: b
  integer(psb_ipk_), intent(out)                    :: info
  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='reallocate_nz'
  logical, parameter :: debug=.false.
end subroutine psb_d_cyy_mold

subroutine  psb_d_cyy_allocate_mnnz(m,n,a,nz) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_allocate_mnnz
  implicit none 
  integer(psb_ipk_), intent(in) :: m,n
  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(in), optional :: nz
  integer(psb_ipk_) :: err_act, info, nz_
  character(len=20)  :: name='allocate_mnz'
  logical, parameter :: debug=.false.
end subroutine psb_d_cyy_allocate_mnnz

subroutine psb_d_cyy_csgetptn(imin,imax,a,nz,ia,ja,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_csgetptn
  implicit none

  class(psb_d_cyy_sparse_mat), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer(psb_ipk_) :: nzin_, jmin_, jmax_, err_act, i
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

end subroutine psb_d_cyy_csgetptn

subroutine psb_d_cyy_csgetrow(imin,imax,a,nz,ia,ja,val,info,&
     & jmin,jmax,iren,append,nzin,rscale,cscale)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_csgetrow
  implicit none

  class(psb_d_cyy_sparse_mat), intent(in) :: a
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_), intent(out)                 :: nz
  integer(psb_ipk_), allocatable, intent(inout)  :: ia(:), ja(:)
  real(psb_dpk_), allocatable,  intent(inout)    :: val(:)
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax, nzin
  logical, intent(in), optional        :: rscale,cscale

  logical :: append_, rscale_, cscale_ 
  integer(psb_ipk_) :: nzin_, jmin_, jmax_, err_act, i
  character(len=20)  :: name='csget'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
end subroutine psb_d_cyy_csgetrow

subroutine psb_d_cyy_csgetblk(imin,imax,a,b,info,&
     & jmin,jmax,iren,append,rscale,cscale)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_csgetblk
  implicit none

  class(psb_d_cyy_sparse_mat), intent(in) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(in)                  :: imin,imax
  integer(psb_ipk_),intent(out)                  :: info
  logical, intent(in), optional        :: append
  integer(psb_ipk_), intent(in), optional        :: iren(:)
  integer(psb_ipk_), intent(in), optional        :: jmin,jmax
  logical, intent(in), optional        :: rscale,cscale
  integer(psb_ipk_) :: err_act, nzin, nzout
  character(len=20)  :: name='csget'
  logical :: append_
  logical, parameter :: debug=.false.

  info = psb_success_
end subroutine psb_d_cyy_csgetblk

subroutine psb_d_cyy_csput(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_csput
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)      :: val(:)
  integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info
  integer(psb_ipk_), intent(in), optional   :: gtl(:)

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_cyy_csput'
  logical, parameter :: debug=.false.
  integer(psb_ipk_) :: nza, i,j,k, nzl, isza, int_err(5)
  info = psb_success_

end subroutine psb_d_cyy_csput

subroutine psb_d_cyy_reinit(a,clear)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_reinit
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout) :: a   
  logical, intent(in), optional :: clear

  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='reinit'
  logical  :: clear_
  logical, parameter :: debug=.false.

  info = psb_success_
end subroutine psb_d_cyy_reinit

subroutine  psb_d_cyy_trim(a)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_trim
  implicit none 
  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  integer(psb_ipk_) :: err_act, info, nz, m 
  character(len=20)  :: name='trim'
  logical, parameter :: debug=.false.

  info = psb_success_
end subroutine psb_d_cyy_trim

subroutine psb_d_cyy_print(iout,a,iv,eirs,eics,head,ivr,ivc)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_print
  implicit none 

  integer(psb_ipk_), intent(in)               :: iout
  class(psb_d_cyy_sparse_mat), intent(in) :: a   
  integer(psb_ipk_), intent(in), optional     :: iv(:)
  integer(psb_ipk_), intent(in), optional     :: eirs,eics
  character(len=*), optional        :: head
  integer(psb_ipk_), intent(in), optional     :: ivr(:), ivc(:)

  integer(psb_ipk_) :: err_act
  character(len=20)  :: name='d_cyy_print'
  logical, parameter :: debug=.false.

  integer(psb_ipk_) :: irs,ics,i,j, nmx, ni, nr, nc, nz
end subroutine psb_d_cyy_print

subroutine psb_d_cp_cyy_from_coo(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cp_cyy_from_coo
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer(psb_ipk_), intent(out)                        :: info

  type(psb_d_coo_sparse_mat)   :: tmp
  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

   info = psb_success_
end subroutine psb_d_cp_cyy_from_coo

subroutine psb_d_cp_cyy_to_coo(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cp_cyy_to_coo
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(in)  :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                      :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, nc,i,j,irw, idl,err_act
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
end subroutine psb_d_cp_cyy_to_coo

subroutine psb_d_mv_cyy_to_coo(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_mv_cyy_to_coo
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout)   :: b
  integer(psb_ipk_), intent(out)                        :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, nc,i,j,irw, idl,err_act
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
end subroutine psb_d_mv_cyy_to_coo

subroutine psb_d_mv_cyy_from_coo(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_mv_cyy_from_coo
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                        :: info

  integer(psb_ipk_), allocatable :: itemp(:)
  !locals
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
end subroutine psb_d_mv_cyy_from_coo

subroutine psb_d_mv_cyy_to_fmt(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_mv_cyy_to_fmt
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(inout)  :: b
  integer(psb_ipk_), intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

end subroutine psb_d_mv_cyy_to_fmt

subroutine psb_d_cp_cyy_to_fmt(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cp_cyy_to_fmt
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(in)   :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                       :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
end subroutine psb_d_cp_cyy_to_fmt

subroutine psb_d_mv_cyy_from_fmt(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_mv_cyy_from_fmt
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout)  :: a
  class(psb_d_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)                         :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nza, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_

end subroutine psb_d_mv_cyy_from_fmt

subroutine psb_d_cp_cyy_from_fmt(a,b,info) 
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cp_cyy_from_fmt
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  class(psb_d_base_sparse_mat), intent(in)   :: b
  integer(psb_ipk_), intent(out)                        :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  logical             :: rwshr_
  integer(psb_ipk_) :: nz, nr, i,j,irw, idl,err_act, nc
  integer(psb_ipk_), Parameter  :: maxtry=8
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)   :: name

  info = psb_success_
end subroutine psb_d_cp_cyy_from_fmt

subroutine psb_d_cyy_cp_from(a,b)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_cp_from
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout) :: a
  type(psb_d_cyy_sparse_mat), intent(in)   :: b
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='cp_from'
  logical, parameter :: debug=.false.

  info = psb_success_

end subroutine psb_d_cyy_cp_from

subroutine psb_d_cyy_mv_from(a,b)
  use psb_base_mod
  use psb_d_cyy_mat_mod, psb_protect_name => psb_d_cyy_mv_from
  implicit none 

  class(psb_d_cyy_sparse_mat), intent(inout)  :: a
  type(psb_d_cyy_sparse_mat), intent(inout) :: b
  integer(psb_ipk_) :: err_act, info
  character(len=20)  :: name='mv_from'
  logical, parameter :: debug=.false.

  call psb_erractionsave(err_act)
  info = psb_success_
end subroutine psb_d_cyy_mv_from

