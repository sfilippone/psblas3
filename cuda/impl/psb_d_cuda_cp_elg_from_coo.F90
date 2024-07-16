!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
subroutine psb_d_cuda_cp_elg_from_coo(a,b,info) 
  
  use psb_base_mod
  use elldev_mod
  use psb_vectordev_mod
  use psb_d_cuda_elg_mat_mod, psb_protect_name => psb_d_cuda_cp_elg_from_coo
  use psi_ext_util_mod
  use psb_cuda_env_mod
  implicit none 

  class(psb_d_cuda_elg_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer(psb_ipk_), intent(out)             :: info

  !locals
  Integer(Psb_ipk_)   :: nza, nr, i,j,k, idl,err_act, nc, nzm, &
       & ir, ic, ld, ldv, hacksize
  integer(psb_ipk_)   :: debug_level, debug_unit
  character(len=20)   :: name
  type(psb_d_coo_sparse_mat)  :: tmp
  integer(psb_ipk_), allocatable  :: idisp(:)

  info = psb_success_
  hacksize = max(1,psb_cuda_WarpSize())
  if (b%is_dev()) call b%sync()

  if (b%is_by_rows()) then

    call psi_d_count_ell_from_coo(a,b,idisp,ldv,nzm,info,hacksize=hacksize)

    
    if (c_associated(a%deviceMat)) then 
      call freeEllDevice(a%deviceMat)
    endif
    
    nr   = b%get_nrows()
    nc   = b%get_ncols()
    nza  = b%get_nzeros()
    info = FallocEllDevice(a%deviceMat,nr,nzm,nza,nc,spgpu_type_double,1)

    if (info == 0) info = psi_CopyCooToElg(nr,nc,nza, hacksize,ldv,nzm, &
         & a%irn,idisp,b%ja,b%val, a%deviceMat)
    call a%set_dev()
  else
    call b%cp_to_coo(tmp,info)
    call psi_d_count_ell_from_coo(a,tmp,idisp,ldv,nzm,info,hacksize=hacksize)

    
    if (c_associated(a%deviceMat)) then 
      call freeEllDevice(a%deviceMat)
    endif
    
    nr   = b%get_nrows()
    nc   = b%get_ncols()
    nza  = b%get_nzeros()
    info = FallocEllDevice(a%deviceMat,nr,nzm,nza,nc,spgpu_type_double,1)

    if (info == 0) info = psi_CopyCooToElg(nr,nc,nza, hacksize,ldv,nzm, &
         & a%irn,idisp,tmp%ja,tmp%val, a%deviceMat)

    call a%set_dev()
  end if

  if (info /= psb_success_) goto 9999

  return

9999 continue
  info = psb_err_alloc_dealloc_
  return

contains
  
  subroutine psi_d_count_ell_from_coo(a,b,idisp,ldv,nzm,info,hacksize) 

    use psb_base_mod
    use psi_ext_util_mod
    implicit none 

    class(psb_d_ell_sparse_mat), intent(inout) :: a
    class(psb_d_coo_sparse_mat), intent(in)    :: b
    integer(psb_ipk_), allocatable, intent(out) :: idisp(:)
    integer(psb_ipk_), intent(out)             :: info, nzm, ldv
    integer(psb_ipk_), intent(in), optional    :: hacksize

    !locals
    Integer(Psb_ipk_) :: nza, nr, i,j,k, idl,err_act, nc, &
         & ir, ic, hsz_
    real(psb_dpk_) :: t0,t1
    logical, parameter :: timing=.true.


    info = psb_success_

    nr  = b%get_nrows()
    nc  = b%get_ncols()
    nza = b%get_nzeros()

    hsz_ = 1
    if (present(hacksize)) then
      if (hacksize> 1) hsz_ = hacksize
    end if
    ! Make ldv a multiple of hacksize 
    ldv = ((nr+hsz_-1)/hsz_)*hsz_

    ! If it is sorted then we can lessen memory impact 
    a%psb_d_base_sparse_mat = b%psb_d_base_sparse_mat

    ! First compute the number of nonzeros in each row.
    call psb_realloc(nr,a%irn,info)
    if (info == psb_success_) call psb_realloc(nr+1,idisp,info) 
    if (info /= psb_success_) return
    if (timing) t0=psb_wtime()

    a%irn = 0
    do i=1, nza
      ir = b%ia(i)
      a%irn(ir) = a%irn(ir) + 1
    end do
    nzm = 0
    a%nzt = 0
    idisp(1) = 0 
    do i=1,nr
      nzm = max(nzm,a%irn(i))
      a%nzt = a%nzt + a%irn(i)
      idisp(i+1) = a%nzt
    end do

  end subroutine psi_d_count_ell_from_coo
  
end subroutine psb_d_cuda_cp_elg_from_coo
