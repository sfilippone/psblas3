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
  

subroutine psb_z_mv_elg_from_fmt(a,b,info) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use elldev_mod
  use psb_vectordev_mod
  use psb_z_elg_mat_mod, psb_protect_name => psb_z_mv_elg_from_fmt
#else 
  use psb_z_elg_mat_mod
#endif
  implicit none 

  class(psb_z_elg_sparse_mat), intent(inout)  :: a
  class(psb_z_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)              :: info

  !locals
  type(psb_z_coo_sparse_mat) :: tmp
  Integer(Psb_ipk_)          :: nza, nr, i,j,irw, idl,err_act, nc, ld, nzm, m
#ifdef HAVE_SPGPU
  type(elldev_parms) :: gpu_parms
#endif

  info = psb_success_

  if (b%is_dev()) call b%sync()
  select type (b)
  type is (psb_z_coo_sparse_mat) 
    call a%mv_from_coo(b,info)

  class is (psb_z_ell_sparse_mat) 
    nzm = size(b%ja,2)  
    m   = b%get_nrows()
    nc  = b%get_ncols()
    nza = b%get_nzeros()
#ifdef HAVE_SPGPU
    gpu_parms = FgetEllDeviceParams(m,nzm,nza,nc,spgpu_type_double,1)
    ld  = gpu_parms%pitch
    nzm = gpu_parms%maxRowSize
#else
    ld  = m 
#endif
    a%psb_z_base_sparse_mat = b%psb_z_base_sparse_mat
    call move_alloc(b%irn,   a%irn)
    call move_alloc(b%idiag, a%idiag)
    call psb_realloc(ld,nzm,a%ja,info) 
    if (info == 0) then 
      a%ja(1:m,1:nzm) = b%ja(1:m,1:nzm)
      deallocate(b%ja,stat=info) 
    end if
    if (info == 0) call psb_realloc(ld,nzm,a%val,info) 
    if (info == 0) then 
      a%val(1:m,1:nzm) = b%val(1:m,1:nzm)
      deallocate(b%val,stat=info) 
    end if
    a%nzt = nza
    call b%free()
#ifdef HAVE_SPGPU
    call a%to_gpu(info)
#endif

  class default
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select

end subroutine psb_z_mv_elg_from_fmt
