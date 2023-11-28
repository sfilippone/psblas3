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
  

subroutine psb_c_cuda_csrg_to_gpu(a,info,nzrm) 

  use psb_base_mod
#ifdef HAVE_SPGPU
  use cusparse_mod
  use psb_c_cuda_csrg_mat_mod, psb_protect_name => psb_c_cuda_csrg_to_gpu
#else 
  use psb_c_cuda_csrg_mat_mod
#endif
  implicit none 
  class(psb_c_cuda_csrg_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_), intent(in), optional     :: nzrm

  integer(psb_ipk_) :: m, nzm, n, pitch,maxrowsize,nz
  integer(psb_ipk_) :: nzdi,i,j,k,nrz
  integer(psb_ipk_), allocatable :: irpdi(:),jadi(:)
  complex(psb_spk_), allocatable :: valdi(:)

  info = 0

#ifdef HAVE_SPGPU
  if ((.not.allocated(a%val)).or.(.not.allocated(a%ja))) return

  m   = a%get_nrows()
  n   = a%get_ncols()
  nz  = a%get_nzeros()
  if (c_associated(a%deviceMat%Mat)) then 
    info = CSRGDeviceFree(a%deviceMat)
  end if
#if CUDA_SHORT_VERSION <= 10 
  if (a%is_unit()) then 
    !
    ! CUSPARSE has the habit of storing the diagonal and then ignoring,
    ! whereas we do not store it. Hence this adapter code. 
    !    
    nzdi = nz + m
    if (info == 0) info = CSRGDeviceAlloc(a%deviceMat,m,n,nzdi)
    if (info == 0) info = CSRGDeviceSetMatIndexBase(a%deviceMat,cusparse_index_base_one)
    if (info == 0) then 
      if (a%is_unit()) then 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_unit)
      else 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
      end if
    end if
    !!! We are explicitly adding the diagonal 
    !! info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
    if ((info == 0) .and. a%is_triangle()) then 
      info = CSRGDeviceSetMatType(a%deviceMat,cusparse_matrix_type_triangular)
      if ((info == 0).and.a%is_upper()) then 
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_upper)
      else
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_lower)
      end if
    end if
    if (info == 0) allocate(irpdi(m+1),jadi(nzdi),valdi(nzdi),stat=info)
    if (info == 0) then 
      irpdi(1) = 1
      if (a%is_triangle().and.a%is_upper()) then 
        do i=1,m
          j        = irpdi(i) 
          jadi(j)  = i
          valdi(j) = cone
          nrz      = a%irp(i+1)-a%irp(i)
          jadi(j+1:j+nrz)  = a%ja(a%irp(i):a%irp(i+1)-1)
          valdi(j+1:j+nrz) = a%val(a%irp(i):a%irp(i+1)-1)
          irpdi(i+1) = j + nrz + 1
          !          write(0,*) 'Row ',i,' : ',irpdi(i:i+1),':',jadi(j:j+nrz),valdi(j:j+nrz)
        end do
      else
        do i=1,m
          j        = irpdi(i) 
          nrz      = a%irp(i+1)-a%irp(i)
          jadi(j+0:j+nrz-1)  = a%ja(a%irp(i):a%irp(i+1)-1)
          valdi(j+0:j+nrz-1) = a%val(a%irp(i):a%irp(i+1)-1)
          jadi(j+nrz)  = i
          valdi(j+nrz) = cone
          irpdi(i+1)   = j + nrz + 1
          !          write(0,*) 'Row ',i,' : ',irpdi(i:i+1),':',jadi(j:j+nrz),valdi(j:j+nrz)
        end do        
      end if
    end if
    if (info == 0) info = CSRGHost2Device(a%deviceMat,m,n,nzdi,irpdi,jadi,valdi)

  else

    if (info == 0) info = CSRGDeviceAlloc(a%deviceMat,m,n,nz)
    if (info == 0) info = CSRGDeviceSetMatIndexBase(a%deviceMat,cusparse_index_base_one)
    if (info == 0) then 
      if (a%is_unit()) then 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_unit)
      else 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
      end if
    end if
    if ((info == 0) .and. a%is_triangle()) then 
      info = CSRGDeviceSetMatType(a%deviceMat,cusparse_matrix_type_triangular)
      if ((info == 0).and.a%is_upper()) then 
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_upper)
      else
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_lower)
      end if
    end if

    if (info == 0) info = CSRGHost2Device(a%deviceMat,m,n,nz,a%irp,a%ja,a%val)
  endif

  if ((info == 0) .and. a%is_triangle()) then 
    info = CSRGDeviceCsrsmAnalysis(a%deviceMat)
  end if

#elif CUDA_VERSION <  11030
  if (a%is_unit()) then 
    !
    ! CUSPARSE has the habit of storing the diagonal and then ignoring,
    ! whereas we do not store it. Hence this adapter code. 
    !    
    nzdi = nz + m
    if (info == 0) info = CSRGDeviceAlloc(a%deviceMat,m,n,nzdi)
!!$    write(0,*) 'Done deviceAlloc'
    if (info == 0) info = CSRGDeviceSetMatIndexBase(a%deviceMat,cusparse_index_base_zero)
!!$    write(0,*) 'Done SetIndexBase'
    if (info == 0) then 
      if (a%is_unit()) then 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_unit)
      else 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
      end if
    end if
    !!! We are explicitly adding the diagonal 
    !! info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
    if ((info == 0) .and. a%is_triangle()) then 
      info = CSRGDeviceSetMatType(a%deviceMat,cusparse_matrix_type_triangular)
      if ((info == 0).and.a%is_upper()) then 
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_upper)
      else
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_lower)
      end if
    end if
    if (info == 0) allocate(irpdi(m+1),jadi(0:nzdi),valdi(0:nzdi),stat=info)
    if (info == 0) then 
      irpdi(1) = 0
      if (a%is_triangle().and.a%is_upper()) then 
        do i=1,m
          j        = irpdi(i) 
          jadi(j)  = i
          valdi(j) = cone
          nrz      = a%irp(i+1)-a%irp(i)
          jadi(j+1:j+nrz)  = a%ja(a%irp(i):a%irp(i+1)-1)-1
          valdi(j+1:j+nrz) = a%val(a%irp(i):a%irp(i+1)-1)
          irpdi(i+1) = j + nrz + 1
          !          write(0,*) 'Row ',i,' : ',irpdi(i:i+1),':',jadi(j:j+nrz),valdi(j:j+nrz)
        end do
      else
        do i=1,m
          j        = irpdi(i) 
          nrz      = a%irp(i+1)-a%irp(i)
          jadi(j+0:j+nrz-1)  = a%ja(a%irp(i):a%irp(i+1)-1)-1
          valdi(j+0:j+nrz-1) = a%val(a%irp(i):a%irp(i+1)-1)
          jadi(j+nrz)  = i
          valdi(j+nrz) = cone
          irpdi(i+1)   = j + nrz + 1
          !          write(0,*) 'Row ',i,' : ',irpdi(i:i+1),':',jadi(j:j+nrz),valdi(j:j+nrz)
        end do        
      end if
    end if
    if (info == 0) info = CSRGHost2Device(a%deviceMat,m,n,nzdi,irpdi,jadi,valdi)

  else

    if (info == 0) info = CSRGDeviceAlloc(a%deviceMat,m,n,nz)
!!$    write(0,*) 'Done deviceAlloc', info
    if (info == 0) info = CSRGDeviceSetMatIndexBase(a%deviceMat,&
         & cusparse_index_base_zero)
!!$    write(0,*) 'Done setIndexBase', info
    if (info == 0) then 
      if (a%is_unit()) then 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_unit)
      else 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
      end if
    end if
    if ((info == 0) .and. a%is_triangle()) then 
      info = CSRGDeviceSetMatType(a%deviceMat,cusparse_matrix_type_triangular)
      if ((info == 0).and.a%is_upper()) then 
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_upper)
      else
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_lower)
      end if
    end if
    nzdi=a%irp(m+1)-1
    if (info == 0) allocate(irpdi(m+1),jadi(max(nzdi,1)),stat=info)
    if (info == 0) then
      irpdi(1:m+1) = a%irp(1:m+1) -1
      jadi(1:nzdi) = a%ja(1:nzdi) -1
    end if
    if (info == 0) info = CSRGHost2Device(a%deviceMat,m,n,nz,irpdi,jadi,a%val)
!!$    write(0,*) 'Done Host2Device', info
  endif


#else

  if (a%is_unit()) then 
    !
    ! CUSPARSE has the habit of storing the diagonal and then ignoring,
    ! whereas we do not store it. Hence this adapter code. 
    !    
    nzdi = nz + m
    if (info == 0) info = CSRGDeviceAlloc(a%deviceMat,m,n,nzdi)
    if (info == 0) then 
      if (a%is_unit()) then 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_unit)
      else 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
      end if
    end if
    !!! We are explicitly adding the diagonal 
    !! info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
    if ((info == 0) .and. a%is_triangle()) then 
!!$      info = CSRGDeviceSetMatType(a%deviceMat,cusparse_matrix_type_triangular)
      if ((info == 0).and.a%is_upper()) then 
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_upper)
      else
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_lower)
      end if
    end if
    if (info == 0) allocate(irpdi(m+1),jadi(nzdi),valdi(nzdi),stat=info)
    if (info == 0) then 
      irpdi(1) = 1
      if (a%is_triangle().and.a%is_upper()) then 
        do i=1,m
          j        = irpdi(i) 
          jadi(j)  = i
          valdi(j) = cone
          nrz      = a%irp(i+1)-a%irp(i)
          jadi(j+1:j+nrz)  = a%ja(a%irp(i):a%irp(i+1)-1)
          valdi(j+1:j+nrz) = a%val(a%irp(i):a%irp(i+1)-1)
          irpdi(i+1) = j + nrz + 1
          !          write(0,*) 'Row ',i,' : ',irpdi(i:i+1),':',jadi(j:j+nrz),valdi(j:j+nrz)
        end do
      else
        do i=1,m
          j        = irpdi(i) 
          nrz      = a%irp(i+1)-a%irp(i)
          jadi(j+0:j+nrz-1)  = a%ja(a%irp(i):a%irp(i+1)-1)
          valdi(j+0:j+nrz-1) = a%val(a%irp(i):a%irp(i+1)-1)
          jadi(j+nrz)  = i
          valdi(j+nrz) = cone
          irpdi(i+1)   = j + nrz + 1
          !          write(0,*) 'Row ',i,' : ',irpdi(i:i+1),':',jadi(j:j+nrz),valdi(j:j+nrz)
        end do        
      end if
    end if
    if (info == 0) info = CSRGHost2Device(a%deviceMat,m,n,nzdi,irpdi,jadi,valdi)

  else

    if (info == 0) info = CSRGDeviceAlloc(a%deviceMat,m,n,nz)
!!$    if (info == 0) info = CSRGDeviceSetMatIndexBase(a%deviceMat,cusparse_index_base_one)
    if (info == 0) then 
      if (a%is_unit()) then 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_unit)
      else 
        info = CSRGDeviceSetMatDiagType(a%deviceMat,cusparse_diag_type_non_unit)
      end if
    end if
    if ((info == 0) .and. a%is_triangle()) then 
!!$      info = CSRGDeviceSetMatType(a%deviceMat,cusparse_matrix_type_triangular)
      if ((info == 0).and.a%is_upper()) then 
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_upper)
      else
        info = CSRGDeviceSetMatFillMode(a%deviceMat,cusparse_fill_mode_lower)
      end if
    end if

    if (info == 0) info = CSRGHost2Device(a%deviceMat,m,n,nz,a%irp,a%ja,a%val)
  endif

!!$  if ((info == 0) .and. a%is_triangle()) then 
!!$    info = CSRGDeviceCsrsmAnalysis(a%deviceMat)
!!$  end if
 
#endif
  call a%set_sync()

  if (info /= 0) then 
    write(0,*) 'Error in CSRG_TO_GPU ',info
  end if
#endif

end subroutine psb_c_cuda_csrg_to_gpu
