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
  

subroutine psb_d_cp_hlg_from_coo(a,b,info) 
  
  use psb_base_mod
#ifdef HAVE_SPGPU
  use hlldev_mod
  use psb_vectordev_mod
  use psb_gpu_env_mod
  use psb_d_hlg_mat_mod, psb_protect_name => psb_d_cp_hlg_from_coo
#else 
  use psb_d_hlg_mat_mod
#endif
  implicit none 

  class(psb_d_hlg_sparse_mat), intent(inout) :: a
  class(psb_d_coo_sparse_mat), intent(in)    :: b
  integer(psb_ipk_), intent(out)             :: info

  !locals
  type(psb_d_coo_sparse_mat) :: tmp
  integer(psb_ipk_)   :: debug_level, debug_unit, hksz
  integer(psb_ipk_), allocatable  :: idisp(:)
  character(len=20)   :: name='hll_from_coo'
  Integer(Psb_ipk_)   :: nza, nr, i,j,irw, idl,err_act, nc, isz,irs
  integer(psb_ipk_)   :: nzm, ir, ic, k, hk, mxrwl, noffs, kc
  integer(psb_ipk_), allocatable :: irn(:), ja(:), hko(:)
  real(psb_dpk_), allocatable :: val(:)
  logical, parameter :: debug=.false.
  
  info = psb_success_
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
#ifdef HAVE_SPGPU
  hksz = max(1,psb_gpu_WarpSize())
#else
  hksz = psi_get_hksz()
#endif

  if (b%is_by_rows()) then

    nr   = b%get_nrows()
    nc   = b%get_ncols()
    nza  = b%get_nzeros()
    if (debug) write(0,*) 'Copying through GPU',nza
    call  psi_compute_hckoff_from_coo(a,noffs,isz,hksz,idisp,b,info)
    if (info /=0) then 
      write(0,*) ' Error from psi_compute_hckoff:',info, noffs,isz
      return
    end if
    if (debug)write(0,*) ' From psi_compute_hckoff:',noffs,isz,a%hkoffs(1:min(10,noffs+1))

    if (c_associated(a%deviceMat)) then 
      call freeHllDevice(a%deviceMat)
    endif
    info = FallochllDevice(a%deviceMat,hksz,nr,nza,isz,spgpu_type_double,1)
    if (info == 0) info = psi_CopyCooToHlg(nr,nc,nza, hksz,noffs,isz,&
         &  a%irn,a%hkoffs,idisp,b%ja, b%val, a%deviceMat)
    call a%set_dev()
  else
    ! This is to guarantee tmp%is_by_rows()
    call b%cp_to_coo(tmp,info)
    call tmp%fix(info)

    nr   = tmp%get_nrows()
    nc   = tmp%get_ncols()
    nza  = tmp%get_nzeros()
    if (debug) write(0,*) 'Copying through GPU'
    call  psi_compute_hckoff_from_coo(a,noffs,isz,hksz,idisp,tmp,info)
    if (info /=0) then 
      write(0,*) ' Error from psi_compute_hckoff:',info, noffs,isz
      return
    end if
    if (debug)write(0,*) ' From psi_compute_hckoff:',noffs,isz,a%hkoffs(1:min(10,noffs+1))

    if (c_associated(a%deviceMat)) then 
      call freeHllDevice(a%deviceMat)
    endif
    info = FallochllDevice(a%deviceMat,hksz,nr,nza,isz,spgpu_type_double,1)
    if (info == 0) info = psi_CopyCooToHlg(nr,nc,nza, hksz,noffs,isz,&
         &  a%irn,a%hkoffs,idisp,tmp%ja, tmp%val, a%deviceMat)

    call tmp%free()
    call a%set_dev()
  end if
  if (info /= 0) goto 9999

  return

9999 continue
  info = psb_err_alloc_dealloc_
  return

contains 
  subroutine psi_compute_hckoff_from_coo(a,noffs,isz,hksz,idisp,b,info)
    use psb_base_mod
    use psi_ext_util_mod
    implicit none 
    class(psb_d_hll_sparse_mat), intent(inout) :: a
    class(psb_d_coo_sparse_mat), intent(in)    :: b
    integer(psb_ipk_), allocatable, intent(out) :: idisp(:)
    integer(psb_ipk_), intent(in)               :: hksz
    integer(psb_ipk_), intent(out)             :: info, noffs, isz

    !locals
    Integer(Psb_ipk_)   :: nza, nr, i,j,irw, idl,err_act, nc, irs
    integer(psb_ipk_)   :: nzm, ir, ic, k, hk, mxrwl, kc
    logical, parameter :: debug=.false.

    info = 0
    nr   = b%get_nrows()
    nc   = b%get_ncols()
    nza  = b%get_nzeros()

    ! If it is sorted then we can lessen memory impact 
    a%psb_d_base_sparse_mat = b%psb_d_base_sparse_mat
    if (debug) write(0,*) 'Start compute hckoff_from_coo',nr,nc,nza
    ! First compute the number of nonzeros in each row.
    call psb_realloc(nr,a%irn,info) 
    if (info == 0) call psb_realloc(nr+1,idisp,info) 
    if (info /= 0) return
    a%irn = 0
    if (debug) then 
      do i=1, nza
        if ((1<=b%ia(i)).and.(b%ia(i)<= nr)) then 
          a%irn(b%ia(i)) = a%irn(b%ia(i)) + 1
        else
          write(0,*) 'Out of bouds IA ',i,b%ia(i),nr
        end if
      end do
    else
      do i=1, nza
        a%irn(b%ia(i)) = a%irn(b%ia(i)) + 1
      end do
    end if
    a%nzt = nza


    ! Second. Figure out the block offsets. 
    call a%set_hksz(hksz)
    noffs = (nr+hksz-1)/hksz
    call psb_realloc(noffs+1,a%hkoffs,info) 
    if (debug) write(0,*) ' noffsets ',noffs,info
    if (info /= 0) return
    a%hkoffs(1) = 0
    j=1
    idisp(1) = 0 
    do i=1,nr,hksz
      ir    = min(hksz,nr-i+1) 
      mxrwl = a%irn(i)
      idisp(i+1) = idisp(i) + a%irn(i) 
      do k=1,ir-1
        idisp(i+k+1) = idisp(i+k) + a%irn(i+k) 
        mxrwl = max(mxrwl,a%irn(i+k))
      end do
      a%hkoffs(j+1) = a%hkoffs(j) + mxrwl*hksz
      j = j + 1 
    end do

    !
    ! At this point a%hkoffs(noffs+1) contains the allocation
    ! size a%ja a%val. 
    ! 
    isz = a%hkoffs(noffs+1)
!!$    write(*,*) 'End of psi_comput_hckoff ',info
  end subroutine psi_compute_hckoff_from_coo

end subroutine psb_d_cp_hlg_from_coo
