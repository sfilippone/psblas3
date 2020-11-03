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
!    Moved here from AMG-AINV, original copyright below.
!
!
!                       AMG-AINV: Approximate Inverse plugin for
!                             AMG4PSBLAS version 1.0
!
!    (C) Copyright 2020
!
!                        Salvatore Filippone  University of Rome Tor Vergata
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
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
!
!
!                       AMG-AINV: Approximate Inverse plugin for
!                             AMG4PSBLAS version 1.0
!
!    (C) Copyright 2020
!
!                        Salvatore Filippone  University of Rome Tor Vergata
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the AMG4PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AMG4PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
module psb_c_biconjg_mod

  interface psb_sparse_biconjg
    module procedure psb_csparse_biconjg
  end interface


  abstract interface
    subroutine psb_csparse_biconjg_variant(n,a,p,z,w,nzrmax,sp_thresh,info)
      use psb_base_mod, only : psb_c_csr_sparse_mat, psb_c_csc_sparse_mat, &
           &  psb_spk_, psb_ipk_
      !
      implicit none
      integer(psb_ipk_), intent(in)               :: n
      type(psb_c_csr_sparse_mat), intent(in)    :: a
      type(psb_c_csc_sparse_mat), intent(inout) :: z,w
      integer(psb_ipk_), intent(in)               :: nzrmax
      real(psb_spk_), intent(in)                :: sp_thresh
      complex(psb_spk_), intent(out)                :: p(:)
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_csparse_biconjg_variant
  end interface


  procedure(psb_csparse_biconjg_variant) :: psb_csparse_biconjg_llk,&
       & psb_csparse_biconjg_s_llk,  psb_csparse_biconjg_s_ft_llk,&
       &  psb_csparse_biconjg_llk_noth, psb_csparse_biconjg_mlk

#if defined(HAVE_TUMA_SAINV)
  procedure(psb_csparse_biconjg_variant)  ::  psb_csparse_tuma_sainv,&
       & psb_csparse_tuma_lainv
#endif


contains

  subroutine psb_csparse_biconjg(alg,n,acsr,p,z,w,nzrmax,sp_thresh,info)
    use psb_base_mod
    use psb_base_ainv_mod
    integer(psb_ipk_), intent(in)            :: alg,n
    type(psb_c_csr_sparse_mat), intent(in) :: acsr
    type(psb_cspmat_type), intent(out)     :: z, w
    integer(psb_ipk_), intent(in)            :: nzrmax
    real(psb_spk_), intent(in)             :: sp_thresh
    complex(psb_spk_), intent(out)             :: p(:)
    integer(psb_ipk_), intent(out)           :: info

    type(psb_c_csc_sparse_mat)             :: zcsc,wcsc
    integer(psb_ipk_) :: i,j,k,nrm
    integer(psb_ipk_) :: err_act
    character(len=20)  :: name='psb_sparse_biconjg'
    integer(psb_ipk_), parameter :: variant=1


    info = psb_success_
    call psb_erractionsave(err_act)
    if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_; goto 9999
    end if

    if (size(p)<n) then
      write(psb_err_unit,*) 'Size of P wrong'
      info = psb_err_internal_error_
      call psb_errpush(psb_err_internal_error_,name,a_err='Allocate')
      goto 9999
    end if

    select case(alg)
    case (psb_ainv_llk_)
      call psb_csparse_biconjg_llk(n,acsr,p,zcsc,wcsc,nzrmax,sp_thresh,info)
    case (psb_ainv_s_llk_)
      call psb_csparse_biconjg_s_llk(n,acsr,p,zcsc,wcsc,nzrmax,sp_thresh,info)
    case (psb_ainv_mlk_)
      call psb_csparse_biconjg_mlk(n,acsr,p,zcsc,wcsc,nzrmax,sp_thresh,info)
    case (psb_ainv_s_ft_llk_)
      call psb_csparse_biconjg_s_ft_llk(n,acsr,p,zcsc,wcsc,nzrmax,sp_thresh,info)
    case (psb_ainv_llk_noth_)
      call psb_csparse_biconjg_llk_noth(n,acsr,p,zcsc,wcsc,nzrmax,sp_thresh,info)
!#if defined(HAVE_TUMA_SAINV)
!    case (psb_ainv_s_tuma_)
!      call psb_csparse_tuma_sainv(n,acsr,p,zcsc,wcsc,nzrmax,sp_thresh,info)
!    case (psb_ainv_l_tuma_)
!      call psb_csparse_tuma_lainv(n,acsr,p,zcsc,wcsc,nzrmax,sp_thresh,info)
!#endif
    case default
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='Invalid alg')
      goto 9999
    end select

    if (info /= 0) then
      info = psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='sparse_orth')
      goto 9999
    end if

    call z%mv_from(zcsc)
    call z%cscnv(info,type='CSR')
    call w%mv_from(wcsc)
    call w%transp()
    call w%cscnv(info,type='CSR')

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_csparse_biconjg


  subroutine psb_c_spmspv(alpha,a,nx,ix,vx,beta,ny,iy,vy, info)
    !
    !  y = A x  sparse-sparse mode, A in CSC
    !
    use psb_base_mod
    implicit none
    integer(psb_ipk_), intent(in)            :: nx, ix(:)
    complex(psb_spk_), intent(in)              :: alpha, beta, vx(:)
    integer(psb_ipk_), intent(inout)         :: ny, iy(:)
    complex(psb_spk_), intent(inout)           :: vy(:)
    type(psb_c_csc_sparse_mat), intent(in) :: a
    integer(psb_ipk_), intent(out)           :: info

    integer(psb_ipk_) :: i,j,k,m,n, nv, na, iszy
    integer(psb_ipk_), allocatable        :: iv(:)
    complex(psb_spk_), allocatable           :: vv(:)

    info = 0
! !$    write(0,*) 'd_spmspv ',alpha,beta
    if (beta == -cone) then
      do i=1, ny
        vy(i) = -vy(i)
      end do
    else if (beta == czero) then
      do i=1, ny
        vy(i) = czero
      end do
    else if (beta /= cone) then
      do i=1, ny
        vy(i) = vy(i) * beta
      end do
    end if
    if (alpha == czero)  return
    iszy = min(size(iy),size(vy))
    m = a%get_nrows()
    n = a%get_ncols()

    if ((ny > m) .or. (nx > n)) then
      write(0,*) 'Wrong input spmspv rows: ',m,ny,&
           & ' cols: ',n,nx
      info = -4
      return
    end if

    allocate(iv(m), vv(m), stat=info)
    if (info /= 0) then
      write(0,*) 'Allocation error in spmspv'
      info = -999
      return
    endif

    do i = 1, nx
      j  = ix(i)
      ! Access column J of A
      k  = a%icp(j)
      na = a%icp(j+1) - a%icp(j)
      call psb_nspaxpby(nv,iv,vv,&
           & (alpha*vx(i)), na, a%ia(k:k+na-1), a%val(k:k+na-1),&
           & cone, ny, iy, vy, info)

      if (info /= 0) then
        write(0,*) 'Internal error in spmspv from nspaxpby'
        info = -998
        return
      endif
      if (nv > iszy) then
        write(0,*) 'Error in spmspv: out of memory for output'
        info = -997
        return
      endif
      ny = nv
      iy(1:ny) = iv(1:ny)
      vy(1:ny) = vv(1:ny)
    end do
  end subroutine psb_c_spmspv


  subroutine psb_c_spvspm(alpha,a,nx,ix,vx,beta,ny,iy,vy, info)
    !
    !  y = x A  sparse-sparse mode, A in CSR
    !
    use psb_base_mod
    implicit none
    integer(psb_ipk_), intent(in)            :: nx, ix(:)
    complex(psb_spk_), intent(in)              :: alpha, beta, vx(:)
    integer(psb_ipk_), intent(inout)         :: ny, iy(:)
    complex(psb_spk_), intent(inout)           :: vy(:)
    type(psb_c_csr_sparse_mat), intent(in) :: a
    integer(psb_ipk_), intent(out)           :: info

    integer(psb_ipk_)              :: i,j,k,m,n, nv, na, iszy
    integer(psb_ipk_), allocatable :: iv(:)
    complex(psb_spk_), allocatable   :: vv(:)

    info = 0
! !$    write(0,*) 'd_spvspm ',alpha,beta
    if (beta == -cone) then
      do i=1, ny
        vy(i) = -vy(i)
      end do
    else if (beta == czero) then
      do i=1, ny
        vy(i) = czero
      end do
    else if (beta /= cone) then
      do i=1, ny
        vy(i) = vy(i) * beta
      end do
    end if
    if (alpha == czero)  return
    iszy = min(size(iy),size(vy))
    m = a%get_nrows()
    n = a%get_ncols()

    if ((ny > m) .or. (nx > n)) then
      write(0,*) 'Wrong input spmspv rows: ',m,ny,&
           & ' cols: ',n,nx
      info = -4
      return
    end if

    allocate(iv(m), vv(m), stat=info)
    if (info /= 0) then
      write(0,*) 'Allocation error in spmspv'
      info = -999
      return
    endif

    do i = 1, nx
      j  = ix(i)
      ! Access column J of A
      k  = a%irp(j)
      na = a%irp(j+1) - a%irp(j)
      call psb_nspaxpby(nv,iv,vv,&
           & (alpha*vx(i)), na, a%ja(k:k+na-1), a%val(k:k+na-1),&
           & cone, ny, iy, vy, info)

      if (info /= 0) then
        write(0,*) 'Internal error in spmspv from nspaxpby'
        info = -998
        return
      endif
      if (nv > iszy) then
        write(0,*) 'Error in spmspv: out of memory for output'
        info = -997
        return
      endif
      ny = nv
      iy(1:ny) = iv(1:ny)
      vy(1:ny) = vv(1:ny)
    end do
  end subroutine psb_c_spvspm

end module psb_c_biconjg_mod
