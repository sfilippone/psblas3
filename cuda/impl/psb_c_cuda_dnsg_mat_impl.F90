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

subroutine psb_c_cuda_dnsg_vect_mv(alpha,a,x,beta,y,info,trans) 
  use psb_base_mod
  use psb_c_cuda_vect_mod
#ifdef HAVE_SPGPU
  use dnsdev_mod
  use psb_c_vectordev_mod
  use psb_c_cuda_dnsg_mat_mod, psb_protect_name => psb_c_cuda_dnsg_vect_mv
#else
  use psb_c_cuda_dnsg_mat_mod
#endif
  implicit none 
  class(psb_c_cuda_dnsg_sparse_mat), intent(in)    :: a
  complex(psb_spk_), intent(in)                 :: alpha, beta
  class(psb_c_base_vect_type), intent(inout) :: x
  class(psb_c_base_vect_type), intent(inout) :: y
  integer(psb_ipk_), intent(out)             :: info
  character, optional, intent(in)            :: trans
  logical           :: tra
  character         :: trans_
  complex(psb_spk_), allocatable      :: rx(:), ry(:)
  Integer(Psb_ipk_) :: err_act, m, n, k 
  character(len=20) :: name='c_cuda_dnsg_vect_mv'

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

  if (trans_ =='N') then
    m = a%get_nrows()
    n = 1
    k = a%get_ncols()
  else
    m = a%get_ncols()
    n = 1
    k = a%get_nrows()
  end if
  select type (xx => x) 
  type is (psb_c_vect_cuda)
    select type(yy => y) 
    type is (psb_c_vect_cuda)
      if (a%is_host()) call a%sync()
      if (xx%is_host()) call xx%sync()
      if (beta /= czero) then 
        if (yy%is_host()) call yy%sync()
      end if
      info = spmvDnsDevice(trans_,m,n,k,alpha,a%deviceMat,&
           & xx%deviceVect,beta,yy%deviceVect)
      if (info /= 0) then 
        call psb_errpush(psb_err_from_subroutine_ai_,name,&
             & a_err='spmvDnsDevice',i_err=(/info,izero,izero,izero,izero/))
        info = psb_err_from_subroutine_ai_
        goto 9999
      end if
      call yy%set_dev()
    class default
      if (a%is_dev()) call a%sync()
      rx = xx%get_vect()
      ry = y%get_vect()
      call a%spmm(alpha,rx,beta,ry,info)
      call y%bld(ry)
    end select
  class default
    if (a%is_dev()) call a%sync()
    rx = x%get_vect()
    ry = y%get_vect()
    call a%spmm(alpha,rx,beta,ry,info)
    call y%bld(ry)
  end select
  
  
  if (info /= 0) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
 
end subroutine psb_c_cuda_dnsg_vect_mv


subroutine psb_c_cuda_dnsg_mold(a,b,info) 
  use psb_base_mod
  use psb_c_cuda_vect_mod
#ifdef HAVE_SPGPU
  use dnsdev_mod
  use psb_c_vectordev_mod
  use psb_c_cuda_dnsg_mat_mod, psb_protect_name => psb_c_cuda_dnsg_mold
#else
  use psb_c_cuda_dnsg_mat_mod
#endif
  implicit none 
  class(psb_c_cuda_dnsg_sparse_mat), intent(in)                  :: a
  class(psb_c_base_sparse_mat), intent(inout), allocatable :: b
  integer(psb_ipk_), intent(out)                           :: info
  Integer(Psb_ipk_)  :: err_act
  character(len=20)  :: name='dnsg_mold'
  logical, parameter :: debug=.false.

  call psb_get_erraction(err_act)
  
  info = 0 
  if (allocated(b)) then 
    call b%free()
    deallocate(b,stat=info)
  end if
  if (info == 0) allocate(psb_c_cuda_dnsg_sparse_mat :: b, stat=info)

  if (info /= psb_success_) then 
    info = psb_err_alloc_dealloc_ 
    call psb_errpush(info, name)
    goto 9999
  end if
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_cuda_dnsg_mold


!!$
!!$  interface 
!!$    subroutine psb_c_cuda_dnsg_inner_vect_sv(alpha,a,x,beta,y,info,trans) 
!!$      import :: psb_ipk_, psb_c_cuda_dnsg_sparse_mat, psb_spk_,  psb_c_base_vect_type
!!$      class(psb_c_cuda_dnsg_sparse_mat), intent(in)    :: a
!!$      complex(psb_spk_), intent(in)                 :: alpha, beta
!!$      class(psb_c_base_vect_type), intent(inout) :: x, y
!!$      integer(psb_ipk_), intent(out)             :: info
!!$      character, optional, intent(in)            :: trans
!!$    end subroutine psb_c_cuda_dnsg_inner_vect_sv
!!$  end interface

!!$  interface
!!$    subroutine  psb_c_cuda_dnsg_reallocate_nz(nz,a) 
!!$      import :: psb_c_cuda_dnsg_sparse_mat, psb_ipk_
!!$      integer(psb_ipk_), intent(in)              :: nz
!!$      class(psb_c_cuda_dnsg_sparse_mat), intent(inout) :: a
!!$    end subroutine psb_c_cuda_dnsg_reallocate_nz
!!$  end interface
!!$
!!$  interface
!!$    subroutine  psb_c_cuda_dnsg_allocate_mnnz(m,n,a,nz) 
!!$      import :: psb_c_cuda_dnsg_sparse_mat, psb_ipk_
!!$      integer(psb_ipk_), intent(in)              :: m,n
!!$      class(psb_c_cuda_dnsg_sparse_mat), intent(inout) :: a
!!$      integer(psb_ipk_), intent(in), optional    :: nz
!!$    end subroutine psb_c_cuda_dnsg_allocate_mnnz
!!$  end interface


subroutine psb_c_cuda_dnsg_to_gpu(a,info) 
  use psb_base_mod
  use psb_c_cuda_vect_mod
#ifdef HAVE_SPGPU
  use dnsdev_mod
  use psb_c_vectordev_mod
  use psb_c_cuda_dnsg_mat_mod, psb_protect_name => psb_c_cuda_dnsg_to_gpu
#else
  use psb_c_cuda_dnsg_mat_mod
#endif
  class(psb_c_cuda_dnsg_sparse_mat), intent(inout) :: a
  integer(psb_ipk_), intent(out)             :: info
  Integer(Psb_ipk_) :: err_act, pitch, lda
  logical, parameter :: debug=.false.
  character(len=20) :: name='c_cuda_dnsg_to_gpu'
  
  call psb_erractionsave(err_act)
  info = psb_success_
#ifdef HAVE_SPGPU
  if (debug) write(0,*) 'DNS_TO_GPU',size(a%val,1),size(a%val,2)
  info = FallocDnsDevice(a%deviceMat,a%get_nrows(),a%get_ncols(),&
       & spgpu_type_complex_float,1)
  if (info == 0) info = writeDnsDevice(a%deviceMat,a%val,size(a%val,1),size(a%val,2))
  if (debug) write(0,*) 'DNS_TO_GPU: From writeDnsDEvice',info
    
  
#endif 
  if (info /= 0) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_cuda_dnsg_to_gpu



subroutine psb_c_cuda_cp_dnsg_from_coo(a,b,info)
  use psb_base_mod
  use psb_c_cuda_vect_mod
#ifdef HAVE_SPGPU
  use dnsdev_mod
  use psb_c_vectordev_mod
  use psb_c_cuda_dnsg_mat_mod, psb_protect_name => psb_c_cuda_cp_dnsg_from_coo
#else
  use psb_c_cuda_dnsg_mat_mod
#endif
  implicit none 

  class(psb_c_cuda_dnsg_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(in)    :: b
  integer(psb_ipk_), intent(out)             :: info
  Integer(Psb_ipk_) :: err_act
  character(len=20) :: name='c_cuda_dnsg_cp_from_coo'
  integer(psb_ipk_)   :: debug_level, debug_unit
  logical, parameter :: debug=.false.
  type(psb_c_coo_sparse_mat)  :: tmp

  call psb_erractionsave(err_act)
  info = psb_success_
  if (b%is_dev()) call b%sync()

  call a%psb_c_dns_sparse_mat%cp_from_coo(b,info)
  if (debug) write(0,*) 'dnsg_cp_from_coo: dns_cp',info  
  if (info == 0) call a%to_gpu(info)
  if (info /= 0) goto 9999
  
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_cuda_cp_dnsg_from_coo
  
subroutine psb_c_cuda_cp_dnsg_from_fmt(a,b,info)
  use psb_base_mod
  use psb_c_cuda_vect_mod
#ifdef HAVE_SPGPU
  use dnsdev_mod
  use psb_c_vectordev_mod
  use psb_c_cuda_dnsg_mat_mod, psb_protect_name => psb_c_cuda_cp_dnsg_from_fmt
#else
  use psb_c_cuda_dnsg_mat_mod
#endif
  implicit none 

  class(psb_c_cuda_dnsg_sparse_mat), intent(inout) :: a
  class(psb_c_base_sparse_mat), intent(in)   :: b
  integer(psb_ipk_), intent(out)             :: info

  type(psb_c_coo_sparse_mat)  :: tmp
  Integer(Psb_ipk_) :: err_act
  character(len=20) :: name='c_cuda_dnsg_cp_from_fmt'

  call psb_erractionsave(err_act)
  info = psb_success_
  if (b%is_dev()) call b%sync()
 
  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%cp_from_coo(b,info)

!!$  class is (psb_c_ell_sparse_mat) 
!!$    nzm = psb_size(b%ja,2)  
!!$    m   = b%get_nrows()
!!$    nc  = b%get_ncols()
!!$    nza = b%get_nzeros()
!!$#ifdef HAVE_SPGPU
!!$    gpu_parms = FgetEllDeviceParams(m,nzm,nza,nc,spgpu_type_double,1)
!!$    ld  = gpu_parms%pitch
!!$    nzm = gpu_parms%maxRowSize
!!$#else
!!$    ld  = m 
!!$#endif
!!$    a%psb_c_base_sparse_mat = b%psb_c_base_sparse_mat
!!$    if (info == 0) call psb_safe_cpy( b%idiag, a%idiag , info)
!!$    if (info == 0) call psb_safe_cpy( b%irn,   a%irn , info)
!!$    if (info == 0) call psb_safe_cpy( b%ja ,   a%ja  , info)
!!$    if (info == 0) call psb_safe_cpy( b%val,   a%val , info)
!!$    if (info == 0) call psb_realloc(ld,nzm,a%ja,info) 
!!$    if (info == 0) then 
!!$      a%ja(1:m,1:nzm) = b%ja(1:m,1:nzm)
!!$    end if
!!$    if (info == 0) call psb_realloc(ld,nzm,a%val,info) 
!!$    if (info == 0) then 
!!$      a%val(1:m,1:nzm) = b%val(1:m,1:nzm)
!!$    end if
!!$    a%nzt = nza
!!$#ifdef HAVE_SPGPU
!!$    call a%to_gpu(info)
!!$#endif

  class default
  
    call b%cp_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
  
  if (info /= 0) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_c_cuda_cp_dnsg_from_fmt

  

subroutine psb_c_cuda_mv_dnsg_from_coo(a,b,info)
  use psb_base_mod
  use psb_c_cuda_vect_mod
#ifdef HAVE_SPGPU
  use dnsdev_mod
  use psb_c_vectordev_mod
  use psb_c_cuda_dnsg_mat_mod, psb_protect_name => psb_c_cuda_mv_dnsg_from_coo
#else
  use psb_c_cuda_dnsg_mat_mod
#endif
  implicit none 
  
  class(psb_c_cuda_dnsg_sparse_mat), intent(inout) :: a
  class(psb_c_coo_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)             :: info

  Integer(Psb_ipk_) :: err_act
  logical, parameter :: debug=.false.
  character(len=20) :: name='c_cuda_dnsg_mv_from_coo'

  call psb_erractionsave(err_act)
  info = psb_success_

    if (.not.b%is_by_rows()) call b%fix(info)
  if (info /= psb_success_) return
  if (b%is_dev()) call b%sync()
  call a%cp_from_coo(b,info)
  if (debug) write(0,*) 'dnsg_mv_from_coo: cp_from_coo:',info
  call b%free()
  if (info /= 0) goto 9999 
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
  
end subroutine psb_c_cuda_mv_dnsg_from_coo

  
subroutine psb_c_cuda_mv_dnsg_from_fmt(a,b,info)
  use psb_base_mod
  use psb_c_cuda_vect_mod
#ifdef HAVE_SPGPU
  use dnsdev_mod
  use psb_c_vectordev_mod
  use psb_c_cuda_dnsg_mat_mod, psb_protect_name => psb_c_cuda_mv_dnsg_from_fmt
#else
  use psb_c_cuda_dnsg_mat_mod
#endif
  implicit none 
  class(psb_c_cuda_dnsg_sparse_mat), intent(inout)  :: a
  class(psb_c_base_sparse_mat), intent(inout) :: b
  integer(psb_ipk_), intent(out)              :: info


  type(psb_c_coo_sparse_mat)  :: tmp
  Integer(Psb_ipk_) :: err_act
  character(len=20) :: name='c_cuda_dnsg_cp_from_fmt'

  call psb_erractionsave(err_act)
  info = psb_success_
  if (b%is_dev()) call b%sync()
 
  select type (b)
  type is (psb_c_coo_sparse_mat) 
    call a%mv_from_coo(b,info)

!!$  class is (psb_c_ell_sparse_mat) 
!!$    nzm = psb_size(b%ja,2)  
!!$    m   = b%get_nrows()
!!$    nc  = b%get_ncols()
!!$    nza = b%get_nzeros()
!!$#ifdef HAVE_SPGPU
!!$    gpu_parms = FgetEllDeviceParams(m,nzm,nza,nc,spgpu_type_double,1)
!!$    ld  = gpu_parms%pitch
!!$    nzm = gpu_parms%maxRowSize
!!$#else
!!$    ld  = m 
!!$#endif
!!$    a%psb_c_base_sparse_mat = b%psb_c_base_sparse_mat
!!$    if (info == 0) call psb_safe_cpy( b%idiag, a%idiag , info)
!!$    if (info == 0) call psb_safe_cpy( b%irn,   a%irn , info)
!!$    if (info == 0) call psb_safe_cpy( b%ja ,   a%ja  , info)
!!$    if (info == 0) call psb_safe_cpy( b%val,   a%val , info)
!!$    if (info == 0) call psb_realloc(ld,nzm,a%ja,info) 
!!$    if (info == 0) then 
!!$      a%ja(1:m,1:nzm) = b%ja(1:m,1:nzm)
!!$    end if
!!$    if (info == 0) call psb_realloc(ld,nzm,a%val,info) 
!!$    if (info == 0) then 
!!$      a%val(1:m,1:nzm) = b%val(1:m,1:nzm)
!!$    end if
!!$    a%nzt = nza
!!$#ifdef HAVE_SPGPU
!!$    call a%to_gpu(info)
!!$#endif

  class default
  
    call b%mv_to_coo(tmp,info)
    if (info == psb_success_) call a%mv_from_coo(tmp,info)
  end select
  
  if (info /= 0) goto 9999
  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return
  
  
end subroutine psb_c_cuda_mv_dnsg_from_fmt
