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
  

subroutine psb_z_hll_csmv(alpha,a,x,beta,y,info,trans) 
  
  use psb_base_mod
  use psb_z_hll_mat_mod, psb_protect_name => psb_z_hll_csmv
  implicit none 
  class(psb_z_hll_sparse_mat), intent(in) :: a
  complex(psb_dpk_), intent(in)             :: alpha, beta, x(:)
  complex(psb_dpk_), intent(inout)          :: y(:)
  integer(psb_ipk_), intent(out)          :: info
  character, optional, intent(in)         :: trans

  character          :: trans_
  integer(psb_ipk_)  :: i,j,k,m,n, nnz, ir, jc, ic, hksz, hkpnt, mxrwl, mmhk
  logical            :: tra, ctra
  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='z_hll_csmv'
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
  tra  = (psb_toupper(trans_) == 'T')
  ctra = (psb_toupper(trans_) == 'C')

  if (tra.or.ctra) then 

    m = a%get_ncols()
    n = a%get_nrows()
    if (size(x,1)<n) then 
      info = 36
      call psb_errpush(info,name,i_err=(/3*ione,n/))
      goto 9999
    end if

    if (size(y,1)<m) then 
      info = 36
      call psb_errpush(info,name,i_err=(/5*ione,m/))
      goto 9999
    end if

    if (beta == zzero) then
      do i = 1, m
        y(i) = zzero
      enddo
    else
      do  i = 1, m
        y(i) = beta*y(i)
      end do
    endif

    hksz = a%get_hksz()
    j=1
    do i=1,n,hksz
      ir    = min(hksz,n-i+1) 
      mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
      hkpnt = a%hkoffs(j) + 1
      call psb_z_hll_csmv_inner(i,ir,mxrwl,a%irn(i),&
           & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
           & a%is_triangle(),a%is_unit(),&
           & x,zone,y,tra,ctra,info) 
      if (info /= psb_success_) goto 9999
      j = j + 1 
    end do


  else if (.not.(tra.or.ctra)) then 

    n    = a%get_ncols()
    m    = a%get_nrows()
    hksz = a%get_hksz()

    if (size(x,1)<n) then 
      info = 36
      call psb_errpush(info,name,i_err=(/3*ione,n/))
      goto 9999
    end if

    if (size(y,1)<m) then 
      info = 36
      call psb_errpush(info,name,i_err=(/5*ione,m/))
      goto 9999
    end if


    if (psi_get_hll_vector()) then 

      hksz = a%get_hksz()
      j    = 1
      mmhk = (m/hksz) * hksz
      if (mmhk > 0) then 
        select case(hksz)
        case(4)
          !$omp parallel do private(i, j,ir,mxrwl, hkpnt)
          do i=1,mmhk,hksz
            j = ((i-1)/hksz)+1
            ir    = hksz
            mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
            if (mxrwl>0) then 
              hkpnt = a%hkoffs(j) + 1
              if (info ==  psb_success_) &
                   & call psb_z_hll_csmv_notra_4(i,mxrwl,a%irn(i),&
                   & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
                   & a%is_triangle(),a%is_unit(),&
                   & x,beta,y,info) 
            end if
            j = j + 1 
          end do
          if (info /= psb_success_) goto 9999

        case(8)
          !$omp parallel do private(i, j,ir,mxrwl, hkpnt)
          do i=1,mmhk,hksz
            j = ((i-1)/hksz)+1
            ir    = hksz
            mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
            if (mxrwl>0) then 
              hkpnt = a%hkoffs(j) + 1
              if (info ==  psb_success_) &
                   &call psb_z_hll_csmv_notra_8(i,mxrwl,a%irn(i),&
                   & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
                   & a%is_triangle(),a%is_unit(),&
                   & x,beta,y,info) 
            end if
            j = j + 1 
          end do
          if (info /= psb_success_) goto 9999

        case(16)
          !$omp parallel do private(i, j,ir,mxrwl, hkpnt)
          do i=1,mmhk,hksz
            j = ((i-1)/hksz)+1
            ir    = hksz
            mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
            if (mxrwl>0) then 
              hkpnt = a%hkoffs(j) + 1
              if (info ==  psb_success_) &
                   & call psb_z_hll_csmv_notra_16(i,mxrwl,a%irn(i),&
                   & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
                   & a%is_triangle(),a%is_unit(),&
                   & x,beta,y,info) 
            end if
            j = j + 1 
          end do
          if (info /= psb_success_) goto 9999

        case(24)
          !$omp parallel do private(i, j,ir,mxrwl, hkpnt)
          do i=1,mmhk,hksz
            j = ((i-1)/hksz)+1
            ir    = hksz
            mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
            if (mxrwl>0) then 
              hkpnt = a%hkoffs(j) + 1
              if (info ==  psb_success_) &
                   & call psb_z_hll_csmv_notra_24(i,mxrwl,a%irn(i),&
                   & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
                   & a%is_triangle(),a%is_unit(),&
                   & x,beta,y,info) 
            end if
            j = j + 1 
          end do
          if (info /= psb_success_) goto 9999

        case(32)
          !$omp parallel do private(i, j,ir,mxrwl, hkpnt)
          do i=1,mmhk,hksz
            j = ((i-1)/hksz)+1
            ir    = hksz
            mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
            if (mxrwl>0) then 
              hkpnt = a%hkoffs(j) + 1
              if (info ==  psb_success_) &
                   & call psb_z_hll_csmv_notra_32(i,mxrwl,a%irn(i),&
                   & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
                   & a%is_triangle(),a%is_unit(),&
                   & x,beta,y,info) 
            end if
            j = j + 1 
          end do
          if (info /= psb_success_) goto 9999

        case default
          !$omp parallel do private(i, j,ir,mxrwl, hkpnt)
          do i=1,mmhk,hksz
            j = ((i-1)/hksz)+1
            ir    = hksz
            mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
            if (mxrwl>0) then 
              hkpnt = a%hkoffs(j) + 1
              if (info ==  psb_success_) &
                   & call psb_z_hll_csmv_inner(i,ir,mxrwl,a%irn(i),&
                   & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
                   & a%is_triangle(),a%is_unit(),&
                   & x,beta,y,tra,ctra,info) 
            end if
            j = j + 1 
          end do
          if (info /= psb_success_) goto 9999
        end select
      end if
      if (mmhk < m) then
        i     = mmhk+1
        ir    = m-mmhk
        mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
        if (mxrwl>0) then 
          hkpnt = a%hkoffs(j) + 1
          call psb_z_hll_csmv_inner(i,ir,mxrwl,a%irn(i),&
               & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
               & a%is_triangle(),a%is_unit(),&
               & x,beta,y,tra,ctra,info) 
          if (info /= psb_success_) goto 9999
        end if
        j = j + 1 
      end if

    else

      j=1
      !$omp parallel do private(i, j,ir,mxrwl, hkpnt)
      do i=1,m,hksz
        j = ((i-1)/hksz)+1
        ir    = min(hksz,m-i+1) 
        mxrwl = (a%hkoffs(j+1) - a%hkoffs(j))/hksz
        hkpnt = a%hkoffs(j) + 1
        if (info ==  psb_success_) &
             & call psb_z_hll_csmv_inner(i,ir,mxrwl,a%irn(i),&
             & alpha,a%ja(hkpnt),hksz,a%val(hkpnt),hksz,&
             & a%is_triangle(),a%is_unit(),&
             & x,beta,y,tra,ctra,info) 
        j = j + 1 
      end do
      if (info /= psb_success_) goto 9999

    end if
  end if
  
  call psb_erractionrestore(err_act)
  return


9999 call psb_error_handler(err_act)
  return

contains

  subroutine psb_z_hll_csmv_inner(ir,m,n,irn,alpha,ja,ldj,val,ldv,&
       & is_triangle,is_unit, x,beta,y,tra,ctra,info) 
    integer(psb_ipk_), intent(in)    :: ir,m,n,ldj,ldv,ja(ldj,*),irn(*)
    complex(psb_dpk_), intent(in)      :: alpha, beta, x(*),val(ldv,*)
    complex(psb_dpk_), intent(inout)   :: y(*)
    logical, intent(in)              :: is_triangle,is_unit,tra,ctra
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_) :: i,j,k, m4, jc
    complex(psb_dpk_)   :: acc(4), tmp

    info = psb_success_
    if (tra) then 

      if (beta == zone) then 
        do i=1,m
          do j=1, irn(i)
            jc = ja(i,j)
            y(jc) = y(jc) + alpha*val(i,j)*x(ir+i-1)
          end do
        end do
      else
        info = -10

      end if

    else if (ctra) then 

      if (beta == zone) then 
        do i=1,m
          do j=1, irn(i)
            jc = ja(i,j)
            y(jc) = y(jc) + alpha*conjg(val(i,j))*x(ir+i-1)
          end do
        end do
      else
        info = -10

      end if

    else if (.not.(tra.or.ctra)) then 

      if (alpha == zzero) then 
        if (beta == zzero) then 
          do i=1,m
            y(ir+i-1) = zzero
          end do
        else
          do i=1,m
            y(ir+i-1) =  beta*y(ir+i-1) 
          end do
        end if

      else
        if (beta == zzero) then 
          do i=1,m
            tmp = zzero
            do j=1, irn(i)
              tmp = tmp + val(i,j)*x(ja(i,j))
            end do
            y(ir+i-1) = alpha*tmp 
          end do
        else
          do i=1,m
            tmp = zzero
            do j=1, irn(i)
              tmp = tmp + val(i,j)*x(ja(i,j))
            end do
            y(ir+i-1) = alpha*tmp + beta*y(ir+i-1)
          end do
        endif
      end if
    end if

    if (is_unit) then 
      do i=1, min(m,n)
        y(i) = y(i) + alpha*x(i)
      end do
    end if

  end subroutine psb_z_hll_csmv_inner

  subroutine psb_z_hll_csmv_notra_8(ir,n,irn,alpha,ja,ldj,val,ldv,&
       & is_triangle,is_unit, x,beta,y,info)
    use psb_base_mod, only : psb_ipk_, psb_dpk_, zzero, psb_success_
    implicit none 
    integer(psb_ipk_), intent(in)    :: ir,n,ldj,ldv,ja(ldj,*),irn(*)
    complex(psb_dpk_), intent(in)      :: alpha, beta, x(*),val(ldv,*)
    complex(psb_dpk_), intent(inout)   :: y(*)
    logical, intent(in)              :: is_triangle,is_unit
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_), parameter :: m=8
    integer(psb_ipk_) :: i,j,k, m4, jc
    complex(psb_dpk_)   :: acc(4), tmp(m)

    info = psb_success_


    tmp(:) = zzero
    if (alpha /= zzero) then 
      do j=1, maxval(irn(1:8))
        tmp(1:8) = tmp(1:8) + val(1:8,j)*x(ja(1:8,j))
      end do
    end if
    if (beta == zzero) then 
      y(ir:ir+8-1) = alpha*tmp(1:8) 
    else
      y(ir:ir+8-1) = alpha*tmp(1:8) + beta*y(ir:ir+8-1)
    end if


    if (is_unit) then 
      do i=1, min(8,n)
        y(ir+i-1) = y(ir+i-1) + alpha*x(ir+i-1)
      end do
    end if

  end subroutine psb_z_hll_csmv_notra_8

  subroutine psb_z_hll_csmv_notra_24(ir,n,irn,alpha,ja,ldj,val,ldv,&
       & is_triangle,is_unit, x,beta,y,info)
    use psb_base_mod, only : psb_ipk_, psb_dpk_, zzero, psb_success_
    implicit none 
    integer(psb_ipk_), intent(in)    :: ir,n,ldj,ldv,ja(ldj,*),irn(*)
    complex(psb_dpk_), intent(in)      :: alpha, beta, x(*),val(ldv,*)
    complex(psb_dpk_), intent(inout)   :: y(*)
    logical, intent(in)              :: is_triangle,is_unit
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_), parameter :: m=24
    integer(psb_ipk_) :: i,j,k, m4, jc
    complex(psb_dpk_)   :: acc(4), tmp(m)

    info = psb_success_


    tmp(:) = zzero
    if (alpha /= zzero) then 
      do j=1, maxval(irn(1:24))
        tmp(1:24) = tmp(1:24) + val(1:24,j)*x(ja(1:24,j))
      end do
    end if
    if (beta == zzero) then 
      y(ir:ir+24-1) = alpha*tmp(1:24) 
    else
      y(ir:ir+24-1) = alpha*tmp(1:24) + beta*y(ir:ir+24-1)
    end if


    if (is_unit) then 
      do i=1, min(24,n)
        y(ir+i-1) = y(ir+i-1) + alpha*x(ir+i-1)
      end do
    end if

  end subroutine psb_z_hll_csmv_notra_24

  subroutine psb_z_hll_csmv_notra_16(ir,n,irn,alpha,ja,ldj,val,ldv,&
       & is_triangle,is_unit, x,beta,y,info)
    use psb_base_mod, only : psb_ipk_, psb_dpk_, zzero, psb_success_
    implicit none 
    integer(psb_ipk_), intent(in)    :: ir,n,ldj,ldv,ja(ldj,*),irn(*)
    complex(psb_dpk_), intent(in)      :: alpha, beta, x(*),val(ldv,*)
    complex(psb_dpk_), intent(inout)   :: y(*)
    logical, intent(in)              :: is_triangle,is_unit
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_), parameter :: m=16
    integer(psb_ipk_) :: i,j,k, m4, jc
    complex(psb_dpk_)   :: acc(4), tmp(m)

    info = psb_success_


    tmp(:) = zzero
    if (alpha /= zzero) then 
      do j=1, maxval(irn(1:16))
        tmp(1:16) = tmp(1:16) + val(1:16,j)*x(ja(1:16,j))
      end do
    end if
    if (beta == zzero) then 
      y(ir:ir+16-1) = alpha*tmp(1:16) 
    else
      y(ir:ir+16-1) = alpha*tmp(1:16) + beta*y(ir:ir+16-1)
    end if


    if (is_unit) then 
      do i=1, min(16,n)
        y(ir+i-1) = y(ir+i-1) + alpha*x(ir+i-1)
      end do
    end if

  end subroutine psb_z_hll_csmv_notra_16

  subroutine psb_z_hll_csmv_notra_32(ir,n,irn,alpha,ja,ldj,val,ldv,&
       & is_triangle,is_unit, x,beta,y,info) 
    use psb_base_mod, only : psb_ipk_, psb_dpk_, zzero, psb_success_
    implicit none 
    integer(psb_ipk_), intent(in)    :: ir,n,ldj,ldv,ja(ldj,*),irn(*)
    complex(psb_dpk_), intent(in)      :: alpha, beta, x(*),val(ldv,*)
    complex(psb_dpk_), intent(inout)   :: y(*)
    logical, intent(in)              :: is_triangle,is_unit
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_), parameter :: m=32
    integer(psb_ipk_) :: i,j,k, m4, jc
    complex(psb_dpk_)   :: acc(4), tmp(m)

    info = psb_success_


    tmp(:) = zzero
    if (alpha /= zzero) then 
      do j=1, maxval(irn(1:32))
        tmp(1:32) = tmp(1:32) + val(1:32,j)*x(ja(1:32,j))
      end do
    end if
    if (beta == zzero) then 
      y(ir:ir+32-1) = alpha*tmp(1:32) 
    else
      y(ir:ir+32-1) = alpha*tmp(1:32) + beta*y(ir:ir+32-1)
    end if


    if (is_unit) then 
      do i=1, min(32,n)
        y(ir+i-1) = y(ir+i-1) + alpha*x(ir+i-1)
      end do
    end if

  end subroutine psb_z_hll_csmv_notra_32

  subroutine psb_z_hll_csmv_notra_4(ir,n,irn,alpha,ja,ldj,val,ldv,&
       & is_triangle,is_unit, x,beta,y,info)
    use psb_base_mod, only : psb_ipk_, psb_dpk_, zzero, psb_success_
    implicit none 
    integer(psb_ipk_), intent(in)    :: ir,n,ldj,ldv,ja(ldj,*),irn(*)
    complex(psb_dpk_), intent(in)      :: alpha, beta, x(*),val(ldv,*)
    complex(psb_dpk_), intent(inout)   :: y(*)
    logical, intent(in)              :: is_triangle,is_unit
    integer(psb_ipk_), intent(out)   :: info

    integer(psb_ipk_), parameter :: m=4
    integer(psb_ipk_) :: i,j,k, m4, jc
    complex(psb_dpk_)   :: acc(4), tmp(m)

    info = psb_success_


    tmp(:) = zzero
    if (alpha /= zzero) then 
      do j=1, maxval(irn(1:4))
        tmp(1:4) = tmp(1:4) + val(1:4,j)*x(ja(1:4,j))
      end do
    end if
    if (beta == zzero) then 
      y(ir:ir+4-1) = alpha*tmp(1:4) 
    else
      y(ir:ir+4-1) = alpha*tmp(1:4) + beta*y(ir:ir+4-1)
    end if


    if (is_unit) then 
      do i=1, min(4,n)
        y(ir+i-1) = y(ir+i-1) + alpha*x(ir+i-1)
      end do
    end if

  end subroutine psb_z_hll_csmv_notra_4

end subroutine psb_z_hll_csmv
