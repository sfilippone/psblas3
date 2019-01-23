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
!    
! File:  psb_srwextd.f90 
! Subroutine: 
! Arguments:
!
! We have a problem here: 1. How to handle well all the formats? 
!                         2. What should we do with rowscale? Does it only 
!                            apply when a%fmt()='COO' ?????? 
!
!
subroutine psb_srwextd(nr,a,info,b,rowscale)
  use psb_mat_mod
  use psb_serial_mod, psb_protect_name => psb_srwextd
  implicit none

  ! Extend matrix A up to NR rows with empty ones (i.e.: all zeroes)
  integer(psb_ipk_), intent(in)                 :: nr
  type(psb_sspmat_type), intent(inout)        :: a
  integer(psb_ipk_),intent(out)                 :: info
  type(psb_sspmat_type), intent(in), optional :: b
  logical,intent(in), optional                 :: rowscale

  integer(psb_ipk_) :: i,j,ja,jb,err_act,nza,nzb
  character(len=20)                 :: name, ch_err
  type(psb_s_coo_sparse_mat)        :: actmp
  logical  rowscale_ 

  name='psb_srwextd'
  info  = psb_success_
  call psb_erractionsave(err_act)

  if (nr > a%get_nrows()) then 
    select type(aa=> a%a) 
    type is (psb_s_csr_sparse_mat)
      if (present(b)) then 
        call psb_rwextd(nr,aa,info,b%a,rowscale)
      else
        call psb_rwextd(nr,aa,info,rowscale=rowscale)
      end if
    type is (psb_s_coo_sparse_mat) 
      if (present(b)) then 
        call psb_rwextd(nr,aa,info,b%a,rowscale=rowscale)
      else
        call psb_rwextd(nr,aa,info,rowscale=rowscale)
      end if
    class default
      call aa%mv_to_coo(actmp,info)
      if (info == psb_success_) then 
        if (present(b)) then 
          call psb_rwextd(nr,actmp,info,b%a,rowscale=rowscale)
        else
          call psb_rwextd(nr,actmp,info,rowscale=rowscale)
        end if
      end if
      if (info == psb_success_) call aa%mv_from_coo(actmp,info)
    end select
  end if
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_srwextd
subroutine psb_sbase_rwextd(nr,a,info,b,rowscale)
  use psb_mat_mod
  use psb_serial_mod, psb_protect_name => psb_sbase_rwextd
  implicit none

  ! Extend matrix A up to NR rows with empty ones (i.e.: all zeroes)
  integer(psb_ipk_), intent(in)                        :: nr
  class(psb_s_base_sparse_mat), intent(inout)        :: a
  integer(psb_ipk_),intent(out)                        :: info
  class(psb_s_base_sparse_mat), intent(in), optional :: b
  logical,intent(in), optional                       :: rowscale

  integer(psb_ipk_) :: i,j,ja,jb,err_act,nza,nzb, ma, mb, na, nb
  character(len=20)                 :: name, ch_err
  logical  rowscale_ 

  name='psb_sbase_rwextd'
  info  = psb_success_
  call psb_erractionsave(err_act)

  if (present(rowscale)) then 
    rowscale_ = rowscale
  else
    rowscale_ = .true.
  end if

  ma = a%get_nrows()
  na = a%get_ncols()


  select type(a) 
  type is (psb_s_csr_sparse_mat)

    call psb_ensure_size(nr+1,a%irp,info)

    if (present(b)) then 
      mb = b%get_nrows()
      nb = b%get_ncols()
      nzb = b%get_nzeros()

      select type (b) 
      type is (psb_s_csr_sparse_mat)
        call psb_ensure_size(size(a%ja)+nzb,a%ja,info)
        call psb_ensure_size(size(a%val)+nzb,a%val,info)
        do i=1, min(nr-ma,mb)
          a%irp(ma+i+1) =  a%irp(ma+i) + b%irp(i+1) - b%irp(i)
          ja = a%irp(ma+i)
          do jb = b%irp(i), b%irp(i+1)-1
            a%val(ja) = b%val(jb)
            a%ja(ja)  = b%ja(jb)
            ja = ja + 1
          end do
        end do
        do j=i,nr-ma
          a%irp(ma+j+1) = a%irp(ma+j)
        end do
      class default 

        write(psb_err_unit,*) 'Implement SPGETBLK in RWEXTD!!!!!!!'
      end select
      call a%set_ncols(max(na,nb))

    else

      do j=ma+2,nr+1
        a%irp(j) = a%irp(j-1)
      end do

    end if

    call a%set_nrows(nr)


  type is (psb_s_coo_sparse_mat) 
    nza = a%get_nzeros()

    if (present(b)) then 
      mb = b%get_nrows()
      nb = b%get_ncols()
      nzb = b%get_nzeros()
      call a%reallocate(nza+nzb)

      select type(b)
      type is (psb_s_coo_sparse_mat) 

        if (rowscale_) then 
          do j=1,nzb
            if ((ma + b%ia(j)) <= nr) then 
              nza = nza + 1
              a%ia(nza)  = ma + b%ia(j)
              a%ja(nza)  = b%ja(j)
              a%val(nza) = b%val(j)
            end if
          enddo
        else
          do j=1,nzb
            if ((ma + b%ia(j)) <= nr) then 
              nza = nza + 1
              a%ia(nza)  = b%ia(j)
              a%ja(nza)  = b%ja(j)
              a%val(nza) = b%val(j)
            end if
          enddo
        endif
        call a%set_nzeros(nza)

      type is (psb_s_csr_sparse_mat) 
        
        do i=1, min(nr-ma,mb)
          do jb = b%irp(i), b%irp(i+1)-1
            nza = nza + 1 
            a%val(nza) = b%val(jb)
            a%ia(nza)  = ma + i
            a%ja(nza)  = b%ja(jb)
          end do
        end do
        call a%set_nzeros(nza)

      class default 
        write(psb_err_unit,*) 'Implement SPGETBLK in RWEXTD!!!!!!!'

      end select

      call a%set_ncols(max(na,nb))
    endif

    call a%set_nrows(nr)

  class default 
    info = psb_err_unsupported_format_
    ch_err=a%get_fmt()
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999

  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_sbase_rwextd


subroutine psb_lsrwextd(nr,a,info,b,rowscale)
  use psb_mat_mod
  use psb_serial_mod, psb_protect_name => psb_lsrwextd
  implicit none

  ! Extend matrix A up to NR rows with empty ones (i.e.: all zeroes)
  integer(psb_lpk_), intent(in)                  :: nr
  type(psb_lsspmat_type), intent(inout)        :: a
  integer(psb_ipk_),intent(out)                  :: info
  type(psb_lsspmat_type), intent(in), optional :: b
  logical,intent(in), optional                 :: rowscale

  integer(psb_lpk_) :: i,j,ja,jb,nza,nzb
  integer(psb_ipk_) :: err_act
  character(len=20)                 :: name, ch_err
  type(psb_ls_coo_sparse_mat)     :: actmp
  logical  rowscale_ 

  name='psb_lsrwextd'
  info  = psb_success_
  call psb_erractionsave(err_act)

  if (nr > a%get_nrows()) then 
    select type(aa=> a%a) 
    type is (psb_ls_csr_sparse_mat)
      if (present(b)) then 
        call psb_rwextd(nr,aa,info,b%a,rowscale)
      else
        call psb_rwextd(nr,aa,info,rowscale=rowscale)
      end if
    type is (psb_ls_coo_sparse_mat) 
      if (present(b)) then 
        call psb_rwextd(nr,aa,info,b%a,rowscale=rowscale)
      else
        call psb_rwextd(nr,aa,info,rowscale=rowscale)
      end if
    class default
      call aa%mv_to_coo(actmp,info)
      if (info == psb_success_) then 
        if (present(b)) then 
          call psb_rwextd(nr,actmp,info,b%a,rowscale=rowscale)
        else
          call psb_rwextd(nr,actmp,info,rowscale=rowscale)
        end if
      end if
      if (info == psb_success_) call aa%mv_from_coo(actmp,info)
    end select
  end if
  if (info /= psb_success_) goto 9999

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_lsrwextd
subroutine psb_lsbase_rwextd(nr,a,info,b,rowscale)
  use psb_mat_mod
  use psb_serial_mod, psb_protect_name => psb_lsbase_rwextd
  implicit none

  ! Extend matrix A up to NR rows with empty ones (i.e.: all zeroes)
  integer(psb_lpk_), intent(in)                         :: nr
  class(psb_ls_base_sparse_mat), intent(inout)        :: a
  integer(psb_ipk_),intent(out)                         :: info
  class(psb_ls_base_sparse_mat), intent(in), optional :: b
  logical,intent(in), optional                       :: rowscale

  integer(psb_lpk_) :: i,j,ja,jb,nza,nzb, ma, mb, na, nb
  integer(psb_ipk_) :: err_act
  character(len=20) :: name, ch_err
  logical  rowscale_ 

  name='psb_lsbase_rwextd'
  info  = psb_success_
  call psb_erractionsave(err_act)

  if (present(rowscale)) then 
    rowscale_ = rowscale
  else
    rowscale_ = .true.
  end if

  ma = a%get_nrows()
  na = a%get_ncols()


  select type(a) 
  type is (psb_ls_csr_sparse_mat)

    call psb_ensure_size(nr+1,a%irp,info)

    if (present(b)) then 
      mb = b%get_nrows()
      nb = b%get_ncols()
      nzb = b%get_nzeros()

      select type (b) 
      type is (psb_ls_csr_sparse_mat)
        call psb_ensure_size(size(a%ja)+nzb,a%ja,info)
        call psb_ensure_size(size(a%val)+nzb,a%val,info)
        do i=1, min(nr-ma,mb)
          a%irp(ma+i+1) =  a%irp(ma+i) + b%irp(i+1) - b%irp(i)
          ja = a%irp(ma+i)
          do jb = b%irp(i), b%irp(i+1)-1
            a%val(ja) = b%val(jb)
            a%ja(ja)  = b%ja(jb)
            ja = ja + 1
          end do
        end do
        do j=i,nr-ma
          a%irp(ma+j+1) = a%irp(ma+j)
        end do
      class default 

        write(psb_err_unit,*) 'Implement SPGETBLK in RWEXTD!!!!!!!'
      end select
      call a%set_ncols(max(na,nb))

    else

      do j=ma+2,nr+1
        a%irp(j) = a%irp(j-1)
      end do

    end if

    call a%set_nrows(nr)


  type is (psb_ls_coo_sparse_mat) 
    nza = a%get_nzeros()

    if (present(b)) then 
      mb = b%get_nrows()
      nb = b%get_ncols()
      nzb = b%get_nzeros()
      call a%reallocate(nza+nzb)

      select type(b)
      type is (psb_ls_coo_sparse_mat) 

        if (rowscale_) then 
          do j=1,nzb
            if ((ma + b%ia(j)) <= nr) then 
              nza = nza + 1
              a%ia(nza)  = ma + b%ia(j)
              a%ja(nza)  = b%ja(j)
              a%val(nza) = b%val(j)
            end if
          enddo
        else
          do j=1,nzb
            if ((ma + b%ia(j)) <= nr) then 
              nza = nza + 1
              a%ia(nza)  = b%ia(j)
              a%ja(nza)  = b%ja(j)
              a%val(nza) = b%val(j)
            end if
          enddo
        endif
        call a%set_nzeros(nza)

      type is (psb_ls_csr_sparse_mat) 
        
        do i=1, min(nr-ma,mb)
          do jb = b%irp(i), b%irp(i+1)-1
            nza = nza + 1 
            a%val(nza) = b%val(jb)
            a%ia(nza)  = ma + i
            a%ja(nza)  = b%ja(jb)
          end do
        end do
        call a%set_nzeros(nza)

      class default 
        write(psb_err_unit,*) 'Implement SPGETBLK in RWEXTD!!!!!!!'

      end select

      call a%set_ncols(max(na,nb))
    endif

    call a%set_nrows(nr)

  class default 
    info = psb_err_unsupported_format_
    ch_err=a%get_fmt()
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999

  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)

  return

end subroutine psb_lsbase_rwextd
