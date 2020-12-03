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

subroutine psb_c_invk_bld(a,fill1, fill2,lmat,d,umat,desc,info,blck)

  use psb_base_mod
  use psb_c_invk_fact_mod, psb_protect_name =>  psb_c_invk_bld
  use psb_c_ilu_fact_mod
  implicit none

  ! Arguments
  type(psb_cspmat_type), intent(in), target   :: a
  integer(psb_ipk_), intent(in)               :: fill1, fill2
  type(psb_cspmat_type), intent(inout)        :: lmat, umat
  complex(psb_spk_), allocatable                 :: d(:)
  Type(psb_desc_type), Intent(inout)          :: desc
  integer(psb_ipk_), intent(out)              :: info
  type(psb_cspmat_type), intent(in), optional :: blck
  !
  integer(psb_ipk_)      :: i, nztota, err_act, n_row, nrow_a, n_col
  type(psb_cspmat_type)  :: atmp
  complex(psb_spk_), allocatable :: pq(:), pd(:)
  integer(psb_ipk_), allocatable :: uplevs(:)
  integer(psb_ipk_)   :: debug_level, debug_unit
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np,me
  integer(psb_ipk_)   :: nzrmax
  character(len=20)   :: name, ch_err


  info = psb_success_
  name='psb_cinvk_bld'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt        = psb_cd_get_context(desc)
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' start'

  !
  ! Check the memory available to hold the incomplete L and U factors
  ! and allocate it if needed
  !
  nrow_a = a%get_nrows()
  nztota = a%get_nzeros()

  if (present(blck)) then
    nztota = nztota + blck%get_nzeros()
  end if
  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & ': out get_nnzeros',nrow_a,nztota,&
       & a%get_nrows(),a%get_ncols(),a%get_nzeros()


  n_row  = psb_cd_get_local_rows(desc)
  n_col  = psb_cd_get_local_cols(desc)
  allocate(pd(n_row),stat=info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999
  end if


  call lmat%allocate(n_row,n_row,info,nz=nztota)
  if (info == psb_success_) call umat%allocate(n_row,n_row,info,nz=nztota)


  call psb_iluk_fact(fill1,psb_ilu_n_,a,lmat,umat,pd,info,blck=blck)
       !,uplevs=uplevs)
  !call psb_ciluk_fact(fill1,psb_ilu_n_,a,lmat,umat,pd,info,blck=blck)

  if (info == psb_success_) call atmp%allocate(n_row,n_row,info,nz=nztota)
  if(info/=0) then
    info=psb_err_from_subroutine_
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  !
  ! Compute the aprox U^-1  and L^-1
  !
  call psb_sparse_invk(n_row,umat,atmp,fill2,info)
  if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
  if (info == psb_success_) call lmat%transp()
  if (info == psb_success_) call psb_sparse_invk(n_row,lmat,atmp,fill2,info)
  if (info == psb_success_) call psb_move_alloc(atmp,lmat,info)
  if (info == psb_success_) call lmat%transp()
  ! Done. Hopefully....

  if (info /= psb_success_) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='invt')
    goto 9999
  end if

  call psb_move_alloc(pd,d,info)
  call lmat%set_asb()
  call lmat%trim()
  call umat%set_asb()
  call umat%trim()

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),' end'


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_c_invk_bld

subroutine psb_csparse_invk(n,a,z,fill_in,info,inlevs)

  use psb_base_mod
  use psb_c_invk_fact_mod, psb_protect_name => psb_csparse_invk

  integer(psb_ipk_), intent(in)           :: n
  type(psb_cspmat_type), intent(in)       :: a
  type(psb_cspmat_type), intent(inout)    :: z
  integer(psb_ipk_), intent(in)           :: fill_in
  integer(psb_ipk_), intent(out)          :: info
  integer(psb_ipk_), intent(in), optional :: inlevs(:)
  !
  integer(psb_ipk_) :: i,j,k, err_act, nz, nzra, nzrz, ipz1,ipz2, nzz, ip1, ip2, l2
  integer(psb_ipk_), allocatable        :: ia(:), ja(:), iz(:), jz(:)
  complex(psb_spk_), allocatable :: zw(:), val(:), valz(:)
  integer(psb_ipk_), allocatable        :: uplevs(:), rowlevs(:), idxs(:)
  complex(psb_spk_), allocatable :: row(:)
  type(psb_c_coo_sparse_mat)  :: trw
  type(psb_c_csr_sparse_mat)  :: acsr, zcsr
  integer(psb_ipk_)           :: ktrw, nidx
  type(psb_i_heap)            :: heap

  real(psb_dpk_)     :: alpha
  character(len=20)  :: name='psb_sp_invk'

  info = psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if

  if (.not.(a%is_triangle().and.a%is_unit().and.a%is_upper())) then
    write(psb_err_unit,*) 'Wrong A '
    info = psb_err_internal_error_
    call psb_errpush(psb_err_internal_error_,name,a_err='wrong A')
    goto 9999
  end if
  call a%cp_to(acsr)
  call trw%allocate(izero,izero,ione)
  if (info == psb_success_) allocate(zw(n),iz(n),valz(n),&
       & row(n),rowlevs(n),stat=info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999
  end if

  allocate(uplevs(acsr%get_nzeros()),stat=info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='Allocate')
    goto 9999
  end if
  uplevs(:)  = 0
  row(:)     = czero
  rowlevs(:) = -(n+1)

  call zcsr%allocate(n,n,n*fill_in)
  call zcsr%set_triangle()
  call zcsr%set_unit(.false.)
  call zcsr%set_upper()
  call psb_ensure_size(n+1, idxs,  info)


  !
  !
  zcsr%irp(1)  = 1
  nzz          = 0

  l2 = 0
  outer: do i = 1, n-1
    ! ZW = e_i
    call psb_invk_copyin(i,n,acsr,ione,n,row,rowlevs,heap,ktrw,trw,info,&
         & sign=-sone,inlevs=inlevs)
    row(i)     = cone
    rowlevs(i) = 0

    ! Update loop
    call psb_invk_inv(fill_in,i,row,rowlevs,heap,&
         & acsr%ja,acsr%irp,acsr%val,uplevs,nidx,idxs,info)

    call psb_invk_copyout(fill_in,i,n,row,rowlevs,nidx,idxs,&
         & l2,zcsr%ja,zcsr%irp,zcsr%val,info)

    nzz = l2
  end do outer
  if (info /= psb_success_) then
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='mainloop')
    goto 9999
  end if
  ipz1 = nzz+1
  call psb_ensure_size(ipz1,zcsr%val,info)
  call psb_ensure_size(ipz1,zcsr%ja,info)
  zcsr%val(ipz1) = cone
  zcsr%ja(ipz1)  = n
  zcsr%irp(n+1)  = ipz1+1
  call zcsr%set_sorted()
  call z%mv_from(zcsr)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_csparse_invk

subroutine psb_c_invk_copyin(i,m,a,jmin,jmax,row,rowlevs,heap,ktrw,trw,info,sign,inlevs)

  use psb_base_mod
  use psb_c_invk_fact_mod, psb_protect_name => psb_c_invk_copyin

  implicit none

  ! Arguments
  type(psb_c_csr_sparse_mat), intent(in)    :: a
  type(psb_c_coo_sparse_mat), intent(inout) :: trw
  integer(psb_ipk_), intent(in)             :: i,m,jmin,jmax
  integer(psb_ipk_), intent(inout)          :: ktrw,info
  integer(psb_ipk_), intent(inout)          :: rowlevs(:)
  complex(psb_spk_), intent(inout)             :: row(:)
  type(psb_i_heap), intent(inout)           :: heap
  real(psb_spk_), optional, intent(in)      :: sign
  integer(psb_ipk_), intent(in), optional   :: inlevs(:)

  ! Local variables
  integer(psb_ipk_)             :: k,j,irb,err_act, nz
  integer(psb_ipk_), parameter  :: nrb=16
  real(psb_spk_)                :: sign_
  character(len=20), parameter  :: name='invk_copyin'
  character(len=20)             :: ch_err

  info = psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if

  call heap%init(info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_init_heap')
    goto 9999
  end if

  if (present(sign)) then
    sign_ = sign
  else
    sign_ = sone
  end if


  !
  ! Take a fast shortcut if the matrix is stored in CSR format
  !
  if (present(inlevs)) then
    do j = a%irp(i), a%irp(i+1) - 1
      k          = a%ja(j)
      if ((jmin<=k).and.(k<=jmax)) then
        row(k)     = sign_ * a%val(j)
        rowlevs(k) = inlevs(j)
        call heap%insert(k,info)
      end if
    end do
  else
    do j = a%irp(i), a%irp(i+1) - 1
      k          = a%ja(j)
      if ((jmin<=k).and.(k<=jmax)) then
        row(k)     = sign_ * a%val(j)
        rowlevs(k) = 0
        call heap%insert(k,info)
      end if
    end do
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_c_invk_copyin


subroutine psb_c_invk_copyout(fill_in,i,m,row,rowlevs,nidx,idxs,&
     &  l2,uia1,uia2,uaspk,info)

  use psb_base_mod
  use psb_c_invk_fact_mod, psb_protect_name => psb_c_invk_copyout

  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)                 :: fill_in, i, m, nidx
  integer(psb_ipk_), intent(inout)              :: l2, info
  integer(psb_ipk_), intent(inout)              :: rowlevs(:), idxs(:)
  integer(psb_ipk_), allocatable, intent(inout) :: uia1(:), uia2(:)
  complex(psb_spk_), allocatable, intent(inout)    :: uaspk(:)
  complex(psb_spk_), intent(inout)                 :: row(:)

  ! Local variables
  integer(psb_ipk_)              :: j,isz,err_act,int_err(5),idxp
  character(len=20), parameter  :: name='psb_ciluk_factint'
  character(len=20)             :: ch_err

  info = psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if

  do idxp=1,nidx

    j = idxs(idxp)

    if (j>=i) then
      !
      ! Copy the upper part of the row
      !
      if (rowlevs(j) <= fill_in) then
        l2     = l2 + 1
        if (size(uaspk) < l2) then
          !
          ! Figure out a good reallocation size!
          !
          isz  = max(int(1.2*l2),l2+100)
          call psb_realloc(isz,uaspk,info)
          if (info == psb_success_) call psb_realloc(isz,uia1,info)
          if (info /= psb_success_) then
            info=psb_err_from_subroutine_
            call psb_errpush(info,name,a_err='Allocate')
            goto 9999
          end if
        end if
        uia1(l2)   = j
        uaspk(l2)  = row(j)
      end if
      !
      ! Re-initialize row(j) and rowlevs(j)
      !
      row(j)     = czero
      rowlevs(j) = -(m+1)
    end if
  end do

  uia2(i+1) = l2 + 1

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_c_invk_copyout

subroutine psb_cinvk_inv(fill_in,i,row,rowlevs,heap,ja,irp,val,uplevs,&
     & nidx,idxs,info)

  use psb_base_mod
  use psb_c_invk_fact_mod, psb_protect_name => psb_cinvk_inv

  implicit none

  ! Arguments
  type(psb_i_heap), intent(inout)               :: heap
  integer(psb_ipk_), intent(in)                 :: i, fill_in
  integer(psb_ipk_), intent(inout)              :: nidx,info
  integer(psb_ipk_), intent(inout)              :: rowlevs(:)
  integer(psb_ipk_), allocatable, intent(inout) :: idxs(:)
  integer(psb_ipk_), intent(in)                 :: ja(:),irp(:),uplevs(:)
  complex(psb_spk_), intent(in)                    :: val(:)
  complex(psb_spk_), intent(inout)                 :: row(:)

  ! Local variables
  integer(psb_ipk_) :: k,j,lrwk,jj,lastk, iret
  real(psb_dpk_)    :: rwk


  info = psb_success_

  call psb_ensure_size(200, idxs,  info)
  if (info /= psb_success_) return
  nidx    = 1
  idxs(1) = i
  lastk   = i

  !
  ! Do while there are indices to be processed
  !
  do
    ! Beware: (iret < 0) means that the heap is empty, not an error.
    call heap%get_first(k,iret)
    if (iret < 0) then
!!$        write(psb_err_unit,*) 'IINVK: ',i,' returning at ',lastk
      return
    end if

    !
    ! Just in case an index has been put on the heap more than once.
    !
    if (k == lastk) cycle

    lastk = k
    nidx = nidx + 1
    if (nidx>size(idxs)) then
      call psb_realloc(nidx+psb_heap_resize,idxs,info)
      if (info /= psb_success_) return
    end if
    idxs(nidx) = k

    if ((row(k) /= czero).and.(rowlevs(k) <= fill_in)) then
      !
      ! Note: since U is scaled while copying it out (see iluk_copyout),
      ! we can use rwk in the update below
      !
      rwk    = row(k)
      lrwk   = rowlevs(k)

      do jj=irp(k),irp(k+1)-1
        j = ja(jj)
        if (j<=k) then
          info = -i
          return
        endif
        !
        ! Insert the index into the heap for further processing.
        ! The fill levels are initialized to a negative value. If we find
        ! one, it means that it is an as yet untouched index, so we need
        ! to insert it; otherwise it is already on the heap, there is no
        ! need to insert it more than once.
        !
        if (rowlevs(j)<0) then
          call heap%insert(j,info)
          if (info /= psb_success_) return
          rowlevs(j) = abs(rowlevs(j))
        end if
        !
        ! Update row(j) and the corresponding fill level
        !
        row(j)     = row(j) - rwk * val(jj)
        rowlevs(j) = min(rowlevs(j),lrwk+uplevs(jj)+1)
      end do

    end if
  end do

end subroutine psb_cinvk_inv
