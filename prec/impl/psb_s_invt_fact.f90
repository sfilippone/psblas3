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
!    Moved here from AMG4PSBLAS, original copyright below.
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

subroutine psb_s_invt_bld(a,fillin,invfill,thresh,invthresh,&
     & lmat,d,umat,desc,info,blck)

  use psb_base_mod
  use psb_s_invt_fact_mod, psb_protect_name =>  psb_s_invt_bld
  use psb_s_ilu_fact_mod
  implicit none

  ! Arguments
  type(psb_sspmat_type), intent(inout), target   :: a
  integer(psb_ipk_), intent(in)               :: fillin,invfill
  real(psb_spk_), intent(in)                :: thresh
  real(psb_spk_), intent(in)                  :: invthresh
  type(psb_sspmat_type), intent(inout)        :: lmat, umat
  real(psb_spk_), allocatable                 :: d(:)
  Type(psb_desc_type), Intent(inout)          :: desc
  integer(psb_ipk_), intent(out)              :: info
  type(psb_sspmat_type), intent(in), optional :: blck
  !
  integer(psb_ipk_) :: i, nztota, err_act, n_row, nrow_a, n_col
  type(psb_sspmat_type)       :: atmp
  real(psb_spk_), allocatable :: pq(:), pd(:), w(:)
  integer(psb_ipk_) :: debug_level, debug_unit
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_)   :: np,me
  integer(psb_ipk_)   :: nzrmax
  real(psb_spk_)    :: sp_thresh
  character(len=20) :: name, ch_err, fname


  info = psb_success_
  name='psb_sinvt_fact'
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
  allocate(pd(n_row),w(n_row),stat=info)
  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999
  end if

  nzrmax    = fillin
  sp_thresh = thresh

  call lmat%allocate(n_row,n_row,info,nz=nztota)
  if (info == psb_success_) call umat%allocate(n_row,n_row,info,nz=nztota)

  if (info == 0) call psb_ilut_fact(nzrmax,sp_thresh,&
       & a,lmat,umat,pd,info,blck=blck,iscale=psb_ilu_scale_maxval_)

  if (info == psb_success_) call atmp%allocate(n_row,n_row,info,nz=nztota)
  if(info/=0) then
    info=psb_err_from_subroutine_
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (.false.) then
!!$    if (debug_level >= psb_debug_inner_) then
    write(fname,'(a,i0,a)') 'invt-lo-',me,'.mtx'
    call lmat%print(fname,head="INVTLOW")
    write(fname,'(a,i0,a)') 'invt-up-',me,'.mtx'
    call umat%print(fname,head="INVTUPP")
  end if

  !
  ! Compute the approx U^-1  and L^-1
  !
  nzrmax    = invfill
  call psb_ssparse_invt(n_row,umat,atmp,nzrmax,invthresh,info)
  if (info == psb_success_) call psb_move_alloc(atmp,umat,info)
  if (info == psb_success_) call lmat%transp()
  if (info == psb_success_) call psb_ssparse_invt(n_row,lmat,atmp,nzrmax,invthresh,info)
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
end subroutine psb_s_invt_bld

subroutine psb_ssparse_invt(n,a,z,nzrmax,sp_thresh,info)

  use psb_base_mod
  use psb_s_invt_fact_mod, psb_protect_name => psb_ssparse_invt

  implicit none
  integer(psb_ipk_), intent(in)        :: n
  type(psb_sspmat_type), intent(inout) :: a
  type(psb_sspmat_type), intent(inout) :: z
  integer(psb_ipk_), intent(in)        :: nzrmax
  real(psb_spk_), intent(in)           :: sp_thresh
  integer(psb_ipk_), intent(out)       :: info
  !
  integer(psb_ipk_) :: i,j,k, err_act, nz, nzra, nzrz, ipz1,ipz2, nzz, ip1, ip2, l2
  integer(psb_ipk_), allocatable :: ia(:), ja(:), iz(:),jz(:)
  real(psb_spk_), allocatable    :: zw(:), val(:), valz(:)
  integer(psb_ipk_), allocatable :: uplevs(:), rowlevs(:),idxs(:)
  real(psb_spk_), allocatable :: row(:)
  type(psb_s_coo_sparse_mat)  :: trw
  type(psb_s_csr_sparse_mat)  :: acsr, zcsr
  integer(psb_ipk_)           :: ktrw, nidx, nlw,nup,jmaxup
  type(psb_i_heap)            :: heap
  real(psb_spk_)     :: alpha, nrmi
  character(len=20)  :: name='psb_sp_invt'

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

  call zcsr%allocate(n,n,n*nzrmax)
  call zcsr%set_triangle()
  call zcsr%set_unit(.false.)
  call zcsr%set_upper()
  !
  !
  nzz        = 0
  row(:)     = szero
  rowlevs(:) = 0
  l2         = 0
  zcsr%irp(1) = 1

  outer: do i = 1, n-1
    ! ZW = e_i
    call psb_s_invt_copyin(i,n,acsr,i,ione,n,nlw,nup,jmaxup,nrmi,row,&
         & heap,rowlevs,ktrw,trw,info,sign=-sone)
    if (info /= 0) exit
    row(i) = sone
    ! Adjust norm
    if (nrmi < sone) then
      nrmi = sqrt(sone + nrmi**2)
    else
      nrmi = nrmi*sqrt(sone+sone/(nrmi**2))
    end if

    call psb_invt_inv(sp_thresh,i,nrmi,row,heap,rowlevs,&
         & acsr%ja,acsr%irp,acsr%val,nidx,idxs,info)
    if (info /= 0) exit
!!$    write(0,*) 'Calling copyout ',nzrmax,nlw,nup,nidx,l2
    call psb_s_invt_copyout(nzrmax,sp_thresh,i,n,nlw,nup,jmaxup,nrmi,row,&
         & nidx,idxs,l2,zcsr%ja,zcsr%irp,zcsr%val,info)
    if (info /= 0) exit
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
  zcsr%val(ipz1) = sone
  zcsr%ja(ipz1)  = n
  zcsr%irp(n+1)  = ipz1+1

  call z%mv_from(zcsr)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_ssparse_invt

subroutine psb_s_invt_copyin(i,m,a,jd,jmin,jmax,nlw,nup,jmaxup,nrmi,row,heap,&
     & irwt,ktrw,trw,info,sign)
  use psb_base_mod
  use psb_d_invt_fact_mod, psb_protect_name => psb_d_invt_copyin
  implicit none
  type(psb_s_csr_sparse_mat), intent(in)    :: a
  type(psb_s_coo_sparse_mat), intent(inout) :: trw
  integer(psb_ipk_), intent(in)             :: i, m,jmin,jmax,jd
  integer(psb_ipk_), intent(inout)          :: ktrw,nlw,nup,jmaxup,info
  integer(psb_ipk_), intent(inout)          :: irwt(:)
  real(psb_spk_), intent(inout)          :: nrmi
  real(psb_spk_), intent(inout)            :: row(:)
  type(psb_i_heap), intent(inout)           :: heap
  real(psb_spk_), intent(in), optional    :: sign
  !
  integer(psb_ipk_)            :: k,j,irb,kin,nz, err_act
  integer(psb_ipk_), parameter :: nrb=16
  real(psb_dpk_)               :: dmaxup, sign_
  real(psb_dpk_), external     :: dnrm2
  character(len=20), parameter :: name='invt_copyin'

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
  sign_ = sone
  if (present(sign)) sign_ = sign
  !
  ! nrmi is the norm of the current sparse row (for the time being,
  ! we use the 2-norm).
  ! NOTE: the 2-norm below includes also elements that are outside
  ! [jmin:jmax] strictly. Is this really important? TO BE CHECKED.
  !

  nlw    = 0
  nup    = 0
  jmaxup = 0
  dmaxup = szero
  nrmi   = szero

  do j = a%irp(i), a%irp(i+1) - 1
    k = a%ja(j)
    if ((jmin<=k).and.(k<=jmax)) then
      row(k)     = sign_ * a%val(j)
      call heap%insert(k,info)
      irwt(k) = 1
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_insert_heap')
        goto 9999
      end if
    end if
    if (k<jd) nlw = nlw + 1
    if (k>jd) then
      nup = nup + 1
      if (abs(row(k))>dmaxup) then
        jmaxup = k
        dmaxup = abs(row(k))
      end if
    end if
  end do
  nz   = a%irp(i+1) - a%irp(i)
  nrmi = dnrm2(nz,a%val(a%irp(i):),ione)

  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_s_invt_copyin

subroutine psb_s_invt_copyout(fill_in,thres,i,m,nlw,nup,jmaxup,nrmi,row, &
     & nidx,idxs,l2,ja,irp,val,info)

  use psb_base_mod
  use psb_s_invt_fact_mod, psb_protect_name => psb_s_invt_copyout

  implicit none

  ! Arguments
  integer(psb_ipk_), intent(in)                 :: fill_in,i,m,nidx,nlw,nup,jmaxup
  integer(psb_ipk_), intent(in)                 :: idxs(:)
  integer(psb_ipk_), intent(inout)              :: l2, info
  integer(psb_ipk_), allocatable, intent(inout) :: ja(:),irp(:)
  real(psb_spk_), intent(in)                    :: thres,nrmi
  real(psb_spk_),allocatable, intent(inout)     :: val(:)
  real(psb_spk_), intent(inout)                 :: row(:)

  ! Local variables
  real(psb_dpk_),allocatable     :: xw(:)
  integer(psb_ipk_), allocatable :: xwid(:), indx(:)
  real(psb_dpk_)                 :: witem, wmin
  integer(psb_ipk_)              :: widx
  integer(psb_ipk_)              :: k,isz,err_act,int_err(5),idxp, nz
  type(psb_d_idx_heap)           :: heap
  character(len=20), parameter   :: name='invt_copyout'
  character(len=20)              :: ch_err
  logical                        :: fndmaxup

  info = psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_; goto 9999
  end if

  !
  ! Here we need to apply also the dropping rule base on the fill-in.
  ! We do it by putting into a heap the elements that are not dropped
  ! by using the 2-norm rule, and then copying them out.
  !
  ! The heap goes down on the entry absolute value, so the first item
  ! is the largest absolute value.
  !
!!$  write(0,*) 'invt_copyout ',nidx,nup+fill_in
  call heap%init(info,dir=psb_asort_down_)

  if (info == psb_success_) allocate(xwid(nidx),xw(nidx),indx(nidx),stat=info)
  if (info /= psb_success_) then
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/3*nidx/),&
         & a_err='real(psb_dpk_)')
    goto 9999
  end if

  !
  ! First the lower part
  !

  nz   = 0
  idxp = 0

  do

    idxp = idxp + 1
    if (idxp > nidx) exit
    if (idxs(idxp) >= i) exit
    widx      = idxs(idxp)
    witem     = row(widx)
    !
    ! Dropping rule based on the 2-norm
    !
    if (abs(witem) < thres*nrmi) cycle

    nz       = nz + 1
    xw(nz)   = witem
    xwid(nz) = widx
    call heap%insert(witem,widx,info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_insert_heap')
      goto 9999
    end if
  end do

  if (nz > 1) then
    write(psb_err_unit,*) 'Warning: lower triangle from invt???? '
  end if


  if (idxp <= size(idxs)) then
    if (idxs(idxp) < i) then
      do
        idxp = idxp + 1
        if (idxp > nidx) exit
        if (idxs(idxp) >= i) exit
      end do
    end if
  end if
  idxp = idxp - 1
  nz   = 0
  wmin=HUGE(wmin)
  if (.false.) then
    do

      idxp = idxp + 1
      if (idxp > nidx) exit
      widx      = idxs(idxp)
      if (widx < i) then
        write(psb_err_unit,*) 'Warning: lower triangle in upper copy',widx,i,idxp,idxs(idxp)
        cycle
      end if
      if (widx > m) then
        cycle
      end if
      witem     = row(widx)
      !
      ! Dropping rule based on the 2-norm. But keep the jmaxup-th entry anyway.
      !
      if ((widx /= jmaxup) .and. (widx /= i) .and. (abs(witem) < thres*nrmi)) then
        cycle
      end if
      if ((widx/=jmaxup).and.(nz > nup+fill_in)) then
        if (abs(witem) < wmin) cycle
      endif
      wmin = min(abs(witem),wmin)
      nz       = nz + 1
      xw(nz)   = witem
      xwid(nz) = widx
      call heap%insert(witem,widx,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_insert_heap')
        goto 9999
      end if

    end do

    !
    ! Now we have to take out the first nup-fill_in entries. But make sure
    ! we include entry jmaxup.
    !
    if (nz <= nup+fill_in) then
      !
      ! Just copy everything from xw
      !
      fndmaxup=.true.
    else
      fndmaxup = .false.
      nz = nup+fill_in
      do k=1,nz
        call heap%get_first(witem,widx,info)
        xw(k)   = witem
        xwid(k) = widx
        if (widx == jmaxup) fndmaxup=.true.
      end do
    end if
    if ((i<jmaxup).and.(jmaxup<=m)) then
      if (.not.fndmaxup) then
        !
        ! Include entry jmaxup, if it is not already there.
        ! Put it in the place of the smallest coefficient.
        !
        xw(nz)   = row(jmaxup)
        xwid(nz) = jmaxup
      endif
    end if

  else if (.true.) then

    do

      idxp = idxp + 1
      if (idxp > nidx) exit
      widx      = idxs(idxp)
      if (widx < i) then
        write(psb_err_unit,*) 'Warning: lower triangle in upper copy',widx,i,idxp,idxs(idxp)
        cycle
      end if
      if (widx > m) then
        cycle
      end if
      witem     = row(widx)
      !
      ! Dropping rule based on the 2-norm. But keep the jmaxup-th entry anyway.
      !
      if ((widx /= i) .and. (abs(witem) < thres*nrmi)) then
        cycle
      end if
      if (nz > nup+fill_in) then
        if (abs(witem) < wmin) cycle
      endif
      wmin = min(abs(witem),wmin)
      nz       = nz + 1
      xw(nz)   = witem
      xwid(nz) = widx
      call heap%insert(witem,widx,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_insert_heap')
        goto 9999
      end if

    end do

    !
    ! Now we have to take out the first nup-fill_in entries. But make sure
    ! we include entry jmaxup.
    !
    if (nz >  nup+fill_in) then
      nz = nup+fill_in
      do k=1,nz
        call heap%get_first(witem,widx,info)
        xw(k)   = witem
        xwid(k) = widx
      end do
    end if
  end if

  !
  ! Now we put things back into ascending column order
  !
  call psb_msort(xwid(1:nz),indx(1:nz),dir=psb_sort_up_)

  !
  ! Copy out the upper part of the row
  !
  do k=1,nz
    l2     = l2 + 1
    if (size(val) < l2) then
      !
      ! Figure out a good reallocation size!
      !
      isz  = max(int(1.2*l2),l2+100)
      call psb_realloc(isz,val,info)
      if (info == psb_success_) call psb_realloc(isz,ja,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='Allocate')
        goto 9999
      end if
    end if
    ja(l2)   = xwid(k)
    val(l2)  = xw(indx(k))
  end do

  !
  ! Set row to zero
  !
  do idxp=1,nidx
    row(idxs(idxp)) = szero
  end do

  irp(i+1) = l2 + 1

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_s_invt_copyout

subroutine psb_s_invt_inv(thres,i,nrmi,row,heap,irwt,ja,irp,val,nidx,idxs,info)

  use psb_base_mod
  use psb_s_invt_fact_mod, psb_protect_name => psb_s_invt_inv

  implicit none

  ! Arguments
  type(psb_i_heap), intent(inout)     :: heap
  integer(psb_ipk_), intent(in)       :: i
  integer(psb_ipk_), intent(inout)    :: nidx,info
  integer(psb_ipk_), intent(inout)    :: irwt(:)
  real(psb_spk_), intent(in)          :: thres,nrmi
  integer(psb_ipk_), allocatable, intent(inout) :: idxs(:)
  integer(psb_ipk_), intent(in)       :: ja(:),irp(:)
  real(psb_spk_), intent(in)          :: val(:)
  real(psb_spk_), intent(inout)       :: row(:)

  ! Local Variables
  integer(psb_ipk_) :: k,j,jj,lastk,iret
  real(psb_dpk_)    :: rwk, alpha

  info  = psb_success_

  call psb_ensure_size(200, idxs,  info)
  if (info /= psb_success_) return
  nidx    = 1
  idxs(1) = i
  lastk   = i
  irwt(i) = 1
!!$  write(0,*) 'Drop Threshold ',thres*nrmi
  !
  ! Do while there are indices to be processed
  !
  do

    call heap%get_first(k,iret)
    if (iret < 0)  exit

    !
    ! An index may have been put on the heap more than once.
    ! Should not happen, but just in case.
    !
    if (k == lastk)  cycle
    lastk = k

    !
    ! Dropping rule based on the threshold: compare the absolute
    ! value of each updated entry of row with thres * 2-norm of row.
    !
    rwk    = row(k)

    if (abs(rwk) < thres*nrmi) then
      !
      ! Drop the entry.
      !
      row(k)  = dzero
      irwt(k) = 0
      cycle
    else
      !
      ! Note: since U is scaled while copying it out (see ilut_copyout),
      ! we can use rwk in the update below.
      !
      do jj=irp(k),irp(k+1)-1
        j = ja(jj)
        if (j<=k) then
          info = -i
          return
        endif
        !
        ! Update row(j) and, if it is not to be discarded, insert
        ! its index into the heap for further processing.
        !
        row(j)     = row(j) - rwk * val(jj)
        if (irwt(j) == 0) then
          if (abs(row(j)) < thres*nrmi) then
            !
            ! Drop the entry.
            !
            row(j)  = dzero
          else
            !
            ! Do the insertion.
            !
            call heap%insert(j,info)
            if (info /= psb_success_) return
            irwt(j) = 1
          end if
        end if
      end do
    end if

    !
    ! If we get here it is an index we need to keep on copyout.
    !

    nidx       = nidx + 1
    call psb_ensure_size(nidx,idxs,info,addsz=psb_heap_resize)
    if (info /= psb_success_) return
    idxs(nidx) = k
    irwt(k)    = 0
  end do

  irwt(i) = 0
end subroutine psb_s_invt_inv
