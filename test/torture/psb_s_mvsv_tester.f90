module psb_s_mvsv_tester
contains
  subroutine s_usmv_2_n_ap3_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=3
    real*4 :: beta=1
    ! 1 1
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 2/)
    real*4 :: VA(3)=(/1, 1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/9, 6/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is ok"
  end subroutine s_usmv_2_n_ap3_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_t_ap3_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=3
    real*4 :: beta=1
    ! 1 0
    ! 1 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 1/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/9, 3/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is ok"
  end subroutine s_usmv_2_t_ap3_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_c_ap3_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=3
    real*4 :: beta=1
    ! 1 2
    ! 0 6

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 2/)
    real*4 :: VA(3)=(/1, 2, 6/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/6, 27/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is ok"
  end subroutine s_usmv_2_c_ap3_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=3
    real*4 :: beta=0
    ! 1 2
    ! 0 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 1/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 2/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/9, 0/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine s_usmv_2_n_ap3_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=3
    real*4 :: beta=0
    ! 1 3
    ! 2 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 1/)
    real*4 :: VA(3)=(/1, 3, 2/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/9, 9/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine s_usmv_2_t_ap3_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=3
    real*4 :: beta=0
    ! 1 0
    ! 1 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 1/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/6, 0/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine s_usmv_2_c_ap3_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_n_ap1_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=1
    real*4 :: beta=1
    ! 1 0
    ! 0 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=1
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(1)=(/1/)
    integer(psb_ipk_) :: JA(1)=(/1/)
    real*4 :: VA(1)=(/1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/4, 3/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is ok"
  end subroutine s_usmv_2_n_ap1_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_t_ap1_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=1
    real*4 :: beta=1
    ! 1 0
    ! 1 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 1/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/5, 3/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is ok"
  end subroutine s_usmv_2_t_ap1_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_c_ap1_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=1
    real*4 :: beta=1
    ! 1 2
    ! 5 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    real*4 :: VA(4)=(/1, 2, 5, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/9, 6/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is ok"
  end subroutine s_usmv_2_c_ap1_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=1
    real*4 :: beta=0
    ! 1 1
    ! 2 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 1/)
    real*4 :: VA(3)=(/1, 1, 2/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/2, 2/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine s_usmv_2_n_ap1_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=1
    real*4 :: beta=0
    ! 1 3
    ! 1 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    real*4 :: VA(4)=(/1, 3, 1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/2, 4/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine s_usmv_2_t_ap1_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=1
    real*4 :: beta=0
    ! 1 0
    ! 2 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    real*4 :: VA(3)=(/1, 2, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/3, 1/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine s_usmv_2_c_ap1_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_n_am1_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-1
    real*4 :: beta=1
    ! 1 3
    ! 0 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 1/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 3/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-1, 3/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is ok"
  end subroutine s_usmv_2_n_am1_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_t_am1_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-1
    real*4 :: beta=1
    ! 1 1
    ! 0 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 1/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/2, 2/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is ok"
  end subroutine s_usmv_2_t_am1_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_c_am1_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-1
    real*4 :: beta=1
    ! 1 0
    ! 1 2

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    real*4 :: VA(3)=(/1, 1, 2/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/1, 1/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is ok"
  end subroutine s_usmv_2_c_am1_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-1
    real*4 :: beta=0
    ! 1 0
    ! 1 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    real*4 :: VA(3)=(/1, 1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-1, -2/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine s_usmv_2_n_am1_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-1
    real*4 :: beta=0
    ! 1 4
    ! 3 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    real*4 :: VA(4)=(/1, 4, 3, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-4, -5/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine s_usmv_2_t_am1_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-1
    real*4 :: beta=0
    ! 1 1
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 2/)
    real*4 :: VA(3)=(/1, 1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-1, -2/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine s_usmv_2_c_am1_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_n_am3_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-3
    real*4 :: beta=1
    ! 1 3
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 2/)
    real*4 :: VA(3)=(/1, 3, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-9, 0/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is ok"
  end subroutine s_usmv_2_n_am3_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_t_am3_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-3
    real*4 :: beta=1
    ! 1 4
    ! 1 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 1/)
    real*4 :: VA(3)=(/1, 4, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-3, -9/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is ok"
  end subroutine s_usmv_2_t_am3_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_c_am3_bp1_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-3
    real*4 :: beta=1
    ! 1 1
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 2/)
    real*4 :: VA(3)=(/1, 1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/0, -3/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is ok"
  end subroutine s_usmv_2_c_am3_bp1_ix1_iy1
  ! 

  subroutine s_usmv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-3
    real*4 :: beta=0
    ! 1 0
    ! 2 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 1/)
    real*4 :: VA(2)=(/1, 2/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-3, -6/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine s_usmv_2_n_am3_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-3
    real*4 :: beta=0
    ! 1 0
    ! 0 0

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=1
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(1)=(/1/)
    integer(psb_ipk_) :: JA(1)=(/1/)
    real*4 :: VA(1)=(/1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-3, 0/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine s_usmv_2_t_am3_bm0_ix1_iy1
  ! 

  subroutine s_usmv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    real*4 :: alpha=-3
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/-3, -3/)! reference cy after 
    real*4 :: bcy(2)=(/3, 3/)! reference bcy before 
    real*4 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spmm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spmm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine s_usmv_2_c_am3_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_n_ap3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=3
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/3, 3/)! reference x 
    real*4 :: cy(2)=(/9, 9/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)
    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine s_ussv_2_n_ap3_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_t_ap3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=3
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/3, 3/)! reference x 
    real*4 :: cy(2)=(/9, 9/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine s_ussv_2_t_ap3_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_c_ap3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=3
    real*4 :: beta=0
    ! 1 0
    ! 1 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    real*4 :: VA(3)=(/1, 1, 1/)
    real*4 :: x(2)=(/6, 3/)! reference x 
    real*4 :: cy(2)=(/9, 9/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,i,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine s_ussv_2_c_ap3_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_n_ap1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=1
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/1, 1/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine s_ussv_2_n_ap1_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_t_ap1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=1
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/1, 1/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine s_ussv_2_t_ap1_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_c_ap1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=1
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/1, 1/)! reference x 
    real*4 :: cy(2)=(/1, 1/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine s_ussv_2_c_ap1_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_n_am1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=-1
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/-1, -1/)! reference x 
    real*4 :: cy(2)=(/1, 1/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine s_ussv_2_n_am1_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_t_am1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=-1
    real*4 :: beta=0
    ! 1 0
    ! 3 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    real*4 :: VA(3)=(/1, 3, 1/)
    real*4 :: x(2)=(/-4, -1/)! reference x 
    real*4 :: cy(2)=(/1, 1/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine s_ussv_2_t_am1_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_c_am1_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=-1
    real*4 :: beta=0
    ! 1 0
    ! 2 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    real*4 :: VA(3)=(/1, 2, 1/)
    real*4 :: x(2)=(/-3, -1/)! reference x 
    real*4 :: cy(2)=(/1, 1/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine s_ussv_2_c_am1_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_n_am3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=-3
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/-3, -3/)! reference x 
    real*4 :: cy(2)=(/9, 9/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine s_ussv_2_n_am3_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_t_am3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=-3
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/-3, -3/)! reference x 
    real*4 :: cy(2)=(/9, 9/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine s_ussv_2_t_am3_bm0_ix1_iy1
  ! 

  subroutine s_ussv_2_c_am3_bm0_ix1_iy1(res,afmt,ictxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_sspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    integer(psb_ipk_) :: ictxt, iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    real*4 :: alpha=-3
    real*4 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    real*4 :: VA(2)=(/1, 1/)
    real*4 :: x(2)=(/-3, -3/)! reference x 
    real*4 :: cy(2)=(/9, 9/)! reference cy after 
    real*4 :: bcy(2)=(/0, 0/)! reference bcy before 
    real*4 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ictxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ictxt)
    call psb_cdall(ictxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ictxt)
    call psb_spins(nnz,IA,JA,VA,a,desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_cdasb(desc_a,info)
    if (info /= psb_success_)goto 9996
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)
    if(info.ne.0)print *,"matrix assembly failed"
    if(info.ne.0)goto 9996

    call psb_spsm(alpha,A,x,beta,y,desc_a,info,transa)
    if(info.ne.0)print *,"psb_spsm failed"
    if(info.ne.0)goto 9996
    do i=1,2
      if(y(i) /= cy(i))print*,"results mismatch:",y,"instead of",cy
      if(y(i) /= cy(i))info=-1
      if(y(i) /= cy(i))goto 9996
    enddo
9996 continue
    if(info /= psb_success_)res=res+1
    call psb_spfree(a,desc_a,info)
    if (info /= psb_success_)goto 9997
9997 continue
    if(info /= psb_success_)res=res+1
    call psb_cdfree(desc_a,info)
    if (info /= psb_success_)goto 9998
9998 continue
    if(info /= psb_success_)res=res+1
9999 continue
    if(info /= psb_success_)res=res+1
    if(res /= 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on s matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine s_ussv_2_c_am3_bm0_ix1_iy1
  ! 
end module psb_s_mvsv_tester
