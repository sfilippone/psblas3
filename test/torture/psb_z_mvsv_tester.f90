module psb_z_mvsv_tester
contains

  subroutine z_usmv_2_n_ap3_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=3
    complex*16 :: beta=1
    ! 1+1i 1+0i
    ! 5+1i 1+1i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    complex*16 :: VA(4)=(/(1.e0,1.e0), (1.e0,0.e0), (5.e0,1.e0), (1,1)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(9.e0,3.e0), (21,6)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=n is ok"
  end subroutine z_usmv_2_n_ap3_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_t_ap3_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=3
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 2+3i 2+2i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (2.e0,3.e0), (2,2)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(12.e0,12.e0), (9,6)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=t is ok"
  end subroutine z_usmv_2_t_ap3_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_c_ap3_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=3
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 2+0i 1+3i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (2.e0,0.e0), (1,3)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(12.e0,-3.e0), (6,-9)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 1 incx=1 incy=1 trans=c is ok"
  end subroutine z_usmv_2_c_ap3_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_n_ap3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=3
    complex*16 :: beta=0
    ! 1+1i 0+0i
    ! 0+0i 0+2i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    complex*16 :: VA(2)=(/(1.e0,1.e0), (0,2)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(3.e0,3.e0), (0,6)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine z_usmv_2_n_ap3_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_t_ap3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=3
    complex*16 :: beta=0
    ! 1+1i 0+1i
    ! 1+0i 3+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    complex*16 :: VA(4)=(/(1.e0,1.e0), (0.e0,1.e0), (1.e0,0.e0), (3,0)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(6.e0,3.e0), (9,3)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine z_usmv_2_t_ap3_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_c_ap3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=3
    complex*16 :: beta=0
    ! 1+1i 0+0i
    ! 1+3i 0+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 1/)
    complex*16 :: VA(2)=(/(1.e0,1.e0), (1,3)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(6.e0,-12.e0), (0,0)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 3 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine z_usmv_2_c_ap3_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_n_ap1_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=1
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 0+3i 0+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 1/)
    complex*16 :: VA(2)=(/(1.e0,1.e0), (0,3)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(4.e0,1.e0), (3,3)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=n is ok"
  end subroutine z_usmv_2_n_ap1_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_t_ap1_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=1
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 0+1i 1+3i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (0.e0,1.e0), (1,3)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(4.e0,2.e0), (4,3)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=t is ok"
  end subroutine z_usmv_2_t_ap1_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_c_ap1_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=1
    complex*16 :: beta=1
    ! 1+1i 1+3i
    ! 0+0i 0+2i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (1.e0,3.e0), (0,2)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(4.e0,-1.e0), (4,-5)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 1 incx=1 incy=1 trans=c is ok"
  end subroutine z_usmv_2_c_ap1_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_n_ap1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=1
    complex*16 :: beta=0
    ! 1+1i 3+2i
    ! 0+0i 0+4i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (3.e0,2.e0), (0,4)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(4.e0,3.e0), (0,4)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine z_usmv_2_n_ap1_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_t_ap1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=1
    complex*16 :: beta=0
    ! 1+1i 0+0i
    ! 0+4i 1+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (0.e0,4.e0), (1,0)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(1.e0,5.e0), (1,0)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine z_usmv_2_t_ap1_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_c_ap1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=1
    complex*16 :: beta=0
    ! 1+1i 0+0i
    ! 1+3i 0+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 1/)
    complex*16 :: VA(2)=(/(1.e0,1.e0), (1,3)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(2.e0,-4.e0), (0,0)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha= 1 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine z_usmv_2_c_ap1_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_n_am1_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-1
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 3+2i 1+1i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (3.e0,2.e0), (1,1)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(2.e0,-1.e0), (-1,-3)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=n is ok"
  end subroutine z_usmv_2_n_am1_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_t_am1_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-1
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 0+3i 0+1i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (0.e0,3.e0), (0,1)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(2.e0,-4.e0), (3,-1)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=t is ok"
  end subroutine z_usmv_2_t_am1_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_c_am1_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-1
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 0+4i 1+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (0.e0,4.e0), (1,0)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(2.e0,5.e0), (2,0)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 1 incx=1 incy=1 trans=c is ok"
  end subroutine z_usmv_2_c_am1_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_n_am1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-1
    complex*16 :: beta=0
    ! 1+1i 0+0i
    ! 5+3i 2+2i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (5.e0,3.e0), (2,2)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(-1.e0,-1.e0), (-7,-5)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine z_usmv_2_n_am1_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_t_am1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-1
    complex*16 :: beta=0
    ! 1+1i 1+0i
    ! 0+3i 3+1i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    complex*16 :: VA(4)=(/(1.e0,1.e0), (1.e0,0.e0), (0.e0,3.e0), (3,1)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(-1.e0,-4.e0), (-4,-1)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine z_usmv_2_t_am1_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_c_am1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-1
    complex*16 :: beta=0
    ! 1+1i 2+0i
    ! 1+0i 0+1i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    complex*16 :: VA(4)=(/(1.e0,1.e0), (2.e0,0.e0), (1.e0,0.e0), (0,1)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(-2.e0,1.e0), (-2,1)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-1 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine z_usmv_2_c_am1_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_n_am3_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-3
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 2+3i 0+1i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (2.e0,3.e0), (0,1)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(0.e0,-3.e0), (-3,-12)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=n is ok"
  end subroutine z_usmv_2_n_am3_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_t_am3_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-3
    complex*16 :: beta=1
    ! 1+1i 0+0i
    ! 1+4i 2+4i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (1.e0,4.e0), (2,4)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(-3.e0,-15.e0), (-3,-12)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=t is ok"
  end subroutine z_usmv_2_t_am3_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_c_am3_bp1_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-3
    complex*16 :: beta=1
    ! 1+1i 0+2i
    ! 2+0i 0+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 1, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 2, 1/)
    complex*16 :: VA(3)=(/(1.e0,1.e0), (0.e0,2.e0), (2,0)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(-6.e0,3.e0), (3,6)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 1 incx=1 incy=1 trans=c is ok"
  end subroutine z_usmv_2_c_am3_bp1_ix1_iy1
  ! 

  subroutine z_usmv_2_n_am3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-3
    complex*16 :: beta=0
    ! 1+1i 0+3i
    ! 0+1i 1+1i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    complex*16 :: VA(4)=(/(1.e0,1.e0), (0.e0,3.e0), (0.e0,1.e0), (1,1)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(-3.e0,-12.e0), (-3,-6)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine z_usmv_2_n_am3_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_t_am3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-3
    complex*16 :: beta=0
    ! 1+1i 0+1i
    ! 0+3i 1+5i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=4
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(4)=(/1, 1, 2, 2/)
    integer(psb_ipk_) :: JA(4)=(/1, 2, 1, 2/)
    complex*16 :: VA(4)=(/(1.e0,1.e0), (0.e0,1.e0), (0.e0,3.e0), (1,5)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(-3.e0,-12.e0), (-3,-18)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine z_usmv_2_t_am3_bm0_ix1_iy1
  ! 

  subroutine z_usmv_2_c_am3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    integer(psb_ipk_) :: incy=1
    complex*16 :: alpha=-3
    complex*16 :: beta=0
    ! 1+1i 0+0i
    ! 0+0i 0+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=1
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(1)=(/1/)
    integer(psb_ipk_) :: JA(1)=(/1/)
    complex*16 :: VA(1)=(/(1,1)/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/(-3.e0,3.e0), (0,0)/)! reference cy after 
    complex*16 :: bcy(2)=(/3, 3/)! reference bcy before 
    complex*16 :: y(2)=(/3, 3/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 usmv alpha=-3 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine z_usmv_2_c_am3_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_n_ap3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=3
    complex*16 :: beta=0
    ! 1+0i 0+0i
    ! 0+2i 1+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,0.e0), (0.e0,2.e0), (1,0)/)
    complex*16 :: x(2)=(/(3.e0,0.e0), (3,6)/)! reference x 
    complex*16 :: cy(2)=(/9, 9/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine z_ussv_2_n_ap3_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_t_ap3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=3
    complex*16 :: beta=0
    ! 1+0i 0+0i
    ! 0+1i 1+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,0.e0), (0.e0,1.e0), (1,0)/)
    complex*16 :: x(2)=(/(3.e0,3.e0), (3,0)/)! reference x 
    complex*16 :: cy(2)=(/9, 9/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine z_ussv_2_t_ap3_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_c_ap3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=3
    complex*16 :: beta=0
    ! 1 0
    ! 1 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/1, 1, 1/)
    complex*16 :: x(2)=(/6, 3/)! reference x 
    complex*16 :: cy(2)=(/9, 9/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 3 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine z_ussv_2_c_ap3_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_n_ap1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=1
    complex*16 :: beta=0
    ! 1+0i 0+0i
    ! 1+5i 1+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,0.e0), (1.e0,5.e0), (1,0)/)
    complex*16 :: x(2)=(/(1.e0,0.e0), (2,5)/)! reference x 
    complex*16 :: cy(2)=(/1, 1/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine z_ussv_2_n_ap1_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_t_ap1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=1
    complex*16 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    complex*16 :: VA(2)=(/1, 1/)
    complex*16 :: x(2)=(/1, 1/)! reference x 
    complex*16 :: cy(2)=(/1, 1/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine z_ussv_2_t_ap1_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_c_ap1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=1
    complex*16 :: beta=0
    ! 1 0
    ! 2 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/1, 2, 1/)
    complex*16 :: x(2)=(/3, 1/)! reference x 
    complex*16 :: cy(2)=(/1, 1/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha= 1 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine z_ussv_2_c_ap1_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_n_am1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=-1
    complex*16 :: beta=0
    ! 1 0
    ! 2 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/1, 2, 1/)
    complex*16 :: x(2)=(/-1, -3/)! reference x 
    complex*16 :: cy(2)=(/1, 1/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine z_ussv_2_n_am1_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_t_am1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=-1
    complex*16 :: beta=0
    ! 1 0
    ! 0 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=2
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(2)=(/1, 2/)
    integer(psb_ipk_) :: JA(2)=(/1, 2/)
    complex*16 :: VA(2)=(/1, 1/)
    complex*16 :: x(2)=(/-1, -1/)! reference x 
    complex*16 :: cy(2)=(/1, 1/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine z_ussv_2_t_am1_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_c_am1_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=-1
    complex*16 :: beta=0
    ! 1 0
    ! 2 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/1, 2, 1/)
    complex*16 :: x(2)=(/-3, -1/)! reference x 
    complex*16 :: cy(2)=(/1, 1/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-1 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine z_ussv_2_c_am1_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_n_am3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='n'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=-3
    complex*16 :: beta=0
    ! 1 0
    ! 1 1

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/1, 1, 1/)
    complex*16 :: x(2)=(/-3, -6/)! reference x 
    complex*16 :: cy(2)=(/9, 9/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=n is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=n is ok"
  end subroutine z_ussv_2_n_am3_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_t_am3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='t'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=-3
    complex*16 :: beta=0
    ! 1+0i 0+0i
    ! 1+3i 1+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,0.e0), (1.e0,3.e0), (1,0)/)
    complex*16 :: x(2)=(/(-6.e0,-9.e0), (-3,0)/)! reference x 
    complex*16 :: cy(2)=(/9, 9/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=t is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=t is ok"
  end subroutine z_ussv_2_t_am3_bm0_ix1_iy1
  ! 

  subroutine z_ussv_2_c_am3_bm0_ix1_iy1(res,afmt,ctxt)
    use psb_base_mod  
    implicit none
    character(len=*) :: afmt
    type(psb_zspmat_type) :: a
    type(psb_desc_type)   :: desc_a
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iam=-1, np=-1
    integer(psb_ipk_) :: info=-1

    integer(psb_ipk_) ::res,istat=0,i
    character::transa='c'
    integer(psb_ipk_) :: incx=1
    complex*16 :: alpha=-3
    complex*16 :: beta=0
    ! 1+0i 0+0i
    ! 2+3i 1+0i

    ! declaration of VA,IA,JA 
    integer(psb_ipk_) :: nnz=3
    integer(psb_ipk_) :: m=2
    integer(psb_ipk_) :: k=2
    integer(psb_ipk_) :: IA(3)=(/1, 2, 2/)
    integer(psb_ipk_) :: JA(3)=(/1, 1, 2/)
    complex*16 :: VA(3)=(/(1.e0,0.e0), (2.e0,3.e0), (1,0)/)
    complex*16 :: x(2)=(/(-9.e0,9.e0), (-3,0)/)! reference x 
    complex*16 :: cy(2)=(/9, 9/)! reference cy after 
    complex*16 :: bcy(2)=(/0, 0/)! reference bcy before 
    complex*16 :: y(2)=(/0, 0/)! y 

    y=bcy
    res=0
    call psb_info(ctxt,iam,np)
    if(iam<0)then
      info=-1
      goto 9999
    endif
    call psb_barrier(ctxt)
    call psb_cdall(ctxt,desc_a,info,nl=m)
    if (info /= psb_success_)goto 9996
    call psb_spall(a,desc_a,info,nnz=nnz)
    if (info /= psb_success_)goto 9996
    call a%set_triangle()
    call a%set_lower()
    call a%set_unit(.false.)

    call psb_barrier(ctxt)
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
    if(res /= 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=c is not ok"
    if(res == 0)print*,"on z matrix    2 x    2 blocked   1 x   1 ussv alpha=-3 beta= 0 incx=1 incy=1 trans=c is ok"
  end subroutine z_ussv_2_c_am3_bm0_ix1_iy1
end module psb_z_mvsv_tester
