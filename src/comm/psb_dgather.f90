! File:  psb_dgather.f90
!
! Subroutine: psb_dgatherm
!   This subroutine gathers pieces of a distributed dense matrix into a local one.
!
! Parameters:
!   globx     -  real,dimension(:,:).          The local matrix into which gather the distributed pieces.
!   locx      -  real,dimension(:,:).          The local piece of the ditributed matrix to be gathered.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   iroot     -  integer.                      The process that has to own the global matrix. If -1 all
!                                              the processes will have a copy.
!   iiglobx   -  integer(optional).            The starting row of the global matrix. 
!   ijglobx   -  integer(optional).            The starting column of the global matrix. 
!   iilocx    -  integer(optional).            The starting row of the local piece of matrix. 
!   ijlocx    -  integer(optional).            The starting column of the local piece of matrix.
!   ik        -  integer(optional).            The number of columns to gather. 
!
subroutine  psb_dgatherm(globx, locx, desc_a, info, iroot,&
     & iiglobx, ijglobx, iilocx,ijlocx,ik)
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)    :: locx(:,:)
  real(kind(1.d0)), intent(out)   :: globx(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, intent(in), optional   :: iroot, iiglobx, ijglobx, iilocx, ijlocx, ik


  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, temp(2), root, iiroot, ilocx, iglobx, jlocx,&
       & jglobx, lda_locx, lda_globx, m, lock, globk, maxk, k, jlx, ilx, i, j, idx
  real(kind(1.d0))         :: locmax(2), amax
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_dgatherm'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (nprow == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= 1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
     root = iroot
     if((root.lt.-1).or.(root.gt.nprow)) then
        info=30
        int_err(1:2)=(/5,root/)
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
     end if
  else
     root = -1
  end if
  if (root==-1) then
     iiroot=0
  endif

  if (present(iiglobx)) then
     iglobx = iiglobx
  else
     iglobx = 1
  end if

  if (present(ijglobx)) then
     jglobx = ijglobx
  else
     jglobx = 1
  end if

  if (present(iilocx)) then
     ilocx = iilocx
  else
     ilocx = 1
  end if

  if (present(ijlocx)) then
     jlocx = ijlocx
  else
     jlocx = 1
  end if

  lda_globx = size(globx,1)
  lda_locx  = size(locx, 1)

  m = desc_a%matrix_data(psb_m_)
  n = desc_a%matrix_data(psb_n_)
  
  lock=size(locx,2)-jlocx+1
  globk=size(globx,2)-jglobx+1
  maxk=min(lock,globk)
  
  if(present(ik)) then
     if(ik.gt.maxk) then
        k=maxk
     else
        k=ik
     end if
  else
     k = maxk
  end if

  if (myrow == iiroot) then
     call igebs2d(icontxt, 'all', ' ', 1, 1, k, 1)
  else
     call igebr2d(icontxt, 'all', ' ', 1, 1, k, 1, iiroot, 0)
  end if

  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx,1),iglobx,jglobx,desc_a%matrix_data,info)
  call psb_chkvect(m,n,size(locx,1),ilocx,jlocx,desc_a%matrix_data,info,ilx,jlx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chk(glob)vect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((ilx.ne.1).or.(iglobx.ne.1)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if
  
  globx(:,:)=0.d0

  do j=1,k
     do i=1,desc_a%matrix_data(psb_n_row_)
        idx = desc_a%loc_to_glob(i)
        globx(idx,jglobx+j-1) = locx(i,jlx+j-1)
     end do
     ! adjust overlapped elements
     i=0
     do while (desc_a%ovrlap_elem(i).ne.-1)
        idx=desc_a%ovrlap_elem(i+psb_ovrlp_elem_)
        idx=desc_a%loc_to_glob(idx)
        globx(idx,jglobx+j-1) = globx(idx,jglobx+j-1)/desc_a%ovrlap_elem(i+psb_n_dom_ovr_)
        i=i+2
     end do
  end do

  call dgsum2d(icontxt,'a',' ',m,k,globx(1,jglobx),size(globx,1),root,mycol)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dgatherm






! Subroutine: psb_dgatherv
!   This subroutine gathers pieces of a distributed dense vector into a local one.
!
! Parameters:
!   globx     -  real,dimension(:).            The local vector into which gather the distributed pieces.
!   locx      -  real,dimension(:).            The local piece of the ditributed vector to be gathered.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   iroot     -  integer.                      The process that has to own the global vector. If -1 all
!                                              the processes will have a copy.
!   iiglobx   -  integer(optional).            The starting row of the global vector. 
!   iilocx    -  integer(optional).            The starting row of the local piece of vector. 
!
subroutine  psb_dgatherv(globx, locx, desc_a, info, iroot,&
     & iiglobx, iilocx)
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  real(kind(1.d0)), intent(in)    :: locx(:)
  real(kind(1.d0)), intent(out)   :: globx(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, intent(in), optional   :: iroot, iiglobx, iilocx


  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, n, iix, jjx, temp(2), root, iiroot, ilocx, iglobx, jlocx,&
       & jglobx, lda_locx, lda_globx, lock, maxk, globk, m, k, jlx, ilx, i, j, idx
  real(kind(1.d0))         :: locmax(2), amax
  real(kind(1.d0)),pointer :: tmpx(:)
  character(len=20)        :: name, ch_err

  name='psb_dgatherv'
  info=0
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (nprow == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= 1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
     root = iroot
     if((root.lt.-1).or.(root.gt.nprow)) then
        info=30
        int_err(1:2)=(/5,root/)
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
     end if
  else
     root = -1
  end if
  if (root==-1) then
     iiroot=0
  endif

  jglobx=1
  if (present(iiglobx)) then
     iglobx = iiglobx
  else
     iglobx = 1
  end if

  jlocx=1
  if (present(iilocx)) then
     ilocx = iilocx
  else
     ilocx = 1
  end if

  lda_globx = size(globx)
  lda_locx  = size(locx)

  m = desc_a%matrix_data(psb_m_)
  n = desc_a%matrix_data(psb_n_)
  
  k = 1

  if (myrow == iiroot) then
     call igebs2d(icontxt, 'all', ' ', 1, 1, k, 1)
  else
     call igebr2d(icontxt, 'all', ' ', 1, 1, k, 1, iiroot, 0)
  end if

  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx),iglobx,jglobx,desc_a%matrix_data,info)
  call psb_chkvect(m,n,size(locx),ilocx,jlocx,desc_a%matrix_data,info,ilx,jlx)
  if(info.ne.0) then
     info=4010
     ch_err='psb_chk(glob)vect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((ilx.ne.1).or.(iglobx.ne.1)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if
  
  globx(:)=0.d0

  do i=1,desc_a%matrix_data(psb_n_row_)
     idx = desc_a%loc_to_glob(i)
     globx(idx) = locx(i)
  end do
  ! adjust overlapped elements
  i=0
  do while (desc_a%ovrlap_elem(i).ne.-1)
     idx=desc_a%ovrlap_elem(i+psb_ovrlp_elem_)
     idx=desc_a%loc_to_glob(idx)
     globx(idx) = globx(idx)/desc_a%ovrlap_elem(i+psb_n_dom_ovr_)
     i=i+2
  end do
  
  call dgsum2d(icontxt,'a',' ',m,k,globx,size(globx),root,mycol)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dgatherv
