! File:  psb_dscatter.f90
!
! Subroutine: psb_dscatterm
!   This subroutine scatters a global matrix locally owned by one process
!   into pieces that are local to alle the processes.
!
! Parameters:
!   globx     -  real,dimension(:,:).          The global matrix to scatter.
!   locx      -  real,dimension(:,:).          The local piece of the ditributed matrix.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   iroot     -  integer(optional).            The process that owns the global matrix. If -1 all
!                                              the processes have a copy.
!   iiglobx   -  integer(optional).            The starting row of the global matrix. 
!   ijglobx   -  integer(optional).            The starting column of the global matrix. 
!   iilocx    -  integer(optional).            The starting row of the local piece of matrix. 
!   ijlocx    -  integer(optional).            The starting column of the local piece of matrix.
!   ik        -  integer(optional).            The number of columns to gather. 
!
subroutine  psb_dscatterm(globx, locx, desc_a, info, iroot,&
     & iiglobx, ijglobx, iilocx,ijlocx,ik)

  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use mpi
  implicit none

  real(kind(1.d0)), intent(out)    :: locx(:,:)
  real(kind(1.d0)), intent(in)     :: globx(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  integer, intent(in), optional    :: iroot,iiglobx,&
       & ijglobx,iilocx,ijlocx,ik


  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, m, n, iix, jjx, temp(2), i, j, idx, nrow, iiroot, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, lock, globk, icomm, k, maxk, root, ilx,&
       & jlx, myrank, rootrank, c, pos
  real(kind(1.d0))         :: locmax(2), amax
  real(kind(1.d0)),pointer :: scatterv(:)
  integer, pointer         :: displ(:), l_t_g_all(:), all_dim(:)
  integer                  :: blacs_pnum
  character(len=20)        :: name, ch_err

  name='psb_scatterm'
  if(psb_get_errstatus().ne.0) return 
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

  call blacs_get(icontxt,10,icomm)
  myrank = blacs_pnum(icontxt,myrow,mycol)

  lda_globx = size(globx)
  lda_locx  = size(locx)

  m = desc_a%matrix_data(psb_m_)
  n = desc_a%matrix_data(psb_n_)
  
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

  nrow=desc_a%matrix_data(psb_n_row_)

  if(root.eq.-1) then
     ! extract my chunk
     do j=1,k
        do i=1, nrow
           idx=desc_a%loc_to_glob(i)
           locx(i,jlocx+j-1)=globx(idx,jglobx+j-1)
        end do
     end do
  else
     rootrank = blacs_pnum(icontxt,root,mycol)
  end if

  ! root has to gather size information
  allocate(displ(nprow),all_dim(nprow))
  call mpi_gather(nrow,1,mpi_integer,all_dim,&
       & nprow,mpi_integer,rootrank,icomm,info)
  
  displ(1)=1
  displ(2:)=all_dim(1:nprow-1)+1

  ! root has to gather loc_glob from each process
  if(myrow.eq.root) then
     allocate(l_t_g_all(sum(all_dim)),scatterv(sum(all_dim)))
  end if
     
  call mpi_gatherv(desc_a%loc_to_glob,nrow,&
       & mpi_integer,l_t_g_all,all_dim,&
       & displ,mpi_integer,rootrank,icomm,info)

  
  do c=1, k
     ! prepare vector to scatter
     if(myrow.eq.root) then
        do i=1,nprow
           pos=displ(i)
           do j=1, all_dim(i)
              idx=l_t_g_all(pos+j-1)
              scatterv(pos+j-1)=globx(idx,jglobx+c-1)
           end do
        end do
     end if
     
     ! scatter !!!
     call mpi_scatterv(scatterv,all_dim,displ,&
          & mpi_double_precision,locx(1,jlocx+c-1),nrow,&
          & mpi_double_precision,rootrank,icomm,info)

  end do
  
  deallocate(all_dim, l_t_g_all, displ, scatterv)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dscatterm





! Subroutine: psb_dscatterv
!   This subroutine scatters a global vector locally owned by one process
!   into pieces that are local to alle the processes.
!
! Parameters:
!   globx     -  real,dimension(:).            The global vector to scatter.
!   locx      -  real,dimension(:).            The local piece of the ditributed vector.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   iroot     -  integer(optional).            The process that owns the global vector. If -1 all
!                                              the processes have a copy.
!
subroutine  psb_dscatterv(globx, locx, desc_a, info, iroot)
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use mpi
  implicit none

  real(kind(1.d0)), intent(out)    :: locx(:)
  real(kind(1.d0)), intent(in)     :: globx(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  integer, intent(in), optional    :: iroot


  ! locals
  integer                  :: int_err(5), icontxt, nprow, npcol, myrow, mycol,&
       & err_act, m, n, iix, jjx, temp(2), i, j, idx, nrow, iiroot, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, lock, globk, root, k, maxk, icomm, myrank,&
       & rootrank, c, pos, ilx, jlx
  real(kind(1.d0))         :: locmax(2), amax
  real(kind(1.d0)),pointer :: scatterv(:)
  integer, pointer         :: displ(:), l_t_g_all(:), all_dim(:)
  integer                  :: blacs_pnum
  character(len=20)        :: name, ch_err

  name='psb_scatterv'
  if(psb_get_errstatus().ne.0) return 
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
  
  call blacs_get(icontxt,10,icomm)
  myrank = blacs_pnum(icontxt,myrow,mycol)

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

  nrow=desc_a%matrix_data(psb_n_row_)

  if(root.eq.-1) then
     ! extract my chunk
     do i=1, nrow
        idx=desc_a%loc_to_glob(i)
        locx(i)=globx(idx)
     end do
  else
     rootrank = blacs_pnum(icontxt,root,mycol)
  end if

  ! root has to gather size information
  allocate(displ(nprow),all_dim(nprow))
  call mpi_gather(nrow,1,mpi_integer,all_dim,&
       & nprow,mpi_integer,rootrank,icomm,info)
  
  displ(1)=1
  displ(2:)=all_dim(1:nprow-1)+1

  ! root has to gather loc_glob from each process
  if(myrow.eq.root) then
     allocate(l_t_g_all(sum(all_dim)),scatterv(sum(all_dim)))
  end if
     
  call mpi_gatherv(desc_a%loc_to_glob,nrow,&
       & mpi_integer,l_t_g_all,all_dim,&
       & displ,mpi_integer,rootrank,icomm,info)

  ! prepare vector to scatter
  if(myrow.eq.root) then
     do i=1,nprow
        pos=displ(i)
        do j=1, all_dim(i)
           idx=l_t_g_all(pos+j-1)
           scatterv(pos+j-1)=globx(idx)
        end do
     end do
  end if

  call mpi_scatterv(scatterv,all_dim,displ,&
       & mpi_double_precision,locx,nrow,&
       & mpi_double_precision,rootrank,icomm,info)

  deallocate(all_dim, l_t_g_all, displ, scatterv)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dscatterv
