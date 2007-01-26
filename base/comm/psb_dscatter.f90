!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
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
!   info      -  integer.                      Error code.
!   iroot     -  integer(optional).            The process that owns the global matrix. If -1 all
!                                              the processes have a copy. Default -1.
!
subroutine  psb_dscatterm(globx, locx, desc_a, info, iroot)

  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use mpi
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(out)    :: locx(:,:)
  real(kind(1.d0)), intent(in)     :: globx(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  integer, intent(in), optional    :: iroot


  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, m, n, i, j, idx, nrow, iiroot, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, lock, globk, icomm, k, maxk, root, ilx,&
       & jlx, myrank, rootrank, c, pos
  real(kind(1.d0)), allocatable  :: scatterv(:)
  integer, allocatable           :: displ(:), l_t_g_all(:), all_dim(:)
  character(len=20)        :: name, ch_err

  name='psb_scatterm'
  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
     root = iroot
     if((root.lt.-1).or.(root.gt.np)) then
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
  
  iglobx = 1
  jglobx = 1
  ilocx = 1
  jlocx = 1
  lda_globx = size(globx,1)
  lda_locx  = size(locx, 1)

  m = psb_cd_get_global_rows(desc_a)
  n = psb_cd_get_global_cols(desc_a)
  
  lock=size(locx,2)-jlocx+1
  globk=size(globx,2)-jglobx+1
  maxk=min(lock,globk)
  k = maxk
  call psb_get_mpicomm(ictxt,icomm)
  call psb_get_rank(myrank,ictxt,me)


  lda_globx = size(globx)
  lda_locx  = size(locx)

  m = psb_cd_get_global_rows(desc_a)
  n = psb_cd_get_global_cols(desc_a)
  
  if (me == iiroot) then
     call igebs2d(ictxt, 'all', ' ', 1, 1, k, 1)
  else
     call igebr2d(ictxt, 'all', ' ', 1, 1, k, 1, iiroot, 0)
  end if

  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx),iglobx,jglobx,desc_a,info)
  if (info == 0) call psb_chkvect(m,n,size(locx),ilocx,jlocx,desc_a,info,ilx,jlx)
  if(info /= 0) then
     info=4010
     ch_err='psb_chk(glob)vect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((ilx /= 1).or.(iglobx /= 1)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  nrow=psb_cd_get_local_rows(desc_a)

  if(root == -1) then
     ! extract my chunk
     do j=1,k
        do i=1, nrow
           idx=desc_a%loc_to_glob(i)
           locx(i,jlocx+j-1)=globx(idx,jglobx+j-1)
        end do
     end do
  else
    call psb_get_rank(rootrank,ictxt,root)
  end if

  ! root has to gather size information
  allocate(displ(np),all_dim(np),stat=info)
  if(info /= 0) then
     info=4010
     ch_err='Allocate'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  call mpi_gather(nrow,1,mpi_integer,all_dim,&
       & np,mpi_integer,rootrank,icomm,info)
  
  displ(1)=1
  displ(2:)=all_dim(1:np-1)+1

  ! root has to gather loc_glob from each process
  if(me == root) then
     allocate(l_t_g_all(sum(all_dim)),scatterv(sum(all_dim)),stat=info)
     if(info /= 0) then
       info=4010
       ch_err='Allocate'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
     end if

  end if
     
  call mpi_gatherv(desc_a%loc_to_glob,nrow,&
       & mpi_integer,l_t_g_all,all_dim,&
       & displ,mpi_integer,rootrank,icomm,info)

  
  do c=1, k
     ! prepare vector to scatter
     if(me == root) then
        do i=1,np
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

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return

end subroutine psb_dscatterm




!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  

! Subroutine: psb_dscatterv
!   This subroutine scatters a global vector locally owned by one process
!   into pieces that are local to alle the processes.
!
! Parameters:
!   globx     -  real,dimension(:).            The global vector to scatter.
!   locx      -  real,dimension(:).            The local piece of the ditributed vector.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Error code.
!   iroot     -  integer(optional).            The process that owns the global vector. If -1 all
!                                              the processes have a copy.
!
subroutine  psb_dscatterv(globx, locx, desc_a, info, iroot)
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use mpi
  use psb_penv_mod
  implicit none

  real(kind(1.d0)), intent(out)    :: locx(:)
  real(kind(1.d0)), intent(in)     :: globx(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  integer, intent(in), optional    :: iroot


  ! locals
  integer                  :: int_err(5), ictxt, np, me, &
       & err_act, m, n, i, j, idx, nrow, iiroot, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, root, k, icomm, myrank,&
       & rootrank, pos, ilx, jlx
  real(kind(1.d0)), allocatable  :: scatterv(:)
  integer, allocatable           :: displ(:), l_t_g_all(:), all_dim(:)
  character(len=20)        :: name, ch_err

  name='psb_scatterv'
  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
     root = iroot
     if((root.lt.-1).or.(root.gt.np)) then
        info=30
        int_err(1:2)=(/5,root/)
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
     end if
  else
     root = -1
  end if
  
  call psb_get_mpicomm(ictxt,icomm)
  call psb_get_rank(myrank,ictxt,me)

  lda_globx = size(globx)
  lda_locx  = size(locx)

  m = psb_cd_get_global_rows(desc_a)
  n = psb_cd_get_global_cols(desc_a)
  
  k = 1

  if (me == iiroot) then
     call igebs2d(ictxt, 'all', ' ', 1, 1, k, 1)
  else
     call igebr2d(ictxt, 'all', ' ', 1, 1, k, 1, iiroot, 0)
  end if

  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx),iglobx,jglobx,desc_a,info)
  if (info == 0) &
       & call psb_chkvect(m,n,size(locx),ilocx,jlocx,desc_a,info,ilx,jlx)
  if(info /= 0) then
     info=4010
     ch_err='psb_chk(glob)vect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((ilx /= 1).or.(iglobx /= 1)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if

  nrow=psb_cd_get_local_rows(desc_a)

  if(root == -1) then
     ! extract my chunk
     do i=1, nrow
        idx=desc_a%loc_to_glob(i)
        locx(i)=globx(idx)
     end do
  else
    call psb_get_rank(rootrank,ictxt,root)
  end if

  ! root has to gather size information
  allocate(displ(np),all_dim(np))
  call mpi_gather(nrow,1,mpi_integer,all_dim,&
       & np,mpi_integer,rootrank,icomm,info)
  
  displ(1)=1
  displ(2:)=all_dim(1:np-1)+1

  ! root has to gather loc_glob from each process
  if(me == root) then
     allocate(l_t_g_all(sum(all_dim)),scatterv(sum(all_dim)))
  end if
     
  call mpi_gatherv(desc_a%loc_to_glob,nrow,&
       & mpi_integer,l_t_g_all,all_dim,&
       & displ,mpi_integer,rootrank,icomm,info)

  ! prepare vector to scatter
  if(me == root) then
     do i=1,np
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

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return

end subroutine psb_dscatterv
