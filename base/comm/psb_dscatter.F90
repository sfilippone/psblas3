!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Arguments:
!   globx     -  real,dimension(:,:).          The global matrix to scatter.
!   locx      -  real,dimension(:,:).          The local piece of the ditributed matrix.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Error code.
!   iroot     -  integer(optional).            The process that owns the global matrix. If -1 all
!                                              the processes have a copy. Default -1.
!
subroutine  psb_dscatterm(globx, locx, desc_a, info, iroot)

  use psb_base_mod, psb_protect_name => psb_dscatterm
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  real(psb_dpk_), intent(out)    :: locx(:,:)
  real(psb_dpk_), intent(in)     :: globx(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  integer, intent(in), optional    :: iroot


  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, m, n, i, j, idx, nrow, iiroot, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, lock, globk, icomm, k, maxk, root, ilx,&
       & jlx, myrank, rootrank, c, pos
  real(psb_dpk_), allocatable  :: scatterv(:)
  integer, allocatable         :: displ(:), l_t_g_all(:), all_dim(:), ltg(:)
  character(len=20)        :: name, ch_err

  name='psb_scatterm'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
    root = iroot
    if((root < -1).or.(root > np)) then
      info=psb_err_input_value_invalid_i_
      int_err(1:2)=(/5,root/)
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
  else
    root = -1
  end if
  if (root == -1) then
    iiroot = psb_root_
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

  call psb_bcast(ictxt,k,root=iiroot)

  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx),iglobx,jglobx,desc_a,info)
  if (info == psb_success_) call psb_chkvect(m,n,size(locx),ilocx,jlocx,desc_a,info,ilx,jlx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chk(glob)vect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((ilx /= 1).or.(iglobx /= 1)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  nrow=psb_cd_get_local_rows(desc_a)

  if ((root == -1).or.(np == 1)) then
    ! extract my chunk
    do j=1,k
      do i=1, nrow
        call psb_loc_to_glob(i,idx,desc_a,info)
        locx(i,jlocx+j-1)=globx(idx,jglobx+j-1)
      end do
    end do
  else
    call psb_get_rank(rootrank,ictxt,root)

    ! root has to gather size information
    allocate(displ(np),all_dim(np),ltg(nrow),stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    do i=1, nrow
      ltg(i) = i
    end do
    call psb_loc_to_glob(ltg(1:nrow),desc_a,info) 

    call mpi_gather(nrow,1,mpi_integer,all_dim,&
         & 1,mpi_integer,rootrank,icomm,info)

    if (me == root) then
      displ(1)=0
      do i=2,np
        displ(i)=displ(i-1)+all_dim(i-1)
      end do

    ! root has to gather loc_glob from each process
      allocate(l_t_g_all(sum(all_dim)),scatterv(sum(all_dim)),stat=info)
      if(info /= psb_success_) then
        info=psb_err_from_subroutine_
        ch_err='Allocate'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

    end if

    call mpi_gatherv(ltg,nrow,&
         & mpi_integer,l_t_g_all,all_dim,&
         & displ,mpi_integer,rootrank,icomm,info)


    do c=1, k
      ! prepare vector to scatter
      if(me == root) then
        do i=1,np
          pos=displ(i)
          do j=1, all_dim(i)
            idx=l_t_g_all(pos+j)
            scatterv(pos+j)=globx(idx,jglobx+c-1)
          end do
        end do
      end if

      ! scatter !!!
      call mpi_scatterv(scatterv,all_dim,displ,&
           & mpi_double_precision,locx(1,jlocx+c-1),nrow,&
           & mpi_double_precision,rootrank,icomm,info)

    end do

    if (me == root) deallocate(all_dim, l_t_g_all, displ, scatterv)
  end if

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
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Arguments:
!   globx     -  real,dimension(:).            The global vector to scatter.
!   locx      -  real,dimension(:).            The local piece of the ditributed vector.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Error code.
!   iroot     -  integer(optional).            The process that owns the global vector. If -1 all
!                                              the processes have a copy.
!
subroutine  psb_dscatterv(globx, locx, desc_a, info, iroot)
  use psb_base_mod, psb_protect_name => psb_dscatterv

#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  real(psb_dpk_), intent(out)    :: locx(:)
  real(psb_dpk_), intent(in)     :: globx(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  integer, intent(in), optional    :: iroot


  ! locals
  integer                  :: int_err(5), ictxt, np, me, &
       & err_act, m, n, i, j, idx, nrow, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, root, k, icomm, myrank,&
       & rootrank, pos, ilx, jlx
  real(psb_dpk_), allocatable  :: scatterv(:)
  integer, allocatable         :: displ(:), l_t_g_all(:), all_dim(:), ltg(:)
  character(len=20)        :: name, ch_err
  integer                  :: debug_level, debug_unit

  name='psb_scatterv'
  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  ictxt=psb_cd_get_context(desc_a)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()


  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
     root = iroot
     if((root < -1).or.(root > np)) then
        info=psb_err_input_value_invalid_i_
        int_err(1:2)=(/5,root/)
        call psb_errpush(info,name,i_err=int_err)
        goto 9999
     end if
  else
     root = -1
  end if
  
  call psb_get_mpicomm(ictxt,icomm)
  call psb_get_rank(myrank,ictxt,me)

  iglobx = 1
  jglobx = 1
  ilocx = 1
  jlocx = 1
  lda_globx = size(globx)
  lda_locx  = size(locx)

  m = psb_cd_get_global_rows(desc_a)
  n = psb_cd_get_global_cols(desc_a)

  k = 1
  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx),iglobx,jglobx,desc_a,info)
  if (info == psb_success_) &
       & call psb_chkvect(m,n,size(locx),ilocx,jlocx,desc_a,info,ilx,jlx)
  if(info /= psb_success_) then
     info=psb_err_from_subroutine_
     ch_err='psb_chk(glob)vect'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if ((ilx /= 1).or.(iglobx /= 1)) then
     info=psb_err_ix_n1_iy_n1_unsupported_
     call psb_errpush(info,name)
     goto 9999
  end if

  nrow = psb_cd_get_local_rows(desc_a)

  if ((root == -1).or.(np == 1)) then
    ! extract my chunk
    do i=1, nrow
      call psb_loc_to_glob(i,idx,desc_a,info)
      locx(i)=globx(idx)
    end do
  else
    call psb_get_rank(rootrank,ictxt,root)

    ! root has to gather size information
    allocate(displ(np),all_dim(np),ltg(nrow))
    do i=1, nrow
      ltg(i) = i
    end do
    call psb_loc_to_glob(ltg(1:nrow),desc_a,info) 


    call mpi_gather(nrow,1,mpi_integer,all_dim,&
         & 1,mpi_integer,rootrank,icomm,info)

    if(me == root) then
      displ(1)=0
      do i=2,np
        displ(i)=displ(i-1) + all_dim(i-1)
      end do
      if (debug_level >= psb_debug_inner_) then 
        write(debug_unit,*) me,' ',trim(name),' displ:',displ(1:np), &
             &' dim',all_dim(1:np), sum(all_dim)
      endif
      
      ! root has to gather loc_glob from each process
      allocate(l_t_g_all(sum(all_dim)),scatterv(sum(all_dim)))
    end if

    call mpi_gatherv(ltg,nrow,&
         & mpi_integer,l_t_g_all,all_dim,&
         & displ,mpi_integer,rootrank,icomm,info)

    ! prepare vector to scatter
    if (me == root) then
      do i=1,np
        pos=displ(i)
        do j=1, all_dim(i)
          idx=l_t_g_all(pos+j)
          scatterv(pos+j)=globx(idx)

        end do
      end do
    end if

    call mpi_scatterv(scatterv,all_dim,displ,&
         & mpi_double_precision,locx,nrow,&
         & mpi_double_precision,rootrank,icomm,info)

    if (me == root) deallocate(all_dim, l_t_g_all, displ, scatterv)
  end if

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
