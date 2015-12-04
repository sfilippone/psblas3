!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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
! File:  psb_sscatter.f90
!
! Subroutine: psb_sscatterm
!   This subroutine scatters a global matrix locally owned by one process
!   into pieces that are local to alle the processes.
!
! Arguments:
!   globx     -  real,dimension(:,:).       The global matrix to scatter.
!   locx      -  real,dimension(:,:).       The local piece of the distributed matrix.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Error code.
!   iroot     -  integer(optional).            The process that owns the global matrix. 
!                                              If -1 all the processes have a copy. 
!                                              Default -1
subroutine  psb_sscatterm(globx, locx, desc_a, info, iroot)

  use psb_base_mod, psb_protect_name => psb_sscatterm
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  real(psb_spk_), intent(out)    :: locx(:,:)
  real(psb_spk_), intent(in)     :: globx(:,:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_), intent(in), optional    :: iroot


  ! locals
  integer(psb_mpik_) :: ictxt, np, me, root, iiroot, icomm, myrank, rootrank
  integer(psb_ipk_) :: ierr(5), err_act, m, n, i, j, idx, nrow, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, lock, globk, k, maxk, ilx,&
       & jlx, c, pos
  real(psb_spk_),allocatable  :: scatterv(:)
  integer(psb_ipk_), allocatable           :: displ(:), l_t_g_all(:), all_dim(:), ltg(:)
  character(len=20)        :: name, ch_err

  name='psb_scatterm'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

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
      ierr(1)=5; ierr(2)=root
      call psb_errpush(info,name,i_err=ierr)
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

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()

  lock=size(locx,2)-jlocx+1
  globk=size(globx,2)-jglobx+1
  maxk=min(lock,globk)
  k = maxk
  call psb_get_mpicomm(ictxt,icomm)
  call psb_get_rank(myrank,ictxt,me)


  lda_globx = size(globx)
  lda_locx  = size(locx)

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()

  call psb_bcast(ictxt,k,root=iiroot)

  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,lda_globx,iglobx,jglobx,desc_a,info)
  if (info == psb_success_) call psb_chkvect(m,n,lda_locx,ilocx,jlocx,desc_a,info,ilx,jlx)
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

  nrow=desc_a%get_local_rows()

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

    call mpi_gather(nrow,1,psb_mpi_ipk_integer,all_dim,&
         & 1,psb_mpi_ipk_integer,rootrank,icomm,info)

    if (me == root) then
      displ(1)=0
      do i=2,np
        displ(i)=displ(i-1)+all_dim(i-1)
      end do

      ! root has to gather loc_glob from each process
      allocate(l_t_g_all(sum(all_dim)),scatterv(sum(all_dim)),stat=info)
    else
      !
      ! This is to keep debugging compilers from being upset by 
      ! calling an external MPI function with an unallocated array;
      ! the Fortran side would complain even if the MPI side does
      ! not use the unallocated stuff. 
      ! 
      allocate(l_t_g_all(1),scatterv(1),stat=info)
    end if
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call mpi_gatherv(ltg,nrow,&
         & psb_mpi_ipk_integer,l_t_g_all,all_dim,&
         & displ,psb_mpi_ipk_integer,rootrank,icomm,info)


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
           & psb_mpi_r_spk_,locx(1,jlocx+c-1),nrow,&
           & psb_mpi_r_spk_,rootrank,icomm,info)

    end do

    deallocate(all_dim, l_t_g_all, displ, ltg, scatterv,stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='deallocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  end if

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

    return

end subroutine psb_sscatterm




!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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

! Subroutine: psb_sscatterv
!   This subroutine scatters a global vector locally owned by one process
!   into pieces that are local to alle the processes.
!
! Arguments:
!   globx     -  real,dimension(:).         The global vector to scatter.
!   locx      -  real,dimension(:).         The local piece of the ditributed vector.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Return code
!   iroot     -  integer(optional).            The process that owns the global vector. If -1 all
!                                              the processes have a copy.
!
subroutine  psb_sscatterv(globx, locx, desc_a, info, iroot)
  use psb_base_mod, psb_protect_name => psb_sscatterv
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif

  real(psb_spk_), intent(out)    :: locx(:)
  real(psb_spk_), intent(in)     :: globx(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_), intent(in), optional    :: iroot


  ! locals
  integer(psb_mpik_) :: ictxt, np, me, root, iiroot, icomm, myrank, rootrank
  integer(psb_ipk_) :: ierr(5), err_act, m, n, i, j, idx, nrow, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, k, pos, ilx, jlx
  real(psb_spk_), allocatable  :: scatterv(:)
  integer(psb_ipk_), allocatable            :: displ(:), l_t_g_all(:), all_dim(:), ltg(:)
  character(len=20)        :: name, ch_err
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_scatterv'
  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  ictxt=desc_a%get_context()
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
        ierr(1) = 5; ierr(2)=root
        call psb_errpush(info,name,i_err=ierr)
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

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()

  k = 1
  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,lda_globx,iglobx,jglobx,desc_a,info)
  if (info == psb_success_) &
       & call psb_chkvect(m,n,lda_locx,ilocx,jlocx,desc_a,info,ilx,jlx)
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

  nrow = desc_a%get_local_rows()

  if ((root == -1).or.(np == 1)) then
    ! extract my chunk
    do i=1, nrow
      call psb_loc_to_glob(i,idx,desc_a,info)
      locx(i)=globx(idx)
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


    call mpi_gather(nrow,1,psb_mpi_ipk_integer,all_dim,&
         & 1,psb_mpi_ipk_integer,rootrank,icomm,info)

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
      allocate(l_t_g_all(sum(all_dim)),scatterv(sum(all_dim)),stat=info)
      
    else
      !
      ! This is to keep debugging compilers from being upset by 
      ! calling an external MPI function with an unallocated array;
      ! the Fortran side would complain even if the MPI side does
      ! not use the unallocated stuff. 
      ! 
      allocate(l_t_g_all(1),scatterv(1),stat=info)
    end if
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call mpi_gatherv(ltg,nrow,&
         & psb_mpi_ipk_integer,l_t_g_all,all_dim,&
         & displ,psb_mpi_ipk_integer,rootrank,icomm,info)

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
         & psb_mpi_r_spk_,locx,nrow,&
         & psb_mpi_r_spk_,rootrank,icomm,info)

    deallocate(all_dim, l_t_g_all, displ, ltg, scatterv,stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='deallocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  end if

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

    return

end subroutine psb_sscatterv

!!$ 
!!$              Parallel Sparse BLAS  version 3.1
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012, 2013
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

! Subroutine: psb_sscatterv
!   This subroutine scatters a global vector locally owned by one process
!   into pieces that are local to alle the processes.
!
! Arguments:
!   globx     -  real,dimension(:).         The global vector to scatter.
!   locx      -  real,dimension(:).         The local piece of the ditributed vector.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Return code
!   iroot     -  integer(optional).            The process that owns the global vector. If -1 all
!                                              the processes have a copy.
!
subroutine  psb_sscatter_vect(globx, locx, desc_a, info, iroot, mold)
  use psb_base_mod, psb_protect_name => psb_sscatter_vect
  implicit none
  type(psb_s_vect_type), intent(inout) :: locx
  real(psb_spk_), intent(in)     :: globx(:)
  type(psb_desc_type), intent(in)  :: desc_a
  integer(psb_ipk_), intent(out)             :: info
  integer(psb_ipk_), intent(in), optional    :: iroot
  class(psb_s_base_vect_type), intent(in), optional :: mold
  
  ! locals
  integer(psb_mpik_) :: ictxt, np, me, root, iiroot, icomm, myrank, rootrank
  integer(psb_ipk_) :: ierr(5), err_act, m, n, i, j, idx, nrow, iglobx, jglobx,&
       & ilocx, jlocx, lda_locx, lda_globx, k, pos, ilx, jlx
  real(psb_spk_), allocatable  :: vlocx(:)
  character(len=20)        :: name, ch_err
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_scatter_vect'
  if (psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  ictxt=desc_a%get_context()
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()


  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif
  
  call psb_scatter(globx, vlocx, desc_a, info, iroot)
  if (info /= psb_success_) then
    info = psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_scatterv')
    goto 9999
  endif
  
  call locx%bld(vlocx,mold)
  
  call psb_erractionrestore(err_act)
  return  
  
9999 call psb_error_handler(ictxt,err_act)
  
  return
  
end subroutine psb_sscatter_vect
