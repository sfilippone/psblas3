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
!
! File: psi_fnd_owner.f90
!
! Subroutine: psi_fnd_owner
!   Figure out who owns  global indices. 
! 
! Arguments: 
!    nv       - integer                 Number of indices required on  the calling
!                                       process 
!    idx(:)   - integer                 Required indices on the calling process.
!                                       Note: the indices should be unique!
!    iprc(:)  - integer(psb_ipk_), allocatable    Output: process identifiers for the corresponding
!                                       indices
!    desc_a   - type(psb_desc_type).    The communication descriptor.        
!    info     - integer.                return code.
! 
subroutine psb_indx_map_fnd_owner(idx,iprc,idxmap,info)
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_realloc_mod
  use psb_indx_map_mod, psb_protect_name => psb_indx_map_fnd_owner
#ifdef MPI_MOD
  use mpi
#endif

  implicit none 
#ifdef MPI_H
  include 'mpif.h'
#endif
  integer(psb_ipk_), intent(in) :: idx(:)
  integer(psb_ipk_), allocatable, intent(out) ::  iprc(:)
  class(psb_indx_map), intent(in) :: idxmap
  integer(psb_ipk_), intent(out) :: info


  integer(psb_ipk_), allocatable :: hsz(:),hidx(:),helem(:),hproc(:),&
       & sdsz(:),sdidx(:), rvsz(:), rvidx(:),answers(:,:),idxsrch(:,:)

  integer(psb_ipk_) :: i,n_row,n_col,err_act,ih,icomm,hsize,ip,isz,k,j,&
       & last_ih, last_j, nv
  integer(psb_ipk_) :: ictxt,np,me
  logical, parameter  :: gettime=.false.
  real(psb_dpk_)      :: t0, t1, t2, t3, t4, tamx, tidx
  character(len=20)   :: name

  info = psb_success_
  name = 'psb_indx_map_fnd_owner'
  call psb_erractionsave(err_act)

  ictxt   = idxmap%get_ctxt()
  icomm   = idxmap%get_mpic()
  n_row   = idxmap%get_lr()
  n_col   = idxmap%get_lc()

  
  call psb_info(ictxt, me, np)
  
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.(idxmap%is_valid())) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='invalid idxmap')
    goto 9999      
  end if

  if (gettime) then 
    t0 = psb_wtime()
  end if

  nv = size(idx)
  !
  ! The basic idea is very simple. 
  ! First we collect (to all) all the requests. 
  Allocate(hidx(np+1),hsz(np),&
       & sdsz(0:np-1),sdidx(0:np-1),&
       & rvsz(0:np-1),rvidx(0:np-1),&
       & stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate') 
    goto 9999      
  end if

  hsz       = 0
  hsz(me+1) = nv
  call psb_amx(ictxt,hsz)
  hidx(1)   = 0
  do i=1, np
    hidx(i+1) = hidx(i) + hsz(i)
  end do
  hsize = hidx(np+1)
  Allocate(helem(hsize),hproc(hsize),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if

  if (gettime) then 
    t3 = psb_wtime()
  end if

  call mpi_allgatherv(idx,hsz(me+1),psb_mpi_ipk_integer,&
       & hproc,hsz,hidx,psb_mpi_ipk_integer,&
       & icomm,info)
  if (gettime) then 
    tamx = psb_wtime() - t3
  end if

  ! Second, we figure out locally whether we own the indices (whoever is 
  ! asking for them). 
  if (gettime) then 
    t3 = psb_wtime()
  end if

  call idxmap%g2l(hproc(1:hsize),helem(1:hsize),info,owned=.true.)
  if (gettime) then 
    tidx = psb_wtime()-t3
  end if
  if (info == psb_err_iarray_outside_bounds_) info = psb_success_
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_idx_cnv')
    goto 9999      
  end if

  ! Third: we build the answers for those indices we own,
  ! with a section for each process asking. 
  hidx = hidx +1 
  j    = 0
  do ip = 0, np-1
    sdidx(ip) = j
    sdsz(ip)  = 0
    do i=hidx(ip+1), hidx(ip+1+1)-1
      if ((0 < helem(i)).and. (helem(i) <= n_row)) then 
        j        = j + 1 
        hproc(j) = hproc(i)
        sdsz(ip) = sdsz(ip) + 1
      end if
    end do
  end do

  if (gettime) then 
    t3 = psb_wtime()
  end if

  ! Collect all the answers with alltoallv (need sizes) 
  call mpi_alltoall(sdsz,1,psb_mpi_ipk_integer,rvsz,1,psb_mpi_def_integer,icomm,info)

  isz = sum(rvsz) 

  allocate(answers(isz,2),idxsrch(nv,2),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if
  j = 0
  do ip=0, np-1
    rvidx(ip) = j
    j         = j + rvsz(ip)
  end do
  call mpi_alltoallv(hproc,sdsz,sdidx,psb_mpi_ipk_integer,&
       & answers(:,1),rvsz,rvidx,psb_mpi_ipk_integer,&
       & icomm,info)
  if (gettime) then 
    tamx = psb_wtime() - t3 + tamx
  end if
  j = 1
  do ip = 0,np-1
    do k=1,rvsz(ip)
      answers(j,2) = ip
      j = j + 1 
    end do
  end do
  ! Sort the answers and the requests, so we can
  ! match them efficiently
  call psb_msort(answers(:,1),ix=answers(:,2),&
       & flag=psb_sort_keep_idx_)
  idxsrch(1:nv,1) = idx(1:nv)
  call psb_msort(idxsrch(1:nv,1),ix=idxsrch(1:nv,2))

  ! Now  extract the answers for our local query
  call psb_realloc(nv,iprc,info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='psb_realloc')
    goto 9999      
  end if
  last_ih = -1 
  last_j  = -1
  j = 1
  do i=1, nv 
    ih = idxsrch(i,1)
    if (ih == last_ih) then 
      iprc(idxsrch(i,2)) = answers(last_j,2)
    else

      do
        if (j > size(answers,1)) then 
          ! Last resort attempt.
          j = psb_ibsrch(ih,size(answers,1,kind=psb_ipk_),answers(:,1))
          if (j == -1) then 
            write(psb_err_unit,*) me,'psi_fnd_owner: searching for ',ih, &
                 & 'not found : ',size(answers,1),':',answers(:,1)
            info = psb_err_internal_error_
            call psb_errpush(psb_err_internal_error_,name,a_err='out bounds srch ih') 
            goto 9999      
          end if
        end if
        if (answers(j,1) == ih) exit
        if (answers(j,1) > ih) then 
          k = j 
          j = psb_ibsrch(ih,k,answers(1:k,1))
          if (j == -1) then 
            write(psb_err_unit,*) me,'psi_fnd_owner: searching for ',ih, &
                 & 'not found : ',size(answers,1),':',answers(:,1)
            info = psb_err_internal_error_
            call psb_errpush(psb_err_internal_error_,name,a_err='out bounds srch ih') 
            goto 9999      
          end if
        end if

        j = j + 1 
      end do
      ! Note that the answers here are given in order
      ! of sending process, so we are implicitly getting
      ! the max process index in case of overlap. 
      last_ih = ih 
      do 
        last_j = j 
        iprc(idxsrch(i,2)) = answers(j,2)
        j = j + 1 
        if (j > size(answers,1)) exit
        if (answers(j,1) /= ih) exit
      end do
    end if
  end do

  if (gettime) then 
    call psb_barrier(ictxt)
    t1 = psb_wtime()
    t1 = t1 -t0 - tamx - tidx   
    call psb_amx(ictxt,tamx)
    call psb_amx(ictxt,tidx)
    call psb_amx(ictxt,t1)
    if (me == psb_root_) then 
      write(psb_out_unit,'(" fnd_owner  idx time  : ",es10.4)') tidx
      write(psb_out_unit,'(" fnd_owner  amx time  : ",es10.4)') tamx
      write(psb_out_unit,'(" fnd_owner remainedr  : ",es10.4)') t1 
    endif
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_indx_map_fnd_owner
