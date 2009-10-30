!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
program sf_sample
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use getp
  implicit none

  ! input parameters
  character(len=40) :: kmethd, ptype, mtrx_file, rhs_file

  ! sparse matrices
  type(psb_s_sparse_mat) :: a, aux_a

  ! preconditioner data
  type(psb_sprec_type)  :: prec

  ! dense matrices
  real(psb_spk_), allocatable, target ::  aux_b(:,:), d(:)
  real(psb_spk_), allocatable , save  :: b_col(:), x_col(:), r_col(:), &
       & x_col_glob(:), r_col_glob(:)
  real(psb_spk_), pointer  :: b_col_glob(:)

  ! communications data structure
  type(psb_desc_type):: desc_a

  integer            :: ictxt, iam, np

  ! solver paramters
  integer            :: iter, itmax, ierr, itrace, ircode, ipart,&
       & methd, istopc, irst
  integer(psb_long_int_k_) :: amatsize, precsize, descsize
  real(psb_spk_)   :: err, eps,cond

  character(len=5)   :: afmt
  character(len=20)  :: name
  character(len=2)   :: filefmt
  integer, parameter :: iunit=12
  integer   :: iparm(20)

  ! other variables
  integer            :: i,info,j,m_problem
  integer            :: internal, m,ii,nnzero
  real(psb_spk_) :: t1, t2, tprec, r_amax, b_amax,&
       &scale,resmx,resmxp
  integer :: nrhs, nrow, n_row, dim, nv, ne
  integer, allocatable :: ivg(:), ipv(:)


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='df_sample'
  if(psb_get_errstatus() /= 0) goto 9999
  info=0
  call psb_set_errverbosity(2)
  !
  !  get parameters
  !
  call get_parms(ictxt,mtrx_file,rhs_file,filefmt,kmethd,ptype,&
       & ipart,afmt,istopc,itmax,itrace,irst,eps)

  call psb_barrier(ictxt)
  t1 = psb_wtime()  
  ! read the input matrix to be processed and (possibly) the rhs 
  nrhs = 1

  if (iam==psb_root_) then
    select case(psb_toupper(filefmt)) 
    case('MM') 
      ! For Matrix Market we have an input file for the matrix
      ! and an (optional) second file for the RHS. 
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
      if (info == 0) then 
        if (rhs_file /= 'NONE') then
          call mm_vet_read(aux_b,info,iunit=iunit,filename=rhs_file)
        end if
      end if
      
    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)
      
    case default
      info = -1 
      write(0,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= 0) then
      write(0,*) 'Error while reading input matrix '
      call psb_abort(ictxt)
    end if
    
    m_problem = aux_a%get_nrows()
    call psb_bcast(ictxt,m_problem)
    
    ! At this point aux_b may still be unallocated
    if (psb_size(aux_b,dim=1)==m_problem) then
      ! if any rhs were present, broadcast the first one
      write(0,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
    else
      write(*,'("Generating an rhs...")')
      write(*,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
        call psb_errpush(4000,name)
        goto 9999
      endif

      b_col_glob => aux_b(:,1)
      do i=1, m_problem
        b_col_glob(i) = 1.d0
      enddo
    endif
    call psb_bcast(ictxt,b_col_glob(1:m_problem))
  else
    call psb_bcast(ictxt,m_problem)
    call psb_realloc(m_problem,1,aux_b,ircode)
    if (ircode /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    endif
    b_col_glob =>aux_b(:,1)
    call psb_bcast(ictxt,b_col_glob(1:m_problem)) 
  end if

  ! switch over different partition types
  if (ipart == 0) then 
    call psb_barrier(ictxt)
    if (iam==psb_root_) write(*,'("Partition type: block")')
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
      call part_block(i,m_problem,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,v=ivg)
  else if (ipart == 2) then 
    if (iam==psb_root_) then 
      write(*,'("Partition type: graph")')
      write(*,'(" ")')
      !      write(0,'("Build type: graph")')
      call build_mtpart(aux_a,np)

    endif
    call psb_barrier(ictxt)
    call distr_mtpart(psb_root_,ictxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,v=ivg)
  else 
    if (iam==psb_root_) write(*,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,parts=part_block)
  end if

  call psb_geall(x_col,desc_a,info)
  x_col(:) =0.0
  call psb_geasb(x_col,desc_a,info)
  call psb_geall(r_col,desc_a,info)
  r_col(:) =0.0
  call psb_geasb(r_col,desc_a,info)
  t2 = psb_wtime() - t1


  call psb_amx(ictxt, t2)

  if (iam==psb_root_) then
    write(*,'(" ")')
    write(*,'("Time to read and partition matrix : ",es12.5)')t2
    write(*,'(" ")')
  end if

  ! 

  call psb_precinit(prec,ptype,info)

  ! building the preconditioner
  t1 = psb_wtime()
  call psb_precbld(a,desc_a,prec,info)
  tprec = psb_wtime()-t1
  if (info /= 0) then
    call psb_errpush(4010,name,a_err='psb_precbld')
    goto 9999
  end if


  call psb_amx(ictxt, tprec)

  if(iam==psb_root_) then
    write(*,'("Preconditioner time: ",es12.5)')tprec
    write(*,'(" ")')
  end if
  cond = -sone
  iparm = 0
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_krylov(kmethd,a,prec,b_col,x_col,eps,desc_a,info,& 
       & itmax=itmax,iter=iter,err=err,itrace=itrace,istop=istopc,&
       & irst=irst,cond=cond)     
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1

  call psb_amx(ictxt,t2)
  call psb_geaxpby(sone,b_col,szero,r_col,desc_a,info)
  call psb_spmm(-sone,a,x_col,sone,r_col,desc_a,info)
  call psb_genrm2s(resmx,r_col,desc_a,info)
  call psb_geamaxs(resmxp,r_col,desc_a,info)

  amatsize = psb_sizeof(a)
  descsize = psb_sizeof(desc_a)
  precsize = psb_sizeof(prec)
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  call psb_sum(ictxt,precsize)
  if (iam==psb_root_) then 
    call psb_precdescr(prec)
    write(*,'("Matrix: ",a)')mtrx_file
    write(*,'("Computed solution on ",i8," processors")')np
    write(*,'("Iterations to convergence: ",i6)')iter
    write(*,'("Error estimate on exit   : ",es12.5)') err
    write(*,'("Time to buil prec.       : ",es12.5)')tprec
    write(*,'("Time to solve matrix     : ",es12.5)')t2
    write(*,'("Time per iteration       : ",es12.5)')t2/(iter)
    write(*,'("Total time               : ",es12.5)')t2+tprec
    write(*,'("Residual norm 2          : ",es12.5)')resmx
    write(*,'("Residual norm inf        : ",es12.5)')resmxp
    write(*,*)"Condition number         : ",cond
    write(*,'("Total memory occupation for A:      ",i12)')amatsize
    write(*,'("Total memory occupation for DESC_A: ",i12)')descsize
    write(*,'("Total memory occupation for PREC:   ",i12)')precsize
  end if

  allocate(x_col_glob(m_problem),r_col_glob(m_problem),stat=ierr)
  if (ierr /= 0) then 
    write(0,*) 'allocation error: no data collection'
  else
    call psb_gather(x_col_glob,x_col,desc_a,info,root=psb_root_)
    call psb_gather(r_col_glob,r_col,desc_a,info,root=psb_root_)
    if (iam==psb_root_) then
      write(0,'(" ")')
      write(0,'("Saving x on file")')
      write(20,*) 'matrix: ',mtrx_file
      write(20,*) 'computed solution on ',np,' processors.'
      write(20,*) 'iterations to convergence: ',iter
      write(20,*) 'error estimate (infinity norm) on exit:', &
           & ' ||r||/(||a||||x||+||b||) = ',err
      write(20,'("Residual norm 2          : ",es12.5)')resmx
      write(20,'("Residual norm inf        : ",es12.5)')resmxp
      write(20,'(a8,4(2x,a20))') 'I','X(I)','R(I)','B(I)'
      do i=1,m_problem
        write(20,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
      enddo
    end if
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))


  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call psb_precfree(prec,info)
  call psb_cdfree(desc_a,info)

9999 continue
  if(info /= 0) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

end program sf_sample
  




