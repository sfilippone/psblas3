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
!
program psb_zkrylov_sample
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod
  use psb_zkrylovsubspace_mod
  use getp
  implicit none

  ! input parameters
  character(len=40) :: kmethd, ptype, mtrx_file, rhs_file

  ! sparse matrices
  type(psb_zspmat_type)  :: a
  type(psb_lzspmat_type) :: aux_a
  ! Krylov space
  type(psb_zpolkrylov) :: kryl

  ! dense matrices
  complex(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:)
  complex(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  complex(psb_dpk_), pointer  :: b_col_glob(:)
  type(psb_z_vect_type)    :: b_col, x_col, r_col

  ! communications data structure
  type(psb_desc_type):: desc_a

  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: iam, np
  integer(psb_lpk_) :: lnp
  ! solver paramters
  integer(psb_ipk_) :: iter, itmax, ierr, itrace, ircode,&
       & methd, istopc, irst
  integer(psb_epk_) :: amatsize, descsize
  real(psb_dpk_)   :: err, eps, cond

  character(len=5)   :: afmt
  character(len=20)  :: name, part
  character(len=2)   :: filefmt
  integer(psb_ipk_), parameter :: iunit=12
  integer(psb_ipk_) :: iparm(20)

  ! other variables
  integer(psb_ipk_) :: i,info,j,m_problem, err_act
  integer(psb_ipk_) :: internal, m,ii,nnzero
  real(psb_dpk_) :: t1, t2, tprec
  real(psb_dpk_) :: r_amax, b_amax, scale,resmx,resmxp
  integer(psb_ipk_) :: nrhs, nrow, n_row, dim, ne, nv
  integer(psb_ipk_), allocatable :: ivg(:)
  integer(psb_ipk_), allocatable :: ipv(:)
  character(len=40)  :: fname, fnout


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)

  if (iam < 0) then
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif


  name='psb_zkrylov_sample'
  if(psb_errstatus_fatal()) goto 9999
  info=psb_success_
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
  !
  !  get parameters
  !
  call get_parms(ctxt,mtrx_file,rhs_file,filefmt,kmethd,ptype,&
       & part,afmt,istopc,itmax,itrace,irst,eps)

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  ! read the input matrix to be processed and (possibly) the rhs
  nrhs = 1

  if (iam == psb_root_) then
    select case(psb_toupper(filefmt))
    case('MM')
      ! For Matrix Market we have an input file for the matrix
      ! and an (optional) second file for the RHS.
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)
      if (info == psb_success_) then
        if (rhs_file /= 'NONE') then
          call mm_array_read(aux_b,info,iunit=iunit,filename=rhs_file)
        end if
      end if

    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)

    case default
      info = -1
      write(psb_err_unit,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= psb_success_) then
      write(psb_err_unit,*) 'Error while reading input matrix '
      call psb_abort(ctxt)
    end if

    m_problem = aux_a%get_nrows()
    call psb_bcast(ctxt,m_problem)

    ! At this point aux_b may still be unallocated
    if (size(aux_b,dim=1) == m_problem) then
      ! if any rhs were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
    else
      write(psb_out_unit,'("Generating an rhs...")')
      write(psb_out_unit,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
         call psb_errpush(psb_err_alloc_dealloc_,name)
         goto 9999
      endif

      b_col_glob => aux_b(:,1)
      do i=1, m_problem
         b_col_glob(i) = zone
      enddo
    endif

  else
    call psb_bcast(ctxt,m_problem)

  end if

  ! switch over different partition types
  select case(psb_toupper(part))
  case('BLOCK')
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ctxt,desc_a,info,fmt=afmt,parts=part_block)

  case('GRAPH')
    if (iam == psb_root_) then
      write(psb_out_unit,'("Partition type: graph vector")')
      write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call aux_a%cscnv(info,type='csr')
      lnp = np
      call build_mtpart(aux_a,lnp)

    endif
    call psb_barrier(ctxt)
    call distr_mtpart(psb_root_,ctxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ctxt,desc_a,info,fmt=afmt,vg=ivg)

  case default
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ctxt,desc_a,info,fmt=afmt,parts=part_block)
  end select

  call psb_scatter(b_col_glob,b_col,desc_a,info,root=psb_root_)
  call psb_geall(x_col,desc_a,info)
  call x_col%zero()
  call psb_geasb(x_col,desc_a,info)
  call psb_geall(r_col,desc_a,info)
  call r_col%zero()
  call psb_geasb(r_col,desc_a,info)
  t2 = psb_wtime() - t1


  call psb_amx(ctxt, t2)

  if (iam == psb_root_) then
     write(psb_out_unit,'(" ")')
     write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
     write(psb_out_unit,'(" ")')
  end if

  ! Allocate the Krylov space
  call kryl%allocate(desc_a,itmax,info)
  ! Build
  call kryl%arnoldi(a,desc_a,b_col,istopc)
  ! Write
  call kryl%write(desc_a,mtrx_file)

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  if (iam == psb_root_) then
    write(psb_out_unit,'("Matrix: ",a)')mtrx_file
    write(psb_out_unit,'("Computed space on ",i8," processors")')np
    write(psb_out_unit,'("Time to build space      : ",es12.5)')t2
    write(psb_out_unit,'("Time per iteration       : ",es12.5)')t2/(iter)
    write(psb_out_unit,'("Total memory occupation for A:      ",i12)')amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize
    write(psb_out_unit,'("Storage format for A              : ",a)')&
         &  a%get_fmt()
    write(psb_out_unit,'("Storage format for DESC_A         : ",a)')&
         &  desc_a%get_fmt()
  end if

  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call kryl%free(desc_a,info)
  call psb_cdfree(desc_a,info)
  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)

  stop
end program psb_zkrylov_sample
