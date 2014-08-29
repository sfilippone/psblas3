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
program power_file_real
  use psb_base_mod
  use psb_util_mod
  use psb_eigen_mod
  implicit none

  ! input parameters
  character(len=100) :: kmethd, ptype, mtrx_file, rhs_file

  ! sparse matrices
  type(psb_dspmat_type) :: a, aux_a

  ! dense matrices
  real(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:)
  real(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob(:)
  type(psb_d_vect_type)    :: b_col, x_col, r_col


  ! communications data structure
  type(psb_desc_type):: desc_a

  integer(psb_ipk_) :: ictxt, iam, np

  ! solver paramters
  integer(psb_ipk_) :: iter, itmax, ierr, itrace, ircode, ipart,&
       & methd, istopc, irst, nr
  integer(psb_long_int_k_) :: amatsize, descsize, annz, nbytes
  real(psb_dpk_)   :: err, eps,cond

  character(len=5)   :: afmt
  character(len=20)  :: name
  character(len=2)   :: filefmt
  integer(psb_ipk_), parameter :: iunit=12
  integer(psb_ipk_)  :: times=0
  integer(psb_ipk_) :: iparm(20)

  ! other variables
  integer(psb_ipk_) :: i,info,j,m_problem
  integer(psb_ipk_) :: internal, m,ii,nnzero
  real(psb_dpk_) :: t1, t2, r_amax, b_amax,&
       &scale,resmx,resmxp, flops, bdwdth
  real(psb_dpk_) :: tt1, tt2, tflops
  real(psb_dpk_) :: lambda, lambda2
  real (psb_dpk_) :: norm, precisione
  integer(psb_ipk_) :: nrhs, nrow, n_row, dim, nv, ne
  integer(psb_ipk_), allocatable :: ivg(:), ipv(:)


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='power_file'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(2)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    !write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    !write(*,*) 'This is the ',trim(name),' sample program'
   write(psb_out_unit,'("debut")')
    read(psb_inp_unit,*) mtrx_file
    read(psb_inp_unit,*) filefmt
    read(psb_inp_unit,*) ipart
    read(psb_inp_unit,*) precisione
    !write (psb_out_unit, '("The precision of the power method is ",F30.20)') precisione     
  end if
  call psb_bcast(ictxt,mtrx_file)
  call psb_bcast(ictxt,filefmt)
  call psb_bcast(ictxt,ipart)
  call psb_bcast(ictxt,precisione)
  rhs_file = 'NONE'
  afmt     = 'CSR'
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
      if (info == psb_success_) then 
        if (rhs_file /= 'NONE') then
          call mm_array_read(aux_b,info,iunit=iunit,filename=rhs_file)
        end if
      end if
      
    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)

    case ('AD')
        call psb_d_adj_read(aux_a,mtrx_file,1,iunit,info)
      
    case default
      info = -1 
      write(psb_err_unit,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= psb_success_) then
      write(psb_err_unit,*) 'Error while reading input matrix '
      call psb_abort(ictxt)
    end if
    
    m_problem = aux_a%get_nrows()
    call psb_bcast(ictxt,m_problem)
    
    ! At this point aux_b may still be unallocated
    if (psb_size(aux_b,dim=1)==m_problem) then
      ! if any rhs were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
    else
      !write(psb_out_unit,'("Generating an rhs...")')
      !write(psb_out_unit,'(" ")')
      call psb_realloc(m_problem,1,aux_b,ircode)
      if (ircode /= 0) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
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
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    endif
    b_col_glob =>aux_b(:,1)
    call psb_bcast(ictxt,b_col_glob(1:m_problem)) 

  end if

  ! switch over different partition types
  if (ipart == 0) then 
    call psb_barrier(ictxt)
    !if (iam==psb_root_) write(psb_out_unit,'("Partition type: block")')
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
      call part_block(i,m_problem,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a,ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,v=ivg)
    
  else if (ipart == 2) then 
    if (iam==psb_root_) then 
      !write(psb_out_unit,'("Partition type: graph")')
      !write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call build_mtpart(aux_a,np)

    endif
    call psb_barrier(ictxt)
    call distr_mtpart(psb_root_,ictxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,v=ivg)

  else 
    !if (iam==psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,parts=part_block)
  end if

  t2 = psb_wtime() - t1
  call psb_amx(ictxt, t2)

  if (iam==psb_root_) then
    !write(psb_out_unit,'(" ")')
    !write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
    !write(psb_out_unit,'(" ")')
  end if

  call psb_barrier(ictxt)
  call psb_d_power_vect(a,x_col,lambda,precisione,iter,500,desc_a,info)

  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  call psb_amx(ictxt,t2) 

  nr       = desc_a%get_global_rows() 
  annz     = a%get_nzeros()
  amatsize = psb_sizeof(a)
  descsize = psb_sizeof(desc_a)
  call psb_sum(ictxt,annz)
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)
  
  if (iam==psb_root_) then 
    flops = 2.d0*times*annz
    tflops=flops
   ! write(psb_out_unit,'("Matrix: ",a)') mtrx_file
   ! write(psb_out_unit,'("Test on                          : ",i20," processors")') np
   ! write(psb_out_unit,'("Size of matrix                   : ",i20,"           ")') nr
   ! write(psb_out_unit,'("Number of nonzeros               : ",i20,"           ")') annz
    !write(psb_out_unit,'("Memory occupation                : ",i20,"           ")') amatsize
   ! write(psb_out_unit,'("Number of flops (",i0," iters)        : ",F20.0,"           ")') times,flops
   ! write(psb_out_unit,'("Time for ",i0," iteration (s)         : ",F20.3)')times, t2
   ! write(psb_out_unit,'("Time per iteration    (ms)         : ",F20.3)') t2*1.d3/(1.d0*times)
    open(15, FILE="resultats.dat", position = 'append',ACTION="WRITE")
    write (15, '(a,F20.6,i20)')mtrx_file,lambda,iter
   ! write (psb_out_unit, '("entropy vect propre :",g20.4)')entropy(x_col) 
    close(15)
  end if

  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call psb_cdfree(desc_a,info)

9999 continue
  if(info /= 0) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

contains
  subroutine adj_read (a,filename,iunit,desc_a,info)
    type(psb_dspmat_type), intent (inout) :: a
    character(len=100) :: filename
    integer (psb_ipk_) :: iunit
    type(psb_desc_type):: desc_a
    integer (psb_ipk_) :: info

    integer(psb_ipk_) :: i,nnzero,nrows
    integer (psb_ipk_) ::  iError
    !integer(psb_ipk_), allocatable :: ia(:),ja(:)
    !real (psb_dpk_), allocatable :: val(:)
    type(psb_d_coo_sparse_mat) :: acoo

    open(iunit, FILE=filename, STATUS="OLD", ACTION="READ")
    read(iunit, *) nrows , nnzero

    call acoo%allocate(nrows,nrows,nnzero)
    do i = 1,nnzero
        read(iunit, *) acoo%ia(i),acoo%ja(i)
        acoo%ia(i)=acoo%ia(i)+1
        acoo%ja(i)=acoo%ja(i)+1
        acoo%val(i)=1.0
    end do
    close(UNIT=iunit)
    !call psb_spall(a,desc_a,info,nnzero)
    !call psb_spins(nnzero, ia, ja, val, a, desc_a, info)
    call acoo%set_nzeros(nnzero)
    call acoo%fix(info)
    call a%mv_from(acoo)
    call a%cscnv(info,type='csr')
  end subroutine adj_read


  real(psb_dpk_) function entropy(vect)
        type(psb_d_vect_type), intent(inout)::vect

        real(psb_dpk_), allocatable :: v(:)
        integer(psb_ipk_) :: nrows

        nrows=vect%get_nrows()
        allocate(v(nrows))
        v=vect%get_vect()
        entropy=0
        do i=1,nrows
                entropy=entropy+v(i)*v(i)*log(v(i)*v(i))
        end do
        return
  end function entropy
end program power_file_real




