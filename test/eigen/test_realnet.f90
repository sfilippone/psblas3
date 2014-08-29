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
program maxA_and_lapl_extremums
  use psb_base_mod
  use psb_util_mod
  use psb_eigen_mod
  implicit none

  ! input parameters
  character(len=100) ::  mtrx_file, rhs_file
  character(len=40) :: kmethd, ptype
  real (psb_dpk_) :: precisione,sigma
  integer(psb_ipk_) :: dim_H,to_read

  ! sparse matrices
  type(psb_dspmat_type) :: a, b, aux_a

  ! dense matrices
  real(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:)
  real(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  real(psb_dpk_), pointer  :: b_col_glob(:)
  type(psb_d_vect_type)    :: b_col, x_col, r_col
  complex(psb_dpk_),allocatable :: eig(:),eigmin(:)


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
  integer(psb_ipk_) :: nnzero
  real(psb_dpk_) :: t1, t2, t3, t4
  real(psb_dpk_) :: lambda, lambda2, norm
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
    write(*,*) 'This is the ',trim(name),' sample program'
    read(psb_inp_unit,*) to_read
    read(psb_inp_unit,*) mtrx_file
    read(psb_inp_unit,*) filefmt
    read(psb_inp_unit,*) ipart
    read(psb_inp_unit,*) precisione
    read(psb_inp_unit,*) dim_H
    read(psb_inp_unit,*) sigma
    read(psb_inp_unit,*) ptype
    read(psb_inp_unit,*) kmethd
  end if
  call psb_bcast(ictxt,to_read)
  call psb_bcast(ictxt,mtrx_file)
  call psb_bcast(ictxt,filefmt)
  call psb_bcast(ictxt,ipart)
  call psb_bcast(ictxt,precisione)
  call psb_bcast(ictxt,dim_H)
  call psb_bcast(ictxt,sigma)
  call psb_bcast(ictxt,ptype)
  call psb_bcast(ictxt,kmethd)
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
        call special_adj_read(aux_a,mtrx_file,to_read,iunit,info)
      
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
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
      call part_block(i,m_problem,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a,ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,v=ivg)
    
  else if (ipart == 2) then 
    if (iam==psb_root_) then 
      call build_mtpart(aux_a,np)

    endif
    call psb_barrier(ictxt)
    call distr_mtpart(psb_root_,ictxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,v=ivg)

  else 
    call psb_matdist(aux_a, a,  ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,parts=part_block)
  end if

  call psb_barrier(ictxt)
  call psb_d_power_vect(a,x_col,lambda,precisione,iter,500,desc_a,info)
  
  t2 = psb_wtime() - t1
  t1 = psb_wtime()
  call psb_d_laplacian(a,b)
  call psb_d_arnoldi(b,dim_H,eig,desc_a,ictxt,info)
  
  t3 = psb_wtime()-t1
  t1 = psb_wtime()
  call psb_d_shift_invert(b,dim_H,sigma,eigmin,ptype,kmethd,desc_a,ictxt,iam,info)
  
  call psb_barrier(ictxt)
  t4 = psb_wtime() - t1
  
  t1=t2+t3+t4
  if (iam==psb_root_) then 
    open (15, FILE="plot_data.dat", position = 'append',ACTION="WRITE")
    write (15,'(i20,F20.6,F20.6,F20.6,F20.6,F20.4)')aux_a%get_nzeros(),lambda,real(eig(dim_H)),real(eig(dim_H-1)),&
                &real(eigmin(dim_H-1)),t1
    close(15)
  end if

  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call psb_spfree(b, desc_a,info)
  call psb_spfree(aux_a, desc_a,info)
  call psb_cdfree(desc_a,info)

9999 continue
  if(info /= 0) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

contains
subroutine special_adj_read (a,filename,nnzero,iunit,info)
    use psb_base_mod
    use psb_eigen_mod, psb_protect_name => psb_d_adj_read

    type(psb_dspmat_type), intent (inout) :: a
    character(len=40) , intent(in):: filename
    integer (psb_ipk_),intent(in) :: iunit
    integer (psb_ipk_),intent(out) :: info

    ! reads a file containing the coordinate description of an adjacence
    ! sparse matrix (values are always 1) and stores it into 'a' 
    !
    ! on entry :
    ! filename : file containing the description of the matrix
    ! nnzero : number of elements you want to read        

    integer(psb_ipk_) :: i,nnzero,nrows
    type(psb_d_coo_sparse_mat) :: acoo

    open(iunit, FILE=filename, STATUS="OLD", ACTION="READ")
    read(iunit, *) nrows

    call acoo%allocate(nrows,nrows,nnzero)
    do i = 1,nnzero
        read(iunit, *) acoo%ia(i),acoo%ja(i)
        acoo%val(i)=1.0
    end do
    close(UNIT=iunit)
    call acoo%set_nzeros(nnzero)
    call acoo%fix(info)
    call a%mv_from(acoo)
    call a%cscnv(info,type='csr')
end subroutine special_adj_read

end program maxA_and_lapl_extremums




