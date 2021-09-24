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
program s_file_spmv
  use psb_base_mod
  use psb_util_mod
  implicit none

  ! input parameters
  character(len=40) :: kmethd, ptype, mtrx_file, rhs_file

  ! sparse matrices
  type(psb_sspmat_type) :: a
  type(psb_lsspmat_type) :: aux_a

  ! dense matrices
  real(psb_spk_), allocatable, target ::  aux_b(:,:), d(:)
  real(psb_spk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  real(psb_spk_), pointer  :: b_col_glob(:)
  type(psb_s_vect_type)    :: b_col, x_col, r_col


  ! communications data structure
  type(psb_desc_type):: desc_a

  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: iam, np

  ! solver paramters
  integer(psb_ipk_) :: iter, itmax, ierr, itrace, ircode, ipart,&
       & methd, istopc, irst, nr
  integer(psb_lpk_) :: lnp
  integer(psb_epk_) :: amatsize, descsize, annz, nbytes
  real(psb_spk_)   :: err, eps,cond

  character(len=5)   :: afmt
  character(len=20)  :: name
  character(len=2)   :: filefmt
  integer(psb_ipk_), parameter :: iunit=12
  integer(psb_ipk_), parameter :: times=20
  integer(psb_ipk_) :: iparm(20)

  ! other variables
  integer(psb_lpk_) :: i,j,m_problem
  integer(psb_ipk_) :: internal, m,ii,nnzero,info
  real(psb_dpk_) :: t1, t2, r_amax, b_amax,&
       &scale,resmx,resmxp, flops, bdwdth
  real(psb_dpk_) :: tt1, tt2, tflops
  integer(psb_ipk_) :: nrhs, nrow, n_row, dim, nv, ne
  integer(psb_ipk_), allocatable :: ivg(:), ipv(:)


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif


  name='s_file_spmv'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_
  call psb_set_errverbosity(2)
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
    read(psb_inp_unit,*) mtrx_file
    read(psb_inp_unit,*) filefmt
    read(psb_inp_unit,*) ipart
  end if
  call psb_bcast(ctxt,mtrx_file)
  call psb_bcast(ctxt,filefmt)
  call psb_bcast(ctxt,ipart)
  rhs_file = 'NONE'
  afmt     = 'CSR'
  call psb_barrier(ctxt)
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
    if (psb_size(aux_b,dim=1)==m_problem) then
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
        b_col_glob(i) = 1.d0
      enddo
    endif

  else

    call psb_bcast(ctxt,m_problem)
    b_col_glob =>aux_b(:,1)

  end if

  ! switch over different partition types
  write(psb_out_unit,'("Number of processors : ",i0)')np
  if (ipart == 0) then 
    call psb_barrier(ctxt)
    if (iam==psb_root_) write(psb_out_unit,'("Partition type: block")')
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
      call part_block(i,m_problem,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a, ctxt,desc_a,info,fmt=afmt,vg=ivg)
    
  else if (ipart == 2) then 
    if (iam==psb_root_) then 
      write(psb_out_unit,'("Partition type: graph")')
      write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call build_mtpart(aux_a,lnp)

    endif
    call psb_barrier(ctxt)
    call distr_mtpart(psb_root_,ctxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ctxt, desc_a,info,fmt=afmt,vg=ivg)

  else 
    if (iam==psb_root_) write(psb_out_unit,'("Partition type: default block")')
    call psb_matdist(aux_a, a,  ctxt, desc_a,info,fmt=afmt,parts=part_block)
  end if

  call psb_geall(x_col,desc_a,info)
  call x_col%set(sone)
  call psb_geasb(x_col,desc_a,info)
  call psb_geall(b_col,desc_a,info)
  call x_col%zero()
  call psb_geasb(b_col,desc_a,info)
  t2 = psb_wtime() - t1


  call psb_amx(ctxt, t2)

  if (iam==psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
    write(psb_out_unit,'(" ")')
  end if


  call psb_barrier(ctxt)
  t1 = psb_wtime()
  do i=1,times 
    call psb_spmm(sone,a,x_col,szero,b_col,desc_a,info,'n')
  end do
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  call psb_amx(ctxt,t2)

  ! FIXME: cache flush needed here
  call psb_barrier(ctxt)
  tt1 = psb_wtime()
  do i=1,times 
    call psb_spmm(sone,a,x_col,szero,b_col,desc_a,info,'t')
  end do
  call psb_barrier(ctxt)
  tt2 = psb_wtime() - tt1
  call psb_amx(ctxt,tt2)

  nr       = desc_a%get_global_rows() 
  annz     = a%get_nzeros()
  amatsize = psb_sizeof(a)
  descsize = psb_sizeof(desc_a)
  call psb_sum(ctxt,annz)
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  
  if (iam==psb_root_) then 
    flops = 2.d0*times*annz
    tflops=flops
    write(psb_out_unit,'("Matrix: ",a)') mtrx_file
    write(psb_out_unit,'("Test on                          : ",i20," processors")') np
    write(psb_out_unit,'("Size of matrix                   : ",i20,"           ")') nr
    write(psb_out_unit,'("Number of nonzeros               : ",i20,"           ")') annz
    write(psb_out_unit,'("Memory occupation                : ",i20,"           ")') amatsize
    write(psb_out_unit,'("Number of flops (",i0," prod)        : ",F20.0,"           ")') times,flops
    flops = flops / (t2)
    tflops = tflops / (tt2)
    write(psb_out_unit,'("Time for ",i0," products (s)         : ",F20.3)')times, t2
    write(psb_out_unit,'("Time per product    (ms)         : ",F20.3)') t2*1.d3/(1.d0*times)
    write(psb_out_unit,'("MFLOPS                           : ",F20.3)') flops/1.d6

    write(psb_out_unit,'("Time for ",i0," products (s) (trans.): ",F20.3)') times,tt2
    write(psb_out_unit,'("Time per product    (ms) (trans.): ",F20.3)') tt2*1.d3/(1.d0*times)
    write(psb_out_unit,'("MFLOPS                   (trans.): ",F20.3)') tflops/1.d6

    !
    ! This computation is valid for CSR
    !
    nbytes = nr*(2*psb_sizeof_sp + psb_sizeof_ip)+ &
         & annz*(psb_sizeof_sp + psb_sizeof_ip)
    bdwdth = times*nbytes/(t2*1.d6)
    write(psb_out_unit,*)
    write(psb_out_unit,'("MBYTES/S                         : ",F20.3)') bdwdth
    bdwdth = times*nbytes/(tt2*1.d6)
    write(psb_out_unit,'("MBYTES/S                  (trans): ",F20.3)') bdwdth
    
  end if

  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call psb_cdfree(desc_a,info)

  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)

  stop

end program s_file_spmv
  




