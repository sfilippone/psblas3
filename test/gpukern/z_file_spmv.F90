!    
!                Parallel Sparse BLAS  GPU plugin
!      (C) Copyright 2013
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
program z_file_spmv
  use psb_base_mod
  use psb_util_mod
  use psb_ext_mod
#ifdef HAVE_GPU
  use psb_gpu_mod
#endif
  use data_input
  implicit none

  ! input parameters
  character(len=200) :: mtrx_file

  ! sparse matrices
  type(psb_zspmat_type) :: a, aux_a, agpu

  ! dense matrices
  complex(psb_dpk_), allocatable, target :: aux_b(:,:), d(:)
  complex(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  complex(psb_dpk_), pointer  :: b_col_glob(:)
  type(psb_z_vect_type)    :: b_col, x_col, r_col
  type(psb_z_vect_type) :: xg, bg, xv, bv
#ifdef HAVE_GPU
  type(psb_z_vect_gpu)  :: vmold
#endif
  complex(psb_dpk_), allocatable :: xc1(:),xc2(:)
  ! communications data structure
  type(psb_desc_type):: desc_a

  type(psb_ctxt_type) :: ctxt
  integer            :: iam, np
  integer(psb_epk_) :: amatsize, agmatsize, precsize, descsize, annz, nbytes
  real(psb_dpk_)    :: damatsize, dgmatsize
  complex(psb_dpk_) :: err, eps

  character(len=5)   :: acfmt, agfmt
  character(len=20)  :: name
  character(len=2)   :: filefmt
  integer, parameter :: iunit=12
  integer, parameter :: times=2000 
  integer, parameter :: ntests=200, ngpu=50, ncnv=20 

  type(psb_z_coo_sparse_mat), target   :: acoo
  type(psb_z_csr_sparse_mat), target   :: acsr
  type(psb_z_ell_sparse_mat), target   :: aell
  type(psb_z_hll_sparse_mat), target   :: ahll
#ifdef HAVE_GPU
  type(psb_z_elg_sparse_mat), target   :: aelg
  type(psb_z_csrg_sparse_mat), target  :: acsrg
  type(psb_z_hybg_sparse_mat), target  :: ahybg
  type(psb_z_hlg_sparse_mat), target   :: ahlg
#endif
  class(psb_z_base_sparse_mat), pointer :: acmold, agmold
  ! other variables
  integer            :: i,info,j,nrt, ns, nr, ipart, ig, nrg
  integer            :: internal, m,ii,nnzero
  real(psb_dpk_) :: t0,t1, t2, tprec, flops
  real(psb_dpk_) :: tt1, tt2, tflops, gt1, gt2,gflops, gtint, bdwdth,&
       & tcnvcsr, tcnvc1, tcnvgpu, tcnvg1
  integer :: nrhs, nrow, n_row, dim, nv, ne
  integer, allocatable :: ivg(:), ipv(:)


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
#ifdef HAVE_GPU
  call psb_gpu_init(ctxt)
#endif
  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif


  name='file_spmv'
  if(psb_get_errstatus() /= 0) goto 9999
  info=psb_success_  
  call psb_set_errverbosity(2)
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
#ifdef HAVE_GPU
  write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice(),' out of', psb_cuda_getDeviceCount()
  write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice(),' is a: ', trim(psb_gpu_DeviceName())  
#endif

  if (iam == 0) then 
    write(*,*) 'Matrix? '
    call read_data(mtrx_file,psb_inp_unit)
    write(*,*) 'file format'
    call read_data(filefmt,psb_inp_unit)
    write(*,*) 'CPU format'
    call read_data(acfmt,psb_inp_unit)
    write(*,*) 'GPU format'
    call read_data(agfmt,psb_inp_unit)
    write(*,*) 'distribution '
    call read_data(ipart,psb_inp_unit)
    write(*,*) 'Read all data, going on'
  end if
  call psb_bcast(ctxt,mtrx_file)
  call psb_bcast(ctxt,filefmt)
  call psb_bcast(ctxt,acfmt)
  call psb_bcast(ctxt,agfmt)
  call psb_bcast(ctxt,ipart)
  call psb_barrier(ctxt)
  t0 = psb_wtime()  
  ! read the input matrix to be processed and (possibly) the rhs 
  nrhs = 1

  if (iam==psb_root_) then
    select case(psb_toupper(filefmt)) 
    case('MM') 
      ! For Matrix Market we have an input file for the matrix
      ! and an (optional) second file for the RHS. 
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mtrx_file)

    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,filename=mtrx_file)

    case default
      info = -1 
      write(psb_err_unit,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= 0) then
      write(psb_err_unit,*) 'Error while reading input matrix '
      call psb_abort(ctxt)
    end if

    !
    ! Always get nnz from original matrix.
    ! Some formats add fill-in and do not keep track
    ! of how many were added. So if the original matrix 
    ! contained some extra zeros, the count of entries
    ! is not recoverable exactly.
    !    
    nrt  = aux_a%get_nrows()
    annz = aux_a%get_nzeros()
    call psb_bcast(ctxt,annz)
    call psb_bcast(ctxt,nrt)

    write(psb_out_unit,'("Generating an rhs...")')
    write(psb_out_unit,'(" ")')
    call psb_realloc(nrt,1,aux_b,info)
    if (info /= 0) then
      call psb_errpush(4000,name)
      goto 9999
    endif

    b_col_glob => aux_b(:,1)
    do i=1, nrt
      b_col_glob(i) = 1.d0
    enddo

  else

    call psb_bcast(ctxt,annz)
    call psb_bcast(ctxt,nrt)

  end if


  select case(psb_toupper(acfmt))
  case('COO')
    acmold => acoo
  case('CSR')
    acmold => acsr
  case('ELL')
    acmold => aell
  case('HLL')
    acmold => ahll
  case default
    write(*,*) 'Unknown format defaulting to CSR'
    acmold => acsr
  end select

#ifdef HAVE_GPU
  select case(psb_toupper(agfmt))
  case('ELG')
    agmold => aelg
  case('HLG')
    agmold => ahlg
  case('CSRG')
    agmold => acsrg
  case('HYBG')
    agmold => ahybg
  case default
    write(*,*) 'Unknown format defaulting to HLG'
    agmold => ahlg
  end select
#endif


  ! switch over different partition types
  if (ipart == 0) then 
    call psb_barrier(ctxt)
    if (iam==psb_root_) write(psb_out_unit,'("Partition type: block")')
    allocate(ivg(nrt),ipv(np))
    do i=1,nrt
      call part_block(i,nrt,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a, ctxt, desc_a,info,v=ivg)
  else if (ipart == 2) then 
    if (iam==psb_root_) then 
      write(psb_out_unit,'("Partition type: graph")')
      write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call build_mtpart(aux_a,np)
    endif
    call psb_barrier(ctxt)
    call distr_mtpart(psb_root_,ctxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ctxt, desc_a,info,v=ivg)
  else 
    if (iam==psb_root_) write(psb_out_unit,'("Partition type default: block")')
    call psb_matdist(aux_a, a,  ctxt,desc_a,info,parts=part_block)
  end if

  call psb_scatter(b_col_glob,bv,desc_a,info,root=psb_root_)

  t2 = psb_wtime() - t0

  call psb_amx(ctxt, t2)

  if (iam==psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
    write(psb_out_unit,'(" ")')
  end if
  call a%cscnv(aux_a,info,mold=acoo)
  tcnvcsr = 0
  tcnvgpu = 0
  nr       = desc_a%get_local_rows()
  nrg      = desc_a%get_global_rows() 
  call psb_geall(x_col,desc_a,info)
  do i=1, nr
    call desc_a%l2g(i,ig,info)
    call psb_geins(ione,(/ig/),(/(zone + (zone*ig)/nrg)/),x_col,desc_a,info)
  end do
  call psb_geasb(x_col,desc_a,info)
  do j=1, ncnv
    call aux_a%cscnv(a,info,mold=acoo)
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call a%cscnv(info,mold=acmold)
    t2 = psb_Wtime() -t1
    call psb_amx(ctxt,t2)
    tcnvcsr = tcnvcsr + t2
    if (j==1) tcnvc1 = t2
    xc1 = x_col%get_vect()
    call xv%bld(xc1)
    call psb_geasb(bv,desc_a,info,scratch=.true.)
    
#ifdef HAVE_GPU
    
    call aux_a%cscnv(agpu,info,mold=acoo)
    call xg%bld(xc1,mold=vmold)
    call psb_geasb(bg,desc_a,info,scratch=.true.,mold=vmold)
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call agpu%cscnv(info,mold=agmold)
    call psb_gpu_DeviceSync()
    t2 = psb_Wtime() -t1
    call psb_amx(ctxt,t2)
    if (j==1) tcnvg1 = t2
    tcnvgpu = tcnvgpu + t2
#endif
  end do

  call psb_barrier(ctxt)
  t1 = psb_wtime()
  do i=1,ntests 
    call psb_spmm(zone,a,xv,zzero,bv,desc_a,info)
  end do
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  call psb_amx(ctxt,t2)

#ifdef HAVE_GPU
  ! FIXME: cache flush needed here
  call psb_barrier(ctxt)
  tt1 = psb_wtime()
  do i=1,ntests 
    call psb_spmm(zone,agpu,xv,zzero,bg,desc_a,info)
    if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
      write(0,*) 'From 1 spmm',info,i,ntests
      call psb_error()
      stop
    end if

  end do
  call psb_gpu_DeviceSync()
  call psb_barrier(ctxt)
  tt2 = psb_wtime() - tt1
  call psb_amx(ctxt,tt2)
  xc1 = bv%get_vect()
  xc2 = bg%get_vect()
  nr       = desc_a%get_local_rows() 
  eps = maxval(abs(xc1(1:nr)-xc2(1:nr)))
  call psb_amx(ctxt,eps)
  if (iam==0) write(*,*) 'Max diff on xGPU',eps

  call xg%sync()
  ! FIXME: cache flush needed here

  call psb_barrier(ctxt)
  gt1 = psb_wtime()
  do i=1,ntests*ngpu
    call psb_spmm(zone,agpu,xg,zzero,bg,desc_a,info)
    if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
      write(0,*) 'From 2 spmm',info,i,ntests
      call psb_error()
      stop
    end if

  end do
  ! For timing purposes we need to make sure all threads
  ! in the device are done. 
  call psb_gpu_DeviceSync()
  call psb_barrier(ctxt)
  gt2 = psb_wtime() - gt1
  call psb_amx(ctxt,gt2)
  call bg%sync()
  xc1 = bv%get_vect()
  xc2 = bg%get_vect()
  call psb_geaxpby(-zone,bg,+zone,bv,desc_a,info)
  eps = psb_geamax(bv,desc_a,info)

  call psb_amx(ctxt,t2)
  nr       = desc_a%get_local_rows() 
  eps = maxval(abs(xc1(1:nr)-xc2(1:nr)))
  call psb_amx(ctxt,eps)
  if (iam==0) write(*,*) 'Max diff on GPU',eps
#endif


  amatsize = a%sizeof()
  agmatsize = agpu%sizeof()
  damatsize = amatsize
  damatsize = damatsize/(1024*1024)
  dgmatsize = agmatsize
  dgmatsize = dgmatsize/(1024*1024)
  descsize = psb_sizeof(desc_a)
  call psb_sum(ctxt,damatsize)
  call psb_sum(ctxt,dgmatsize)
  call psb_sum(ctxt,descsize)

  if (iam == psb_root_) then
    write(psb_out_unit,'("Matrix: ",a)') mtrx_file
    write(psb_out_unit,&
         &'("Test on                          : ",i20," processors")') np
    write(psb_out_unit,&
         &'("Size of matrix                   : ",i20,"           ")') nrt
    write(psb_out_unit,&
         &'("Number of nonzeros               : ",i20,"           ")') annz
    write(psb_out_unit,&
         &'("Memory occupation CPU  (MBytes)  : ",f20.2,"           ")') damatsize
    write(psb_out_unit,&
         &'("Memory occupation GPU  (MBytes)  : ",f20.2,"           ")') dgmatsize
    write(psb_out_unit,&
         &'("Memory occupation CPU  (Bytes)   : ",i24,"           ")') amatsize
    write(psb_out_unit,&
         &'("Memory occupation GPU  (Bytes)   : ",i24,"           ")') agmatsize
    flops  = ntests*(2.d0*annz)
    tflops = flops
    gflops = flops * ngpu
    write(psb_out_unit,'("Storage type for    A: ",a)') a%get_fmt()
#ifdef HAVE_GPU
    write(psb_out_unit,'("Storage type for AGPU: ",a)') agpu%get_fmt()
    write(psb_out_unit,'("Time to convert A from COO to CPU (1): ",F20.9)')&
         & tcnvc1
    write(psb_out_unit,'("Time to convert A from COO to CPU (t): ",F20.9)')&
         & tcnvcsr
    write(psb_out_unit,'("Time to convert A from COO to CPU (a): ",F20.9)')&
         & tcnvcsr/ncnv
    write(psb_out_unit,'("Time to convert A from COO to GPU (1): ",F20.9)')&
         & tcnvg1
    write(psb_out_unit,'("Time to convert A from COO to GPU (t): ",F20.9)')&
         & tcnvgpu
    write(psb_out_unit,'("Time to convert A from COO to GPU (a): ",F20.9)')&
         & tcnvgpu/ncnv

#endif
    write(psb_out_unit,&
         & '("Number of flops (",i0," prod)        : ",F20.0,"           ")') &
         &  ntests,flops

    flops  = flops / (t2)
    tflops = tflops / (tt2)
    gflops = gflops / (gt2)
    write(psb_out_unit,'("Time for ",i6," products (s) (CPU)   : ",F20.3)')&
         &  ntests,t2
    write(psb_out_unit,'("Time per product    (ms)     (CPU)   : ",F20.3)')&
         & t2*1.d3/(1.d0*ntests)
    write(psb_out_unit,'("MFLOPS                       (CPU)   : ",F20.3)')&
         & flops/1.d6
#ifdef HAVE_GPU

    write(psb_out_unit,'("Time for ",i6," products (s) (xGPU)  : ",F20.3)')&
         & ntests, tt2
    write(psb_out_unit,'("Time per product    (ms)     (xGPU)  : ",F20.3)')&
         & tt2*1.d3/(1.d0*ntests)
    write(psb_out_unit,'("MFLOPS                       (xGPU)  : ",F20.3)')&
         & tflops/1.d6

    write(psb_out_unit,'("Time for ",i6," products (s) (GPU)   : ",F20.3)')&
         & ngpu*ntests,gt2
    write(psb_out_unit,'("Time per product    (ms)     (GPU)   : ",F20.3)')&
         & gt2*1.d3/(1.d0*ntests*ngpu)
    write(psb_out_unit,'("MFLOPS                       (GPU)   : ",F20.3)')&
         & gflops/1.d6
#endif
    !
    ! This computation assumes the data movement associated with CSR:
    ! it is minimal in terms of coefficients. Other formats may either move
    ! more data (padding etc.) or less data (if they can save on the indices). 
    !
    nbytes = nr*(2*2*psb_sizeof_dp + psb_sizeof_ip)+&
         & annz*(2*psb_sizeof_dp + psb_sizeof_ip)
    bdwdth = ntests*nbytes/(t2*1.d6)
    write(psb_out_unit,*)
    write(psb_out_unit,'("MBYTES/S                  (CPU)  : ",F20.3)') bdwdth
#ifdef HAVE_GPU
    bdwdth = ngpu*ntests*nbytes/(gt2*1.d6)
    write(psb_out_unit,'("MBYTES/S                  (GPU)  : ",F20.3)') bdwdth
#endif
    write(psb_out_unit,'("Storage type for DESC_A: ",a)') desc_a%indxmap%get_fmt()
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize

  end if

  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_gefree(xv, desc_a,info)
  call psb_gefree(bv, desc_a,info)
  call psb_spfree(a, desc_a,info)
#ifdef HAVE_GPU
  call psb_gefree(xg, desc_a,info)
  call psb_gefree(bg, desc_a,info)
  call psb_spfree(agpu,desc_a,info)
  call psb_gpu_exit()
#endif
  call psb_cdfree(desc_a,info)

  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

end program z_file_spmv
  




