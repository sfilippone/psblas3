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
!         software without specific prior written permission.
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
! File: dpdegenmv.f90
!
! Program: pdegenmv
! This sample program measures the performance of the matrix-vector product.
! The matrix is generated in the same way as for the pdegen test case of
! the main PSBLAS library.
!
!
program pdgenmv
  use psb_base_mod
  use psb_util_mod 
  use psb_ext_mod
  use psb_d_pde3d_mod
#ifdef PSB_OPENACC
  use psb_oacc_mod
#endif
#ifdef PSB_HAVE_CUDA
  use psb_cuda_mod
#endif
  implicit none

  ! input parameters
  character(len=5)  :: acfmt, agfmt
  integer   :: idim
  logical   :: tnd
  ! miscellaneous 
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_) :: t1, t2, tprec, flops, tflops,&
       & tt1, tt2, gt1, gt2, gflops, bdwdth,&
       & tcnvcsr, tcnvc1, tcnvgpu, tcnvg1

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a, agpu, aux_a
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense matrices
  type(psb_d_vect_type), target :: xv, bv, xg, bg  
  real(psb_dpk_), allocatable :: x1(:), x2(:), x0(:)
  ! blacs parameters
  type(psb_ctxt_type) :: ctxt
  integer            :: iam, np

  ! solver parameters
  integer(psb_epk_) :: amatsize, precsize, descsize, annz, nbytes
  real(psb_dpk_)   :: err, eps, tnv, tng,tdot, dnrm2,ddot
  integer, parameter :: ntests=8, ngpu=4, ncnv=3
  type(psb_d_coo_sparse_mat), target   :: acoo
  type(psb_d_csr_sparse_mat), target   :: acsr
  type(psb_d_ell_sparse_mat), target   :: aell
  type(psb_d_hll_sparse_mat), target   :: ahll
  type(psb_d_dia_sparse_mat), target   :: adia
  type(psb_d_hdia_sparse_mat), target   :: ahdia
  type(psb_d_base_vect_type), target   :: dvect
  type(psb_i_base_vect_type), target   :: ivect
  class(psb_d_base_sparse_mat), pointer :: acmold => acsr
  class(psb_d_base_vect_type), pointer :: vmold => dvect
  class(psb_i_base_vect_type), pointer :: imold => ivect
  class(psb_d_base_vect_type), pointer :: vgmold => dvect
  class(psb_i_base_vect_type), pointer :: igmold => ivect
#if defined(PSB_HAVE_CUDA)
  type(psb_d_cuda_hlg_sparse_mat), target   :: ahlg
  type(psb_d_cuda_hdiag_sparse_mat), target :: ahdiag
  type(psb_d_cuda_csrg_sparse_mat), target   :: acsrg
  type(psb_d_cuda_elg_sparse_mat), target    :: aelg
  type(psb_d_vect_cuda), target         :: dvgpu
  type(psb_i_vect_cuda), target         :: ivgpu
#endif
#if defined(PSB_OPENACC)
  type(psb_d_oacc_hll_sparse_mat), target   :: ahlo
  type(psb_d_oacc_csr_sparse_mat), target   :: acsro
  type(psb_d_oacc_ell_sparse_mat), target    :: aelo
  type(psb_d_vect_oacc), target         :: dvoac
  type(psb_i_vect_oacc), target         :: ivoac
#endif
  class(psb_d_base_sparse_mat), pointer :: agmold
  ! other variables
  logical, parameter :: dump=.false.
  integer(psb_ipk_)  :: info, i, j, nr, nrg
  integer(psb_lpk_)  :: ig
  character(len=20)  :: name,ch_err
  character(len=40)  :: fname

  info=psb_success_


  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)

#ifdef PSB_OPENACC
  call psb_oacc_init(ctxt)
#endif
#ifdef PSB_HAVE_CUDA
  call psb_cuda_init(ctxt)
#endif

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ctxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='pdegenmv-oacc'
  !
  ! Hello world
  !
  if (iam == psb_root_) then 
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if
#ifdef PSB_OPENACC
!!$  write(*,*) 'Process ',iam,' running on device: ', psb_oacc_getDevice(),' out of', psb_oacc_getDeviceCount()
!!$  write(*,*) 'Process ',iam,' device ', psb_oacc_getDevice(),' is a: ', trim(psb_oacc_DeviceName())  
#endif
#if defined(PSB_HAVE_CUDA)
  write(*,*) 'Process ',iam,' running on device: ', psb_cuda_getDevice(),' out of', psb_cuda_getDeviceCount()
  write(*,*) 'Process ',iam,' device ', psb_cuda_getDevice(),' is a: ', trim(psb_cuda_DeviceName())  
#endif
  !
  !  get parameters
  !
  call get_parms(ctxt,acfmt,agfmt,idim,tnd)
  call psb_init_timers()
  !
  !  allocate and fill in the coefficient matrix and initial vectors
  !
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  call psb_d_gen_pde3d(ctxt,idim,a,bv,xv,desc_a,'CSR  ',info,partition=3,tnd=tnd)  
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='create_matrix'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (iam == psb_root_) write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(psb_out_unit,'(" ")')

  if (dump) then 
    write(fname,'(a,i3.3,a,i3.3,a,i3.3,a)')  'pde',idim,'-',iam,'-',np,'.mtx'
    call a%print(fname,head='PDEGEN test matrix')
  end if

  select case(psb_toupper(acfmt))
  case('ELL')
    acmold => aell
  case('HLL')
    acmold => ahll
  case('DIA')
    acmold => adia
  case('HDIA')
    acmold => ahdia
  case('CSR')
    acmold => acsr
  case('COO')
    acmold => acoo
  case default
    write(*,*) 'Unknown format defaulting to HLL'
    acmold => ahll
  end select
  call a%cscnv(info,mold=acmold)
  if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
    write(0,*) 'From cscnv ',info
    call psb_error()
    stop
  end if

  select case(psb_toupper(agfmt))
#if defined(PSB_HAVE_CUDA)
  case('ELG')
    agmold => aelg
    vgmold  => dvgpu
    igmold  => ivgpu
  case('HLG')
    call psi_set_hksz(32)
    agmold => ahlg
    vgmold  => dvgpu
    igmold  => ivgpu
  case('HDIAG')
    agmold => ahdiag
    vgmold  => dvgpu
    igmold  => ivgpu
  case('CSRG')
    agmold => acsrg
    vgmold  => dvgpu
    igmold  => ivgpu
#endif
#if defined(PSB_OPENACC)
  case('ELO')
    agmold => aelo
    vgmold  => dvoac
    igmold  => ivoac
  case('HLO')
    call psi_set_hksz(32)
    agmold => ahlo
    vgmold  => dvoac
    igmold  => ivoac
!!$  case('HDIAG')
!!$    amold => ahdiag
  case('CSRO')
    agmold => acsro
    vgmold  => dvoac
    igmold  => ivoac
#endif
  end select
  call a%cscnv(agpu,info,mold=agmold)
  if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
    write(0,*) 'From cscnv ',info
    call psb_error()
    stop
  end if
  call desc_a%cnv(mold=igmold)
  
  call psb_geasb(bg,desc_a,info,scratch=.true.,mold=vgmold)
  call psb_geasb(xg,desc_a,info,scratch=.true.,mold=vgmold)
  nr       = desc_a%get_local_rows()
  nrg      = desc_a%get_global_rows() 
  call psb_geall(x0,desc_a,info)
  do i=1, nr
    call desc_a%l2g(i,ig,info)
    x0(i) = 1.0 + (1.0*ig)/(nrg**2)
  end do
  call a%cscnv(aux_a,info,mold=acoo)
  tcnvcsr = 0
  tcnvgpu = 0
  call psb_geall(x1,desc_a,info)
  do j=1, ncnv
    call aux_a%cscnv(a,info,mold=acoo)
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call a%cscnv(info,mold=acmold)
    t2 = psb_Wtime() -t1
    call psb_amx(ctxt,t2)
    tcnvcsr = tcnvcsr + t2
    if (j==1) tcnvc1 = t2
    call psb_geasb(x1,desc_a,info)
    call xv%bld(x0)
    call psb_geasb(bv,desc_a,info,scratch=.true.)
    
#if defined(PSB_OPENACC)||defined(PSB_HAVE_CUDA)   
    call aux_a%cscnv(agpu,info,mold=acoo)
    call xg%bld(x0,mold=vgmold)
    call psb_geasb(bg,desc_a,info,scratch=.true.,mold=vgmold)
    call psb_barrier(ctxt)
    t1 = psb_wtime()
    call agpu%cscnv(info,mold=agmold)
!!$    call psb_oacc_DeviceSync()
    t2 = psb_Wtime() -t1
    call psb_amx(ctxt,t2)
    if (j==1) tcnvg1 = t2
    tcnvgpu = tcnvgpu + t2
#endif
  end do


  call xv%set(x0)
  call psb_barrier(ctxt)
  t1 = psb_wtime()
  do i=1,ntests 
    call psb_spmm(done,a,xv,dzero,bv,desc_a,info)
  end do
  call psb_barrier(ctxt)
  t2 = psb_wtime() - t1
  call psb_amx(ctxt,t2)

#if defined(PSB_OPENACC)||defined(PSB_HAVE_CUDA)
  call xg%set(x0)

  ! FIXME: cache flush needed here
  x1 = bv%get_vect()
  x2 = bg%get_vect()
  
  call psb_barrier(ctxt)
  tt1 = psb_wtime()
  do i=1,ntests 
    call psb_spmm(done,agpu,xv,dzero,bg,desc_a,info)
    if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
      write(0,*) 'From 1 spmm',info,i,ntests
      call psb_error()
      stop
    end if

  end do
!!$  call psb_oacc_DeviceSync()
  call psb_barrier(ctxt)
  tt2 = psb_wtime() - tt1
  call psb_amx(ctxt,tt2)
  x1 = bv%get_vect()
  x2 = bg%get_vect()
  nr       = desc_a%get_local_rows() 
  eps = maxval(abs(x1(1:nr)-x2(1:nr)))
  call psb_amx(ctxt,eps)
  if (iam==0) write(*,*) 'Max diff on xGPU',eps

  ! FIXME: cache flush needed here
  call xg%set(x0)
  call xg%sync()
  call psb_barrier(ctxt)
  gt1 = psb_wtime()
  do i=1,ntests*ngpu
    call psb_spmm(done,agpu,xg,dzero,bg,desc_a,info)
     ! For timing purposes we need to make sure all threads
    ! in the device are done. 
    if ((info /= 0).or.(psb_get_errstatus()/=0)) then 
      write(0,*) 'From 2 spmm',info,i,ntests
      call psb_error()
      stop
    end if

  end do
!!$  call psb_oacc_DeviceSync()
  call psb_barrier(ctxt)
  gt2 = psb_wtime() - gt1
  call psb_amx(ctxt,gt2)
  call bg%sync()
  x1 = bv%get_vect()
  x2 = bg%get_vect()
  tnv = psb_genrm2(bv,desc_a,info)
  tng = psb_genrm2(bg,desc_a,info)
  tdot = psb_gedot(bg,bg,desc_a,info)
  write(0,*) ' bv ',tnv,' bg ',tng, ' dot ',tdot,eps,&
       & dnrm2(desc_a%get_local_rows(),x2,1),&
       & ddot(desc_a%get_local_rows(),x1,1,x2,1)
  call psb_geaxpby(-done,bg,+done,bv,desc_a,info)
  eps = psb_geamax(bv,desc_a,info)

  call psb_amx(ctxt,t2)
  eps = maxval(abs(x1(1:nr)-x2(1:nr)))
  call psb_amx(ctxt,eps)
  if (iam==0) write(*,*) 'Max diff on GPU',eps
  if (dump) then 
    write(fname,'(a,i3.3,a,i3.3,a)')'XCPU-out-',iam,'-',np,'.mtx'
    call mm_array_write(x1(1:nr),'Local part CPU',info,filename=fname)
    write(fname,'(a,i3.3,a,i3.3,a)')'XGPU-out-',iam,'-',np,'.mtx'
    call mm_array_write(x2(1:nr),'Local part GPU',info,filename=fname)
  end if

#endif
  annz     = a%get_nzeros()
  amatsize = a%sizeof()
  descsize = psb_sizeof(desc_a)
  call psb_sum(ctxt,nr)
  call psb_sum(ctxt,annz)
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)
  
  if (iam == psb_root_) then
    write(psb_out_unit,&
         & '("Matrix: ell1 ",i0)') idim
    write(psb_out_unit,&
         &'("Test on                          : ",i20," processors")') np
    write(psb_out_unit,&
         &'("Size of matrix                   : ",i20,"           ")') nr
    write(psb_out_unit,&
         &'("Number of nonzeros               : ",i20,"           ")') annz
    write(psb_out_unit,&
         &'("Memory occupation                : ",i20,"           ")') amatsize
    flops  = ntests*(2.d0*annz)
    tflops = flops
    gflops = flops * ngpu
    write(psb_out_unit,'("Storage type for    A: ",a)') a%get_fmt()
#if defined(PSB_OPENACC)||defined(PSB_HAVE_CUDA)
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
#if defined(PSB_OPENACC)||defined(PSB_HAVE_CUDA)
    write(psb_out_unit,'("Time for ",i6," products (s) (xGPU)  : ",F20.3)')&
         & ntests, tt2
    write(psb_out_unit,'("Time per product    (ms)     (xGPU)  : ",F20.3)')&
         & tt2*1.d3/(1.d0*ntests)
    write(psb_out_unit,'("MFLOPS                       (xGPU)  : ",F20.3)')&
         & tflops/1.d6

    write(psb_out_unit,'("Time for ",i6," products (s) (GPU.)  : ",F20.3)')&
         & ngpu*ntests,gt2
    write(psb_out_unit,'("Time per product    (ms)     (GPU.)  : ",F20.3)')&
         & gt2*1.d3/(1.d0*ntests*ngpu)
    write(psb_out_unit,'("MFLOPS                       (GPU.)  : ",F20.3)')&
         & gflops/1.d6
#endif
    !
    ! This computation assumes the data movement associated with CSR:
    ! it is minimal in terms of coefficients. Other formats may either move
    ! more data (padding etc.) or less data (if they can save on the indices). 
    !
    nbytes = nr*(2*psb_sizeof_dp + psb_sizeof_ip)+&
         & annz*(psb_sizeof_dp + psb_sizeof_ip)
    bdwdth = ntests*nbytes/(t2*1.d6)
    write(psb_out_unit,*)
    write(psb_out_unit,'("MBYTES/S sust. effective bandwidth  (CPU)  : ",F20.3)') bdwdth
#if defined(PSB_OPENACC)||defined(PSB_HAVE_CUDA)
    bdwdth = ngpu*ntests*nbytes/(gt2*1.d6)
    write(psb_out_unit,'("MBYTES/S sust. effective bandwidth  (GPU)  : ",F20.3)') bdwdth
!!$    bdwdth = psb_oacc_MemoryPeakBandwidth()
    write(psb_out_unit,'("MBYTES/S peak bandwidth             (GPU)  : ",F20.3)') bdwdth
#endif
    write(psb_out_unit,'("Storage type for DESC_A: ",a)') desc_a%indxmap%get_fmt()
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize

  end if
  call psb_print_timers(ctxt)

  !  
  !  cleanup storage and exit
  !
  call psb_gefree(bv,desc_a,info)
  call psb_gefree(xv,desc_a,info)
  call psb_spfree(a,desc_a,info)
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
#ifdef PSB_OPENACC
  call psb_oacc_exit()
#endif
#ifdef PSB_HAVE_CUDA
  call psb_cuda_exit()
#endif
  call psb_exit(ctxt)
  stop

9999 continue
  call psb_error(ctxt)

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ctxt,acfmt,agfmt,idim,tnd)
    type(psb_ctxt_type) :: ctxt
    character(len=*) :: agfmt, acfmt
    integer      :: idim
    logical :: tnd
    integer      :: np, iam
    integer      :: intbuf(10), ip

    call psb_info(ctxt, iam, np)

    if (iam == 0) then
      write(*,*) 'CPU side format?'
      read(psb_inp_unit,*) acfmt
      write(*,*) 'DEVICE side format?'
      read(psb_inp_unit,*) agfmt
      write(*,*) 'Size of discretization cube?'
      read(psb_inp_unit,*) idim
      write(*,*) 'Try comm/comp overlap?'
      read(psb_inp_unit,*) tnd
    endif
    call psb_bcast(ctxt,acfmt)
    call psb_bcast(ctxt,agfmt)
    call psb_bcast(ctxt,idim)
    call psb_bcast(ctxt,tnd)
    
    if (iam == 0) then
      write(psb_out_unit,'("Testing matrix       : ell1")')      
      write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4,"x",i4)')idim,idim,idim
      write(psb_out_unit,'("Number of processors : ",i0)')np
      write(psb_out_unit,'("Data distribution    : BLOCK")')
      write(psb_out_unit,'(" ")')
      write(psb_out_unit,'("Storage formats        ",a)') acfmt,' ',agfmt
      write(psb_out_unit,'("Testing overlap ND     ",l8)') tnd
    end if
    return

  end subroutine get_parms

end program pdgenmv
