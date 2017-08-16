!    
!                Parallel Sparse BLAS  version 3.5.1
!      (C) Copyright 2015
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: ppde.f90
!
program pdgenspmv
  use psb_base_mod
  use psb_util_mod
  implicit none

  ! input parameters
  character(len=20) :: kmethd, ptype
  character(len=5)  :: afmt
  integer(psb_ipk_) :: idim

  ! miscellaneous 
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_) :: t1, t2, tprec, flops, tflops, tflops1, tt1, tt2, bdwdth, th, th1 

  ! sparse matrix and preconditioner
  type(psb_dspmat_type) :: a, ad, ah, ah1
  type(psb_d_csr_sparse_mat) :: acsr
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! dense matrices
  type(psb_d_vect_type)  :: xv,bv, vtst,bvh, xvh
  real(psb_dpk_), allocatable :: tst(:), work(:), temp(:)
  real(psb_dpk_), allocatable :: xvc(:)[:], bvc(:)[:]
  type(event_type), allocatable :: ready[:]
  ! blacs parameters
  integer(psb_ipk_) :: ictxt, iam, np
  !  type(psb_d_csre_sparse_mat) :: acsre

  ! solver parameters
  integer(psb_ipk_) :: iter, itmax,itrace, istopc, irst, nr,nrl,ncl, lwork, nclmx
  integer(psb_ipk_) :: ip, img, nxch, p1,p2
  integer(psb_long_int_k_) :: amatsize, precsize, descsize, d2size, annz, nbytes, ahnnz
  real(psb_dpk_)   :: err, eps, err2
  integer(psb_ipk_) :: times
  integer(psb_ipk_), parameter :: iwarm=2

  logical :: indirect_remote_get
  ! other variables
  integer(psb_ipk_) :: info, i
  character(len=20)  :: name,ch_err
  character(len=40)  :: fname

  info=psb_success_


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='pde90'
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
  call get_parms(ictxt,afmt,idim, times)

  !
  !  allocate and fill in the coefficient matrix, rhs and initial guess 
  !
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  call psb_gen_pde3d(ictxt,idim,a,bv,xv,desc_a,afmt,&
       & a1,a2,a3,b1,b2,b3,c,g,info)  
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_gen_pde3d'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_geasb(bvh,desc_a,info,scratch=.true.)
  if (iam == psb_root_) write(psb_out_unit,'("Overall matrix creation time : ",es12.5)')t2
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
!!$  write(fname,'(a,i3.3,a,i3.3,a,i3.3,a)') 'testmat-',idim,'-',np,'-',iam,'.mtx'
!!$  call a%print(fname,head="psb-testing")
  
  
  !
  ! Coarrays for experiment
  !
  nrl = desc_a%get_local_rows()
  ncl = desc_a%get_local_cols()
  nclmx = ncl
  call psb_amx(ictxt,nclmx)
!!$  write(*,*) iam,'NCLMX',nclmx,ncl
  allocate(xvc(nclmx)[*])
  allocate(bvc(nclmx)[*])
  allocate(ready[*])
  allocate(temp(nclmx))
  xvc(1:ncl) = done *(/(i,i=1,ncl)/)
  bvc(:) = dzero
  call xv%set(xvc(1:ncl))


  do i=1,iwarm
    call psb_spmm(done,a,xv,dzero,bv,desc_a,info,'n')
  end do


  ! FIXME: cache flush needed here
  do i=1,iwarm
    call psb_spmm(done,a,xv,dzero,bv,desc_a,info,'n')
  end do
  call psb_barrier(ictxt)
  tt1 = psb_wtime()
  do i=1,times 
    call psb_spmm(done,a,xv,dzero,bv,desc_a,info,'n')
  end do
  call psb_barrier(ictxt)
  th = psb_wtime() - tt1
  call psb_amx(ictxt,th)

  indirect_remote_get = .false.
  if (.true.) then 
    associate(xchg => desc_a%halo_xch)
      ! FIXME: cache flush needed here
      nxch = size(xchg%prcs_xch)
!!$      if (this_image()==2) write(0,*) this_image(),nxch,nrl,ncl,' Exchanging with ',xchg%prcs_xch+1
      temp = 0.d0
      xvc(nrl+1:) = 0.d0
      do i=1,iwarm
        ! Sync images
        sync images(xchg%prcs_xch+1)
        if (indirect_remote_get) then 
          do ip = 1, nxch
            img = xchg%prcs_xch(ip) + 1
            p1  = xchg%loc_rcv_bnd(ip)
            p2  = xchg%loc_rcv_bnd(ip+1)-1
!!$            if (this_image()==2) write(0,*) this_image(),'Boundaries ',&
!!$                 & p1,p2,' :',xchg%loc_rcv_idx(p1:p2),':',xchg%rmt_rcv_idx(p1:p2)
!!$        xvc(xchg%loc_rcv_idx(p1:p2)) = xvc(xchg%rmt_rcv_idx(p1:p2))[img]
            
            temp(p1:p2) = xvc(xchg%rmt_rcv_idx(p1:p2))[img]
            xvc(xchg%loc_rcv_idx(p1:p2)) = temp(p1:p2)
!!$            if (this_image()==2) write(0,*) this_image(),' :x: ',ip,' : ',&
!!$                 &xvc(xchg%loc_rcv_idx(p1:p2)),' : ',xv%v%v(xchg%loc_rcv_idx(p1:p2))
          end do
!!$        if (this_image()==2) write(0,*) this_image(),' :x: ',&
!!$             &xvc(nrl+1:ncl),' : ',xv%v%v(nrl+1:ncl)
          call  a%csmv(done,xvc,dzero,bvc,info)
        else
          do ip = 1, nxch
            img = xchg%prcs_xch(ip) + 1
            p1  = xchg%loc_snd_bnd(ip)
            p2  = xchg%loc_snd_bnd(ip+1)-1
            !xvc(xchg%loc_rcv_idx(p1:p2)) = xvc(xchg%rmt_rcv_idx(p1:p2))[img]
!!$            write(0,*) this_image(),'Boundaries ',&
!!$                 & p1,p2,' :',xchg%loc_snd_idx(p1:p2),':',xchg%rmt_snd_idx(p1:p2)            
            temp(p1:p2) = xvc(xchg%loc_snd_idx(p1:p2))
            xvc(xchg%rmt_snd_idx(p1:p2))[img] = temp(p1:p2)
          end do
          sync images(xchg%prcs_xch+1)
          call  a%csmv(done,xvc,dzero,bvc,info)
          
        end if
      end do
      call psb_barrier(ictxt)
      tt1 = psb_wtime()
      do i=1,times
        ! Sync images
!        sync images(xchg%prcs_xch+1)
!!$        do ip=1,nxch
!!$          img = xchg%prcs_xch(ip) + 1
!!$          event post(ready[img])
!!$        end do
!!$        event wait(ready, until_count=nxch)
        if (indirect_remote_get) then         
          do ip = 1, nxch
            img = xchg%prcs_xch(ip) + 1
            sync images (img)
            p1  = xchg%loc_rcv_bnd(ip)
            p2  = xchg%loc_rcv_bnd(ip+1)-1
            !xvc(xchg%loc_rcv_idx(p1:p2)) = xvc(xchg%rmt_rcv_idx(p1:p2))[img]         
            temp(p1:p2) = xvc(xchg%rmt_rcv_idx(p1:p2))[img]
            xvc(xchg%loc_rcv_idx(p1:p2)) = temp(p1:p2)
          end do
        else
          do ip = 1, nxch
            img = xchg%prcs_xch(ip) + 1
            p1  = xchg%loc_snd_bnd(ip)
            p2  = xchg%loc_snd_bnd(ip+1)-1
            !xvc(xchg%loc_rcv_idx(p1:p2)) = xvc(xchg%rmt_rcv_idx(p1:p2))[img]         
            temp(p1:p2) = xvc(xchg%loc_snd_idx(p1:p2))
            xvc(xchg%rmt_snd_idx(p1:p2))[img] = temp(p1:p2)
          end do
          sync images(xchg%prcs_xch+1)
        end if
        call  a%csmv(done,xvc,dzero,bvc,info)
      end do
      call psb_barrier(ictxt)
      tt2 = psb_wtime() - tt1
      call psb_amx(ictxt,tt2)
    end associate
    call bvh%set(bvc(1:ncl))
    call xvh%set(xvc(1:ncl))
  endif
!!$  do i=2,2
!!$    sync all
!!$    if (i==this_image()) then 
!!$      write(0,*) this_image(),' Xdiff ',xv%v%v(nrl+1:ncl) -xvc(nrl+1:ncl)
!!$      write(0,*) this_image(),' X1 ',xv%v%v(nrl+1:ncl) 
!!$      write(0,*) this_image(),' X2 ',xvc(nrl+1:ncl)
!!$    end if
!!$  end do
  
  call psb_geaxpby(-done,bv,done,bvh,desc_a,info)
  err = psb_genrm2(bvh,desc_a,info)
  call psb_geaxpby(-done,xv,done,xvh,desc_a,info)
  err2 = psb_genrm2(xvh,desc_a,info)

  nr       = desc_a%get_global_rows() 
  annz     = a%get_nzeros()
  amatsize = a%sizeof()
  descsize = psb_sizeof(desc_a)
  call psb_sum(ictxt,annz)
  call psb_sum(ictxt,ahnnz)
  call psb_sum(ictxt,amatsize)
  call psb_sum(ictxt,descsize)

  if (iam == psb_root_) then
    flops   = 2.d0*times*annz
    tflops  = flops

    write(psb_out_unit,'("Matrix: ell1 ",i0)') idim
    write(psb_out_unit,'("Test on                          : ",i20," processors")') np
    write(psb_out_unit,'("Size of matrix                   : ",i20,"           ")') nr
    write(psb_out_unit,'("Number of nonzeros               : ",i20,"           ")') annz
    write(psb_out_unit,'("Memory occupation                : ",i20,"           ")') amatsize
    write(psb_out_unit,'("Number of flops (",i0," prod)        : ",F20.0,"           ")') times,flops
    flops   = flops / (th)
    tflops  = tflops / (tt2)

    write(psb_out_unit,'("Time for ",i0," products (s)         : ",F20.3)')times, th
    write(psb_out_unit,'("Time per product    (ms)         : ",F20.3)') th*1.d3/(1.d0*times)
    write(psb_out_unit,'("MFLOPS                           : ",F20.3)') flops/1.d6

    write(psb_out_unit,'("Time for ",i0," products (s) nw      : ",F20.3)')times, tt2
    write(psb_out_unit,'("Time per product    (ms)    nw   : ",F20.3)') tt2*1.d3/(1.d0*times)
    write(psb_out_unit,'("MFLOPS                           : ",F20.3)') tflops/1.d6

    write(psb_out_unit,'("Difference                       : ",E20.12)') err
    write(psb_out_unit,'("Difference                       : ",E20.12)') err2

    !
    ! This computation is valid for CSR
    !
    nbytes = nr*(2*psb_sizeof_dp + psb_sizeof_int)+&
         & annz*(psb_sizeof_dp + psb_sizeof_int)
    bdwdth = times*nbytes/(th*1.d6)
    write(psb_out_unit,*)
    write(psb_out_unit,'("MBYTES/S                         : ",F20.3)') bdwdth

    write(psb_out_unit,'("Storage type for DESC_A: ",a)') desc_a%get_fmt()
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize

  end if


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

  call psb_exit(ictxt)
  stop

9999 call psb_error(ictxt)
  
  stop

contains
  !
  ! get iteration parameters from standard input
  !
  subroutine  get_parms(ictxt,afmt,idim, times)
    integer(psb_ipk_) :: ictxt
    character(len=*) :: afmt
    integer(psb_ipk_) :: idim, times
    integer(psb_ipk_) :: np, iam
    integer(psb_ipk_) :: intbuf(10), ip

    call psb_info(ictxt, iam, np)

    if (iam == 0) then
      read(psb_inp_unit,*) afmt
      read(psb_inp_unit,*) idim
      read(psb_inp_unit,*) times
    endif
    call psb_bcast(ictxt,afmt)
    call psb_bcast(ictxt,idim)
    call psb_bcast(ictxt,times)

    if (iam == 0) then
      write(psb_out_unit,'("Testing matrix       : ell1")')      
      write(psb_out_unit,'("Grid dimensions      : ",i4,"x",i4,"x",i4)')idim,idim,idim
      write(psb_out_unit,'("Number of processors : ",i0)')np
      write(psb_out_unit,'("Data distribution    : BLOCK")')
      write(psb_out_unit,'(" ")')
    end if
    return

  end subroutine get_parms
  !
  !  print an error message 
  !  
  subroutine pr_usage(iout)
    integer(psb_ipk_) :: iout
    write(iout,*)'incorrect parameter(s) found'
    write(iout,*)' usage:  pde90 methd prec dim &
         &[istop itmax itrace]'  
    write(iout,*)' where:'
    write(iout,*)'     methd:    cgstab cgs rgmres bicgstabl' 
    write(iout,*)'     prec :    bjac diag none'
    write(iout,*)'     dim       number of points along each axis'
    write(iout,*)'               the size of the resulting linear '
    write(iout,*)'               system is dim**3'
    write(iout,*)'     istop     stopping criterion  1, 2  '
    write(iout,*)'     itmax     maximum number of iterations [500] '
    write(iout,*)'     itrace    <=0  (no tracing, default) or '  
    write(iout,*)'               >= 1 do tracing every itrace'
    write(iout,*)'               iterations ' 
  end subroutine pr_usage

  !
  ! functions parametrizing the differential equation 
  !  
  function b1(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) :: b1
    real(psb_dpk_), intent(in) :: x,y,z
    b1=1.d0/sqrt(3.d0)
  end function b1
  function b2(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  b2
    real(psb_dpk_), intent(in) :: x,y,z
    b2=1.d0/sqrt(3.d0)
  end function b2
  function b3(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  b3
    real(psb_dpk_), intent(in) :: x,y,z      
    b3=1.d0/sqrt(3.d0)
  end function b3
  function c(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  c
    real(psb_dpk_), intent(in) :: x,y,z      
    c=0.d0
  end function c
  function a1(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a1   
    real(psb_dpk_), intent(in) :: x,y,z
    a1=1.d0/80
  end function a1
  function a2(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a2
    real(psb_dpk_), intent(in) :: x,y,z
    a2=1.d0/80
  end function a2
  function a3(x,y,z)
    use psb_base_mod, only : psb_dpk_
    real(psb_dpk_) ::  a3
    real(psb_dpk_), intent(in) :: x,y,z
    a3=1.d0/80
  end function a3
  function g(x,y,z)
    use psb_base_mod, only : psb_dpk_, done
    real(psb_dpk_) ::  g
    real(psb_dpk_), intent(in) :: x,y,z
    g = dzero
    if (x == done) then
      g = done
    else if (x == dzero) then 
      g = exp(y**2-z**2)
    end if
  end function g
end program
