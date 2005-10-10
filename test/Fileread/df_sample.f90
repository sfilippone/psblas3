program df_sample
  use psb_all_mod
  use mat_dist
  use read_mat
  use partgraph
  use getp
  implicit none

  ! input parameters
  character*20 :: cmethd, prec, mtrx_file, rhs_file
  character*80 :: charbuf

  interface 
    !   .....user passed subroutine.....
    subroutine part_block(global_indx,n,np,pv,nv)
      implicit none
      integer, intent(in)  :: global_indx, n, np
      integer, intent(out) :: nv
      integer, intent(out) :: pv(*) 
    end subroutine part_block
  end interface   ! local variables
  interface 
    !   .....user passed subroutine.....
    subroutine part_blk2(global_indx,n,np,pv,nv)
      implicit none
      integer, intent(in)  :: global_indx, n, np
      integer, intent(out) :: nv
      integer, intent(out) :: pv(*) 
    end subroutine part_blk2
  end interface   ! local variables


  character, parameter  :: order='r'
  real(kind(1.d0)), pointer, save :: b_col(:), x_col(:), r_col(:), &
       & b_col_glob(:), x_col_glob(:), r_col_glob(:)

  integer              :: iargc
  real(kind(1.d0)) :: mpi_wtime, t1, t2, tprec, r_amax, b_amax,bb(1,1),&
       &scale,resmx,resmxp
  integer :: nrhs, nrow, n_row, dim, nv
  logical :: amroot
  integer, pointer :: ivg(:), ipv(:)

  external iargc, mpi_wtime

  ! sparse matrices
  type(psb_dspmat_type) :: a, aux_a
  type(psb_dprec_type)  :: pre

  ! dense matrices
  real(kind(1.d0)), pointer ::  aux_b(:,:), vdiag(:), d(:)


  ! communications data structure
  type(psb_desc_type):: desc_a, desc_a_out

  ! blacs parameters
  integer            :: nprow, npcol, ictxt, iam, np, myprow, mypcol

  ! solver paramters
  integer            :: iter, itmax, ierr, itrace, ircode, ipart,&
       & methd, istopc, ml, iprec, novr

  character(len=5)   :: afmt
  character(len=20)  :: name
  real(kind(1.d0))   :: err, eps
  integer   :: iparm(20)

  ! other variables
  integer            :: i,info,j
  integer            :: internal, m,ii,nnzero

  integer            :: m_problem, nproc
  

  ! initialize blacs
  call blacs_pinfo(iam, np)
  call blacs_get(izero, izero, ictxt)

  ! rectangular grid,  np x 1

  call blacs_gridinit(ictxt, order, np, ione)
  call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
  amroot = (myprow==0).and.(mypcol==0)

  name='df_sample'
  info=0
  !
  !  get parameters
  !
  call get_parms(ictxt,mtrx_file,rhs_file,cmethd,prec,&
       & ipart,afmt,istopc,itmax,itrace,novr,iprec,eps)

  call blacs_barrier(ictxt,'a')
  t1 = mpi_wtime()  
  ! read the input matrix to be processed and (possibly) the rhs 
  nrhs = 1
  nproc = nprow 

  if (amroot) then
    nullify(aux_b)
    call readmat(mtrx_file, aux_a, ictxt)

    m_problem = aux_a%m
    call igebs2d(ictxt,'a',' ',1,1,m_problem,1)

    if(rhs_file /= 'NONE') then
       !  reading an rhs
       call read_rhs(rhs_file,aux_b,ictxt)
    end if

    if (associated(aux_b).and.size(aux_b,1)==m_problem) then
      ! if any rhs were present, broadcast the first one
      write(0,'("Ok, got an rhs ")')
      b_col_glob =>aux_b(:,1)
    else
      write(0,'("Generating an rhs ")')
      allocate(aux_b(m_problem,1), stat=ircode)
      if (ircode /= 0) then
         call psb_errpush(4000,name)
         goto 9999
      endif

      b_col_glob => aux_b(:,1)
      do i=1, m_problem
         b_col_glob(i) = 1.d0
      enddo      
    endif
    call dgebs2d(ictxt,'a',' ',m_problem,1,b_col_glob,m_problem) 
  else
    call igebr2d(ictxt,'a',' ',1,1,m_problem,1,0,0)
    allocate(aux_b(m_problem,1), stat=ircode)
    if (ircode /= 0) then
       call psb_errpush(4000,name)
       goto 9999
    endif
    b_col_glob =>aux_b(:,1)
    call dgebr2d(ictxt,'a',' ',m_problem,1,b_col_glob,m_problem,0,0) 
  end if

  ! switch over different partition types
  if (ipart.eq.0) then 
     call blacs_barrier(ictxt,'a')
     write(6,'("Partition type: block")')
     allocate(ivg(m_problem),ipv(np))
     do i=1,m_problem
        call part_block(i,m_problem,np,ipv,nv)
        ivg(i) = ipv(1)
     enddo
     call matdist(aux_a, a, ivg, ictxt, &
          & desc_a,b_col_glob,b_col,info,fmt=afmt)
  else  if (ipart.eq.1) then 
    call blacs_barrier(ictxt,'a')
    write(6,'("Partition type: blk2")')
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
       call part_blk2(i,m_problem,np,ipv,nv)
       ivg(i) = ipv(1)
    enddo
    call matdist(aux_a, a, ivg, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt)
  else if (ipart.eq.2) then 
    write(6,'("Partition type: graph")')
    if (amroot) then 
      write(0,'("Build type: graph")')
      call build_grppart(aux_a%m,aux_a%fida,aux_a%ia1,aux_a%ia2,np)
    endif
    call blacs_barrier(ictxt,'a')
    call distr_grppart(0,0,ictxt)
    call getv_grppart(ivg)
    call matdist(aux_a, a, ivg, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt)
  else 
    write(6,'("Partition type: block")')
    call matdist(aux_a, a, part_block, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt)
  end if
  
  call psb_alloc(m_problem,x_col,desc_a,info)
  x_col(:) =0.0
  call psb_asb(x_col,desc_a,info)
  call psb_alloc(m_problem,r_col,desc_a,info)
  r_col(:) =0.0
  call psb_asb(r_col,desc_a,info)
  t2 = mpi_wtime() - t1
  
  
  call dgamx2d(ictxt, 'a', ' ', ione, ione, t2, ione,&
       & t1, t1, -1, -1, -1)
  
  if (amroot) then
     write(6,'("Time to read and partition matrix : ",es10.4)')t2
  end if
  
  !
  !  prepare the preconditioning matrix. note the availability
  !  of optional parameters
  !

  if (amroot) write(6,'("Preconditioner : ",a)')prec(1:6)
  
  ! zero initial guess.
  select case(iprec)
  case(noprec_)
    call psb_precset(pre,'noprec')
  case(diagsc_)             
    call psb_precset(pre,'diagsc')
  case(bja_)             
    call psb_precset(pre,'ilu')
  case(asm_)             
    call psb_precset(pre,'asm',iv=(/novr,halo_,sum_/))
  case(ash_)             
    call psb_precset(pre,'asm',iv=(/novr,nohalo_,sum_/))
  case(ras_)             
    call psb_precset(pre,'asm',iv=(/novr,halo_,none_/))
  case(rash_)             
    call psb_precset(pre,'asm',iv=(/novr,nohalo_,none_/))
  end select

  ! building the preconditioner
  t1 = mpi_wtime()
  call psb_precbld(a,pre,desc_a,info)
  tprec = mpi_wtime()-t1
  
  
  call dgamx2d(ictxt,'a',' ',ione, ione,tprec,ione,t1,t1,-1,-1,-1)
  
  write(6,'("Preconditioner time: ",es10.4)')tprec
  if (info /= 0) then
    write(0,*) 'error in preconditioner :',info
    call blacs_abort(ictxt,-1)
    stop
  end if
  
  iparm = 0
  call blacs_barrier(ictxt,'all')
  t1 = mpi_wtime()
  if (cmethd.eq.'BICGSTAB') then
    call  psb_bicgstab(a,pre,b_col,x_col,eps,desc_a,info,& 
       & itmax,iter,err,itrace,istop=istopc)     
  else if (cmethd.eq.'BICG') then
    call  psb_bicg(a,pre,b_col,x_col,eps,desc_a,info,& 
       & itmax,iter,err,itrace)     
  else if (cmethd.eq.'CGS') then
    call  psb_cgs(a,pre,b_col,x_col,eps,desc_a,info,& 
       & itmax,iter,err,itrace)     
  else if (cmethd.eq.'BICGSTABL') then
    call  psb_bicgstabl(a,pre,b_col,x_col,eps,desc_a,info,& 
       & itmax,iter,err,ierr,itrace,ml)     
  endif
  call blacs_barrier(ictxt,'all')
  t2 = mpi_wtime() - t1
  write(0,*)'Calling gamx2d'
  call dgamx2d(ictxt,'a',' ',ione, ione,t2,ione,t1,t1,-1,-1,-1)
  write(0,*)'Calling axpby'
  call psb_axpby(1.d0,b_col,0.d0,r_col,desc_a,info)
  write(0,*)'Calling spmm'
  call psb_spmm(-1.d0,a,x_col,1.d0,r_col,desc_a,info)
  write(0,*)'Calling nrm2'
  call psb_nrm2(resmx,r_col,desc_a,info)
  write(0,*)'Calling amax'
  call psb_amax(resmxp,r_col,desc_a,info)

!!$  iter=iparm(5)
!!$  err = rparm(2)
  if (amroot) then 
!    call psb_prec_descr(6,pre)
    write(6,'("Matrix: ",a)')mtrx_file
    write(6,'("Computed solution on ",i4," processors")')nprow
    write(6,'("Iterations to convergence: ",i)')iter
    write(6,'("Error indicator on exit: ",f7.2)')err
    write(6,'("Time to buil prec.   : ",es10.4)')tprec
    write(6,'("Time to solve matrix : ",es10.4)')t2
    write(6,'("Time per iteration   : ",es10.4)')t2/(iter)
    write(6,'("Total time           : ",es10.4)')t2+tprec
    write(6,'("Residual norm 2   = ",f7.2)')resmx
    write(6,'("Residual norm inf = ",f7.2)')resmxp
  end if

  allocate(x_col_glob(m_problem),r_col_glob(m_problem),stat=ierr)
  if (ierr.ne.0) then 
    write(0,*) 'allocation error: no data collection'
  else
    call psb_gather(x_col_glob,x_col,desc_a,info,iroot=0)
    call psb_gather(r_col_glob,r_col,desc_a,info,iroot=0)
    if (amroot) then
      write(0,*) 'Saving x on file'
      write(20,*) 'matrix: ',mtrx_file
      write(20,*) 'computed solution on ',nprow,' processors.'
      write(20,*) 'iterations to convergence: ',iter
      write(20,*) 'error indicator (infinity norm) on exit:', &
           & ' ||r||/(||a||||x||+||b||) = ',err
      write(20,*) 'max residual = ',resmx, resmxp
      do i=1,m_problem
        write(20,998) i,x_col_glob(i),r_col_glob(i),b_col_glob(i)
      enddo
    end if
  end if
998 format(i8,4(2x,g20.14))
993 format(i6,4(1x,e12.6))

  
  call psb_free(b_col, desc_a,info)
  call psb_free(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call psb_precfree(pre,info)
  call psb_dscfree(desc_a,info)

9999 continue
  if(info /= 0) then
     call psb_error(ictxt)
     call blacs_gridexit(ictxt)
     call blacs_exit(0)
  else
     call blacs_gridexit(ictxt)
     call blacs_exit(0)
  end if

  stop
  
end program df_sample
  




