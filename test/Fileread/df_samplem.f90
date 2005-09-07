program df_samplem
  use f90sparse
  use mat_dist
  use read_mat
  use partgraph
  implicit none

  ! input parameters
  character*20 :: cmethd, prec, mtrx_file, rhs_file
  character*80 :: charbuf
  integer      :: inparms(20)
  double precision ddot
  external ddot
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
  end interface   ! Local variables


  integer, parameter    :: izero=0, ione=1
  character, parameter  :: order='r'
  real(kind(1.d0)), pointer, save :: b_col(:), x_col(:), r_col(:), &
       & b_col_glob(:), x_col_glob(:), r_col_glob(:), b_glob(:,:), &
       &z(:), q(:),z1(:), xm(:,:), ym(:,:)
  integer              :: iargc, check_descr, convert_descr
  real(kind(1.d0)), parameter :: dzero = 0.d0, one = 1.d0
  real(kind(1.d0)) :: mpi_wtime, t1, t2, tprec, r_amax, b_amax,bb(1,1),&
       &lambda,scale,resmx,resmxp
  integer :: nrhs, nrow, nx1, nx2, n_row, dim,iread,ip,io,no,nmats,imat,irenum
  logical :: amroot
  external iargc, mpi_wtime
  integer bsze,overlap, nn, nprecs, novrs
  common/part/bsze,overlap
  integer, pointer :: work(:), precs(:), ovrs(:)
  ! sparse matrices
  type(d_spmat) :: a, aux_a, h
  type(d_prec) :: pre
!!$  type(d_precn) :: aprc
  ! dense matrices
  real(kind(1.d0)), pointer ::  aux_b(:,:) , aux1(:), aux2(:), vdiag(:), &
       &  aux_g(:,:), aux_x(:,:), d(:)

  ! communications data structure
  type(desc_type)    :: desc_a, desc_a_out

  ! blacs parameters
  integer            :: nprow, npcol, ictxt, iam, np, myprow, mypcol

  ! solver paramters
  integer            :: iter, itmax, ierr, itrace, ircode, ipart,&
       & methd, istopc, ml
  integer, pointer :: ierrv(:)
  character(len=5)   :: afmt
  real(kind(1.d0))   :: err, eps
  integer   iparm(20)
  real(kind(1.d0)) rparm(20)

  ! other variables
  integer            :: i,info,j
  integer            :: internal, m,ii,nnzero

  ! common area
  integer m_problem, nproc


  allocate(ierrv(6))
  ! initialize blacs
  call blacs_pinfo(iam, np)
  call blacs_get(izero, izero, ictxt)

  ! rectangular Grid,  Np x 1

  call blacs_gridinit(ictxt, order, np, ione)
  call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
  amroot = (myprow==0).and.(mypcol==0)

  !
  !  Get parameters from file
  !
  if (amroot) then 
    read(*,*) cmethd
    do i = 1, len(cmethd)
      inparms(i) = iachar(cmethd(i:i))
    end do
    call igebs2d(ictxt,'all',' ',20,1,inparms,20)

    read(*,*) afmt

    do i = 1, len(afmt)
      inparms(i) = iachar(afmt(i:i))
    end do
    call igebs2d(ictxt,'all',' ',20,1,inparms,20)
    read(*,*) ipart
    read(*,*) itmax
    read(*,*) itrace
    read(*,*) istopc
    read(*,*) irenum
    read(*,*) nprecs
    inparms(1) = ipart
    inparms(2) = itmax
    inparms(3) = itrace
    inparms(4) = irenum
    inparms(5) = nprecs
    inparms(6) = istopc
    call igebs2d(ictxt,'all',' ',6,1,inparms,20)    
!!$    write(0,*) 'Sent nprecs: ',nprecs
    allocate(precs(nprecs))
    read(*,*) (precs(i),i=1,nprecs)
    call igebs2d(ictxt,'all',' ',nprecs,1,precs,nprecs)
    read(*,*) novrs
    call igebs2d(ictxt,'all',' ',1,1,novrs,1)
!!$    write(0,*) 'Sent novrs: ',novrs    
    allocate(ovrs(novrs))
    read(*,*) (ovrs(i),i=1,novrs)
    call igebs2d(ictxt,'all',' ',novrs,1,ovrs,novrs)
    read(*,*) eps
    call dgebs2d(ictxt,'all',' ',1,1,eps,1)    
    read(*,*) nmats
    call igebs2d(ictxt,'all',' ',1,1,nmats,1)    


  else
    call igebr2d(ictxt,'a',' ',20,1,inparms,20,0,0)
    do i = 1, len(cmethd)
      cmethd(i:i) = achar(inparms(i))
    end do
    call igebr2d(ictxt,'a',' ',20,1,inparms,20,0,0)
    do i = 1, len(afmt) 
      afmt(i:i) = achar(inparms(i))
    end do

    call igebr2d(ictxt,'all',' ',6,1,inparms,20,0,0)
    ipart  = inparms(1)
    itmax  = inparms(2)
    itrace = inparms(3)
    irenum = inparms(4)
    nprecs = inparms(5)
    istopc = inparms(6)
!!$    write(0,*) 'Recvd nprecs: ',nprecs
    allocate(precs(nprecs))
    call igebr2d(ictxt,'all',' ',nprecs,1,precs,nprecs,0,0)
    call igebr2d(ictxt,'all',' ',1,1,novrs,1,0,0)
!!$    write(0,*) 'Recvd novrs: ',novrs    
    allocate(ovrs(novrs))   
    call igebr2d(ictxt,'all',' ',novrs,1,ovrs,novrs,0,0)
    call dgebr2d(ictxt,'all',' ',1,1,eps,1,0,0)
    call igebr2d(ictxt,'all',' ',1,1,nmats,1,0,0)
  endif


  do imat=1,nmats

    if (amroot) then 
      read(*,*) mtrx_file, rhs_file
!!$      write(0,*) 'Read mtx rhs : "',mtrx_file,'" "',rhs_file,'"'
!!$      do i = 1, len(mtrx_file)
!!$        inparms(i) = iachar(mtrx_file(i:i))
!!$      end do
!!$      ! broadcast parameters to all processors
!!$      call igebs2d(ictxt,'all',' ',20,1,inparms,20)
!!$      do i = 1, len(rhs_file)
!!$        inparms(i) = iachar(rhs_file(i:i))
!!$      end do
!!$      ! broadcast parameters to all processors
!!$      call igebs2d(ictxt,'all',' ',20,1,inparms,20)
!!$      write(0,*) 'Sent mtx rhs : "',mtrx_file,'" "',rhs_file,'"'
!!$    else
!!$      call igebr2d(ictxt,'a',' ',20,1,inparms,20,0,0)
!!$      do i = 1, len(mtrx_file)
!!$        mtrx_file(i:i) = achar(inparms(i))
!!$      end do
!!$      call igebr2d(ictxt,'a',' ',20,1,inparms,20,0,0)
!!$      do i = 1, len(rhs_file)
!!$        rhs_file(i:i) = achar(inparms(i))
!!$      end do
!!$      write(0,*) 'Recvd mtx rhs : "',mtrx_file,'" "',rhs_file,'"'
    endif

    call blacs_barrier(ictxt,'A')
    t1 = mpi_wtime()  
    ! read the input matrix to be processed and (possibly) the RHS 
    nrhs = 1
    nproc = nprow 

    if (amroot) then
      nullify(aux_b)
      call readmat(mtrx_file, aux_a, ictxt)
!!$      write(0,*) 'from readmat:  ',aux_a%fida,aux_a%m,':',&
!!$           &aux_a%ia2(aux_a%m+1)-1,':',aux_a%ia1(1:10)
      m_problem = aux_a%m
      call igebs2d(ictxt,'a',' ',1,1,m_problem,1)

      if (rhs_file /= 'none') then
        !  reading an rhs
        call read_rhs(rhs_file,aux_b,ictxt)
      end if

      if (associated(aux_b).and.size(aux_b,1)==m_problem) then
        ! if any rhs were present, broadcast the first one
!!$        write(0,*) 'ok, got an rhs ',aux_b(m_problem,1)
        b_col_glob =>aux_b(:,1)
      else
!!$        write(0,*) 'inventing an rhs '
        allocate(aux_b(m_problem,1), stat=ircode)
        if (ircode /= 0) then
          write(0,*) 'memory allocation failure in df_sample'
          call blacs_abort(ictxt,-1)
          stop
        endif
        b_col_glob =>aux_b(:,1)
        do i=1, m_problem
          b_col_glob(i) = 1.d0
        enddo
      endif
      call dgebs2d(ictxt,'a',' ',m_problem,1,b_col_glob,m_problem) 
    else
      call igebr2d(ictxt,'a',' ',1,1,m_problem,1,0,0)
!!$      write(0,*) 'Receiving AUX_B'
      allocate(aux_b(m_problem,1), stat=ircode)
      if (ircode /= 0) then
        write(0,*) 'Memory allocation failure in DF_SAMPLE'
        call blacs_abort(ictxt,-1)
        stop
      endif
      b_col_glob =>aux_b(:,1)
      call dgebr2d(ictxt,'A',' ',m_problem,1,b_col_glob,m_problem,0,0) 
    end if

    ! switch over different partition types
    if (ipart.eq.0) then 
      call blacs_barrier(ictxt,'A')
!!$      write(0,*) 'Partition type: BLOCK'
      call matdist(aux_a, a, part_block, ictxt, &
           & desc_a,b_col_glob,b_col,fmt=afmt)
    else  if (ipart.eq.1) then 
      call blacs_barrier(ictxt,'A')
!!$      write(0,*) 'partition type: BLK2'
      call matdist(aux_a, a, part_blk2, ictxt, &
           & desc_a,b_col_glob,b_col,fmt=afmt)
    else if (ipart.eq.2) then 
!!$      write(0,*) 'partition type: GRAPH'
      if (amroot) then 
!!$      write(0,*) 'Call BUILD',size(aux_a%ia1),size(aux_a%ia2),np
!!$        write(0,*) 'Build type: GRAPH ',aux_a%fida,&
!!$             &aux_a%m
        call build_grppart(aux_a%m,aux_a%fida,aux_a%ia1,aux_a%ia2,np)
      endif

      call distr_grppart(0,0,ictxt)

      call matdist(aux_a, a, part_graph, ictxt, &
           & desc_a,b_col_glob,b_col,fmt=afmt)
    else 
!!$      write(6,*) 'Partition type: BLOCK'
      call matdist(aux_a, a, part_block, ictxt, &
           & desc_a,b_col_glob,b_col,fmt=afmt)
    end if

    call f90_psdsall(m_problem,x_col,ierrv,desc_a)
    x_col(:) =0.0
    call f90_psdsasb(x_col,ierrv,desc_a)
    call f90_psdsall(m_problem,r_col,ierrv,desc_a)
    r_col(:) =0.0
    call f90_psdsasb(r_col,ierrv,desc_a)
    t2 = mpi_wtime() - t1


    dim=size(a%aspk)


    call dgamx2d(ictxt, 'a', ' ', ione, ione, t2, ione,&
         & t1, t1, -1, -1, -1)

    if (amroot) then
      write(6,*) 'time to Read and Partition Matrix : ',T2
    END IF

    !
    !  Prepare the preconditioning matrix. Note the availability
    !  of optional parameters
    !
    do ip=1,nprecs
      pre%prec=precs(ip)


      if (pre%prec>2) then 
        no=novrs
      else
        no=1
      endif
      do io=1, no
        if (pre%prec>2) then
          pre%n_ovr=ovrs(io)
        else
          pre%n_ovr=0
        endif
        pre%irenum = irenum
!!$        if (amroot) write(0,*) 'Preconditioner : ',&
!!$             &PRE%PREC,pre%n_ovr

!!$  do i=1,a%m
!!$     do j=a%ia2(i),a%ia2(i+1)-1
!!$        write(0,*)'a ',i,a%ia1(j),a%aspk(j)
!!$     end do
!!$  end do
!!$
!!$  write(0,*)'halo_index',desc_a%halo_index(:)
!!$  write(0,*)'ovrlap_index',desc_a%ovrlap_index(:)
!!$  write(0,*)'ovrlap_elem',desc_a%ovrlap_elem(:)

        ! Zero initial guess.
        call f90_psaxpby(0.d0,b_col,0.d0,x_col,desc_a)
        call blacs_barrier(ictxt,'All')
        t1 = mpi_wtime()
        call preconditioner(a,pre,desc_a,info)!,'f')
        tprec = mpi_wtime()-t1
        
        write(0,*) myprow,' Preconditioner Time :',TPREC,' ',&
               &pre%prec

        call DGAMX2D(ICTXT,'A',' ',IONE, IONE,TPREC,IONE,T1,T1,-1,-1,-1)
        if (amroot) then
          write(0,*) 'Preconditioner Time :',TPREC,' ',&
               &pre%prec
        endif
        if (info /= 0) then
          write(0,*) 'error in preconditioner :',info
          call blacs_abort(ictxt,-1)
          stop
        end if
        if (pre%prec>=ras2lv_) then 
          write(*,*) myprow, 'Aggregation checks: ',pre%na_f1,pre%nn_f1,pre%na_tot
          if (amroot) write(*,*) 'Output local aggregates ',pre%nlaggr(:)
        end if
        iparm = 0
        rparm = 0.d0   
        call blacs_barrier(ictxt,'all')
        t1 = mpi_wtime()
        if (cmethd.eq.'BICGSTAB') Then
          call  f90_bicgstab(a,pre,b_col,x_col,eps,desc_a,& 
               & itmax,iter,err,ierr,itrace,istop=istopc)     
!!$  ELSE IF (CMETHD.EQ.'BICG') Then
!!$    CALL  F90_BICG(A,IPREC,L,U,VDIAG,B_COL,X_COL,EPS,DESC_A,& 
!!$       & ITMAX,ITER,ERR,IERR,ITRACE)     
!!$  ELSE IF (CMETHD.EQ.'CGS') Then
!!$    CALL  F90_CGS(A,IPREC,L,U,VDIAG,B_COL,X_COL,EPS,DESC_A,& 
!!$       & ITMAX,ITER,ERR,IERR,ITRACE)     
!!$  ELSE IF (CMETHD.EQ.'BICGSTABL') Then
!!$    CALL  F90_BICGSTABL(A,IPREC,L,U,VDIAG,B_COL,X_COL,EPS,DESC_A,& 
!!$       & ITMAX,ITER,ERR,IERR,ITRACE,ML)     
        endif
        call blacs_barrier(ictxt,'all')
        t2 = mpi_wtime() - t1
        call dgamx2d(ictxt,'a',' ',ione, ione,t2,ione,t1,t1,-1,-1,-1)
        call f90_psaxpby(1.d0,b_col,0.d0,r_col,desc_A)
        call f90_psspmm(-1.d0,a,x_col,1.d0,r_col,desc_a)
        call f90_nrm2(resmx,r_col,desc_a)
        call f90_amax(resmxp,r_col,desc_a)

        if (amroot) then 
          write(6,*) 'methd iprec istopc   : ',pre%prec, pre%n_ovr, istopc
!!$    write(6,*) 'Number of iterations : ',iter
!!$    WRITE(6,*) 'Error on exit        : ',ERR
          write(6,*) 'Matrix: ',mtrx_file
          write(6,*) 'Computed solution on ',NPROW,' processors.'
          write(6,*) 'Iterations to convergence: ',iter
          write(6,*) 'Error indicator on exit:',err
          write(6,*) 'Time to Buil Prec.   : ',TPREC
          write(6,*) 'Time to Solve Matrix : ',T2
          WRITE(6,*) 'Time per iteration   : ',T2/(ITER)
          write(6,*) 'Total Time           : ',T2+TPREC
          write(6,*) 'Residual norm 2   = ',resmx
          write(6,*) 'Residual norm inf = ',resmxp
          write(6,*) 
          write(6,*) 
          
          write(8,'(a18,3(1x,i2),1x,i5,5(1x,g9.4))') mtrx_file,nprow,pre%prec,pre%n_ovr,&
               & iter,tprec,t2,t2+tprec,resmx,resmxp
        END IF
!!$        write(0,*) 'Done matrix ',imat,ip,io
        call blacs_barrier(ictxt,'all')
        call f90_psprecfree(pre,info)
!!$        write(0,*) 'Done precfree'
        call blacs_barrier(ictxt,'all')
      end do
    end do
    deallocate(aux_b)
    if (amroot) call spfree(aux_a,info)
    call f90_psdsfree(b_col, desc_a)
    call f90_psdsfree(x_col, desc_a)
    call f90_psdsfree(r_col, desc_a)
    call f90_psspfree(a, desc_a)
    call f90_psdscfree(desc_a,info)

  end do

  call blacs_gridexit(ictxt)
  call blacs_exit(0)

end program df_samplem
  




