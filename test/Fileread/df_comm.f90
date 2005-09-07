program df_comm
  use f90sparse
  use mat_dist
  use read_mat
  use partgraph
  use comminfo
  implicit none

  ! input parameters
  character(len=20) :: cmethd, prec
  character(len=40) :: mtrx_file
  character(len=20) :: mtrx_name, out_file, tmpstr
  character(len=200) :: charbuf
  integer      :: inparms(20)
  double precision ddot
  external ddot


  integer, parameter    :: izero=0, ione=1
  character, parameter  :: order='r'
  real(kind(1.d0)), pointer, save :: b_col(:), x_col(:), r_col(:), &
       & b_col_glob(:), x_col_glob(:), r_col_glob(:), b_glob(:,:), &
       &z(:), q(:),z1(:), xm(:,:), ym(:,:)
  integer              :: iargc, check_descr, convert_descr
  real(kind(1.d0)), parameter :: dzero = 0.d0, one = 1.d0
  real(kind(1.d0)) :: mpi_wtime, t1, t2, t3, tprec, tslv, ttot, &
       &r_amax, b_amax,bb(1,1), lambda,scale,resmx,resmxp, omega
  integer :: nrhs, nrow, nx1, nx2, n_row, n_col, dim,iread,ip,io,no,nmats,&
       & imat,irenum, igsmth, matop
  logical :: amroot
  external iargc, mpi_wtime
  integer bsze,overlap, nn, nprecs, novrs
  common/part/bsze,overlap
  integer, pointer :: work(:), precs(:), ovrs(:), comm_info(:,:)
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
       & methd, istopc, irst
  integer, pointer :: ierrv(:), ivg(:)
  character(len=5)   :: afmt
  real(kind(1.d0))   :: err, eps
  integer   iparm(20)
  real(kind(1.d0)) rparm(20)

  ! other variables
  integer            :: i,info,j, ntryslv
  integer            :: internal, m,ii,nnzero, itryslv
  integer, parameter :: ncsw=4, ntcs=4
  real(kind(1.d0)), pointer :: tsthal(:,:)
  real(kind(1.d0))   :: tswmpi(ntcs,ncsw),tswsyn(ntcs,ncsw),tswsdrv(ntcs,ncsw)
  ! common area
  integer m_problem, nproc, nc


  allocate(ierrv(6))
  ! initialize blacs
  call blacs_pinfo(iam, np)
  call blacs_get(izero, izero, ictxt)

  ! rectangular Grid,  Np x 1

  call blacs_gridinit(ictxt, order, np, ione)
  call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
  amroot = (myprow==0).and.(mypcol==0)
  afmt='CSR'
  !
  !  Get parameters from file
  !
  if (amroot) then 
     read(*,*) ipart
     call igebs2d(ictxt,'all',' ',1,1,ipart,1)    
     read(*,*) nmats
     call igebs2d(ictxt,'all',' ',1,1,nmats,1)    


  else
     call igebr2d(ictxt,'all',' ',1,1,ipart,1,0,0)
     call igebr2d(ictxt,'all',' ',1,1,nmats,1,0,0)
  endif

  do imat=1,nmats

     if (amroot) then 
!!$      read(*,*) mtrx_file,rhs_file 
        read(*,'(a)') charbuf
        charbuf=adjustl(charbuf)
        i=index(charbuf," ")
        mtrx_file=charbuf(1:i-1)
        i=index(mtrx_file,"/",back=.true.)
        mtrx_name=trim(mtrx_file(i+1:))
        write(0,*) 'Read mtx rhs : "', mtrx_file
        i=index(mtrx_file,'.mtx')
        write(tmpstr,'("_",i1,"_",i2.2,".out")')ipart,np
        out_file=mtrx_file(1:i-1)//trim(tmpstr)
        write(0,*)' out file is ', out_file
        open(2,file=out_file,action='write')
        write(2,'("Matrix:  ",a20)')mtrx_name
        write(2,'("Number of processes:  ",i2)')np
        write(2,'("Partitioning:  ",i1)')ipart
        write(2,'(" ")')
        write(2,'(" ")')
        close(2)
     endif

     call blacs_barrier(ictxt,'A')
     t1 = mpi_wtime()  
     ! read the input matrix to be processed and (possibly) the RHS 
     nrhs = 1
     nproc = nprow 

     if (amroot) then
        nullify(aux_b)
        call readmat(mtrx_file, aux_a, ictxt)

        m_problem = aux_a%m
        call igebs2d(ictxt,'a',' ',1,1,m_problem,1)
        
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
        
        call dgebs2d(ictxt,'a',' ',m_problem,1,b_col_glob,m_problem) 
     else
        call igebr2d(ictxt,'a',' ',1,1,m_problem,1,0,0)

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
        write(0,*) 'Partition type: BLOCK'
        allocate(ivg(m_problem))
        call bld_partblock(m_problem,np,ivg)
        call matdist(aux_a, a, ivg, ictxt, &
             & desc_a,b_col_glob,b_col,fmt=afmt)
     else  if (ipart.eq.1) then 
        call blacs_barrier(ictxt,'A')
        write(0,*) 'partition type: BLK2'
        allocate(ivg(m_problem))
        call bld_partblk2(m_problem,np,ivg)
        call matdist(aux_a, a, ivg, ictxt, &
             & desc_a,b_col_glob,b_col,fmt=afmt)
     else if (ipart.eq.2) then 
        write(0,*) 'partition type: GRAPH'
        if (amroot) then 
           call build_grppart(aux_a%m,aux_a%fida,aux_a%ia1,aux_a%ia2,np)
        endif
        call distr_grppart(0,0,ictxt)
        if (.false.) then 
           call matdist(aux_a, a, part_graph, ictxt, &
                & desc_a,b_col_glob,b_col,fmt=afmt)
        else
           call getv_grppart(ivg)
           call matdist(aux_a, a, ivg, ictxt, &
                & desc_a,b_col_glob,b_col,fmt=afmt)
        endif
     end if

     if(amroot) then
        allocate(comm_info(np,np))
        comm_info(:,:)=0
     end if
     call blacs_barrier(ictxt,'all')
     write(0,'("Getting communication info")')
     call get_comminfo(ictxt,desc_a,comm_info)
     if(amroot) then
        open(2,file=out_file,action='write',position='append')
        write(2,'("Exchange table:")')
        do i=1,np
           write(2,*)'Row ',i,' : ',comm_info(i,:)
        end do
        write(2,'(" ")')
        write(2,'(" ")')
        close(2)
     end if


     n_row = desc_a%matrix_data(n_row_)
     n_col = desc_a%matrix_data(n_col_)
     write(0,'("Allocating vectors")')
     call f90_psdsall(m_problem,ncsw,tsthal,ierrv,desc_a)
     forall (i=1:n_row)
        forall (j=1:ncsw)
           tsthal(i,j) = j * desc_a%loc_to_glob(i)
        end forall
     end forall
     tsthal(n_row+1:,:) = -1.d0

     tswmpi  = 1.d200
     tswsyn  = 1.d200
     tswsdrv = 1.d200

     write(0,'("Cycling")')
     do nc=1, ncsw
        do i=1, ntcs
           tsthal(n_row+1:,:) = -1.d0
           t1=mpi_wtime()
           call f90_pshalo(tsthal(:,1:nc),desc_a,trans='T',mode=SWAP_MPI)
           t2=mpi_wtime()-t1
           call dgamx2d(ictxt,'a',' ',ione, ione,t2,ione,t1,t1,-1,-1,-1)
           tswmpi(i,nc) = t2
           !! Add correctness check
           tsthal(n_row+1:,:) = -1.d0
           t1=mpi_wtime()
           call f90_pshalo(tsthal(:,1:nc),desc_a,trans='T',mode=SWAP_SYNC)
           t2=mpi_wtime()-t1
           call dgamx2d(ictxt,'a',' ',ione, ione,t2,ione,t1,t1,-1,-1,-1)
           tswsyn(i,nc) = t2
           !! Add correctness check
           tsthal(n_row+1:,:) = -1.d0
           t1=mpi_wtime()
           call f90_pshalo(tsthal(:,1:nc),desc_a,trans='T',mode=IOR(SWAP_SEND,SWAP_RECV))
           t2=mpi_wtime()-t1
           call dgamx2d(ictxt,'a',' ',ione, ione,t2,ione,t1,t1,-1,-1,-1)
           tswsdrv(i,nc) = t2
           !! Add correctness check
        end do
     end do

     if (amroot) then 
        open(2,file=out_file,action='write',position='append')
        do nc=1, ncsw
           write(*,'(a18,1x,a4,1(1x,i2),1x,i5,5(1x,g9.4))') mtrx_name,'MPI',&
                & nprow,nc,minval(tswmpi(:,nc)),maxval(tswmpi(:,nc)),&
                & sum(tswmpi(:,nc))/ntcs
           write(2,'(a18,1x,a4,1(1x,i2),1x,i5,5(1x,g9.4))') mtrx_name,'MPI',&
                & nprow,nc,minval(tswmpi(:,nc)),maxval(tswmpi(:,nc)),&
                & sum(tswmpi(:,nc))/ntcs
           write(*,'(a18,1x,a4,1(1x,i2),1x,i5,5(1x,g9.4))') mtrx_name,'SYNC',&
                & nprow,nc,minval(tswsyn(:,nc)),maxval(tswsyn(:,nc)),&
                & sum(tswsyn(:,nc))/ntcs
           write(2,'(a18,1x,a4,1(1x,i2),1x,i5,5(1x,g9.4))') mtrx_name,'SYNC',&
                & nprow,nc,minval(tswsyn(:,nc)),maxval(tswsyn(:,nc)),&
                & sum(tswsyn(:,nc))/ntcs
           write(*,'(a18,1x,a4,1(1x,i2),1x,i5,5(1x,g9.4))') mtrx_name,'SDRV',&
                & nprow,nc,minval(tswsdrv(:,nc)),maxval(tswsdrv(:,nc)),&
                & sum(tswsdrv(:,nc))/ntcs
           write(2,'(a18,1x,a4,1(1x,i2),1x,i5,5(1x,g9.4))') mtrx_name,'SDRV',&
                & nprow,nc,minval(tswsdrv(:,nc)),maxval(tswsdrv(:,nc)),&
                & sum(tswsdrv(:,nc))/ntcs
        end do
        close(2)
     end if
     call f90_psdsfree(tsthal, desc_a)
     call f90_psdsfree(b_col, desc_a)
     call f90_psspfree(a, desc_a)
     call f90_psdscfree(desc_a,info)

  end do
  deallocate(ovrs,precs,ierrv, stat=info)
  write(0,*) 'Info from deallocate ',info
  call blacs_gridexit(ictxt)
  call blacs_exit(0)

end program df_comm




