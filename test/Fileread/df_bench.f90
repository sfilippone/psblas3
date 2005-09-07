program df_bench
  use f90sparse
  use mat_dist
  use read_mat
  use partgraph
  use errormod
  implicit none

  ! input parameters
  character(len=20) :: cmethd, prec
  character(len=40) :: mtrx_file, rhs_file
  character(len=20) :: mtrx_name, name, ch_err
  character(len=10) :: ptype
  character(len=200) :: charbuf
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
  real(kind(1.d0)) :: mpi_wtime, t1, t2, t3, tprec, tslv, ttot, &
       &r_amax, b_amax,bb(1,1), lambda,scale,resmx,resmxp, omega
  integer :: nrhs, nrow, nx1, nx2, n_row, n_col, dim,iread,ip,io,no,nmats,&
       & imat,irenum, igsmth, matop, jacswp
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
       & methd, istopc, irst,nv
  integer, pointer :: ivg(:), ipv(:)
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


  ! initialize blacs
  call blacs_pinfo(iam, np)
  call blacs_get(izero, izero, ictxt)

  ! rectangular Grid,  Np x 1

  call blacs_gridinit(ictxt, order, np, ione)
  call blacs_gridinfo(ictxt, nprow, npcol, myprow, mypcol)
  amroot = (myprow==0).and.(mypcol==0)


  info=0
  name='df_bench'
  call psb_set_errverbosity(2)
  call psb_set_erraction(0)
  !  call psb_erractionsave(err_act)
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
     read(*,*) irst
     read(*,*) irenum
     read(*,*) ntryslv
     read(*,*) nprecs
     inparms(1) = ipart
     inparms(2) = itmax
     inparms(3) = itrace
     inparms(4) = irenum
     inparms(5) = ntryslv
     inparms(6) = nprecs
     inparms(7) = istopc
     inparms(8) = irst
     call igebs2d(ictxt,'all',' ',8,1,inparms,20)    
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
    read(*,*) omega
    call dgebs2d(ictxt,'all',' ',1,1,omega,1)    
    read(*,*) igsmth
    call igebs2d(ictxt,'all',' ',1,1,igsmth,1)    
    read(*,*) matop
    call igebs2d(ictxt,'all',' ',1,1,matop,1)    
    read(*,*) jacswp
    call igebs2d(ictxt,'all',' ',1,1,jacswp,1)    
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

     call igebr2d(ictxt,'all',' ',8,1,inparms,20,0,0)
     ipart  = inparms(1)
     itmax  = inparms(2)
     itrace = inparms(3)
     irenum = inparms(4)
     ntryslv= inparms(5)
     nprecs = inparms(6)
     istopc = inparms(7)
     irst   = inparms(8)
!!$    write(0,*) 'Recvd nprecs: ',nprecs
     allocate(precs(nprecs))
     call igebr2d(ictxt,'all',' ',nprecs,1,precs,nprecs,0,0)
     call igebr2d(ictxt,'all',' ',1,1,novrs,1,0,0)
!!$    write(0,*) 'Recvd novrs: ',novrs    
    allocate(ovrs(novrs))   
    call igebr2d(ictxt,'all',' ',novrs,1,ovrs,novrs,0,0)
    call dgebr2d(ictxt,'all',' ',1,1,eps,1,0,0)
    call dgebr2d(ictxt,'all',' ',1,1,omega,1,0,0)
    call igebr2d(ictxt,'all',' ',1,1,igsmth,1,0,0)
    call igebr2d(ictxt,'all',' ',1,1,matop,1,0,0)
    call igebr2d(ictxt,'all',' ',1,1,jacswp,1,0,0)
    call igebr2d(ictxt,'all',' ',1,1,nmats,1,0,0)
  endif

  do imat=1,nmats

     if (amroot) then 
!!$      read(*,*) mtrx_file,rhs_file 
        read(*,'(a)') charbuf
        charbuf=adjustl(charbuf)
        i=index(charbuf," ")
        mtrx_file=charbuf(1:i-1)
        rhs_file=adjustl(charbuf(i+1:))
        i=index(mtrx_file,"/",back=.true.)
        mtrx_name=trim(mtrx_file(i+1:))
        write(0,*) 'Read mtx rhs : "',&
             & mtrx_file,'" "',rhs_file,'" "',mtrx_name,'"'
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
        if (.true.) then 
          allocate(ivg(m_problem),ipv(np))
          do i=1,m_problem
            call part_block(i,m_problem,np,ipv,nv)
            ivg(i) = ipv(1)
          enddo
          call matdist(aux_a, a, ivg, ictxt, &
               & desc_a,b_col_glob,b_col,info,fmt=afmt)
        else
          call matdist(aux_a, a, part_block, ictxt, &
               & desc_a,b_col_glob,b_col,info,fmt=afmt)
        endif
     else  if (ipart.eq.1) then 
        call blacs_barrier(ictxt,'A')
        call matdist(aux_a, a, part_blk2, ictxt, &
             & desc_a,b_col_glob,b_col,info,fmt=afmt)
     else if (ipart.eq.2) then 
        if (amroot) then 
           call build_grppart(aux_a%m,aux_a%fida,aux_a%ia1,aux_a%ia2,np)
        endif
        call distr_grppart(0,0,ictxt)
        if (.false.) then 
           call matdist(aux_a, a, part_graph, ictxt, &
                & desc_a,b_col_glob,b_col,info,fmt=afmt)
        else
           call getv_grppart(ivg)
           call matdist(aux_a, a, ivg, ictxt, &
                & desc_a,b_col_glob,b_col,info,fmt=afmt)
        endif
     else 
        call matdist(aux_a, a, part_block, ictxt, &
             & desc_a,b_col_glob,b_col,info,fmt=afmt)
     end if

     if(info /= 0) then
        info=4010
        ch_err='matdist'
        goto 9999
     end if

     call f90_psdsall(m_problem,x_col,desc_a,info)
     if(info /= 0) then
        info=4010
        ch_err='f90_psdsall'
        goto 9999
     end if
     x_col(:) =0.0
     call f90_psdsasb(x_col,desc_a,info)
     call f90_psdsasb(b_col,desc_a,info)
     if(info /= 0) then
        info=4010
        ch_err='f90_psdsasb'
        goto 9999
     end if

     call f90_psdsall(m_problem,r_col,desc_a,info)
     if(info /= 0) then
        info=4010
        ch_err='f90_psdsall'
        goto 9999
     end if
     r_col(:) =0.0
     call f90_psdsasb(r_col,desc_a,info)
     if(info /= 0) then
        info=4010
        ch_err='f90_psdsasb'
        goto 9999
     end if
     t2 = mpi_wtime() - t1


     dim=size(a%aspk)

     call dgamx2d(ictxt, 'a', ' ', ione, ione, t2, ione,&
          & t1, t1, -1, -1, -1)

     if (amroot) then
        write(6,*) 'time to Read and Partition Matrix : ',T2
     END IF
!!$    call blacs_barrier(ictxt,'all')
!!$    do i=0, nprow-1
!!$      if (myprow==i) then 
!!$        write(6,*) 'Main descriptor for process ',i,' ',mtrx_file
!!$        call descprt(6,desc_a,short=.true.)
!!$      endif
!!$      call blacs_barrier(ictxt,'all')
!!$    enddo


     !
     !  Prepare the preconditioning matrix. Note the availability
     !  of optional parameters
     !

     do ip=1,nprecs

        pre%prec=precs(ip)
        if (precs(ip) > 2) then 
           no=novrs
        else
           no=1
        endif

        do io=1, no

           ttot = 1.d300

           do itryslv=1,ntryslv
              ! Zero initial guess.
              select case(precs(ip))
              case(noprec_)
                 ptype='noprec'
                 call psb_precset(pre,ptype)
              case(diagsc_)             
                 ptype='diagsc'
                 call psb_precset(pre,ptype)
              case(bja_)             
                 ptype='bja'
                 call psb_precset(pre,ptype)
              case(asm_)             
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,sum_/))
              case(ash_)             
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),nohalo_,sum_/))
              case(ras_)             
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
              case(50+ras_)
                 pre%prec = ras_
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_,f_slu_/))
              case(rash_)             
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),nohalo_,none_/))
              case(ras2lv_) 
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 call psb_precset(pre,ptype,&
                      &iv=(/add_ml_prec_,glb_aggr_,pre_smooth_,igsmth,matop/),rs=0.d0)
              case(ras2lvm_) 
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 call psb_precset(pre,ptype,&
                      & iv=(/mult_ml_prec_,glb_aggr_,pre_smooth_,igsmth,matop/),rs=0.d0)
              case(lv2mras_) 
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 call psb_precset(pre,ptype,&
                      & iv=(/mult_ml_prec_,glb_aggr_,post_smooth_,igsmth,matop/),rs=0.d0)
              case(50+lv2mras_) 
                 pre%prec = lv2mras_
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 call psb_precset(pre,ptype,&
                      & iv=(/mult_ml_prec_,glb_aggr_,post_smooth_,igsmth,matop,f_slu_/),&
                      & rs=0.d0)
              case(lv2smth_) 
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,'ml',&
                        & iv=(/mult_ml_prec_,glb_aggr_,post_smooth_,igsmth,matop/),&
                        & rs=omega)
                 else
                   call psb_precset(pre,'ml',&
                        & iv=(/mult_ml_prec_,glb_aggr_,post_smooth_,igsmth,matop/))
                 endif
              case(50+lv2smth_) 
                 pre%prec = lv2smth_
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,ptype,&
                        & iv=(/mult_ml_prec_,glb_aggr_,post_smooth_,igsmth,matop,f_slu_/),&
                        & rs=omega)
                 else
                   call psb_precset(pre,ptype,&
                        & iv=(/mult_ml_prec_,glb_aggr_,post_smooth_,igsmth,matop,f_slu_/))
                 endif
              case(lv2lsm_) 
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,ptype,&
                        & iv=(/mult_ml_prec_,loc_aggr_,post_smooth_,igsmth,matop/),&
                        & rs=omega)
                 else
                   call psb_precset(pre,ptype,&
                        & iv=(/mult_ml_prec_,loc_aggr_,post_smooth_,igsmth,matop/))
                 endif
              case(50+lv2lsm_) 
                 pre%prec = lv2lsm_
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,ptype,&
                        & iv=(/mult_ml_prec_,loc_aggr_,post_smooth_,igsmth,matop,f_slu_/),&
                        & rs=omega)
                 else
                   call psb_precset(pre,ptype,&
                        & iv=(/mult_ml_prec_,loc_aggr_,post_smooth_,igsmth,matop,f_slu_/))
                 endif
              case(sl2sm_) 
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,sum_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,ptype,&
                      & iv=(/add_ml_prec_,loc_aggr_,post_smooth_,igsmth,matop/),rs=omega)
                 else
                   call psb_precset(pre,ptype,&
                      & iv=(/add_ml_prec_,loc_aggr_,post_smooth_,igsmth,matop/))
                 endif
              case(new_loc_smth_) 
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,ptype,&
                        & iv=(/new_ml_prec_,new_loc_aggr_,post_smooth_,1,&
                        & matop,f_ilu_n_,jacswp/), rs=omega)
                 else
                   call psb_precset(pre,ptype,&
                        & iv=(/new_ml_prec_,new_loc_aggr_,post_smooth_,1,&
                        & matop,f_ilu_n_,jacswp/))
                 endif
               case(50+new_loc_smth_) 
                 pre%prec = new_loc_smth_
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,ptype,&
                        & iv=(/new_ml_prec_,new_loc_aggr_,post_smooth_,1,&
                        & matop,f_slu_,jacswp/), rs=omega)
                 else
                   call psb_precset(pre,ptype,&
                        & iv=(/new_ml_prec_,new_loc_aggr_,post_smooth_,1,&
                        & matop,f_slu_,jacswp/))
                 endif
              case(new_glb_smth_) 
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,ptype,&
                        & iv=(/new_ml_prec_,new_glb_aggr_,post_smooth_,1,&
                        & matop,f_ilu_n_,1/), rs=omega)
                 else
                   call psb_precset(pre,ptype,&
                        & iv=(/new_ml_prec_,new_glb_aggr_,post_smooth_,1,&
                        & matop,f_ilu_n_,1/))
                 endif
              case(50+new_glb_smth_) 
                 pre%prec = new_glb_smth_
                 ptype='asm'
                 call psb_precset(pre,ptype,iv=(/ovrs(io),halo_,none_/))
                 ptype='ml'
                 if (omega>0.0d0) then 
                   call psb_precset(pre,ptype,&
                        & iv=(/new_ml_prec_,new_glb_aggr_,post_smooth_,1,&
                        & matop,f_slu_,1/), rs=omega)
                 else
                   call psb_precset(pre,ptype,&
                        & iv=(/new_ml_prec_,new_glb_aggr_,post_smooth_,1,&
                        & matop,f_slu_,1/))
                 endif
              case default
                 write(0,*) 'Unknown iprec, defaulting to BJA'
                 ptype='bja'
                 call psb_precset(pre,ptype)
              end select

              call blacs_barrier(ictxt,'All')
              call f90_psaxpby(0.d0,b_col,0.d0,x_col,desc_a,info)
              if(info /= 0) then
                 info=4010
                 ch_err='f90_psaxpby'
                 goto 9999
              end if


              call blacs_barrier(ictxt,'All')
              t1 = mpi_wtime()
              call psb_precbld(a,pre,desc_a,info)!,'f')
              t2 = mpi_wtime()-t1

              if(info /= 0) then
                 info=4010
                 ch_err='psb_precbld'
                 goto 9999
              end if

              call dgamx2d(ictxt,'a',' ',ione, ione,t2,ione,t1,t1,-1,-1,-1)

              if (info /= 0) then
                 write(0,*) 'error in preconditioner :',info
                 call blacs_abort(ictxt,-1)
                 stop
              end if

              iparm = 0
              rparm = 0.d0   

              call blacs_barrier(ictxt,'all')
              t1 = mpi_wtime()
              if (cmethd.eq.'CG') Then
                 call  f90_cg(a,pre,b_col,x_col,eps,desc_a,info,& 
                      & itmax,iter,err,itrace,istop=istopc)     
                 if(info /= 0) then
                    info=4010
                    ch_err='f90_cg'
                    goto 9999
                 end if
              else if (cmethd.eq.'BICGSTAB') Then
                 call  f90_bicgstab(a,pre,b_col,x_col,eps,desc_a,info,& 
                      & itmax,iter,err,itrace,istop=istopc)     
                 if(info /= 0) then
                    info=4010
                    ch_err='f90_bicgstab'
                    goto 9999
                 end if
              ELSE IF (CMETHD.EQ.'BICG') Then
                 call  f90_bicg(a,pre,b_col,x_col,eps,desc_a,info,& 
                      & itmax,iter,err,itrace,istop=istopc)     
                 if(info /= 0) then
                    info=4010
                    ch_err='f90_bicg'
                    goto 9999
                 end if
              ELSE IF (CMETHD.EQ.'CGS') Then
                 call  f90_cgs(a,pre,b_col,x_col,eps,desc_a,info,& 
                      & itmax,iter,err,itrace,istop=istopc)     
                 if(info /= 0) then
                    info=4010
                    ch_err='f90_cg'
                    goto 9999
                 end if
              ELSE IF (CMETHD.EQ.'GMRES') Then
                 call  f90_rgmres(a,pre,b_col,x_col,eps,desc_a,info,& 
                      & itmax,iter,err,itrace,irst=irst,istop=istopc)     
                 if(info /= 0) then
                    info=4010
                    ch_err='f90_rgmres'
                    goto 9999
                 end if
              else
                 write(*,*) 'Unknown method : "',cmethd,'"'
              ENDIF
              call blacs_barrier(ictxt,'all')
              t3 = mpi_wtime() - t1
              call dgamx2d(ictxt,'a',' ',ione, ione,t3,ione,t1,t1,-1,-1,-1)

              if (amroot) then
                 write(10,'(a18,3(1x,i2),1x,i5,5(1x,g9.4))') mtrx_name,nprow,&
                      & precs(ip),pre%baseprec%iprcparm(n_ovr_),iter,t2,t3,t2+t3
              endif
              if (itryslv < ntryslv) then 
                 call psb_precfree(pre,info)
                 if(info /= 0) then
                    info=4010
                    ch_err='psb_precfree'
                    goto 9999
                 end if
              else
                 if (amroot) call prec_descr(6,pre)
              end if
              if ((t2+t3)<ttot) then
                 tprec = t2
                 tslv = t3
                 ttot = t2+t3
              end if

           end do


           call f90_psaxpby(one,b_col,dzero,r_col,desc_A,info)
           call f90_psspmm(-one,a,x_col,one,r_col,desc_a,info)
           call f90_nrm2(resmx,r_col,desc_a,info)
           call f90_amax(resmxp,r_col,desc_a,info)
           if(info /= 0) then
              info=4011
              goto 9999
           end if

           if (amroot) then 
              write(6,*) ' iprec istopc   : ',precs(ip), istopc, ovrs(io)
!!$!    write(6,*) 'Number of iterations : ',iter
!!$!    WRITE(6,*) 'Error on exit        : ',ERR
              write(6,*) 'Matrix: ',mtrx_name
              write(6,*) 'Computed solution on ',NPROW,' processors.'
              write(6,*) 'Iterations to convergence: ',iter
              write(6,*) 'Error indicator on exit:',err
              write(6,*) 'Time to Buil Prec.   : ',TPREC
              write(6,*) 'Time to Solve Matrix : ',TSLV
              WRITE(6,*) 'Time per iteration   : ',TSLV/(ITER)
              write(6,*) 'Total Time           : ',ttot
              write(6,*) 'Residual norm 2   = ',resmx
              write(6,*) 'Residual norm inf = ',resmxp
              write(6,*) 
              write(6,*) 

              write(8,'(a18,3(1x,i2),1x,i5,5(1x,g9.4))') mtrx_name,nprow,precs(ip),&
                   & pre%baseprec%iprcparm(n_ovr_),&
                   & iter,tprec,tslv,ttot,resmx,resmxp
              if (associated(pre%smthprec)) then 
                 write(11,'(a18,9(1x,i2),1(1x,g9.4),1x,i5,5(1x,g9.4))') &
                      & mtrx_name,nprow,&
                      & pre%baseprec%iprcparm(p_type_),pre%baseprec%iprcparm(n_ovr_),&
                      & pre%baseprec%iprcparm(restr_),pre%baseprec%iprcparm(prol_),&
                      & pre%smthprec%iprcparm(ml_type_),pre%smthprec%iprcparm(aggr_alg_),&
                      & pre%smthprec%iprcparm(smth_pos_),pre%smthprec%iprcparm(glb_smth_),&
                      & pre%smthprec%dprcparm(smooth_omega_),&
                      & iter,tprec,tslv,ttot,resmx,resmxp
              else
                 write(11,'(a18,9(1x,i2),1(1x,g9.4),1x,i5,5(1x,g9.4))') &
                      & mtrx_name,nprow,&
                      & pre%baseprec%iprcparm(p_type_),pre%baseprec%iprcparm(n_ovr_),&
                      & pre%baseprec%iprcparm(prol_),pre%baseprec%iprcparm(restr_),&
                      & 0, 0, 0,&
                      & 0,  0.0,&
                      & iter,tprec,tslv,ttot,resmx,resmxp
              end if
           end if
           call psb_precfree(pre,info)
           if(info /= 0) then
              info=4011
              goto 9999
           end if
        end do
     end do
     deallocate(aux_b,stat=info)
     if (amroot) call spfree(aux_a,info)
!!$
     call f90_psdsfree(b_col, desc_a,info)
     call f90_psdsfree(x_col, desc_a,info)
     call f90_psdsfree(r_col, desc_a,info)
     call f90_psspfree(a, desc_a,info)
     call f90_psdscfree(desc_a,info)
     if(info /= 0) then
        info=4011
        goto 9999
     end if

  end do
  deallocate(ovrs,precs,stat=info)
  write(0,*) 'Info from deallocate ',info

9999 continue
  if(info /= 0) then
     call psb_errpush(info,name,a_err=ch_err)
     call psb_error(ictxt)
     call blacs_gridexit(ictxt)
     call blacs_exit(0)
  else
     call blacs_gridexit(ictxt)
     call blacs_exit(0)
  end if

end program df_bench
  




