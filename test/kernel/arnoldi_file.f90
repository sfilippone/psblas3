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
program arnoldi_file
  use psb_base_mod
  use psb_util_mod
  implicit none

  ! input parameters
  character(len=40) :: kmethd, ptype, mtrx_file, rhs_file

  ! sparse matrices
  type(psb_zspmat_type) :: a, aux_a

  ! dense matrices
  complex(psb_dpk_), allocatable, target ::  aux_b(:,:), d(:)
  complex(psb_dpk_), allocatable , save  :: x_col_glob(:), r_col_glob(:)
  complex(psb_dpk_), allocatable, target ::  H(:,:),eig(:),work(:),Z(:,:)
  integer, allocatable :: indexes(:)
  type(psb_z_vect_type), allocatable, target ::  V(:)
  complex(psb_dpk_), pointer  :: b_col_glob(:)
  type(psb_z_vect_type)    :: b_col, x_col, r_col


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
  integer(psb_ipk_) :: internal, m,ii,nnzero, dim_H, alloc_stat
  real(psb_dpk_) :: t1, t2, r_amax, b_amax,&
       &scale,resmx,resmxp, flops, bdwdth
  real(psb_dpk_) :: tt1, tt2, tflops
  real (psb_dpk_) :: norm
  complex (psb_dpk_) :: dotprod 
  integer(psb_ipk_) :: nrhs, nrow, n_row, nv, ne
  integer(psb_ipk_), allocatable :: ivg(:), ipv(:)


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then 
    ! This should not happen, but just in case
    call psb_exit(ictxt)
    stop
  endif


  name='arnoldi_file'
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
    read(psb_inp_unit,*) dim_H
    write (psb_out_unit, '("Hessenberg matrix dim is ",i20)') dim_H     
  end if
  call psb_bcast(ictxt,mtrx_file)
  call psb_bcast(ictxt,filefmt)
  call psb_bcast(ictxt,ipart)
  call psb_bcast(ictxt,dim_H)
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
      ! For Harwell-Boeig we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_b,filename=mtrx_file)
      
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
    if (iam==psb_root_) write(psb_out_unit,'("Partition type: block")')
    allocate(ivg(m_problem),ipv(np))
    do i=1,m_problem
      call part_block(i,m_problem,np,ipv,nv)
      ivg(i) = ipv(1)
    enddo
    call psb_matdist(aux_a, a,ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,v=ivg)
    
  else if (ipart == 2) then 
    if (iam==psb_root_) then 
      write(psb_out_unit,'("Partition type: graph")')
      write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call build_mtpart(aux_a,np)

    endif
    call psb_barrier(ictxt)
    call distr_mtpart(psb_root_,ictxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,v=ivg)

  else 
    if (iam==psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ictxt, &
         & desc_a,b_col_glob,b_col,info,fmt=afmt,parts=part_block)
  end if

  allocate(H(dim_H,dim_H),stat = alloc_stat)
  do i=1,dim_h
        do j=1,dim_H
                H(i,j)=zzero
        enddo
  enddo
  allocate(V(dim_H+1),stat = alloc_stat)
  do i=1,dim_H+1
	  call psb_geall(V(i),desc_a,info)
	  call psb_geasb(V(i),desc_a,info)
  enddo
  call V(1)%set(zone)
  t2 = psb_wtime() - t1
  call psb_amx(ictxt, t2)

  if (iam==psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
    write(psb_out_unit,'(" ")')
  end if


  call psb_barrier(ictxt)
  t1 = psb_wtime() 
  
  norm = psb_genrm2(V(1),desc_a,info)
  if (iam==psb_root_) write(psb_out_unit,'("norma iniziale    :"F20.3)')norm
!  H(2,1)=dcmplx(norm,0.0)
  H(2,1)=norm*zone
  norm = 1/norm
  !normalisation of V(1)  
  call psb_geaxpby(zzero,V(1),norm*zone,V(1),desc_a, info)
  do i=2,dim_H+1
          call psb_spmm(zone,a,V(i-1),zzero,V(i),desc_a,info,'n') !we do V(i)=a*V(i-1)
          ! Gram-Schmitt's reorthogonalisation
          do j=1,i-1 
                dotprod= psb_gedot(V(i),V(j),desc_a,info) ! dotprod = (V(i) dot V(j))
                call psb_geaxpby(-dotprod,V(j),zone,V(i),desc_a, info) !V(i)=V(i)-V(j)*dotprod
                if(iam==psb_root_) write(*,'("dotprod : "i5, g20.4,g20.4)') i,real(dotprod),aimag(dotprod)
                H(j,i-1)=dotprod
          end do
        norm = psb_genrm2(V(i),desc_a,info)
        if (iam==psb_root_) then
                write(psb_out_unit,'("norma finale    :"i20,F20.3)')i,norm
                write(psb_out_unit,'("")')
        end if
        if (i .ne. dim_H+1) then
                H(i,i-1)=cmplx(norm,0.0)
                !H(i,i-1)=norm*zone
        endif
        norm=1/norm
        call psb_geaxpby(zzero,V(i),norm*zone,V(i),desc_a, info)        
  enddo
  do i=2,dim_H
        if (iam==psb_root_) write(*,'("basse diagonale de H : "i5, g20.4,g20.4)') i,real(H(i,i-1)),aimag(H(i,i-1))
  enddo

  write(psb_out_unit,'("")')

  allocate(eig(dim_H),work(dim_h),stat = info)
  allocate(Z(dim_H,dim_H),stat = alloc_stat)
  call ZHSEQR('E','N',dim_H,1,dim_H,H,dim_H,eig,Z,dim_H,work,3*dim_H,info)
       
  !sort H's eigenvalues
  allocate(indexes(1:dim_H))
 call psb_qsort(eig,indexes,psb_alsort_up_,psb_sort_ovw_idx_)
  
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
    write(psb_out_unit,'("Matrix: ",a)') mtrx_file
    write(psb_out_unit,'("Test on                          : ",i20," processors")') np
    write(psb_out_unit,'("Size of matrix                   : ",i20,"           ")') nr
    write(psb_out_unit,'("Number of nonzeros               : ",i20,"           ")') annz
    
    write(*,'("valeurs propres de H : ")')
   OPEN(unit=10, file=mtrx_file(1:3)//'.txt')
    do i=dim_H/3,dim_H
       write(psb_out_unit,'(g20.4,g20.4)')real(eig(i)),aimag(eig(i))
       write(10,'(g20.4,g20.4)')real(eig(i)),aimag(eig(i))
    enddo
    CLOSE(10)
   
    tflops = tflops / (tt2)
   !
    ! This computation is valid for CSR
    !
    nbytes = nr*(2*psb_sizeof_sp + psb_sizeof_int)+ &
         & annz*(psb_sizeof_sp + psb_sizeof_int)
    bdwdth = times*nbytes/(t2*1.d6)
    bdwdth = times*nbytes/(tt2*1.d6)
    
  end if

  call psb_gefree(b_col, desc_a,info)
  call psb_gefree(x_col, desc_a,info)
  call psb_spfree(a, desc_a,info)
  call psb_cdfree(desc_a,info)
  do i=1,dim_H
        call psb_gefree(V(i),desc_a,info)
  end do
  DEALLOCATE (H)
  DEALLOCATE (eig,V)
  DEALLOCATE (work)
  DEALLOCATE (Z)
9999 continue
  if(info /= 0) then
    call psb_error(ictxt)
  end if
  call psb_exit(ictxt)
  stop

end program arnoldi_file



