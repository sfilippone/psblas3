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
subroutine psi_i_extract_dep_list(ictxt,is_bld,is_upd,desc_str,dep_list,&
     & length_dl,dl_lda,mode,info)

  !    internal routine
  !    == = ============= 
  !   
  !    _____called by psi_crea_halo and psi_crea_ovrlap ______
  !
  ! purpose
  ! == = ====
  !   process root (pid=0) extracts for each process "k" the ordered list of process
  !   to which "k" must communicate. this list with its order is extracted from
  !   desc_str list
  ! 
  !   
  !  input
  ! == = ====
  !  desc_data :integer array
  !  explanation:
  !  name		 explanation
  !  ------------------ -------------------------------------------------------
  !  desc_data	 array of integer that contains some local and global
  !		 information of matrix.
  !
  !
  !  now we explain each of the above vectors.
  !
  !  let a be a generic sparse matrix. we denote with matdata_a the matrix_data
  !  array for matrix a.
  !  data stored in matrix_data array are:
  !
  !  notation        stored in		     explanation
  !  --------------- ---------------------- -------------------------------------
  !  dec_type        matdata_a[psb_dec_type_]   decomposition type
  !  m 	           matdata_a[m_]          total number of equations
  !  n 	           matdata_a[n_]          total number of variables
  !  n_row           matdata_a[psb_n_row_]      number of local equations
  !  n_col           matdata_a[psb_n_col_]      number of local variables
  !  psb_ctxt_a          matdata_a[ctxt_]       the blacs context handle, indicating
  !	     	                          the global context of the operation
  !					  on the matrix.
  !					  the context itself is global.
  !  desc_str integer array
  !  explanation:
  !  let desc_str_p be the array desc_str for local process.
  !! this is composed of variable dimension blocks for each process to 
  !  communicate to.
  !  each block contain indexes of local halo elements to exchange with other 
  !  process.
  !  let p be the pointer to the first element of a block in desc_str_p.
  !  this block is stored in desc_str_p as :
  !
  !  notation        stored in		          explanation
  !  --------------- --------------------------- -----------------------------------
  !  process_id      desc_str_p[p+psb_proc_id_]      identifier of process which exchange
  !						  data with.
  !  n_elements_recv desc_str_p[p+n_elem_recv_]  number of elements to receive.
  !  elements_recv   desc_str_p[p+elem_recv_+i]  indexes of local elements to
  !					          receive. these are stored in the
  !					          array from location p+elem_recv_ to
  !					          location p+elem_recv_+
  !						  desc_str_p[p+n_elem_recv_]-1.
  !  if desc_data(psb_dec_type_) == 0 
  !  then also will be:
  !  n_elements_send desc_str_p[p+n_elem_send_]  number of elements to send.
  !  elements_send   desc_str_p[p+elem_send_+i]  indexes of local elements to
  !					          send. these are stored in the
  !					          array from location p+elem_send_ to
  !					          location p+elem_send_+
  !						  desc_str_p[p+n_elem_send_]-1.
  !  list is ended by -1 value
  !  
  !  np     integer (global input)
  !         number of grid process.
  !
  !  mode   integer (global input)
  !         if mode =0 then will be inserted also duplicate element in 
  !         a same dependence list
  !         if mode =1 then not will be inserted duplicate element in
  !          a same dependence list
  !  output
  !  == = ==
  !  only for root (pid=0) process:
  !  dep_list  integer array(dl_lda,0:np)
  !            dependence list dep_list(*,i) is the list of process identifiers to which process i
  !            must communicate with. this list with its order is extracted from
  !            desc_str list.
  !   length_dl  integer array(0:np)
  !             length_dl(i) is the length of dep_list(*,i) list
  use psi_mod, psb_protect_name => psi_i_extract_dep_list
#ifdef MPI_MOD
  use mpi
#endif
  use psb_penv_mod
  use psb_const_mod
  use psb_error_mod
  use psb_desc_mod
  use psb_sort_mod
  use psb_timers_mod
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  !     ....scalar parameters...
  logical,  intent(in)           :: is_bld, is_upd
  integer(psb_ipk_), intent(in)  :: ictxt,mode
  integer(psb_ipk_), intent(out) :: dl_lda
  integer(psb_ipk_), intent(in)  :: desc_str(*)
  integer(psb_ipk_), allocatable, intent(out) :: dep_list(:,:),length_dl(:)
  integer(psb_ipk_), intent(out) :: info
  !     .....local arrays....
  integer(psb_ipk_) :: int_err(5)
  integer(psb_ipk_), allocatable :: itmp(:)

  !     .....local scalars...
  integer(psb_ipk_) :: i,pointer_dep_list,proc,j,err_act
  integer(psb_ipk_) :: err
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_mpk_) :: iictxt, icomm, me, np, minfo
  logical, parameter :: dist_symm_list=.false., print_dl=.false., profile=.true.
  logical, parameter  :: do_timings=.false.
  integer(psb_ipk_), save  :: idx_phase1=-1, idx_phase2=-1, idx_phase3=-1
  character  name*20
  name='psi_extrct_dl'

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  iictxt = ictxt 
  info = psb_success_
  if ((do_timings).and.(idx_phase1==-1))       &
       & idx_phase1 = psb_get_timer_idx("PSI_XTR_DL: phase1 ")
  if ((do_timings).and.(idx_phase2==-1))       &
       & idx_phase2 = psb_get_timer_idx("PSI_XTR_DL: phase2")
!!$  if ((do_timings).and.(idx_phase3==-1))       &
!!$       & idx_phase3 = psb_get_timer_idx("PSI_XTR_DL: phase3")

  call psb_info(iictxt,me,np)
  if (do_timings) call psb_tic(idx_phase1)

  allocate(itmp(2*np+1),length_dl(0:np),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    goto 9999
  end if
  do i=0,np 
    length_dl(i) = 0
    itmp(i+1)    = -1
  enddo
  i=1
  if (debug_level >= psb_debug_inner_)&
       & write(debug_unit,*) me,' ',trim(name),': start ',info

  pointer_dep_list=1
  if (is_bld) then 
    do while (desc_str(i) /= -1)
      if (debug_level >= psb_debug_inner_)&
           & write(debug_unit,*) me,' ',trim(name),' : looping ',i,&
           &    desc_str(i),desc_str(i+1),desc_str(i+2)

      !        ...with different decomposition type we have different
      !           structure of indices  lists............................
      if ((desc_str(i+1) /= 0).or.(desc_str(i+2) /= 0)) then
        !           ..if number of element to be exchanged !=0
        proc=desc_str(i)
        if ((proc < 0).or.(proc >= np)) then
          if (debug_level >= psb_debug_inner_)&
               & write(debug_unit,*) me,' ',trim(name),': error ',i,desc_str(i)
          info = 9999
          int_err(1) = i
          int_err(2) = desc_str(i)
          goto 9999
        endif
        !            if((me == 1).and.(proc == 3))write(psb_err_unit,*)'found 3'
        if (mode == 1) then
          !              ...search if already exist proc 
          !                 in itmp(*)...  
          j=1
          do while ((j < pointer_dep_list).and.&
               & (itmp(j) /= proc))
            j=j+1
          enddo

          if (j == pointer_dep_list) then
            !                 ...if not found.....
            itmp(pointer_dep_list)=proc
            pointer_dep_list=pointer_dep_list+1
          endif
        else if (mode == 0) then
          itmp(pointer_dep_list)=proc
          pointer_dep_list=pointer_dep_list+1
        endif
      endif
      i=i+desc_str(i+1)+2
    enddo

  else if (is_upd) then

    do while (desc_str(i) /= -1)
      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name),': looping ',i,desc_str(i)

      !        ...with different decomposition type we have different
      !           structure of indices  lists............................
      if (desc_str(i+1) /= 0) then

        proc=desc_str(i)
        !        ..if number of element to be exchanged !=0

        if (mode == 1) then
          !              ...search if already exist proc....                 
          j=1
          do while ((j < pointer_dep_list).and.&
               & (itmp(j) /= proc))
            j=j+1
          enddo
          if (j == pointer_dep_list) then
            !                 ...if not found.....
            itmp(pointer_dep_list)=proc
            pointer_dep_list=pointer_dep_list+1
          endif
        else if (mode == 0) then
          itmp(pointer_dep_list)=proc
          pointer_dep_list=pointer_dep_list+1
        endif
      endif
      i=i+desc_str(i+1)+2
    enddo
  else
    info = 2020
    goto 9999
  endif

  length_dl(me)=pointer_dep_list-1
  if (do_timings) call psb_toc(idx_phase1)
  
  if (do_timings) call psb_tic(idx_phase2)
  if (dist_symm_list) then 
    call psb_realloc(length_dl(me),itmp,info)
    call psi_symm_dep_list(itmp,ictxt,info)       
    dl_lda = max(size(itmp),1)
    call psb_max(iictxt, dl_lda)
    
    if (debug_level >= psb_debug_inner_) &
         & write(debug_unit,*) me,' ',trim(name),': Dep_list length ',length_dl(me),dl_lda
    call psb_realloc(dl_lda,itmp,info)
    !  dl_lda = min(np,2*dl_lda)
    allocate(dep_list(dl_lda,0:np),stat=info)
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999      
    end if

    call psb_sum(iictxt,length_dl(0:np))
    icomm = psb_get_mpi_comm(iictxt)
    call mpi_allgather(itmp,dl_lda,psb_mpi_ipk_,&
         & dep_list,dl_lda,psb_mpi_ipk_,icomm,minfo)
    info = minfo  
    if (info == 0) deallocate(itmp,stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_dealloc_
      goto 9999
    endif
  else
    block
      integer(psb_ipk_), allocatable :: list1(:,:), ldl2(:), list2(:,:)
      integer(psb_ipk_) :: i,j,ip,dlsym, ldu, mdl, l1, l2

      dl_lda = max(length_dl(me),1)
      call psb_max(iictxt, dl_lda)
      if (debug_level >= psb_debug_inner_) &
           & write(debug_unit,*) me,' ',trim(name),': Dep_list length ',length_dl(me),dl_lda
      allocate(dep_list(dl_lda,0:np),stat=info)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999      
      end if
      call psb_sum(iictxt,length_dl(0:np))
      icomm = psb_get_mpi_comm(iictxt)
      call mpi_allgather(itmp,dl_lda,psb_mpi_ipk_,&
           & dep_list,dl_lda,psb_mpi_ipk_,icomm,minfo)
      info = minfo  
      if (info /= psb_success_) then 
        info=psb_err_alloc_dealloc_
        goto 9999
      endif
      allocate(ldl2(0:np),stat=info)
      ldl2 = 0 
      do j=0, np-1
        do i=1,length_dl(j)
          ip = dep_list(i,j)
          ldl2(ip) = ldl2(ip) + 1
        end do
      end do
      dlsym = maxval(ldl2)
      allocate(list2(dlsym,0:np),stat=info)
      call move_alloc(dep_list,list1)
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999      
      end if
      ldl2 = 0 
      do j=0, np-1
        do i=1,length_dl(j)
          ip       = list1(i,j)
          ldl2(ip) = ldl2(ip) + 1
          list2(ldl2(ip),ip) = j
        end do
      end do
      mdl = 0 
      do j = 0, np-1
        l1 = length_dl(j)
        l2 = ldl2(j) 
        itmp(1:l1)       = list1(1:l1,j)
        itmp(l1+1:l1+l2) = list2(1:l2,j)
        ldu = l1 + l2
        !if (me == 0) write(0,*) 'Iter ',j,':',l1,l2,':',itmp(1:l1),':',itmp(l1+1:l1+l2)
        call psb_msort_unique(itmp(1:l1+l2),ldu)
        mdl = max(mdl, ldu)
        !if (me == 0) write(0,*) 'Iter ',j,':',ldu,':',itmp(1:ldu)
      end do
      dl_lda = mdl
      allocate(dep_list(dl_lda,0:np),stat=info)
      dep_list = -1
      if (info /= psb_success_) then 
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999      
      end if
      do j = 0, np-1
        l1 = length_dl(j)
        l2 = ldl2(j) 
        itmp(1:l1)       = list1(1:l1,j)
        itmp(l1+1:l1+l2) = list2(1:l2,j)
        ldu = l1 + l2 
        call psb_msort_unique(itmp(1:l1+l2),ldu)
        length_dl(j) = ldu
        dep_list(1:ldu,j) = itmp(1:ldu)
      end do
      
    end block
  endif
  if (print_dl) then
    if (me == 0) then
      write(0,*) ' Dep_list '
      do i=0,np-1
        j = length_dl(i) 
        write(0,*) 'Proc ',i,':',dep_list(1:j,i)
      end do
      flush(0)
    end if
    call psb_barrier(ictxt)
  end if
  if (do_timings) call psb_toc(idx_phase2)
  if ((profile).and.(me==0))  then
    block
      integer(psb_ipk_) :: dlmax, dlavg
      dlmax = maxval(length_dl(:))
      dlavg = (sum(length_dl(:))+np-1)/np
      if (dlmax>0) write(0,*) 'Dependency list : max:',dlmax,&
           & '  avg:',dlavg, ((dlmax>np/2).or.((dlavg>=np/4).and.(np>128)))
      
    end block
  end if

  call psb_erractionrestore(err_act)
  return


9999 continue

  call psb_errpush(info,name,i_err=int_err)
  call psb_error_handler(err_act)

  return

end subroutine psi_i_extract_dep_list
