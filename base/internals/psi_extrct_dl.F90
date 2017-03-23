!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
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
subroutine psi_extract_dep_list(ictxt,is_bld,is_upd,desc_str,dep_list,&
     & length_dl,np,dl_lda,mode,info)

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
  use psi_mod, psb_protect_name => psi_extract_dep_list
#ifdef MPI_MOD
  use mpi
#endif
  use psb_penv_mod
  use psb_const_mod
  use psb_error_mod
  use psb_desc_mod
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  !     ....scalar parameters...
  logical :: is_bld, is_upd
  integer(psb_ipk_) :: ictxt
  integer(psb_ipk_) :: np,dl_lda,mode, info

  !     ....array parameters....
  integer(psb_ipk_) ::  desc_str(*),dep_list(dl_lda,0:np),length_dl(0:np)
  integer(psb_ipk_), allocatable :: itmp(:)
  !     .....local arrays....
  integer(psb_ipk_) :: int_err(5)

  !     .....local scalars...
  integer(psb_ipk_) :: i,pointer_dep_list,proc,j,err_act
  integer(psb_ipk_) :: err
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_mpik_) :: iictxt, icomm, me, npr, dl_mpi, minfo
  character  name*20
  name='psi_extrct_dl'

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  iictxt = ictxt 
  info = psb_success_

  call psb_info(iictxt,me,npr)
  do i=0,np 
    length_dl(i) = 0
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
        if ((proc < 0).or.(proc >= npr)) then
          if (debug_level >= psb_debug_inner_)&
               & write(debug_unit,*) me,' ',trim(name),': error ',i,desc_str(i)
          info = 9999
          int_err(1) = i
          int_err(2) = desc_str(i)
          goto 998
        endif
        !            if((me == 1).and.(proc == 3))write(psb_err_unit,*)'found 3'
        if (mode == 1) then
          !              ...search if already exist proc 
          !                 in dep_list(*,me)...  
          j=1
          do while ((j < pointer_dep_list).and.&
               & (dep_list(j,me) /= proc))
            j=j+1
          enddo

          if (j == pointer_dep_list) then
            !                 ...if not found.....
            dep_list(pointer_dep_list,me)=proc
            pointer_dep_list=pointer_dep_list+1
          endif
        else if (mode == 0) then
          if (pointer_dep_list > dl_lda) then
            info = psb_err_alloc_dealloc_
            goto 998
          endif
          dep_list(pointer_dep_list,me)=proc
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
               & (dep_list(j,me) /= proc))
            j=j+1
          enddo
          if (j == pointer_dep_list) then
            !                 ...if not found.....
            if (pointer_dep_list > dl_lda) then
              info = psb_err_alloc_dealloc_
              goto 998
            endif
            dep_list(pointer_dep_list,me)=proc
            pointer_dep_list=pointer_dep_list+1
          endif
        else if (mode == 0) then
          if (pointer_dep_list > dl_lda) then
            info = psb_err_alloc_dealloc_
            goto 998
          endif
          dep_list(pointer_dep_list,me)=proc
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

  !     ... check for errors...
998 continue 
  if (debug_level >= psb_debug_inner_)&
       & write(debug_unit,*) me,' ',trim(name),': info ',info
  err = info

  if (err /= 0) goto 9999

  call psb_sum(iictxt,length_dl(0:np))
  call psb_get_mpicomm(iictxt,icomm )
  allocate(itmp(dl_lda),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    goto 9999
  endif
  itmp(1:dl_lda) = dep_list(1:dl_lda,me)
  dl_mpi = dl_lda
  call mpi_allgather(itmp,dl_mpi,psb_mpi_ipk_integer,&
       & dep_list,dl_mpi,psb_mpi_ipk_integer,icomm,minfo)
  info = minfo
  if (info == 0) deallocate(itmp,stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return


9999 continue

  call psb_errpush(info,name,i_err=int_err)
  call psb_error_handler(err_act)

  return

end subroutine psi_extract_dep_list
