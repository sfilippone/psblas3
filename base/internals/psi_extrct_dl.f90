!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
subroutine psi_extract_dep_list(desc_data,desc_str,dep_list,&
     & length_dl,np,dl_lda,mode,info)

  !    internal routine
  !    ================ 
  !   
  !    _____called by psi_crea_halo and psi_crea_ovrlap ______
  !
  ! purpose
  ! =======
  !   process root (pid=0) extracts for each process "k" the ordered list of process
  !   to which "k" must communicate. this list with its order is extracted from
  !   desc_str list
  ! 
  !   
  !  input
  ! =======
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
  !  =====
  !  only for root (pid=0) process:
  !  dep_list  integer array(dl_lda,0:np)
  !            dependence list dep_list(*,i) is the list of process identifiers to which process i
  !            must communicate with. this list with its order is extracted from
  !            desc_str list.
  !   length_dl  integer array(0:np)
  !             length_dl(i) is the length of dep_list(*,i) list
  use mpi
  use psb_penv_mod
  use psb_const_mod
  use psb_error_mod
  use psb_descriptor_type
  implicit none
  !     ....scalar parameters...
  integer np,dl_lda,mode, info

  !     ....array parameters....
  integer ::  desc_str(*),desc_data(*),dep_list(dl_lda,0:np),length_dl(0:np)
  integer, pointer :: itmp(:)
  !     .....local arrays....
  integer int_err(5)
  double precision real_err(5)

  !     .....local scalars...
  integer i,nprow,npcol,me,mycol,pointer_dep_list,proc,j,err_act
  integer ictxt, err, icomm
  logical, parameter :: debug=.false.
  character  name*20
  name='psi_extrct_dl'

  call psb_erractionsave(err_act)

  info = 0
  ictxt = desc_data(psb_ctxt_)


  call psb_info(ictxt,me,nprow)
  do i=0,np 
    length_dl(i) = 0
  enddo
  i=1
  if (debug) write(0,*) 'extract: info ',info,desc_data(psb_dec_type_)
  pointer_dep_list=1
  if (psb_is_bld_dec(desc_data(psb_dec_type_))) then 
    do while (desc_str(i) /= -1)
      if (debug) write(0,*) me,' extract: looping ',i,&
           &    desc_str(i),desc_str(i+1),desc_str(i+2)

      !        ...with different decomposition type we have different
      !           structure of indices  lists............................
      if ((desc_str(i+1) /= 0).or.(desc_str(i+2) /= 0)) then
        !           ..if number of element to be exchanged !=0
        proc=desc_str(i)
        if ((proc < 0).or.(proc.ge.nprow)) then
          if (debug) write(0,*) 'extract error ',i,desc_str(i)
          info = 9999
          int_err(1) = i
          int_err(2) = desc_str(i)
          goto 998
        endif
        !            if((me == 1).and.(proc == 3))write(0,*)'found 3'
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
          if (pointer_dep_list.gt.dl_lda) then
            info = 4000
            goto 998
          endif
          dep_list(pointer_dep_list,me)=proc
          pointer_dep_list=pointer_dep_list+1
        endif
      endif
      i=i+desc_str(i+1)+2
    enddo
  else if (psb_is_upd_dec(desc_data(psb_dec_type_))) then
    do while (desc_str(i) /= -1)
      if (debug) write(0,*) 'extract: looping ',i,desc_str(i)

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
            if (pointer_dep_list.gt.dl_lda) then
              info = 4000
              goto 998
            endif
            dep_list(pointer_dep_list,me)=proc
            pointer_dep_list=pointer_dep_list+1
          endif
        else if (mode == 0) then
          if (pointer_dep_list.gt.dl_lda) then
            info = 4000
            goto 998
          endif
          dep_list(pointer_dep_list,me)=proc
          pointer_dep_list=pointer_dep_list+1
        endif
      endif
      i=i+desc_str(i+1)+2
    enddo
  else
    write(0,*) 'invalid dec_type',desc_data(psb_dec_type_)
    info = 2020
    goto 9999
  endif

  length_dl(me)=pointer_dep_list-1

  !     ... check for errors...
998 continue 
  if (debug) write(0,*) 'extract: info ',info
  err = info

  if (err /= 0) goto 9999

  call psb_sum(ictxt,length_dl(0:np))
  call psb_get_mpicomm(ictxt,icomm )
  allocate(itmp(dl_lda),stat=info)
  if (info /= 0) goto 9999
  itmp(1:dl_lda) = dep_list(1:dl_lda,me)
  call mpi_allgather(itmp,dl_lda,mpi_integer,&
       & dep_list,dl_lda,mpi_integer,icomm,info)
  deallocate(itmp)

  call psb_erractionrestore(err_act)
  return


9999 continue

  call psb_errpush(info,name,i_err=int_err)
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error()
  endif
  return

end subroutine psi_extract_dep_list
