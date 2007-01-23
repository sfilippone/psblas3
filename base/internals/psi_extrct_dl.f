C
C             Parallel Sparse BLAS  v2.0
C   (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
C                      Alfredo Buttari        University of Rome Tor Vergata
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions
C are met:
C   1. Redistributions of source code must retain the above copyright
C      notice, this list of conditions and the following disclaimer.
C   2. Redistributions in binary form must reproduce the above copyright
C      notice, this list of conditions, and the following disclaimer in the
C      documentation and/or other materials provided with the distribution.
C   3. The name of the PSBLAS group or the names of its contributors may
C      not be used to endorse or promote products derived from this
C      software without specific written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
C TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
C PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
C BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
C CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
C SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
C INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
C CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
C ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
C POSSIBILITY OF SUCH DAMAGE.
C
C 
      subroutine psi_extract_dep_list(desc_data,
     +  desc_str,dep_list,
     +  length_dl,np,dl_lda,mode,info)

c    internal routine
c    ================ 
c   
c    _____called by psi_crea_halo and psi_crea_ovrlap ______
c
c purpose
c =======
c   process root (pid=0) extracts for each process "k" the ordered list of process
c   to which "k" must communicate. this list with its order is extracted from
c   desc_str list
c
c   
c  input
c =======
c  desc_data :integer array
c  explanation:
c  name		 explanation
c  ------------------ -------------------------------------------------------
c  desc_data	 array of integer that contains some local and global
c		 information of matrix.
c
c
c  now we explain each of the above vectors.
c
c  let a be a generic sparse matrix. we denote with matdata_a the matrix_data
c  array for matrix a.
c  data stored in matrix_data array are:
c
c  notation        stored in		     explanation
c  --------------- ---------------------- -------------------------------------
c  dec_type        matdata_a[psb_dec_type_]   decomposition type
c  m 	           matdata_a[m_]          total number of equations
c  n 	           matdata_a[n_]          total number of variables
c  n_row           matdata_a[psb_n_row_]      number of local equations
c  n_col           matdata_a[psb_n_col_]      number of local variables
c  psb_ctxt_a          matdata_a[ctxt_]       the blacs context handle, indicating
c	     	                          the global context of the operation
c					  on the matrix.
c					  the context itself is global.
c  desc_str integer array
c  explanation:
c  let desc_str_p be the array desc_str for local process.
c  this is composed of variable dimension blocks for each process to 
c  communicate to.
c  each block contain indexes of local halo elements to exchange with other 
c  process.
c  let p be the pointer to the first element of a block in desc_str_p.
c  this block is stored in desc_str_p as :
c
c  notation        stored in		          explanation
c  --------------- --------------------------- -----------------------------------
c  process_id      desc_str_p[p+psb_proc_id_]      identifier of process which exchange
c						  data with.
c  n_elements_recv desc_str_p[p+n_elem_recv_]  number of elements to receive.
c  elements_recv   desc_str_p[p+elem_recv_+i]  indexes of local elements to
c					          receive. these are stored in the
c					          array from location p+elem_recv_ to
c					          location p+elem_recv_+
c						  desc_str_p[p+n_elem_recv_]-1.
c  if desc_data(psb_dec_type_) == 0 
c  then also will be:
c  n_elements_send desc_str_p[p+n_elem_send_]  number of elements to send.
c  elements_send   desc_str_p[p+elem_send_+i]  indexes of local elements to
c					          send. these are stored in the
c					          array from location p+elem_send_ to
c					          location p+elem_send_+
c						  desc_str_p[p+n_elem_send_]-1.
c  list is ended by -1 value
c  
c  np     integer (global input)
c         number of grid process.
c
c  mode   integer (global input)
c         if mode =0 then will be inserted also duplicate element in 
c         a same dependence list
c         if mode =1 then not will be inserted duplicate element in
c          a same dependence list
c  output
c  =====
c  only for root (pid=0) process:
c  dep_list  integer array(dl_lda,0:np)
c            dependence list dep_list(*,i) is the list of process identifiers to which process i
c            must communicate with. this list with its order is extracted from
c            desc_str list.
c   length_dl  integer array(0:np)
c             length_dl(i) is the length of dep_list(*,i) list
      use psb_penv_mod
      use psb_const_mod
      use psb_error_mod
      use psb_descriptor_type
      implicit none
      include 'mpif.h'
c     ....scalar parameters...
      integer np,dl_lda,mode, info

c     ....array parameters....
      integer desc_str(*),desc_data(*),
     +  dep_list(dl_lda,0:np),length_dl(0:np)
      integer, pointer :: itmp(:)
c     .....local arrays....
      integer int_err(5)
      double precision real_err(5)

c     .....local scalars...
      integer i,nprow,npcol,me,mycol,pointer_dep_list,proc,j,err_act
      integer ictxt, err, icomm
      logical debug
      parameter (debug=.false.)
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
      if (debug) write(0,*) 'extract: info ',info,
     +  desc_data(psb_dec_type_)
      pointer_dep_list=1
c$$$      if (desc_data(psb_dec_type_).eq.psb_desc_bld_) then
      if (psb_is_bld_dec(desc_data(psb_dec_type_))) then 
        do while (desc_str(i).ne.-1)
          if (debug) write(0,*) me,' extract: looping ',i,
     +      desc_str(i),desc_str(i+1),desc_str(i+2)
          
c        ...with different decomposition type we have different
c           structure of indices  lists............................
          if ((desc_str(i+1).ne.0).or.(desc_str(i+2).ne.0)) then
c           ..if number of element to be exchanged !=0
            proc=desc_str(i)
            if ((proc.lt.0).or.(proc.ge.nprow)) then
              if (debug) write(0,*) 'extract error ',i,desc_str(i)
              info = 9999
              int_err(1) = i
              int_err(2) = desc_str(i)
              goto 998
            endif
!            if((me.eq.1).and.(proc.eq.3))write(0,*)'found 3'
            if (mode.eq.1) then
c              ...search if already exist proc 
c                 in dep_list(*,me)...  
              j=1
              do while ((j.lt.pointer_dep_list).and.
     +          (dep_list(j,me).ne.proc))
                j=j+1
              enddo
              
              if (j.eq.pointer_dep_list) then
c                 ...if not found.....
                dep_list(pointer_dep_list,me)=proc
                pointer_dep_list=pointer_dep_list+1
              endif
            else if (mode.eq.0) then
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
c$$$      else if (desc_data(psb_dec_type_).eq.psb_desc_upd_) then
      else if (psb_is_upd_dec(desc_data(psb_dec_type_))) then
        do while (desc_str(i).ne.-1)
          if (debug) write(0,*) 'extract: looping ',i,desc_str(i)
          
c        ...with different decomposition type we have different
c           structure of indices  lists............................
          if (desc_str(i+1).ne.0) then
            
            proc=desc_str(i)
c        ..if number of element to be exchanged !=0
            
            if (mode.eq.1) then
c              ...search if already exist proc....                 
              j=1
              do while ((j.lt.pointer_dep_list).and.
     +          (dep_list(j,me).ne.proc))
                j=j+1
              enddo
              if (j.eq.pointer_dep_list) then
c                 ...if not found.....
                if (pointer_dep_list.gt.dl_lda) then
                  info = 4000
                  goto 998
                endif
                dep_list(pointer_dep_list,me)=proc
                pointer_dep_list=pointer_dep_list+1
              endif
            else if (mode.eq.0) then
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

c     ... check for errors...
 998  continue 
      if (debug) write(0,*) 'extract: info ',info
      err = info

      if (err.ne.0) goto 9999

      call psb_sum(ictxt,length_dl(0:np))
      call psb_get_mpicomm(ictxt,icomm )
      allocate(itmp(dl_lda),stat=info)
      if (info /= 0) goto 9999
      itmp(1:dl_lda) = dep_list(1:dl_lda,me)
      call mpi_allgather(itmp,dl_lda,mpi_integer,
     +  dep_list,dl_lda,mpi_integer,icomm,info)
      deallocate(itmp)

      call psb_erractionrestore(err_act)
      return


 9999 continue

      call psb_errpush(info,name,i_err=int_err)
      call psb_erractionrestore(err_act)
      if (err_act.eq.psb_act_ret_) then
        return
      else
        call psb_error()
      endif
      return

      end
