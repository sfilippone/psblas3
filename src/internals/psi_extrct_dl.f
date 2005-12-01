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

      implicit none
      include 'psb_const.fh'
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
      integer icontxt, err, icomm
      logical debug
      parameter (debug=.false.)
      character  name*20
      name='psi_extrct_dl'
      call fcpsb_get_erraction(err_act)

      info = 0
      icontxt = desc_data(psb_ctxt_)


      call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)
      do i=0,np 
        length_dl(i) = 0
      enddo
      i=1
      if (debug) write(0,*) 'extract: info ',info,
     +     desc_data(psb_dec_type_)
      pointer_dep_list=1
      if (desc_data(psb_dec_type_).eq.psb_desc_bld_) then
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
              info = 3999
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
      else if (desc_data(psb_dec_type_).eq.psb_desc_upd_) then
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
      endif

      length_dl(me)=pointer_dep_list-1

c     ... check for errors...
 998  continue 
      if (debug) write(0,*) 'extract: info ',info
      err = info
c$$$      call igamx2d(icontxt, all, topdef, ione, ione, err, ione, 
c$$$     +  i, i, -ione ,-ione,-ione)

      if (err.ne.0) goto 9999

      if (.true.) then 
        call igsum2d(icontxt,'all',' ',np+1,1,length_dl,np+1,-1,-1)
        call blacs_get(icontxt,10,icomm )
        allocate(itmp(dl_lda))
        itmp(1:dl_lda) = dep_list(1:dl_lda,me)
        call mpi_allgather(itmp,dl_lda,mpi_integer,
     +    dep_list,dl_lda,mpi_integer,icomm,info)
        deallocate(itmp)

      else
        
        if (me.eq.psb_root_) then
          do proc=0,np-1
            if (proc.ne.psb_root_) then
              if (debug) write(0,*) 'receiving from: ',proc
c              ...receive from proc length of its dependence list....
              call igerv2d(icontxt,1,1,length_dl(proc),1,
     +          proc,mycol)

c              ...receive from proc its dependence list....
              call igerv2d(icontxt,length_dl(proc),1,
     +          dep_list(1,proc),length_dl(proc),proc,mycol)

            endif
          enddo
        else if (me.ne.psb_root_) then
c        ...send to root dependence list length.....
          if (debug) write(0,*) 'sending to: ',me,psb_root_
          call igesd2d(icontxt,1,1,length_dl(me),1,psb_root_,mycol)
          if (debug) write(0,*) 'sending to: ',me,psb_root_
c        ...send to root dependence list....
          call igesd2d(icontxt,length_dl(me),1,dep_list(1,me),
     +      length_dl(me),psb_root_,mycol)

        endif
      end if 

      return

 9999 continue
      call fcpsb_errpush(info,name,int_err)
      if(err_act.eq.act_abort) then
         call fcpsb_perror(icontxt)
      endif
      return

      end
