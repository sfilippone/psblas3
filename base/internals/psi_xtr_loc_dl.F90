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
subroutine psi_i_xtr_loc_dl(ictxt,is_bld,is_upd,desc_str,loc_dl,length_dl,info)

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
  !  output
  !  == = ==
  !  loc_dl    integer array(:)
  !            dependence list of current process
  !           
  use psi_mod, psb_protect_name => psi_i_xtr_loc_dl
#ifdef MPI_MOD
  use mpi
#endif
  use psb_penv_mod
  use psb_const_mod
  use psb_error_mod
  use psb_desc_mod
  use psb_sort_mod
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  !     ....scalar parameters...
  logical,  intent(in)            :: is_bld, is_upd
  type(psb_ctxt_type), intent(in) :: ictxt
  integer(psb_ipk_), intent(in)   :: desc_str(:)
  integer(psb_ipk_), allocatable, intent(out) :: loc_dl(:), length_dl(:)
  integer(psb_ipk_), intent(out)  :: info
  !     .....local arrays....
  integer(psb_ipk_) :: int_err(5)

  !     .....local scalars...
  integer(psb_ipk_) :: i,pdl,proc,j,err_act, ldu
  integer(psb_ipk_) :: err
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_ipk_) :: me, np
  logical, parameter :: dist_symm_list=.true., print_dl=.false., profile=.true.
  character  name*20
  name='psi_extrct_dl'

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  info = psb_success_

  call psb_info(ictxt,me,np)
  pdl = size(desc_str)
  allocate(loc_dl(pdl+1),length_dl(0:np),stat=info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    goto 9999
  end if

  if (debug_level >= psb_debug_inner_)&
       & write(debug_unit,*) me,' ',trim(name),': start ',info

  loc_dl = -1
  i      = 1
  pdl    = 0
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
        pdl=pdl+1
        loc_dl(pdl)=proc
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
        pdl=pdl+1
        loc_dl(pdl)=proc
      endif
      i=i+desc_str(i+1)+2
    enddo
    
  else
    info = 2020
    goto 9999
  endif
  call psb_msort_unique(loc_dl(1:pdl),ldu)
  pdl = ldu
  call psb_realloc(pdl,loc_dl,info)
  call psi_symm_dep_list(loc_dl,ictxt,info)       
  pdl = size(loc_dl)
  length_dl     = 0
  length_dl(me) = pdl
  call psb_sum(ictxt,length_dl)
  
  call psb_erractionrestore(err_act)
  return


9999 continue

  call psb_errpush(info,name,i_err=int_err)
  call psb_error_handler(err_act)

  return

end subroutine psi_i_xtr_loc_dl
