subroutine psi_crea_index(desc_a,index_in,index_out,glob_idx,info)
  
  use psb_realloc_mod
  use psb_descriptor_type
  use psb_error_mod
  implicit none

  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info
  integer, intent(in)              :: index_in(:)
  integer, pointer                 :: index_out(:)
  logical                          :: glob_idx
  
!         ....local scalars...      
  integer    :: me,npcol,mycol,nprow,i,j,k,&
       & mode, int_err(5), err, err_act, np,&
       & dl_lda, icontxt
!         ...parameters...
  integer, pointer     :: dep_list(:,:), length_dl(:)
  integer,parameter    :: root=0,no_comm=-1
  logical,parameter    :: debug=.false.
  character(len=20)    :: name, ch_err

  interface
     subroutine psi_compute_size(desc_data,&
          & index_in, dl_lda, info)
       integer  :: info, dl_lda
       integer  :: desc_data(:), index_in(:)
     end subroutine psi_compute_size
  end interface

  interface
     subroutine psi_sort_dl(dep_list,l_dep_list,np,info)
       integer :: np,dep_list(:,:), l_dep_list(:), info
     end subroutine psi_sort_dl
  end interface

  interface
     subroutine psi_desc_index(desc_data,index_in,dep_list,&
          & length_dl,loc_to_glob,glob_to_loc,desc_index,&
          & isglob_in,info)
       integer :: desc_data(:),index_in(:),dep_list(:)
       integer :: loc_to_glob(:),glob_to_loc(:)
       integer,pointer :: desc_index(:)
       integer :: length_dl, info
       logical :: isglob_in
     end subroutine psi_desc_index
  end interface

  info = 0
  name='psi_crea_index'
  call psb_erractionsave(err_act)

  icontxt = desc_a%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,np,npcol,me,mycol)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  else if (npcol /= 1) then
    info = 2030
    int_err(1) = npcol
    call psb_errpush(info,name)
    goto 9999
  endif

  ! allocate dependency list
  call psi_compute_size(desc_a%matrix_data, index_in, dl_lda, info)
  allocate(dep_list(max(1,dl_lda),0:np-1),length_dl(0:np-1))
  ! ...extract dependence list (ordered list of identifer process
  !    which every process must communcate with...
  if (debug) write(*,*) 'crea_halo: calling extract_dep_list'
  mode = 1
  call psi_extract_dep_list(desc_a%matrix_data,index_in,&
       & dep_list,length_dl,np,dl_lda,mode,info)
  if(info /= 0) then
     call psb_errpush(4010,name,a_err='extrct_dl')
     goto 9999
  end if

  if (debug) write(*,*) 'crea_index: from extract_dep_list',&
       &     me,length_dl(0),index_in(1), ':',dep_list(:length_dl(me),me)
  ! ...now process root contains dependence list of all processes...
  if (debug) write(*,*) 'crea_halo: root sorting dep list'

  ! ....i must order communication in in halo
  call psi_dl_check(dep_list,dl_lda,np,length_dl)

  ! ....now i can sort dependence list......
  call psi_sort_dl(dep_list,length_dl,np,info)
  if(info.ne.0) then
     call psb_errpush(4010,name,a_err='psi_sort_dl')
     goto 9999
  end if

  ! ...create desc_halo array.....
  if(debug) write(0,*)'in psi_crea_index calling psi_desc_index',&
       & size(index_out)
  call psi_desc_index(desc_a%matrix_data,index_in,dep_list(1:,me),&
       & length_dl(me),desc_a%loc_to_glob,desc_a%glob_to_loc,&
       & index_out,glob_idx,info)
  if(info.ne.0) then
     call psb_errpush(4010,name,a_err='psi_desc_index')
     goto 9999
  end if
  
  deallocate(dep_list,length_dl)
  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return
end subroutine psi_crea_index
