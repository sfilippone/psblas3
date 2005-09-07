subroutine psi_sort_dl(dep_list,l_dep_list,np,info)
  !
  !     interface between former sort_dep_list subroutine
  !     and new srtlist
  !
  use psb_error_mod
  implicit none

  integer :: np,dep_list(:,:), l_dep_list(:)
  integer :: idg, iupd, idgp, iedges, iidx, iich,ndgmx, isz, err_act
  integer :: i, info
  integer, pointer   :: work(:)
  logical, parameter :: debug=.false.
  character(len=20)        :: name, ch_err
  
  name='psi_sort_dl'
  info=0
  call psb_erractionsave(err_act)
  
  info = 0
  ndgmx = 0
  do i=1,np
     ndgmx = ndgmx + l_dep_list(i)
     if (debug) write(0,*) i,l_dep_list(i)
  enddo
  idg = 1
  iupd = idg+np
  idgp = iupd+np
  iedges = idgp + ndgmx
  iidx = iedges + 2*ndgmx
  iich = iidx + ndgmx
  isz = iich + ndgmx
  if (debug)write(0,*) 'psi_sort_dl: ndgmx ',ndgmx,isz
  
  allocate(work(isz))
  ! call srtlist(dep_list, dl_lda, l_dep_list, np, info)
  call srtlist(dep_list,size(dep_list,1),l_dep_list,np,work(idg),&
       & work(idgp),work(iupd),work(iedges),work(iidx),work(iich),info)

  if (info .ne. 0) then
     call psb_errpush(4010,name,a_err='srtlist')
     goto 9999
  endif
  
  deallocate(work)
  call psb_erractionrestore(err_act)
  return  
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return
end subroutine psi_sort_dl


      
