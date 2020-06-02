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
subroutine psi_i_sort_dl(dep_list,l_dep_list,np,info)
  !
  !     interface between former sort_dep_list subroutine
  !     and new srtlist
  !
  use psi_mod, psb_protect_name => psi_i_sort_dl
  use psb_const_mod
  use psb_error_mod
  implicit none

  integer(psb_ipk_) :: np,dep_list(:,:), l_dep_list(:)
  integer(psb_ipk_) :: idg, iupd, idgp, iedges, iidx, iich,ndgmx, isz, err_act
  integer(psb_ipk_) :: i, info
  integer(psb_ipk_), allocatable :: work(:)
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)    :: name
  
  name='psi_sort_dl'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  
  info = psb_success_
  ndgmx = 0
  do i=1,np
     ndgmx = ndgmx + l_dep_list(i)
     if (debug_level >= psb_debug_inner_)&
          & write(debug_unit,*) name,': ',i,l_dep_list(i)
  enddo
  idg = 1
  iupd = idg+np
  idgp = iupd+np
  iedges = idgp + ndgmx
  iidx = iedges + 2*ndgmx
  iich = iidx + ndgmx
  isz = iich + ndgmx
  if (debug_level >= psb_debug_inner_)&
       & write(debug_unit,*) name,': ndgmx ',ndgmx,isz
  allocate(work(isz))
  ! call srtlist(dep_list, dl_lda, l_dep_list, np, info)
  call srtlist(dep_list,size(dep_list,1,kind=psb_ipk_),l_dep_list,np,work(idg),&
       & work(idgp),work(iupd),work(iedges),work(iidx),work(iich),info)
  if (info  /=  psb_success_) then
     call psb_errpush(psb_err_from_subroutine_,name,a_err='srtlist')
     goto 9999
  endif
  
  deallocate(work)
  call psb_erractionrestore(err_act)
  return  
  
9999 call psb_error_handler(err_act)

  return

end subroutine psi_i_sort_dl

!**********************************************************************
!                                                                     *
!               The communication step among processors at each       *
!               matrix-vector  product is a variable all-to-all       *
!               collective communication that we reimplement          *
!               in terms of point-to-point communications.            *
!               The data in input is a list of dependencies:          *
!               for each node a list of all the nodes it has to       *
!               communicate with. The lists are guaranteed to be      *
!               symmetric, i.e. for each pair (I,J) there is a        *
!               pair (J,I). The idea is to organize the ordering      *
!               so   that at each communication step as many          *
!               processors as possible are communicating at the       *
!               same time, i.e. a step is defined by the fact         *
!               that all edges (I,J) in it have no common node.       *
!                                                                     *
!               Formulation of the problem is:                        *
!                 Given an undirected graph (forest):                 *
!                 Find the shortest series of steps to cancel all     *
!                 graph edges, where at each step all edges belonging *
!                 to a matching in the graph are canceled.            *
!                                                                     *
!               An obvious lower bound to the optimum number of steps *
!               is the largest degree of any node in the graph.       *
!                                                                     *
!               The algorithm proceeds as follows:                    *
!               1. Build a list of all edges, e.g. copy the           *
!                  dependencies lists keeping only (I,J) with I<J     *
!               2. Compute an auxiliary vector with the degree of     *
!                  each node of the graph.                            *
!               3. While there are edges in the graph do:             *
!               4. Weight the edges with the sum of the degrees       *
!                  of their nodes and sort them into descending order *
!               5. Scan the list of edges; if neither node of the     *
!                  edge has been marked yet, cancel the edge and mark *
!                  the two nodes                                      *
!               6. If no edge was chosen but the graph is nonempty    *
!                  raise an error condition                           *
!               7. Queue the edges in the matchin to the output       *
!                  sequence;                                          *
!               8. Decrease by 1 the degree of all marked nodes,      *
!                  then clear all marks                               *
!               9. Cycle to 3.                                        *
!              10. For each node: scan the edge sequence; if an       *
!                  edge has the node as an endpoint, queue the other  *
!                  node in the dependency list for the current one    *
!                                                                     *
!**********************************************************************
subroutine psi_i_csr_sort_dl(dl_ptr,c_dep_list,l_dep_list,ictxt,info)
  use psi_mod, psb_protect_name => psi_i_csr_sort_dl
  use psb_const_mod
  use psb_error_mod
  use psb_sort_mod
  implicit none
  
  integer(psb_ipk_), intent(inout) :: c_dep_list(:), dl_ptr(0:), l_dep_list(0:)
  integer(psb_ipk_), intent(in)    :: ictxt
  integer(psb_ipk_), intent(out)   :: info
  ! Local variables
  integer(psb_ipk_), allocatable  ::  dg(:), dgp(:),&
       &  idx(:), upd(:), edges(:,:), ich(:)
  integer(psb_ipk_) ::  i, j, nedges, ip1, ip2, nch, ip, iedge,&
       &  i1, ix, ist, iswap(2)
  logical :: internal_error 
  integer(psb_ipk_) :: me, np

  info = 0
  call psb_info(ictxt,me,np)
  nedges = size(c_dep_list)
  
  allocate(dg(0:np-1),dgp(nedges),edges(2,nedges),upd(0:np-1),&
       & idx(nedges),ich(nedges),stat = info)
  
  if (info /= 0) then
    info = -9
    return
  end if
  !
  !  1. Compute an auxiliary vector with the degree of 
  !     each node of the graph.
  dg(0:np-1) = l_dep_list(0:np-1)
  !
  !  2. Build a list of all edges, e.g. copy the      
  !     dependencies lists keeping only (I,J) with I<J
  !
  nedges  = 0
  do i = 0, np-1
    do j = dl_ptr(i),dl_ptr(i+1) - 1
      ip = c_dep_list(j)
      if (i<=ip) then
        nedges = nedges + 1
        edges(1,nedges) = i
        edges(2,nedges) = ip
      end if
    end do
  end do

  !
  !   3. Loop over all edges
  !   
  ist = 1 
  do while (ist <= nedges)
    !
    !  4. Weight the edges with the sum of the degrees       
    !     of their nodes and sort them into descending order 
    upd(:) = 0
    do i = ist, nedges
      dgp(i) = (dg(edges(1,i)) + dg(edges(2,i)))
    end do
    call psb_msort(dgp(ist:nedges),ix=idx(ist:nedges),dir=psb_sort_down_)
    
    !  5. Scan the list of edges; if neither node of the     
    !     edge has been marked yet, take out the edge and mark 
    !     the two nodes                                      
    i1 = ist         
    nch = 0
    do i = ist, nedges
      ix  = idx(i)+ist-1
      ip1 = edges(1,ix)
      ip2 = edges(2,ix)
      if ((upd(ip1)==0).and.(upd(ip2)==0)) then 
        upd(ip1) = -1
        upd(ip2) = -1
        nch      = nch + 1
        ich(nch) = ix
      end if
    end do
    !
    !  6. If no edge was chosen but the graph is nonempty 
    !     raise an error condition                        
    if (nch == 0) then
      write(psb_err_unit,*)&
           & 'srtlist ?????? impossible error !!!!!?????',&
           &  nedges,ist
      do i=ist, nedges
        ix = idx(i)+ist-1
        write(psb_err_unit,*)&
             & 'SRTLIST: Edge:',ix,edges(1,ix),&
             & edges(2,ix),dgp(ix)
      end do
      info = psb_err_input_value_invalid_i_
      return
    end if
    !
    ! 7. Queue the edges in the matching to the output  
    !    sequence; decrease by 1 the degree of all marked
    !    nodes, then clear all marks
    !    
    call psb_msort(ich(1:nch))
    do i=1, nch
      iswap(1:2)        = edges(1:2,ist)
      edges(1:2,ist)    = edges(1:2,ich(i))
      edges(1:2,ich(i)) = iswap(1:2)
      ist               = ist + 1 
    end do
    do i=0, np-1
      dg(i) = dg(i) + upd(i) 
    end do
  end do
  internal_error = .false.
  do i=0, np-1
    if (dg(i) /= 0) then
      internal_error = .true.
      if (me == 0) write(psb_err_unit,*)&
           &  'csr_SRTLIST Error on exit:',i,dg(i)
    end if
    dg(i) = 0
  end do
  if (internal_error .and. (me==0)) then
    write(0,*) 'Error on srt_list. Input:'
    do i = 0, np-1
      write(0,*) 'Proc: ',i,' list: '
      write(0,*) c_dep_list(dl_ptr(i):dl_ptr(i+1) - 1)
    end do
  end if
  !
  ! 10. Scan the edge sequence;
  !     for each edge, take each one of its
  !     endpoints and queue the other
  !     node in the endpoint dependency list
  !     
  do j=1,nedges
    i     = edges(1,j)
    ix    = dl_ptr(i)
    c_dep_list(ix+dg(i)) = edges(2,j)
    dg(i) = dg(i)+1
    
    i     = edges(2,j)
    ix    = dl_ptr(i)
    c_dep_list(ix+dg(i)) = edges(1,j)
    dg(i) = dg(i)+1
    !
    ! If there are any self loops, adjust for error condition
    ! check
    !
    if (edges(1,j) == edges(2,j)) dg(i) = dg(i) -1 
  end do

  do i=0, np-1
    if (dg(i) /= l_dep_list(i)) then 
      if (me == 0) write(psb_err_unit,*) &
           &   'SRTLIST Mismatch on output',i,dg(i),l_dep_list(i)
    end if
  end do
  
end subroutine psi_i_csr_sort_dl
