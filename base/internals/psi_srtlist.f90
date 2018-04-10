!
!             Parallel Sparse BLAS  version 3.5
!   (C) Copyright 2006, 2010, 2015, 2017
!       Salvatore Filippone     
!       Alfredo Buttari        CNRS-IRIT, Toulouse
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!   1. Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!   2. Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions, and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!   3. The name of the PSBLAS group or the names of its contributors may
!      not be used to endorse or promote products derived from this
!      software without specific written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
! 
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
subroutine srtlist(dep_list,dl_lda,ldl,np,dg,dgp,upd, edges,idx,ich,info)
  use psb_serial_mod
  implicit none
  integer(psb_ipk_) ::  np, dl_lda, info
  integer(psb_ipk_) ::  dep_list(dl_lda,*), ldl(*),dg(*), dgp(*),&
       &  idx(*), upd(*),edges(2,*),ich(*)
  integer(psb_ipk_) ::  i,j, nedges,ip1,ip2,nch,ip,iedge,&
       &  i1,ix,ist,iswap(2)
  integer(psb_ipk_) :: no_comm
  parameter (no_comm=-1)


  if (np .lt. 0) then 
    info = 1
    return
  endif

  !
  !  dg contains number of communications 
  !
  do i=1, np
    dg(i)=ldl(i)
  enddo


  nedges = 0
  do i=1, np
    do j=1, dg(i) 
      ip = dep_list(j,i) + 1            
      if (ip.gt.i)  nedges = nedges + 1
    enddo
  enddo

  iedge = 0 
  do i=1, np
    do j=1, dg(i) 
      ip = dep_list(j,i) + 1            
      if (ip.gt.i) then 
        iedge = iedge + 1
        edges(1,iedge) = i
        edges(2,iedge) = ip
      endif
    enddo
  enddo

  ist = 1
  do while (ist.le.nedges)         

    do i=1, np
      upd(i) = 0      
    enddo
    do i=ist, nedges
      dgp(i) = -(dg(edges(1,i))+dg(edges(2,i)))
    enddo

    call psb_msort(dgp(ist:nedges),ix=idx(ist:nedges))
    i1 = ist         
    nch = 0
    do i = ist, nedges
      ix = idx(i)+ist-1
      ip1 = edges(1,ix)
      ip2 = edges(2,ix)
      if ((upd(ip1).eq.0).and.(upd(ip2).eq.0)) then 
        upd(ip1) = -1
        upd(ip2) = -1
        nch      = nch + 1
        ich(nch) = ix
      endif
    enddo
    if (nch.eq.0) then
      write(psb_err_unit,*)&
           & 'srtlist ?????? impossible error !!!!!?????',&
           &  nedges,ist
      do i=ist, nedges
        ix = idx(i)+ist-1
        write(psb_err_unit,*)&
             & 'SRTLIST: Edge:',ix,edges(1,ix),&
             & edges(2,ix),dgp(ix)
      enddo
      info = psb_err_input_value_invalid_i_
      return
    endif
    call psb_msort(ich(1:nch))
    do i=1, nch
      iswap(1)        = edges(1,ist)
      iswap(2)        = edges(2,ist)
      edges(1,ist)    = edges(1,ich(i))
      edges(2,ist)    = edges(2,ich(i))
      edges(1,ich(i)) = iswap(1)
      edges(2,ich(i)) = iswap(2)
      ist = ist + 1
    enddo
    do i=1, np 
      dg(i) = dg(i) + upd(i) 
    enddo
  enddo

  do i=1, np
    if (dg(i).ne.0) then 
      write(psb_err_unit,*)&
           &  'SRTLIST Error on exit:',i,dg(i)
    endif
    dg(i) = 0
  enddo
  do j=1,nedges
    i = edges(1,j) 
    dg(i) = dg(i)+1
    dep_list(dg(i),i) = edges(2,j)-1
    i = edges(2,j)
    dg(i) = dg(i)+1
    dep_list(dg(i),i) = edges(1,j)-1
  enddo
  do i=1, np
    if (dg(i).ne.ldl(i)) then 
      write(psb_err_unit,*) &
           &   'SRTLIST Mismatch on output',i,dg(i),ldl(i)
    endif
  enddo


  return           
end subroutine srtlist


