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
!
! Purpose: 
!  Provide a set of subroutines to define a data distribution based on 
!  a graph partitioning routine.
! 
!  Subroutines:
!  
!  BUILD_GRPPART(A,NPARTS): This subroutine will be called by the root
!    process to build define the data distribuition mapping. 
!      Input parameters:
!        TYPE(D_SPMAT) :: A   The input matrix. The coefficients are
!                             ignored; only the structure is used.
!        INTEGER       :: NPARTS  How many parts we are requiring to the 
!                                 partition utility
! 
!  DISTR_GRPPART(RROOT,CROOT,ICTXT): This subroutine will be called by
!      all processes to distribute the information computed by the root
!      process, to be used subsequently.
!
!
!  PART_GRAPH : The subroutine to be passed to PSBLAS sparse library;
!      uses information prepared by the previous two subroutines.
!
module psb_graphpart_mod
  public part_graph, build_grppart, distr_grppart,&
       & getv_grppart, build_usrpart, free_part
  private 
  integer, allocatable, save :: graph_vect(:)
  
contains
  
  subroutine part_graph(global_indx,n,np,pv,nv)
    
    integer, intent(in)  :: global_indx, n, np
    integer, intent(out) :: nv
    integer, intent(out) :: pv(*)
    
    IF (.not.allocated(graph_vect)) then
       write(0,*) 'Fatal error in PART_GRAPH: vector GRAPH_VECT ',&
	    & 'not initialized'
       return
    endif
    if ((global_indx<1).or.(global_indx > size(graph_vect))) then       
       write(0,*) 'Fatal error in PART_GRAPH: index GLOBAL_INDX ',&
	    & 'outside GRAPH_VECT bounds',global_indx,size(graph_vect)
       return
    endif
    nv = 1
    pv(nv) = graph_vect(global_indx)
    return
  end subroutine part_graph


  subroutine distr_grppart(root, ictxt)
    use psb_base_mod
    integer    :: root, ictxt
    integer    :: n, me, np

    call psb_info(ictxt,me,np)

    if (.not.((root>=0).and.(root<np))) then 
      write(0,*) 'Fatal error in DISTR_GRPPART: invalid ROOT  ',&
           & 'coordinates '
      call psb_abort(ictxt)
      return
    endif

    if (me == root) then 
      if (.not.allocated(graph_vect)) then
        write(0,*) 'Fatal error in DISTR_GRPPART: vector GRAPH_VECT ',&
             & 'not initialized'
        call psb_abort(ictxt)
        return
      endif
      n = size(graph_vect)
      call psb_bcast(ictxt,n,root=root)
    else 
      call psb_bcast(ictxt,n,root=root)

      allocate(graph_vect(n),stat=info)
      if (info /= 0) then
        write(0,*) 'Fatal error in DISTR_GRPPART: memory allocation ',&
             & ' failure.'
        return
      endif
    endif
    call psb_bcast(ictxt,graph_vect(1:n),root=root)

    return

  end subroutine distr_grppart
  
  subroutine  getv_grppart(ivg)
    integer, allocatable, intent(out)  :: ivg(:)
    if (allocated(graph_vect)) then 
      allocate(ivg(size(graph_vect)))
      ivg(:) = graph_vect(:)
    end if
  end subroutine getv_grppart
  

  subroutine build_grppart(n,fida,ia1,ia2,nparts)
    use psb_base_mod
    integer       :: nparts
    integer       :: ia1(:), ia2(:)
    integer       :: n, i, ib, ii,numflag,nedc,wgflag
    character(len=5)     :: fida
    integer, parameter :: nb=512
    real(kind(1.d0)), parameter :: seed=12345.d0
    real(kind(1.d0)) :: XV(NB)
    integer          :: iopt(10),idummy(2),jdummy(2)
    interface 
      subroutine METIS_PartGraphRecursive(n,ixadj,iadj,ivwg,iajw,&
           & wgflag,numflag,nparts,iopt,nedc,part)
        integer :: n,wgflag,numflag,nparts,nedc
        integer :: ixadj(*),iadj(*),ivwg(*),iajw(*),iopt(*),part(*)
      end subroutine METIS_PartGraphRecursive
    end interface    
    
    allocate(graph_vect(n),stat=info)
    
    if (info /= 0) then
       write(0,*) 'Fatal error in BUILD_GRPPART: memory allocation ',&
	    & ' failure.'
       return
    endif
    if (nparts.gt.1) then
      if (toupper(fida).eq.'CSR') then 
        iopt(1) = 0
        numflag  = 1
        wgflag   = 0

        call METIS_PartGraphRecursive(n,ia2,ia1,idummy,jdummy,&
             & wgflag,numflag,nparts,iopt,nedc,graph_vect)

        do i=1, n
          graph_vect(i) = graph_vect(i) - 1
        enddo
      else
        write(0,*) 'Fatal error in BUILD_GRPPART: matrix format ',&
             & ' failure. ', FIDA
        return
      endif
    else
      do i=1, n
        graph_vect(i) = 0
      enddo
    endif
    
    return

  end subroutine build_grppart 

  subroutine build_usrpart(n,v,nparts)
    integer       :: nparts
    integer       :: v(:)
    integer       :: n, i, ib, ii,numflag,nedc,wgflag

    if ((n<=0) .or. (nparts<1)) then 
      write(0,*) 'Invalid input to BUILD_USRPART ',n,nparts
      return
    endif
    
    allocate(graph_vect(n),stat=info)
    
    if (info /= 0) then
       write(0,*) 'Fatal error in BUILD_USRPART: memory allocation ',&
	    & ' failure.'
       return
    endif

    do i=1, n
      if ((0<=v(i)).and.(v(i)<nparts)) then 
        graph_vect(i) = v(i)
      else
        write(0,*) 'Invalid V input to BUILD_USRPART',i,v(i),nparts
      endif
    end do
        
    return

  end subroutine build_usrpart

  subroutine free_part(info)
    integer :: info
    
    deallocate(graph_vect,stat=info)
    return
  end subroutine free_part    

end module psb_graphpart_mod

