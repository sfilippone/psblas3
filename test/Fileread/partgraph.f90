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
MODULE PARTGRAPH
  public part_graph, build_grppart, distr_grppart,&
       & getv_grppart, build_usrpart, free_part
  private 
  INTEGER, POINTER, SAVE :: GRAPH_VECT(:)
  
CONTAINS
  
  SUBROUTINE PART_GRAPH(GLOBAL_INDX,N,NP,PV,NV)
    
    INTEGER, INTENT(IN)  :: GLOBAL_INDX, N, NP
    INTEGER, INTENT(OUT) :: NV
    INTEGER, INTENT(OUT) :: PV(*)
    
    IF (.NOT.ASSOCIATED(GRAPH_VECT)) THEN
       WRITE(0,*) 'Fatal error in PART_GRAPH: vector GRAPH_VECT ',&
	    & 'not initialized'
       RETURN
    ENDIF
    IF ((GLOBAL_INDX<1).OR.(GLOBAL_INDX > SIZE(GRAPH_VECT))) THEN       
       WRITE(0,*) 'Fatal error in PART_GRAPH: index GLOBAL_INDX ',&
	    & 'outside GRAPH_VECT bounds',global_indx,size(graph_vect)
       RETURN
    ENDIF
    NV = 1
    PV(NV) = GRAPH_VECT(GLOBAL_INDX)
    RETURN
  END SUBROUTINE PART_GRAPH


  SUBROUTINE DISTR_GRPPART(RROOT, CROOT, ICTXT)
    INTEGER    :: RROOT, CROOT, ICTXT
    INTEGER    :: N, MER, MEC, NPR, NPC
    
    CALL BLACS_GRIDINFO(ICTXT,NPR,NPC,MER,MEC)
    
    IF (.NOT.((RROOT>=0).AND.(RROOT<NPR).AND.&
	 & (CROOT>=0).AND.(CROOT<NPC))) THEN 
       WRITE(0,*) 'Fatal error in DISTR_GRPPART: invalid ROOT  ',&
	    & 'coordinates '
       CALL BLACS_ABORT(ICTXT,-1)
       RETURN
    ENDIF

    IF ((MER == RROOT) .AND.(MEC == CROOT)) THEN 
       IF (.NOT.ASSOCIATED(GRAPH_VECT)) THEN
	  WRITE(0,*) 'Fatal error in DISTR_GRPPART: vector GRAPH_VECT ',&
	       & 'not initialized'
	  CALL BLACS_ABORT(ICTXT,-1)
	  RETURN
       ENDIF
       N = SIZE(GRAPH_VECT)
       CALL IGEBS2D(ICTXT,'All',' ',1,1,N,1)
       CALL IGEBS2D(ICTXT,'All',' ',N,1,GRAPH_VECT,N)
    ELSE 
       CALL IGEBR2D(ICTXT,'All',' ',1,1,N,1,RROOT,CROOT)
!!$       IF (ASSOCIATED(GRAPH_VECT)) THEN
!!$	  DEALLOCATE(GRAPH_VECT)
!!$       ENDIF
       ALLOCATE(GRAPH_VECT(N),STAT=INFO)
       IF (INFO /= 0) THEN
	  WRITE(0,*) 'Fatal error in DISTR_GRPPART: memory allocation ',&
	       & ' failure.'
	  RETURN
       ENDIF       
       CALL IGEBR2D(ICTXT,'All',' ',N,1,GRAPH_VECT,N,RROOT,CROOT)
    ENDIF

    RETURN
    
  END SUBROUTINE DISTR_GRPPART
  
  subroutine  getv_grppart(ivg)
    integer, pointer :: ivg(:)
    if (associated(graph_vect)) then 
      allocate(ivg(size(graph_vect)))
      ivg(:) = graph_vect(:)
    else
      ivg => null()
    end if
  end subroutine getv_grppart
  

  SUBROUTINE BUILD_GRPPART(N,FIDA,IA1,IA2,NPARTS)
    INTEGER       :: NPARTS
    INTEGER       :: IA1(:), IA2(:)
    INTEGER       :: N, I, IB, II,numflag,nedc,wgflag
    CHARACTER(LEN=5)     :: FIDA
    INTEGER, PARAMETER :: NB=512
    REAL(KIND(1.D0)), PARAMETER :: SEED=12345.D0
    REAL(KIND(1.D0)) :: XV(NB)
    integer          :: iopt(10),idummy(2),jdummy(2)
    interface 
      subroutine METIS_PartGraphRecursive(n,ixadj,iadj,ivwg,iajw,&
           & wgflag,numflag,nparts,iopt,nedc,part)
        integer :: n,wgflag,numflag,nparts,nedc
        integer :: ixadj(*),iadj(*),ivwg(*),iajw(*),iopt(*),part(*)
      end subroutine METIS_PartGraphRecursive
    end interface    
    

!!$    IF (ASSOCIATED(GRAPH_VECT)) THEN
!!$       DEALLOCATE(GRAPH_VECT)
!!$    ENDIF
    
    ALLOCATE(GRAPH_VECT(N),STAT=INFO)
    
    IF (INFO /= 0) THEN
       WRITE(0,*) 'Fatal error in BUILD_GRPPART: memory allocation ',&
	    & ' failure.'
       RETURN
    ENDIF
    IF (NPARTS.GT.1) THEN
      IF (FIDA.EQ.'CSR') THEN 
        iopt(1) = 0
        numflag  = 1
        wgflag   = 0
!!$
!!$        write(0,*)'CSR structure ', size(ia2),size(ia1),&
!!$             & ia2(n+1),minval(ia1(1:ia2(n+1)-1)),maxval(ia1(1:ia2(n+1)-1))
        call METIS_PartGraphRecursive(n,ia2,ia1,idummy,jdummy,&
             & wgflag,numflag,nparts,iopt,nedc,graph_vect)
!!$        write(0,*)'Edge cut from Metis ',nedc
        DO I=1, N
          GRAPH_VECT(I) = GRAPH_VECT(I) - 1
        ENDDO
      ELSE
        WRITE(0,*) 'Fatal error in BUILD_GRPPART: matrix format ',&
             & ' failure. ', FIDA
        RETURN
      ENDIF
    ELSE
      DO I=1, N
        GRAPH_VECT(I) = 0
      ENDDO
    ENDIF
    
    RETURN

  END SUBROUTINE BUILD_GRPPART 

  SUBROUTINE BUILD_USRPART(N,V,NPARTS)
    INTEGER       :: NPARTS
    INTEGER       :: V(:)
    INTEGER       :: N, I, IB, II,numflag,nedc,wgflag
    CHARACTER(LEN=5)     :: FIDA

    if ((n<=0) .or. (nparts<1)) then 
      write(0,*) 'Invalid input to BUILD_USRPART ',n,nparts
      return
    endif


!!$    IF (ASSOCIATED(GRAPH_VECT)) THEN
!!$       DEALLOCATE(GRAPH_VECT)
!!$    ENDIF
    
    ALLOCATE(GRAPH_VECT(N),STAT=INFO)
    
    IF (INFO /= 0) THEN
       WRITE(0,*) 'Fatal error in BUILD_USRPART: memory allocation ',&
	    & ' failure.'
       RETURN
    ENDIF

    do i=1, n
      if ((0<=v(i)).and.(v(i)<nparts)) then 
        graph_vect(i) = v(i)
      else
        write(0,*) 'Invalid V input to BUILD_USRPART',i,v(i),nparts
      endif
    end do
        
    RETURN

  END SUBROUTINE BUILD_USRPART

  subroutine free_part(info)
    integer :: info
    
    deallocate(graph_vect,stat=info)
    return
  end subroutine free_part    

END MODULE PARTGRAPH

