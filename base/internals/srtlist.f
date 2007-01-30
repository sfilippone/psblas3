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
***********************************************************************
*                                                                     *
*               The communication step among processors at each       *
*               matrix-vector  product is a variable all-to-all       *
*               collective communication that we reimplement          *
*               in terms of point-to-point communications.            *
*               The data in input is a list of dependencies:          *
*               for each node a list of all the nodes it has to       *
*               communicate with. The lists are guaranteed to be      *
*               symmetric, i.e. for each pair (I,J) there is a        *
*               pair (J,I). The idea is to organize the ordering      *
*               so   that at each communication step as many          *
*               processors as possible are communicating at the       *
*               same time, i.e. a step is defined by the fact         *
*               that all edges (I,J) in it have no common node.       *
*                                                                     *
*               Formulation of the problem is:                        *
*                 Given an undirected graph (forest):                 *
*                 Find the shortest series of steps to cancel all     *
*                 graph edges, where at each step all edges belonging *
*                 to a matching in the graph are canceled.            *
*                                                                     *
*               An obvious lower bound to the optimum number of steps *
*               is the largest degree of any node in the graph.       *
*                                                                     *
*               The algorithm proceeds as follows:                    *
*               1. Build a list of all edges, e.g. copy the           *
*                  dependencies lists keeping only (I,J) with I<J     *
*               2. Compute an auxiliary vector with the degree of     *
*                  each node of the graph.                            *
*               3. While there are edges in the graph do:             *
*               4. Weight the edges with the sum of the degrees       *
*                  of their nodes and sort them into descending order *
*               5. Scan the list of edges; if neither node of the     *
*                  edge has been marked yet, cancel the edge and mark *
*                  the two nodes                                      *
*               6. If no edge was chosen but the graph is nonempty    *
*                  raise an error condition                           *
*               7. Queue the edges in the matchin to the output       *
*                  sequence;                                          *
*               8. Decrease by 1 the degree of all marked nodes,      *
*                  then clear all marks                               *
*               9. Cycle to 3.                                        *
*              10. For each node: scan the edge sequence; if an       *
*                  edge has the node as an endpoint, queue the other  *
*                  node in the dependency list for the current one    *
*                                                                     *
***********************************************************************
      SUBROUTINE SRTLIST(DEP_LIST,DL_LDA,LDL,NP,dg,dgp,upd,
     +  edges,idx,ich,INFO)
      IMPLICIT NONE
      INTEGER  NP, DL_LDA, INFO
      INTEGER  DEP_LIST(DL_LDA,*), LDL(*),DG(*), DGP(*), IDX(*),
     +  UPD(*),EDGES(2,*),ICH(*)
      INTEGER  I,J, NEDGES,IP1,IP2,NCH,IP,IEDGE,I1,IX,IST,ISWAP(2)
      INTEGER NO_COMM
      PARAMETER (NO_COMM=-1)
      
      
      IF (NP .LT. 0) THEN 
        INFO = 1
        RETURN
      ENDIF
      
C
C  dg contains number of communications 
C
      DO I=1, NP
        DG(I)=LDL(I)
      ENDDO

      NEDGES = 0
      DO I=1, NP
        DO J=1, DG(I) 
          IP = DEP_LIST(J,I) + 1            
c$$$            write(0,*) 'SRTLIST Input :',i,ip
          IF (IP.GT.I)
     +      NEDGES = NEDGES + 1
        ENDDO
      ENDDO
      
      IEDGE = 0 
      DO I=1, NP
        DO J=1, DG(I) 
          IP = DEP_LIST(J,I) + 1            
          IF (IP.GT.I) THEN 
            IEDGE = IEDGE + 1
            EDGES(1,IEDGE) = I
            EDGES(2,IEDGE) = IP
          ENDIF
        ENDDO
      ENDDO
      
      IST = 1
      DO WHILE (IST.LE.NEDGES)         
        
        DO I=1, NP
          UPD(I) = 0      
        ENDDO
        DO I=IST, NEDGES
          DGP(I) = -(DG(EDGES(1,I))+DG(EDGES(2,I)))
        ENDDO

        CALL ISRX(NEDGES-IST+1,DGP(IST),IDX(IST))
        I1 = IST         
        NCH = 0
        DO I = IST, NEDGES
          IX = IDX(I)+IST-1
          IP1 = EDGES(1,IX)
          IP2 = EDGES(2,IX)
          IF ((UPD(IP1).eq.0).AND.(UPD(IP2).eq.0)) THEN 
            UPD(IP1) = -1
            UPD(IP2) = -1
            NCH      = NCH + 1
            ICH(NCH) = IX
          ENDIF
        ENDDO
        IF (NCH.eq.0) THEN
          write(0,*) 'SRTLIST ?????? Impossible error !!!!!?????',
     +      nedges,ist
          do i=ist, nedges
            IX = IDX(I)+IST-1
            write(0,*) 'SRTLIST: Edge:',ix,edges(1,ix),
     +        edges(2,ix),dgp(ix)
          enddo
          info = 30
          return
        ENDIF
        CALL ISR(NCH,ICH)
        DO I=1, NCH
          ISWAP(1)        = EDGES(1,IST)
          ISWAP(2)        = EDGES(2,IST)
          EDGES(1,IST)    = EDGES(1,ICH(I))
          EDGES(2,IST)    = EDGES(2,ICH(I))
          EDGES(1,ICH(I)) = ISWAP(1)
          EDGES(2,ICH(I)) = ISWAP(2)
          IST = IST + 1
        ENDDO    
        DO I=1, NP 
          DG(I) = DG(I) + UPD(I) 
        ENDDO
      ENDDO
      
      DO I=1, NP
        IF (DG(I).NE.0) THEN 
          WRITE(0,*) 'SRTLIST Error on exit:',i,dg(i)
        ENDIF
        DG(I) = 0
      ENDDO
      DO J=1,NEDGES
        I = EDGES(1,J) 
        DG(I) = DG(I)+1
        DEP_LIST(DG(I),I) = EDGES(2,J)-1
        I = EDGES(2,J)
        DG(I) = DG(I)+1
        DEP_LIST(DG(I),I) = EDGES(1,J)-1
      ENDDO
      DO I=1, NP
        IF (DG(I).NE.LDL(I)) THEN 
          WRITE(0,*) 'SRTLIST Mismatch on output',i,dg(i),ldl(i)
        ENDIF
      ENDDO
      
c$$$      WRITE(0,*) 'Output communication:',t2-t1
c$$$      do i=1,np
c$$$         do j=1,ldl(i)
c$$$            write(0,*)'SRTLIST', i,dep_list(j,i)+1
c$$$         enddo
c$$$      enddo

      RETURN           
      END


