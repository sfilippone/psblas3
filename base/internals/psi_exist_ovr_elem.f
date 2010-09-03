C
C             Parallel Sparse BLAS  version 3.0
C   (C) Copyright 2006, 2007, 2008, 2009, 2010
C                      Salvatore Filippone    University of Rome Tor Vergata
C                      Alfredo Buttari        CNRS-IRIT, Toulouse
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
      INTEGER FUNCTION PSI_EXIST_OVR_ELEM(OVR_ELEM,
     +   DIM_LIST,ELEM_SEARCHED)

C    PURPOSE:
C    == = ====
C
C    If ELEM_SEARCHED exist in the list OVR_ELEM returns its position in
C    the list, else returns -1
C
C
C     INPUT
C     == = ===
C     OVRLAP_ELEMENT_D.: Contains for all overlap points belonging to 
C                        the current process:
C                          1. overlap point index
C                          2. Number of domains sharing that overlap point
C                        the end is marked by a -1...............................
C
C    DIM_LIST..........: Dimension of list OVRLAP_ELEMENT_D
C
C    ELEM_SEARCHED.....:point's  Local index identifier to be searched.

      IMPLICIT NONE

C     ....Scalars parameters....
      INTEGER DIM_LIST,ELEM_SEARCHED
C     ...Array Parameters....
      INTEGER OVR_ELEM(DIM_LIST,*)
      
C     ...Local Scalars....
      INTEGER I

      I=1
      DO WHILE ((I.LE.DIM_LIST).AND.(OVR_ELEM(I,1).NE.ELEM_SEARCHED))
         I=I+1
      ENDDO
      IF ((I.LE.DIM_LIST).AND.(OVR_ELEM(I,1).EQ.ELEM_SEARCHED)) THEN
         PSI_EXIST_OVR_ELEM=I
      ELSE
         PSI_EXIST_OVR_ELEM=-1
      ENDIF
      END
         
