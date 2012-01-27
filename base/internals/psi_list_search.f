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
      INTEGER FUNCTION PSI_LIST_SEARCH(LIST,LENGHT_LIST,ELEM)
      use psb_const_mod
C      !RETURNS POSITION OF ELEM IN A ARRAY LIST
C      !OF LENGHT LENGHT_LIST, IF THIS ELEMENT NOT EXISTS
C      !RETURNS -1
      INTEGER(psb_ipk_) :: LIST(*)
      INTEGER(psb_ipk_) :: LENGHT_LIST
      INTEGER(psb_ipk_) :: ELEM
      
      INTEGER(psb_ipk_) :: I

      I=1
      DO WHILE ((I.LE.LENGHT_LIST).AND.(LIST(I).NE.ELEM))
         I=I+1
      ENDDO
      IF (I.LE.LENGHT_LIST)  THEN 
         IF (LIST(I).EQ.ELEM) THEN
            PSI_LIST_SEARCH=I
         ELSE
            PSI_LIST_SEARCH=-1
         ENDIF
      ELSE
         PSI_LIST_SEARCH=-1
      ENDIF
      END
         
