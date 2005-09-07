      INTEGER FUNCTION PSI_EXIST_OVR_ELEM(OVR_ELEM,
     +   DIM_LIST,ELEM_SEARCHED)

C    PURPOSE:
C    =======
C
C    If ELEM_SEARCHED exist in the list OVR_ELEM returns its position in
C    the list, else returns -1
C
C
C     INPUT
C     ======
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
C     ...Array Parameters....
      INTEGER OVR_ELEM(*)

C     ....Scalars parameters....
      INTEGER DIM_LIST,ELEM_SEARCHED
      
C     ...Local Scalars....
      INTEGER I

      I=1
      DO WHILE ((I.LE.DIM_LIST).AND.(OVR_ELEM(I).NE.ELEM_SEARCHED))
         I=I+2
      ENDDO
      IF ((I.LE.DIM_LIST).AND.(OVR_ELEM(I).EQ.ELEM_SEARCHED)) THEN
         PSI_EXIST_OVR_ELEM=I
      ELSE
         PSI_EXIST_OVR_ELEM=-1
      ENDIF
      END
         
