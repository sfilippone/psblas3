*     SUBROUTINE DESYM(NROW,A,JA,IA,AS,JAS,IAS,IAW,NNZERO)                  *
*                                                                           *
*     Purpose                                                               *
*     =======                                                               *
*         Utility routine to convert from symmetric storage                 *       
*      to full format (CSR mode).                                           * 
*                                                                           *
*     Parameter                                                             *
*     =========                                                             *
*     INPUT=                                                                *
*                                                                           *
*     SYMBOLIC NAME: NROW                                                   *
*     POSITION:      Parameter No.1                                         *
*     ATTRIBUTES:    INTEGER                                                *
*     VALUES:        NROW>0                                                 *
*     DESCRIPTION:   On entry NROW specifies the number of rows of the      *
*                    input sparse matrix. The number of column of the input *
*                    sparse matrix mest be the same.                        *
*                    Unchanged on exit.                                     *
*                                                                           *
*     SYMBOLIC NAME: A                                                      *
*     POSITION:      Parameter No.2                                         *
*     ATTRIBUTES:    DOUBLE PRECISION ARRAY of Dimension (NNZERO)           *
*     VALUES:                                                               *
*     DESCRIPTION:   A specifies the values of the input  sparse matrix.    *
*                    This matrix  is stored in CSR mode                     * 
*                    Unchanged on exit.                                     *
*                                                                           *
*     SYMBOLIC NAME: JA                                                     *
*     POSITION:      Parameter No. 3                                        *
*     ATTRIBUTES:    INTEGER  ARRAY(IA(NNZERO))                             *
*     VALUES:         >  0                                                  *
*     DESCRIPTION:   Column indices stored by rows refered to the input     *
*                    sparse matrix.                                         *
*                    Unchanged on exit.                                     *
*                                                                           *
*     SYMBOLIC NAME: IA                                                     *
*     POSITION:      Parameter No. 4                                        *
*     ATTRIBUTES:    INTEGER ARRAY(NROW+1)                                  *
*     VALUES:        >0; increasing.                                        *
*     DESCRIPTION:   Row pointer array: it contains the starting            *
*                    position of each row of A in array JA.                 *
*                    Unchanged on exit.                                     *
*                                                                           *
*     SYMBOLIC NAME: IAW                                                    *
*     POSITION:      Parameter No. 7                                        *
*     ATTRIBUTES:    INTEGER ARRAY of Dimension (NROW+1)                    *
*     VALUES:        >0;                                                    *
*     DESCRIPTION:   Work Area.                                             *
*                                                                           *
*     SYMBOLIC NAME: WORK                                                   *
*     POSITION:      Parameter No. 8                                        *
*     ATTRIBUTES:    REAL*8  ARRAY of Dimension (NROW+1)                    *
*     VALUES:        >0;                                                    *
*     DESCRIPTION:   Work Area.                                             *
*                                                                           *
*     SYMBOLIC NAME: NNZERO                                                 *
*     POSITION:      Parameter No. 9                                        *
*     ATTRIBUTES:    INTEGER                                                *
*     VALUES:        >0;                                                    *
*     DESCRIPTION:   On entry contains: the number of the non zero          *
*                    entry of the input matrix.                             *
*                    Unchanged on exit.                                     *
*      OUTPUT==                                                             *
*                                                                           *
*                                                                           *
*     SYMBOLIC NAME: AS                                                     *
*     POSITION:      Parameter No.5                                         *
*     ATTRIBUTES:    DOUBLE PRECISION ARRAY of Dimension (*)                *
*     VALUES:                                                               *
*     DESCRIPTION:   On exit A specifies the values of the output  sparse   *
*                    matrix.                                                *
*                    This matrix  correspondes to A rapresented in FULL-CSR *
*                    mode                                                   * 
*                                                                           *
*     SYMBOLIC NAME: JAS                                                    *
*     POSITION:      Parameter No. 6                                        *
*     ATTRIBUTES:    INTEGER  ARRAY(IAS(NROW+1)-1)                          *
*     VALUES:         >  0                                                  *
*     DESCRIPTION:   Column indices stored by rows refered to the output    *
*                    sparse matrix.                                         *
*                                                                           *
*     SYMBOLIC NAME: IAS                                                    *
*     POSITION:      Parameter No. S                                        *
*     ATTRIBUTES:    INTEGER ARRAY(NROW+1)                                  *
*     VALUES:        >0; increasing.                                        *
*     DESCRIPTION:   Row pointer array: it contains the starting            *
*                    position of each row of AS in array JAS.               *
*****************************************************************************

C 
      SUBROUTINE ZDESYM(NROW,A,JA,IA,AS,JAS,IAS,AUX,WORK,NNZERO,PTR,
     + NZR, VALUE)
      IMPLICIT NONE
C      .. Scalar Arguments ..                                              
      INTEGER NROW,NNZERO,VALUE,INDEX,PTR,NZR
C     .. Array Arguments ..                                                     
      COMPLEX*16 A(*),AS(*),WORK(*)                                 
      INTEGER IA(*),IAS(*),JAS(*),JA(*),AUX(*)                
C     .. Local Scalars ..                                                       
      INTEGER I,IAW1,IAW2,IAWT,J,JPT,K,KPT,LDIM,NZL,JS, IRET, NEL,DIAGEL
C      REAL*8  BUF
C     ..                                                                        

      NEL = 0
      DIAGEL=0
      
      DO I=1, NNZERO
         IF(JA(I).LE.IA(I)) THEN
            NEL = NEL+1
            AS(I)  = A(I)
            JAS(I) = JA(I)
            IAS(I) = IA(I)
            IF(JA(I).NE.IA(I)) THEN !This control avoids malfunctions in the cases 
                                    ! where the matrix is declared symmetric but all
                                    !his elements are explicitly stored
                                    ! see young1c.mtx from "Matrix Market"
               AS(NNZERO+I)  = A(I)
               JAS(NNZERO+I) = IA(I)
               IAS(NNZERO+I) = JA(I)
            ELSE
               DIAGEL = DIAGEL+1
            END IF
         END IF         
      END DO
         
C     .... Order with key IAS ...
            CALL MRGSRT(2*NNZERO,IAS,AUX,IRET)
            IF (IRET.EQ.0) CALL ZREORDVN(2*NNZERO,AS,IAS,JAS,AUX)
C     .... Order with key JAS ...
            
            I    = 1
            J    = I
            DO WHILE (I.LE.(2*NNZERO))
               DO WHILE ((IAS(J).EQ.IAS(I)).AND.
     +              (J.LE.2*NNZERO))
                  J = J+1
               ENDDO
               NZL = J - I
               CALL MRGSRT(NZL,JAS(I),AUX,IRET)
               IF (IRET.EQ.0) CALL ZREORDVN(NZL,AS(I),IAS(I),JAS(I),
     +                                      AUX)
               I = J
            ENDDO
            NZR = NEL*2 - DIAGEL
            PTR = 2*NNZERO-NZR+1

      RETURN                                                                    
                                                                                
      END                                                                       




