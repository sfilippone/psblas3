C     SUBROUTINE DCREL(TRANS,M,N,DESCRA,A,IA1,IA2,IP1,DESCRN,
C                      AN,IAN1,IAN2,IP2,LAN,LIAN1,LIAN2,
C                      IAUX,LIAUX,IERRV)
C
C     Purpose: CSR to ELL format conversion
C     =======
C
C     Parameter:
C     =========
C
C     ...
C     IAN2  -  Vector: first element is max number of columns in matrices
C              ARN,IAN1, elements to M+1 are column index of diagonal
C              in ARN,IAN1 (in future releases)
C     ...
C
C
      SUBROUTINE DCREL(TRANS,M,N,DESCRA,A,IA1,IA2,IP1,DESCRN,
     *                 AN,IAN1,IAN2,IP2,LAN,LIAN1,LIAN2,
     *                 IAUX,LIAUX,IERRV)
      IMPLICIT NONE
C
C     .. Scalar Arguments ..
      INTEGER          LAN, LIAUX, LIAN1, LIAN2, M, N
      CHARACTER        TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION A(*), AN(*)
      INTEGER          IA1(*), IA2(*), IAN1(*), IAN2(*), IP1(*), IP2(*),
     *                 IAUX(LIAUX), IERRV(*)
      CHARACTER        DESCRA*11, DESCRN*11
C     .. Local Scalars ..
      INTEGER          I, J, LWORKR
C     .. External Subroutines ..
      EXTERNAL XSPERR
C     .. Executable Statements ..
C

C
C       Check for argument errors
C
      IF (TRANS.NE.'T' .AND. TRANS.NE.'N') THEN
         CALL XSPERR('TRANS   ',ICHAR(TRANS),1,'DCREL',IERRV)
      ENDIF
      IF (M.LE.0) THEN
         CALL XSPERR('MATDIM  ',M,2,'DCREL',IERRV)
      ENDIF
      IF (N.LE.0) THEN
         CALL XSPERR('MATDIM  ',N,3,'DCREL',IERRV)
      ENDIF
      IF(LIAN2.LT.1) THEN
         LIAN2 = 1
         CALL XSPERR('MATST   ',LIAN2,16,'DCREL',IERRV)
      ENDIF
      IF (TRANS.EQ.'N') THEN
         LWORKR = 0
      ELSE IF (TRANS.EQ.'T') THEN
         LWORKR =  N
      ENDIF
      IF (LIAUX.LT.LWORKR) THEN
         CALL XSPERR('LWORK   ',LIAUX,18,'DCREL',IERRV)
         LIAUX = LWORKR
      ENDIF
      IF (IERRV(1).NE.0)    RETURN


      descrn(1:3) = descra(1:3)
      IP1(1)=0
      IP2(1)=0

      IF(TRANS.EQ.'N') THEN
C
C       Input matrix need not be permuted
C
         IAN2(1)=IA2(2)-IA2(1)
         DO I = 2, M
            IAN2(1) = MAX0(IAN2(1),IA2(I+1)-IA2(I))
         ENDDO

         IF (LAN.LT.(M*IAN2(1))) THEN
            CALL XSPERR('MATST   ',LAN,14,'DCREL',IERRV)
            LAN = M * IAN2(1)
         ENDIF
         IF (LIAN1.LT.(M*IAN2(1))) THEN
            CALL XSPERR('MATST   ',LIAN1,15,'DCREL',IERRV)
            LIAN1 = M * IAN2(1)
         ENDIF
         IF (IERRV(1).NE.0)    RETURN

         DO I=1,M
            DO J=IA2(I),IA2(I+1)-1
               AN(M*(J-IA2(I))+I)=A(J)
               IAN1(M*(J-IA2(I))+I)=IA1(J)
            ENDDO
            DO J=IA2(I+1)-IA2(I)+1,IAN2(1)
               AN(M*(J-1)+I)=0.D0
               IAN1(M*(J-1)+I)=IAN2((J-2)*M+I)
            ENDDO
         ENDDO

      ELSE
C
C       Input matrix has to be permuted
C

         DO J=1,N
            IAUX(I)=0
         ENDDO
         DO I=1,M
            DO J=IA2(I),IA2(I+1)-1
               IAUX(IA1(J))=IAUX(IA1(J))+1
            ENDDO
         ENDDO
         IAN2(1)=IAUX(1)
         DO I = 2, M
            IAN2(1) = MAX0(IAN2(1),IAUX(I))
         ENDDO

         IF (LAN.LT.(N*IAN2(1))) THEN
            CALL XSPERR('MATST   ',LAN,14,'DCREL',IERRV)
            LAN = N * IAN2(1)
         ENDIF
         IF (LIAN1.LT.(N*IAN2(1))) THEN
            CALL XSPERR('MATST   ',LIAN1,15,'DCREL',IERRV)
            LIAN1 = N * IAN2(1)
         ENDIF
         IF (IERRV(1).NE.0)    RETURN

         DO J=1,N
            IAUX(I)=0
         ENDDO
         DO I=1,M
            DO J=IA2(I),IA2(I+1)-1
               IAUX(IA1(J))=IAUX(IA1(J))+1
               AN  (N*(IAUX(IA1(J)))+IA1(J))=A(J)
               IAN1(N*(IAUX(IA1(J)))+IA1(J))=I
            ENDDO
         ENDDO
         DO I=1,N
            DO J=IAUX(I)+1,IAN2(I)
               AN  (N*(J-1)+I)=0.D0
               IAN1(N*(J-1)+I)=IAN1(N*IAUX(I)+I)
            ENDDO
         ENDDO

      ENDIF


      RETURN
      END

