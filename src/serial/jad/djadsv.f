C     This routine compute only the casew with Unit diagonal ...
      SUBROUTINE DJADSV(UNITD,NROW,NG,A,KA,IA,JA,X,Y,IERROR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IA(3,*),JA(*),KA(*),A(*),X(*),Y(*)
      CHARACTER UNITD
      INTEGER   IERROR

      IERROR=0

      IF (UNITD.EQ.'U') THEN
         
         IF (NG.EQ.0) THEN
            DO I = 1, NROW
               Y(I) = X(I)
            ENDDO
         ENDIF

         DO IPG=1,NG
            DO I = IA(1,IPG),IA(1,IPG+1)-1
               Y(I) = X(I)
            END DO
*
*        LOOP ON COLUMNS
*        ---------------
*
            IP2 = IA(2,IPG)
            DO K = IP2, IA(3,IPG)-1
               IPX = IA(1,IPG)
               DO I = JA(K), JA(K+1)-1
                  Y(IPX) = Y(IPX) - A(I)*Y(KA(I))
                  IPX = IPX+1
               ENDDO
            ENDDO
*
*
*       LOOP ON ROWS
*       ---------------
*
            IPX = IA(1,IPG)
            DO K = IA(3,IPG), IA(2,IPG+1)-1
               DO I = JA(K), JA(K+1)-1
                  Y(IPX) = Y(IPX) - A(I)*Y(KA(I))
               ENDDO
               IPX = IPX + 1
            ENDDO

****************************************
         END DO                 !END LOOP ON IPG=1,NG
****************************************
      ELSE
         WRITE(0,*) 'ERROR in DJADSV'
      ENDIF
      RETURN
      END

