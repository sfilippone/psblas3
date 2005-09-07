MODULE GETP

  PUBLIC GET_PARMS
  PUBLIC PR_USAGE

CONTAINS
  !
  ! Get iteration parameters from the command line
  !
  SUBROUTINE  GET_PARMS(ICONTXT,MTRX_FILE,RHS_FILE,CMETHD,PREC,IPART,&
       & AFMT,ISTOPC,ITMAX,ITRACE,ML,IPREC,EPS)
    integer      :: icontxt
    Character*20 :: CMETHD, PREC, MTRX_FILE, RHS_FILE
    Integer      :: IRET, ISTOPC,ITMAX,ITRACE,IPART,IPREC,ML
    Character*40 :: CHARBUF
    real(kind(1.d0)) :: eps
    character    :: afmt*5
    INTEGER      :: IARGC, NPROW, NPCOL, MYPROW, MYPCOL
    EXTERNAL     IARGC
    INTEGER      :: INPARMS(20), IP 
    
    CALL BLACS_GRIDINFO(ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL)
    IF (MYPROW==0) THEN
      ! Read Input Parameters
      READ(*,*) IP
      IF (IP.GE.3) THEN
        READ(*,*) MTRX_FILE
        READ(*,*) RHS_FILE
        READ(*,*) CMETHD
        READ(*,*) PREC
        READ(*,*) AFMT

        ! Convert strings in array
        DO I = 1, LEN(MTRX_FILE)
          INPARMS(I) = IACHAR(MTRX_FILE(I:I))
        END DO
        ! Broadcast parameters to all processors
        CALL IGEBS2D(ICONTXT,'ALL',' ',20,1,INPARMS,20)

        ! Convert strings in array
        DO I = 1, LEN(CMETHD)
          INPARMS(I) = IACHAR(CMETHD(I:I))
        END DO
        ! Broadcast parameters to all processors
        CALL IGEBS2D(ICONTXT,'ALL',' ',20,1,INPARMS,20)

        DO I = 1, LEN(PREC)
          INPARMS(I) = IACHAR(PREC(I:I))
        END DO
        ! Broadcast parameters to all processors
        CALL IGEBS2D(ICONTXT,'ALL',' ',20,1,INPARMS,20)

        DO I = 1, LEN(AFMT)
          INPARMS(I) = IACHAR(AFMT(I:I))
        END DO
        ! Broadcast parameters to all processors
        CALL IGEBS2D(ICONTXT,'ALL',' ',20,1,INPARMS,20)

        READ(*,*) IPART
        IF (IP.GE.6) THEN
          READ(*,*) ISTOPC
        ELSE
          ISTOPC=1        
        ENDIF
        IF (IP.GE.7) THEN
          READ(*,*) ITMAX
        ELSE
          ITMAX=500
        ENDIF
        IF (IP.GE.8) THEN
          READ(*,*) ITRACE
        ELSE
          ITRACE=-1
        ENDIF
        IF (IP.GE.9) THEN
          READ(*,*) IPREC
        ELSE
          IPREC=0
        ENDIF
        IF (IP.GE.10) THEN
          READ(*,*) ML
        ELSE
          ML  = 1
        ENDIF
        IF (IP.GE.11) THEN
          READ(*,*) EPS
        ELSE
          EPS=1.D-6
        ENDIF
        ! Broadcast parameters to all processors    

        INPARMS(1) = IPART
        INPARMS(2) = ISTOPC
        INPARMS(3) = ITMAX
        INPARMS(4) = ITRACE
        INPARMS(5) = IPREC
        INPARMS(6) = ML
        CALL IGEBS2D(ICONTXT,'ALL',' ',6,1,INPARMS,6)
        CALL DGEBS2D(ICONTXT,'ALL',' ',1,1,EPS,1)

        WRITE(6,*)'Solving matrix:  ',mtrx_file
        WRITE(6,*)' with BLOCK data distribution, NP=',NPROW,&
             & ' Preconditioner=',PREC
      else
        CALL PR_USAGE(0)
        CALL BLACS_ABORT(ICONTXT,-1)
        STOP 1
      END IF
    ELSE
      ! Receive Parameters
      CALL IGEBR2D(ICONTXT,'A',' ',20,1,INPARMS,20,0,0)
      DO I = 1, 20
        MTRX_FILE(I:I) = ACHAR(INPARMS(I))
      END DO
      
      CALL IGEBR2D(ICONTXT,'A',' ',20,1,INPARMS,20,0,0)
      DO I = 1, 20
        CMETHD(I:I) = ACHAR(INPARMS(I))
      END DO
      
      CALL IGEBR2D(ICONTXT,'A',' ',20,1,INPARMS,20,0,0)
      DO I = 1, 20
        PREC(I:I) = ACHAR(INPARMS(I))
      END DO
      CALL IGEBR2D(ICONTXT,'A',' ',20,1,INPARMS,20,0,0)
      DO I = 1, 20
        AFMT(I:I) = ACHAR(INPARMS(I))
      END DO
      
      CALL IGEBR2D(ICONTXT,'A',' ',6,1,INPARMS,6,0,0)

      IPART  =  INPARMS(1) 
      ISTOPC =  INPARMS(2) 
      ITMAX  =  INPARMS(3) 
      ITRACE =  INPARMS(4) 
      IPREC  =  INPARMS(5) 
      ML     =  INPARMS(6) 
      CALL DGEBR2D(ICONTXT,'A',' ',1,1,EPS,1,0,0)     
    END IF
    
  END SUBROUTINE GET_PARMS
  SUBROUTINE PR_USAGE(IOUT)
    INTEGER IOUT
    WRITE(IOUT, *) ' Number of parameters is incorrect!'
    WRITE(IOUT, *) ' Use: hb_sample mtrx_file methd prec [ptype &
         &itmax istopc itrace]' 
    WRITE(IOUT, *) ' Where:'
    WRITE(IOUT, *) '     mtrx_file      is stored in HB format'
    WRITE(IOUT, *) '     methd          may be: CGSTAB '
    WRITE(IOUT, *) '     prec           may be: ILU DIAGSC NONE'
    WRITE(IOUT, *) '     ptype          Partition strategy default 0'
    WRITE(IOUT, *) '                    0: BLOCK partition '
    WRITE(IOUT, *) '     itmax          Max iterations [500]        '
    WRITE(IOUT, *) '     istopc         Stopping criterion [1]      '
    WRITE(IOUT, *) '     itrace         0  (no tracing, default) or '
    WRITE(IOUT, *) '                    >= 0 do tracing every ITRACE'
    WRITE(IOUT, *) '                    iterations ' 
  END SUBROUTINE PR_USAGE
END MODULE GETP
