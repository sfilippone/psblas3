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
MODULE GETP

  PUBLIC GET_PARMS
  PUBLIC PR_USAGE

CONTAINS
  !
  ! Get iteration parameters from the command line
  !
  SUBROUTINE  GET_PARMS(ICONTXT,MTRX_FILE,RHS_FILE,CMETHD,IPART,&
       & AFMT,ISTOPC,ITMAX,ITRACE,NOVR,IPREC,EPS)
    integer      :: icontxt
    Character*40 :: CMETHD, MTRX_FILE, RHS_FILE
    Integer      :: IRET, ISTOPC,ITMAX,ITRACE,IPART,IPREC,NOVR
    Character*40 :: CHARBUF
    real(kind(1.d0)) :: eps
    character    :: afmt*5
    INTEGER      :: IARGC, NPROW, NPCOL, MYPROW, MYPCOL
    EXTERNAL     IARGC
    INTEGER      :: INPARMS(40), IP 
    
    CALL BLACS_GRIDINFO(ICONTXT, NPROW, NPCOL, MYPROW, MYPCOL)
    IF (MYPROW==0) THEN
      ! Read Input Parameters
      READ(*,*) IP
      IF (IP.GE.3) THEN
        READ(*,*) MTRX_FILE
        READ(*,*) RHS_FILE
        READ(*,*) CMETHD
        READ(*,*) AFMT

        ! Convert strings in array
        DO I = 1, LEN(MTRX_FILE)
          INPARMS(I) = IACHAR(MTRX_FILE(I:I))
        END DO
        ! Broadcast parameters to all processors
        CALL IGEBS2D(ICONTXT,'ALL',' ',40,1,INPARMS,40)

        ! Convert strings in array
        DO I = 1, LEN(CMETHD)
          INPARMS(I) = IACHAR(CMETHD(I:I))
        END DO
        ! Broadcast parameters to all processors
        CALL IGEBS2D(ICONTXT,'ALL',' ',40,1,INPARMS,40)

        DO I = 1, LEN(AFMT)
          INPARMS(I) = IACHAR(AFMT(I:I))
        END DO
        ! Broadcast parameters to all processors
        CALL IGEBS2D(ICONTXT,'ALL',' ',40,1,INPARMS,40)

        READ(*,*) IPART
        IF (IP.GE.5) THEN
          READ(*,*) ISTOPC
        ELSE
          ISTOPC=1        
        ENDIF
        IF (IP.GE.6) THEN
          READ(*,*) ITMAX
        ELSE
          ITMAX=500
        ENDIF
        IF (IP.GE.7) THEN
          READ(*,*) ITRACE
        ELSE
          ITRACE=-1
        ENDIF
        IF (IP.GE.8) THEN
          READ(*,*) IPREC
        ELSE
          IPREC=0
        ENDIF
        IF (IP.GE.9) THEN
          READ(*,*) NOVR
        ELSE
          NOVR  = 1
        ENDIF
        IF (IP.GE.10) THEN
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
        INPARMS(6) = NOVR
        CALL IGEBS2D(ICONTXT,'ALL',' ',6,1,INPARMS,6)
        CALL DGEBS2D(ICONTXT,'ALL',' ',1,1,EPS,1)

          write(*,'("Solving matrix       : ",a40)')mtrx_file      
          write(*,'("Number of processors : ",i3)')nprow
          write(*,'("Data distribution    : ",i2)')ipart
          write(*,'("Preconditioner       : ",i2)')iprec
          if(iprec.gt.2) write(*,'("Overlapping levels   : ",i2)')novr
          write(*,'("Iterative method     : ",a40)')cmethd
          write(*,'("Storage format       : ",a3)')afmt(1:3)
          write(*,'(" ")')
      else
        CALL PR_USAGE(0)
        CALL BLACS_ABORT(ICONTXT,-1)
        STOP 1
      END IF
    ELSE
      ! Receive Parameters
      CALL IGEBR2D(ICONTXT,'A',' ',40,1,INPARMS,40,0,0)
      DO I = 1, 40
        MTRX_FILE(I:I) = ACHAR(INPARMS(I))
      END DO
      
      CALL IGEBR2D(ICONTXT,'A',' ',40,1,INPARMS,40,0,0)
      DO I = 1, 40
        CMETHD(I:I) = ACHAR(INPARMS(I))
      END DO

      CALL IGEBR2D(ICONTXT,'A',' ',40,1,INPARMS,40,0,0)
      DO I = 1, LEN(AFMT)
        AFMT(I:I) = ACHAR(INPARMS(I))
      END DO
      
      CALL IGEBR2D(ICONTXT,'A',' ',6,1,INPARMS,6,0,0)

      IPART  =  INPARMS(1) 
      ISTOPC =  INPARMS(2) 
      ITMAX  =  INPARMS(3) 
      ITRACE =  INPARMS(4) 
      IPREC  =  INPARMS(5) 
      NOVR     =  INPARMS(6) 
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
    WRITE(IOUT, *) '     ptype          Partition strategy default 0'
    WRITE(IOUT, *) '                    0: BLOCK partition '
    WRITE(IOUT, *) '     itmax          Max iterations [500]        '
    WRITE(IOUT, *) '     istopc         Stopping criterion [1]      '
    WRITE(IOUT, *) '     itrace         0  (no tracing, default) or '
    WRITE(IOUT, *) '                    >= 0 do tracing every ITRACE'
    WRITE(IOUT, *) '                    iterations ' 
  END SUBROUTINE PR_USAGE
END MODULE GETP
