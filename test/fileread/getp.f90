!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
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
Module getp
  interface get_parms
    module procedure get_dparms, get_sparms
  end interface

contains
  !
  ! Get iteration parameters from the command line
  !
  subroutine  get_dparms(ictxt,mtrx_file,rhs_file,filefmt,kmethd,ptype,ipart,&
       & afmt,istopc,itmax,itrace,irst,eps)
    use psb_sparse_mod
    integer      :: ictxt
    character(len=2)  :: filefmt
    character(len=40) :: kmethd, mtrx_file, rhs_file, ptype
    integer      :: iret, istopc,itmax,itrace,ipart,irst
    character(len=40) :: charbuf
    real(psb_dpk_) :: eps
    character    :: afmt*5
    integer      :: np, iam
    integer      :: inparms(40), ip 

    call psb_info(ictxt,iam,np)
    if (iam==0) then
      ! Read Input Parameters
      read(*,*) ip
      if (ip >= 5) then
        read(*,*) mtrx_file
        read(*,*) rhs_file
        read(*,*) filefmt
        read(*,*) kmethd
        read(*,*) ptype
        read(*,*) afmt


        call psb_bcast(ictxt,mtrx_file)
        call psb_bcast(ictxt,rhs_file)
        call psb_bcast(ictxt,filefmt)
        call psb_bcast(ictxt,kmethd)
        call psb_bcast(ictxt,ptype)
        call psb_bcast(ictxt,afmt)

        read(*,*) ipart
        if (ip >= 7) then
          read(*,*) istopc
        else
          istopc=1        
        endif
        if (ip >= 8) then
          read(*,*) itmax
        else
          itmax=500
        endif
        if (ip >= 9) then
          read(*,*) itrace
        else
          itrace=-1
        endif
        if (ip >= 10) then
          read(*,*) irst
        else
          irst  = 1
        endif
        if (ip >= 11) then
          read(*,*) eps
        else
          eps=1.d-6
        endif
        inparms(1) = ipart
        inparms(2) = istopc
        inparms(3) = itmax
        inparms(4) = itrace
        inparms(5) = irst
        call psb_bcast(ictxt,inparms(1:5))
        call psb_bcast(ictxt,eps)

        write(*,'("Solving matrix       : ",a)')  mtrx_file      
        write(*,'("Number of processors : ",i3)') np
        write(*,'("Data distribution    : ",i2)') ipart
        write(*,'("Iterative method     : ",a)')  kmethd
        write(*,'("Preconditioner       : ",a)')  ptype
        write(*,'("Restart parameter    : ",i2)') irst
        write(*,'("Storage format       : ",a)')  afmt(1:3)
        write(*,'(" ")')
      else
        write(0,*) 'Wrong format for input file'
        call psb_abort(ictxt)
        stop 1
      end if
    else
      ! Receive Parameters
      call psb_bcast(ictxt,mtrx_file)
      call psb_bcast(ictxt,rhs_file)
      call psb_bcast(ictxt,filefmt)
      call psb_bcast(ictxt,kmethd)
      call psb_bcast(ictxt,ptype)
      call psb_bcast(ictxt,afmt)

      call psb_bcast(ictxt,inparms(1:5))
      ipart  =  inparms(1) 
      istopc =  inparms(2) 
      itmax  =  inparms(3) 
      itrace =  inparms(4) 
      irst   =  inparms(5) 
      call psb_bcast(ictxt,eps)

    end if

  end subroutine get_dparms
  
  subroutine  get_sparms(ictxt,mtrx_file,rhs_file,filefmt,kmethd,ptype,ipart,&
       & afmt,istopc,itmax,itrace,irst,eps)
    use psb_sparse_mod
    integer      :: ictxt
    character(len=2)  :: filefmt
    character(len=40) :: kmethd, mtrx_file, rhs_file, ptype
    integer      :: iret, istopc,itmax,itrace,ipart,irst
    character(len=40) :: charbuf
    real(psb_spk_) :: eps
    character    :: afmt*5
    integer      :: np, iam
    integer      :: inparms(40), ip 

    call psb_info(ictxt,iam,np)
    if (iam==0) then
      ! Read Input Parameters
      read(*,*) ip
      if (ip >= 5) then
        read(*,*) mtrx_file
        read(*,*) rhs_file
        read(*,*) filefmt
        read(*,*) kmethd
        read(*,*) ptype
        read(*,*) afmt


        call psb_bcast(ictxt,mtrx_file)
        call psb_bcast(ictxt,rhs_file)
        call psb_bcast(ictxt,filefmt)
        call psb_bcast(ictxt,kmethd)
        call psb_bcast(ictxt,ptype)
        call psb_bcast(ictxt,afmt)

        read(*,*) ipart
        if (ip >= 7) then
          read(*,*) istopc
        else
          istopc=1        
        endif
        if (ip >= 8) then
          read(*,*) itmax
        else
          itmax=500
        endif
        if (ip >= 9) then
          read(*,*) itrace
        else
          itrace=-1
        endif
        if (ip >= 10) then
          read(*,*) irst
        else
          irst  = 1
        endif
        if (ip >= 11) then
          read(*,*) eps
        else
          eps=1.d-6
        endif
        inparms(1) = ipart
        inparms(2) = istopc
        inparms(3) = itmax
        inparms(4) = itrace
        inparms(5) = irst
        call psb_bcast(ictxt,inparms(1:5))
        call psb_bcast(ictxt,eps)

        write(*,'("Solving matrix       : ",a)')  mtrx_file      
        write(*,'("Number of processors : ",i3)') np
        write(*,'("Data distribution    : ",i2)') ipart
        write(*,'("Iterative method     : ",a)')  kmethd
        write(*,'("Preconditioner       : ",a)')  ptype
        write(*,'("Restart parameter    : ",i2)') irst
        write(*,'("Storage format       : ",a)')  afmt(1:3)
        write(*,'(" ")')
      else
        write(0,*) 'Wrong format for input file'
        call psb_abort(ictxt)
        stop 1
      end if
    else
      ! Receive Parameters
      call psb_bcast(ictxt,mtrx_file)
      call psb_bcast(ictxt,rhs_file)
      call psb_bcast(ictxt,filefmt)
      call psb_bcast(ictxt,kmethd)
      call psb_bcast(ictxt,ptype)
      call psb_bcast(ictxt,afmt)

      call psb_bcast(ictxt,inparms(1:5))
      ipart  =  inparms(1) 
      istopc =  inparms(2) 
      itmax  =  inparms(3) 
      itrace =  inparms(4) 
      irst   =  inparms(5) 
      call psb_bcast(ictxt,eps)

    end if

  end subroutine get_sparms

end module getp
