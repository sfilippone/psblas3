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
module psb_realloc_mod
  use psb_const_mod
  implicit none
  
  Interface psb_realloc
    module procedure psb_dreallocate1i
    module procedure psb_dreallocate2i
    module procedure psb_dreallocate2i1d
    module procedure psb_dreallocate1d
    module procedure psb_dreallocated2
    module procedure psb_dreallocatei2
    module procedure psb_dreallocate2i1z
    module procedure psb_dreallocate1z
    module procedure psb_dreallocatez2
  end Interface

Contains

  Subroutine psb_dreallocate1i(len,rrax,info,pad)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    Integer,pointer :: rrax(:)
    integer         :: info
    integer, optional, intent(in) :: pad
    ! ...Local Variables
    Integer,Pointer :: tmp(:)
    Integer :: dim, err_act, err,i
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_dreallocate1i'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info=0
    if (debug) write(0,*) 'reallocate I',len
    if (associated(rrax)) then 
      dim=size(rrax)
      If (dim /= len) Then
        Allocate(tmp(len),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        if (.true.) then 
          do i=1, min(len,dim)
            tmp(i)=rrax(i)
          end do
        else
          tmp(1:min(len,dim))=rrax(1:min(len,dim))
        end if
        deallocate(rrax,stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        rrax=>tmp
      end if
    else
      dim = 0
      allocate(rrax(len),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(dim+1:len) = pad
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return


  End Subroutine psb_dreallocate1i





  Subroutine psb_dreallocate1d(len,rrax,info,pad)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    Real(kind(1.d0)),pointer :: rrax(:)
    integer :: info
    real(kind(1.d0)), optional, intent(in) :: pad

    ! ...Local Variables
    Real(kind(1.d0)),Pointer :: tmp(:)
    Integer :: dim,err_act,err,i, m
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_dreallocate1d'
    call psb_erractionsave(err_act)
    info = 0 
    if (debug) write(0,*) 'reallocate D',len

    if (associated(rrax)) then 
      dim=size(rrax)

      If (dim /= len) Then
        Allocate(tmp(len),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        m = min(dim,len)
!!$        write(0,*) 'DA: copying ',min(len,dim)
        if (.true.) then 
          do i=1,m
            tmp(i) = rrax(i)
          end do
        else
          tmp(1:m) = rrax(1:m)
        end if
!!$        write(0,*) 'DA: copying done ',m
        Deallocate(rrax,stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        rrax=>tmp
      End If
    else
      dim = 0
      Allocate(rrax(len),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(dim+1:len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocate1d


  Subroutine psb_dreallocate1z(len,rrax,info,pad)
    use psb_error_mod

    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    complex(kind(1.d0)),pointer :: rrax(:)
    integer :: info
    complex(kind(1.d0)), optional, intent(in) :: pad

    ! ...Local Variables
    complex(kind(1.d0)),Pointer :: tmp(:)
    Integer :: dim,err_act,err,i, m
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_dreallocate1z'
    call psb_erractionsave(err_act)
    info = 0
    if (debug) write(0,*) 'reallocate Z',len    

    if (associated(rrax)) then 
      dim=size(rrax)

      If (dim /= len) Then
        Allocate(tmp(len),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        m = min(dim,len)
!!$        write(0,*) 'DA: copying ',min(len,dim)
        if (.true.) then 
          do i=1,m
            tmp(i) = rrax(i)
          end do
        else
          tmp(1:m) = rrax(1:m)
        end if
!!$        write(0,*) 'DA: copying done ',m
        Deallocate(rrax,stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        rrax=>tmp
      End If
    else
      dim = 0
      Allocate(rrax(len),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(dim+1:len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocate1z



  Subroutine psb_dreallocated2(len1,len2,rrax,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len1,len2
    Real(kind(1.d0)),pointer :: rrax(:,:)
    integer :: info
    real(kind(1.d0)), optional, intent(in) :: pad

    ! ...Local Variables
    Real(kind(1.d0)),Pointer :: tmp(:,:)
    Integer :: dim,err_act,err,i, m, dim2
    character(len=20)  :: name

    name='psb_dreallocated2'
    call psb_erractionsave(err_act)
    info = 0 

    if (associated(rrax)) then 
      dim=size(rrax,1)
      dim2=size(rrax,2)

      If (dim /= len1) Then
        Allocate(tmp(len1,len2),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        m = min(dim,len1)
!!$        write(0,*) 'DA: copying ',min(len,dim)
        if (.true.) then 
          do i=1,m
            tmp(i,1:min(len2,dim2)) = rrax(i,1:min(len2,dim2))
          end do
        else
          tmp(1:m,1:min(len2,dim2)) = rrax(1:m,1:min(len2,dim2))
        end if
!!$        write(0,*) 'DA: copying done ',m
        Deallocate(rrax,stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        rrax=>tmp
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(len1,len2),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(dim+1:len1,:) = pad
      rrax(:,dim2+1:len2) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocated2

  Subroutine psb_dreallocatez2(len1,len2,rrax,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len1,len2
    complex(kind(1.d0)),pointer :: rrax(:,:)
    integer :: info
    complex(kind(1.d0)), optional, intent(in) :: pad

    ! ...Local Variables
    complex(kind(1.d0)),Pointer :: tmp(:,:)
    Integer :: dim,err_act,err,i, m, dim2
    character(len=20)  :: name

    name='psb_dreallocatez2'
    call psb_erractionsave(err_act)
    info = 0 

    if (associated(rrax)) then 
      dim=size(rrax,1)
      dim2=size(rrax,2)

      If (dim /= len1) Then
        Allocate(tmp(len1,len2),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        m = min(dim,len1)
!!$        write(0,*) 'DA: copying ',min(len,dim)
        if (.true.) then 
          do i=1,m
            tmp(i,1:min(len2,dim2)) = rrax(i,1:min(len2,dim2))
          end do
        else
          tmp(1:m,1:min(len2,dim2)) = rrax(1:m,1:min(len2,dim2))
        end if
!!$        write(0,*) 'DA: copying done ',m
        Deallocate(rrax,stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        rrax=>tmp
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(len1,len2),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(dim+1:len1,:) = pad
      rrax(:,dim2+1:len2) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocatez2

  Subroutine psb_dreallocatei2(len1,len2,rrax,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len1,len2
    integer,pointer :: rrax(:,:)
    integer :: info
    integer, optional, intent(in) :: pad

    ! ...Local Variables
    integer,Pointer :: tmp(:,:)
    Integer :: dim,err_act,err,i, m, dim2
    character(len=20)  :: name

    name='psb_dreallocatei2'
    call psb_erractionsave(err_act)
    info = 0 

    if (associated(rrax)) then 
      dim=size(rrax,1)
      dim2=size(rrax,2)

      If (dim /= len1) Then
        Allocate(tmp(len1,len2),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        m = min(dim,len1)
!!$        write(0,*) 'DA: copying ',min(len,dim)
        if (.true.) then 
          do i=1,m
            tmp(i,1:min(len2,dim2)) = rrax(i,1:min(len2,dim2))
          end do
        else
          tmp(1:m,1:min(len2,dim2)) = rrax(1:m,1:min(len2,dim2))
        end if
!!$        write(0,*) 'DA: copying done ',m
        Deallocate(rrax,stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        rrax=>tmp
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(len1,len2),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(dim+1:len1,:) = pad
      rrax(:,dim2+1:len2) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocatei2


  Subroutine psb_dreallocate2i(len,rrax,y,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  

    Integer,Intent(in) :: len
    Integer,pointer :: rrax(:),y(:)
    integer :: info
    integer, optional, intent(in) :: pad
    character(len=20)  :: name
    integer :: err_act, err

    name='psb_dreallocate2i'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info=0
    call psb_dreallocate1i(len,rrax,info,pad=pad)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1i(len,y,info,pad=pad)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return

  End Subroutine psb_dreallocate2i




  Subroutine psb_dreallocate2i1d(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    Integer,pointer :: rrax(:),y(:)
    Real(Kind(1.d0)),pointer :: z(:)
    integer :: info
    character(len=20)  :: name
    integer :: err_act, err

    name='psb_dreallocate2i1d'
    call psb_erractionsave(err_act)


    info = 0
    call psb_dreallocate1i(len,rrax,info)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1i(len,y,info)    
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1d(len,z,info)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return
  End Subroutine psb_dreallocate2i1d



  Subroutine psb_dreallocate2i1z(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    Integer,pointer :: rrax(:),y(:)
    complex(Kind(1.d0)),pointer :: z(:)
    integer :: info
    character(len=20)  :: name
    integer :: err_act, err

    name='psb_dreallocate2i1d'
    call psb_erractionsave(err_act)


    info = 0
    call psb_dreallocate1i(len,rrax,info)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1i(len,y,info)    
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_dreallocate1z(len,z,info)
    if (info /= 0) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act.eq.act_ret) then
      return
    else
      call psb_error()
    end if
    return
  End Subroutine psb_dreallocate2i1z


end module psb_realloc_mod
