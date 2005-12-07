module psb_realloc_mod
  implicit none

  Interface psb_realloc
    module procedure psb_dreallocate1i
    module procedure psb_dreallocate2i
    module procedure psb_dreallocate2i1d
    module procedure psb_dreallocate1d
    module procedure psb_dreallocated2
  end Interface

  Interface psb_realloc1it
    module procedure  psb_dreallocate1it
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

    name='psb_dreallocate1i'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info=0
    if (associated(rrax)) then 
      dim=size(rrax)
      If (dim /= len) Then
        Allocate(tmp(len),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
!!$        write(0,*) 'IA: copying ',len,dim
        if (.true.) then 
          do i=1, min(len,dim)
            tmp(i)=rrax(i)
          end do
        else
          tmp(1:min(len,dim))=rrax(1:min(len,dim))
        end if
!!$        write(0,*) 'IA: copying done'
        Deallocate(rrax,stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        rrax=>tmp
      End If
    else
!!$      write(0,*) 'IA: allocating ',len
      allocate(rrax(len),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
!!$      write(0,*) 'IA: padding'
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

    name='psb_dreallocate1d'
    call psb_erractionsave(err_act)

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



  Subroutine psb_dreallocated2(len1,len2,rrax,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len1,len2
    Real(kind(1.d0)),pointer :: rrax(:,:)
    integer :: info
    real(kind(1.d0)), optional, intent(in) :: pad

    ! ...Local Variables
    Real(kind(1.d0)),Pointer :: tmp(:,:)
    Integer :: dim,err_act,err,i, m
    character(len=20)  :: name

    name='psb_dreallocated2'
    call psb_erractionsave(err_act)

    if (associated(rrax)) then 
      dim=size(rrax,1)

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
            tmp(i,:) = rrax(i,:)
          end do
        else
          tmp(1:m,:) = rrax(1:m,:)
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
      Allocate(rrax(len1,len2),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(dim+1:len1,:) = pad
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


  Subroutine psb_dreallocate1it(len,rrax,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  
    Integer,Intent(in) :: len
    Integer,pointer :: rrax(:)
    integer         :: info
    integer, optional, intent(in) :: pad
    ! ...Local Variables
    Integer,Pointer :: tmp(:)
    Integer :: dim,err_act,err
    character(len=20)  :: name

    name='psb_dreallocate1it'
    call psb_erractionsave(err_act)

    if(psb_get_errstatus().ne.0) return 
    info=0
    if (associated(rrax)) then 
      dim=size(rrax)
      If (dim /= len) Then
        Allocate(tmp(len),stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
!!$        write(0,*) 'IA: copying ',min(len,dim)
        tmp(1:min(len,dim))=rrax(1:min(len,dim))
!!$        write(0,*) 'IA: copying done'
        Deallocate(rrax,stat=info)
        if (info /= 0) then
          err=4000
          call psb_errpush(err,name)
          goto 9999
        end if
        rrax=>tmp
      End If
    else
      allocate(rrax(len),stat=info)
      if (info /= 0) then
        err=4000
        call psb_errpush(err,name)
        goto 9999
      end if
    endif
    if (present(pad)) then 
!!$      write(0,*) 'IA: padding'
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

  End Subroutine psb_dreallocate1it


end module psb_realloc_mod
