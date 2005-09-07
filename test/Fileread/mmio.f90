module mmio
  use typesp
  public mm_mat_read, mm_mat_write
  interface mm_mat_read
    module procedure dmm_mat_read, zmm_mat_read
  end interface
  interface mm_mat_write
    module procedure dmm_mat_write
  end interface
  private desym,zdesym

contains

  subroutine dmm_mat_read(a, iret, iunit, filename)   
    use typesp
    implicit none
    type(d_spmat), intent(out)  :: a
    integer, intent(out)        :: iret
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character      :: mmheader*15, fmt*15, object*10, type*10, sym*15
    character(1024)      :: line
    integer        :: indcrd,  ptrcrd, totcrd,&
         & valcrd, rhscrd, nrow, ncol, nnzero, neltvl, nrhs, nrhsix
    real(kind(1.0d0)), pointer  :: as_loc(:), dwork(:)
    integer, pointer            :: ia1_loc(:), ia2_loc(:), iwork(:), tmp(:), aux(:)
    integer                     :: ircode, i,iel,ptr,nzr,infile,&
         & j, liwork, ldwork, root, nprow, npcol, myprow, mypcol
    logical, parameter :: debug=.false.

    iret = 0

    if (present(filename)) then
      if (filename=='-') then 
        infile=5
      else
        if (present(iunit)) then 
          infile=iunit
        else
          infile=99
        endif
        open(infile,file=filename, status='OLD', err=901, action='READ')
      endif
    else 
      if (present(iunit)) then 
        infile=iunit
      else
        infile=5
      endif
    endif

    read(infile,fmt=*,end=902) mmheader, object, fmt, type, sym
    call lowerc(object,1,10)
    call lowerc(fmt,1,15)

    if ( (object .ne. 'matrix').or.(fmt.ne.'coordinate')) then
      write(0,*) 'READ_MATRIX: input file type not yet supported'
      iret=909
      return
    end if
    if (debug) write(*,*) mmheader,':', object, ':',fmt,':', type,':', sym

    do 
      read(infile,fmt='(a)') line
      if (line(1:1) /= '%')  exit
    end do
    if (debug) write(*,*) 'Line on input : "',line,'"'
    read(line,fmt=*) nrow,ncol,nnzero
    if (debug) write(*,*) 'Out: ',nrow,ncol,nnzero
    a%m    = nrow
    a%k    = ncol
    a%fida = 'CSR'
    a%descra='G'
    call lowerc(type,1,10)
    call lowerc(sym,1,15)

    if ((type == 'real').and.(sym == 'general')) then
      allocate(a%aspk(nnzero), a%ia1(nnzero), a%ia2(nrow+1),&
           & a%pl(nrow),a%pr(nrow), tmp(nnzero+1), aux(nnzero+2),stat = ircode)
      if (ircode /= 0)   goto 993
      do i=1,nnzero
        read(infile,fmt=*,end=902) tmp(i),a%ia1(i),a%aspk(i)
      end do

      call mrgsrt(nnzero,tmp,aux,ircode)
      if (ircode.eq.0) call reordvn(nnzero,a%aspk,tmp,a%ia1,aux)
      !     .... Order with key a%ia1 (COLUMN INDEX) ...
      i    = 1
      j    = i
      !     .... order with key tmp (row index) ...
      do 
        if (i > nnzero) exit
        do 
          if (j > nnzero) exit
          if (tmp(j) /= tmp(i)) exit
          j = j+1
          !              if (j.eq.(nnzero+1)) exit
        enddo
        iel = j - i
        call mrgsrt(iel,a%ia1(i),aux,ircode)
        if (ircode == 0) call reordvn(iel,a%aspk(i),tmp(i),&
             &   a%ia1(i), aux)
        i = j
      enddo

      ! convert to csr format     
      iel = 1
      a%ia2(1) = 1
      do i = a%ia2(1), nrow

        do 
          if (tmp(iel) /= i) exit
          iel = iel + 1
          if (iel > nnzero) exit
        enddo
        a%ia2(i+1) = iel
      enddo
      deallocate(aux,tmp)

    else if ((type == 'real').and.(sym == 'symmetric')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read

      allocate(a%aspk(2*nnzero),a%ia1(2*nnzero),&
           & a%ia2(2*nnzero),as_loc(2*nnzero),&
           & ia1_loc(2*nnzero),ia2_loc(2*nnzero),&
           & a%pl(nrow),a%pr(nrow), stat = ircode)

      if (ircode /= 0)   goto 993

      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),a%aspk(i)
      end do

      liwork = 2*nnzero+2
      allocate(iwork(liwork), stat = ircode)
      if (ircode /= 0)   goto 993  
      ! After this call NNZERO contains the actual value for
      ! desymetrized matrix
      call desym(nrow, a%aspk, a%ia2, a%ia1, as_loc, ia2_loc,&
           & ia1_loc, iwork, nnzero, nzr)     
      
      call spreall(a,nzr,ircode)
      if (ircode /= 0)   goto 993
      allocate(tmp(nzr),stat=ircode)
      if (ircode /= 0)   goto 993
      if (.false.) then 
        a%aspk(1:nzr) = as_loc(1:nzr)
        a%ia1(1:nzr)  = ia2_loc(1:nzr)
        tmp(1:nzr)    = ia1_loc(1:nzr)
      else
        write(0,*) 'After DESYM: ',nzr,ia2_loc(1:10)
        do i=1,nzr
          a%aspk(i) = as_loc(i)
          a%ia1(i)  = ia2_loc(i)
          tmp(i)    = ia1_loc(i)
        end do
      endif

      iel = 1
      a%ia2(1) = 1
      do i = 1, nrow
        do 
          if (tmp(iel) /= i) exit
          iel = iel + 1
          if (iel > nzr) exit
        enddo
        a%ia2(i+1) = iel
      enddo

      deallocate(as_loc, ia1_loc, ia2_loc,tmp,iwork)
    else
      write(0,*) 'read_matrix: matrix type not yet supported'
      iret=904
    end if
    if (infile/=5) close(infile)
    return 

    ! open failed
901 iret=901
    write(0,*) 'read_matrix: could not open file ',filename,' for input'
    return
902 iret=902
    write(0,*) 'READ_MATRIX: Unexpected end of file '
    return
993 iret=993
    write(0,*) 'READ_MATRIX: Memory allocation failure'
    return
  end subroutine dmm_mat_read


  subroutine zmm_mat_read(a, iret, iunit, filename)   
    use typesp
    implicit none
    type(z_spmat), intent(out)  :: a
    integer, intent(out)        :: iret
    integer, optional, intent(in)          :: iunit
    character(len=*), optional, intent(in) :: filename
    character      :: mmheader*15, fmt*15, object*10, type*10, sym*15, line*1024
    integer        :: indcrd,  ptrcrd, totcrd,&
         & valcrd, rhscrd, nrow, ncol, nnzero, neltvl, nrhs, nrhsix
    complex(kind(1.0d0)), pointer  :: as_loc(:), dwork(:)
    integer, pointer            :: ia1_loc(:), ia2_loc(:), iwork(:), tmp(:), aux(:)
    integer                     :: ircode, i,iel,ptr,nzr,infile,&
         & j, liwork, ldwork, root, nprow, npcol, myprow, mypcol

    iret = 0

    if (present(filename)) then
      if (filename=='-') then 
        infile=5
      else
        if (present(iunit)) then 
          infile=iunit
        else
          infile=99
        endif
        open(infile,file=filename, status='OLD', err=901, action='READ')
      endif
    else 
      if (present(iunit)) then 
        infile=iunit
      else
        infile=5
      endif
    endif

    read(infile,fmt=*,end=902) mmheader, object, fmt, type, sym
    call lowerc(object,1,10)
    call lowerc(fmt,1,15)

    if ( (object .ne. 'matrix').or.(fmt.ne.'coordinate')) then
      write(0,*) 'READ_MATRIX: input file type not yet supported'
      iret=909
      return
    end if

    do 
      read(infile,fmt='(a)') line
      if (line(1:1) /= '%')  exit
    end do
    read(line,fmt=*) nrow,ncol,nnzero

    a%m    = nrow
    a%k    = ncol
    a%fida = 'CSR'
    call lowerc(type,1,10)
    call lowerc(sym,1,15)

    if ((type == 'complex').and.(sym == 'general')) then

      
      allocate(a%aspk(nnzero), a%ia1(nnzero), a%ia2(nrow+1),&
           & a%pl(nrow),a%pr(nrow), tmp(nnzero+1), aux(nnzero+2),stat = ircode)
      if (ircode /= 0)   goto 993
      do i=1,nnzero
        read(infile,fmt=*,end=902) tmp(i),a%ia1(i),a%aspk(i)
      end do

      call mrgsrt(nnzero,tmp,aux,ircode)
      if (ircode.eq.0) call zreordvn(nnzero,a%aspk,tmp,a%ia1,aux)
      !     .... Order with key a%ia1 (COLUMN INDEX) ...
      i    = 1
      j    = i
      !     .... order with key tmp (row index) ...
      do 
        if (i > nnzero) exit
        do 
          if (j > nnzero) exit
          if (tmp(j) /= tmp(i)) exit
          j = j+1
          !              if (j.eq.(nnzero+1)) exit
        enddo
        iel = j - i
        call mrgsrt(iel,a%ia1(i),aux,ircode)
        if (ircode == 0) call zreordvn(iel,a%aspk(i),tmp(i),&
             &   a%ia1(i), aux)
        i = j
      enddo

      ! convert to csr format     
      iel = 1
      a%ia2(1) = 1
      do i = a%ia2(1), nrow

        do 
          if (iel > nnzero) exit
          if (tmp(iel) /= i) exit
          iel = iel + 1
        enddo
        a%ia2(i+1) = iel
      enddo
      deallocate(aux,tmp)

    else if ((type == 'complex').and.(sym == 'symmetric')) then
      ! we are generally working with non-symmetric matrices, so
      ! we de-symmetrize what we are about to read
      allocate(a%aspk(2*nnzero),a%ia1(2*nnzero),&
           & a%ia2(2*nnzero),as_loc(2*nnzero),&
           & ia1_loc(2*nnzero),ia2_loc(2*nnzero),&
           &a%pl(nrow),a%pr(nrow), stat = ircode)
      if (ircode /= 0)   goto 993
      do i=1,nnzero
        read(infile,fmt=*,end=902) a%ia1(i),a%ia2(i),a%aspk(i)
      end do

      liwork = 2*nnzero+2
      allocate(iwork(liwork), stat = ircode)
      if (ircode /= 0)   goto 993  
      ! After this call NNZERO contains the actual value for
      ! desymetrized matrix
      call zdesym(nrow, a%aspk, a%ia2, a%ia1, as_loc, ia2_loc,&
           & ia1_loc, iwork, nnzero, nzr)     

      deallocate(a%aspk,a%ia1,a%ia2)
      nnzero=nzr
!!$      call spreall(a,nzr,ircode)
      if (ircode /= 0)   goto 993
      allocate(tmp(nzr),stat=ircode)
      if (ircode /= 0)   goto 993

      a%aspk(1:nzr) = as_loc(1:nzr)
      a%ia1(1:nzr)  = ia2_loc(1:nzr)
      tmp(1:nzr)    = ia1_loc(1:nzr)


      iel = 1
      a%ia2(1) = 1
      do i = 1, nrow
        do 
          if (tmp(iel) /= i) exit
          iel = iel + 1
          if (iel > nzr) exit
        enddo
        a%ia2(i+1) = iel
      enddo

      deallocate(as_loc, ia1_loc, ia2_loc,tmp,iwork)
    else
      write(0,*) 'read_matrix: matrix type not yet supported'
      iret=904
    end if
    if (infile/=5) close(infile)
    return 

    ! open failed
901 iret=901
    write(0,*) 'read_matrix: could not open file ',filename,' for input'
    return
902 iret=902
    write(0,*) 'READ_MATRIX: Unexpected end of file '
    return
993 iret=993
    write(0,*) 'READ_MATRIX: Memory allocation failure'
    return
  end subroutine zmm_mat_read



  subroutine dmm_mat_write(a,mtitle,iret,eiout,filename)
    use typesp
    implicit none
    type(d_spmat), intent(in)  :: a
    integer, intent(out)        :: iret
    character(len=*), intent(in) :: mtitle
    integer, optional, intent(in)          :: eiout
    character(len=*), optional, intent(in) :: filename
    integer                     :: iout


    iret = 0

    if (present(filename)) then 
      if (filename=='-') then 
        iout=6
      else
        if (present(eiout)) then 
          iout = eiout
        else
          iout=99
        endif
        open(iout,file=filename, err=901, action='WRITE')
      endif
    else 
      if (present(eiout)) then 
        iout = eiout   
      else
        iout=6
      endif
    endif

    call dcsprt(a%m,a%k,a%fida,a%descra,a%aspk,a%ia1,a%ia2,a%infoa,&
         & mtitle,iout,iret)

    if (iout /= 6) close(iout)


!!$    write(outfile(9:),998) '.xrhs'         
!!$    open (iout,file=outfile,status='replace',err=901)
!!$    write(iout,fmt=997)
!!$    write(iout,fmt=996) mtitle
!!$    write(iout,fmt=995) 'Number of equations  ',nrow
!!$    write(iout,fmt=995) 'Number of iterations to convergence ',iter
!!$    write(iout,fmt=996)
!!$    write(iout,fmt=996)  'index      comp. solution     Right hand side'  
!!$    write(iout,fmt=997)
!!$    do i=1, nrow
!!$      write(iout,993) i,x(i),rhs(i)
!!$993   format(i5,4(1x,e12.6))
!!$    enddo
!!$    close(iout)
!!$    !$$$        call system('gzip -f9 '//outfile)

    return

901 continue 
    iret=901
    write(0,*) 'Error while opening ',filename
    return
  end subroutine dmm_mat_write


!!$  subroutine lowerc(string,pos,len)
!!$    integer pos, len
!!$    character(len=*) string
!!$
!!$    character(len=26), parameter :: lcase='abcdefghijklmnopqrstuvwxyz',&
!!$         &  ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
!!$
!!$    do i=pos,pos+len-1
!!$      k = index(ucase,string(i:i))
!!$      if (k.ne.0) string(i:i) = lcase(k:k)
!!$    enddo
!!$    return
!!$  end subroutine lowerc

  subroutine desym(nrow,a,ja,ia,as,jas,ias,aux,nnzero,nzr)
    implicit none
    !      .. scalar arguments ..                                              
    integer :: nrow,nnzero,value,index,ptr, nzr
    !     .. array arguments ..                                                     
    real(kind(1.d0)) ::  a(*),as(*)              
    integer :: ia(*),ias(*),jas(*),ja(*),aux(*)                
    !     .. local scalars ..                                                       
    integer :: i,iaw1,iaw2,iawt,j,jpt,k,kpt,ldim,nzl,js,iret,nel,diagel
    !     ..                                                                        

    nel = 0
    diagel=0

    do i=1, nnzero
      as(i)  = a(i)
      jas(i) = ja(i)
      ias(i) = ia(i)
      if(ja(i) < ia(i)) then !this control avoids malfunctions in the cases 
        ! where the matrix is declared symmetric but all its elements are
        ! explicitly stored see young1c.mtx from "Matrix Market".
        ! Nominally Matrix Market only stores lower triangle. 
        nel    = nel+1
        as(nnzero+nel)  = a(i)
        jas(nnzero+nel) = ia(i)
        ias(nnzero+nel) = ja(i)
      end if
    end do
    if (nel == 0) then ! Something strange is going on
      write(0,*) 'Warning: DESYM did not copy anything in the upper triangle. '
      write(0,*) '         This feels wrong!!!!! '
    endif
    !     .... order with key ias ...
    nzr = nnzero + nel
    call mrgsrt(nzr,ias,aux,iret)
    if (iret == 0) call reordvn(nzr,as,ias,jas,aux)
    !     .... order with key jas ...

    i    = 1
    j    = i
    do 
      if (i > nzr) exit
      do
        if (j > nzr) exit
        if (ias(j) /= ias(i)) exit
        j = j+1
      enddo
      nzl = j - i
      call mrgsrt(nzl,jas(i),aux,iret)
      if (iret.eq.0) call reordvn(nzl,as(i),ias(i),jas(i),aux)
      i = j

    enddo

    return                                                                    
  end subroutine desym

    subroutine zdesym(nrow,a,ja,ia,as,jas,ias,aux,nnzero,nzr)
      implicit none
      !      .. scalar arguments ..                                              
      integer :: nrow,nnzero,value,index,ptr, nzr
      !     .. array arguments ..                                                     
      complex(kind(1.d0)) ::  a(*),as(*)                               
      integer :: ia(*),ias(*),jas(*),ja(*),aux(*)                
      !     .. local scalars ..                                                       
      integer :: i,iaw1,iaw2,iawt,j,jpt,k,kpt,ldim,nzl,js,iret,nel,diagel
      !     ..                                                                        

      nel = 0
      diagel=0

      do i=1, nnzero
        as(i)  = a(i)
        jas(i) = ja(i)
        ias(i) = ia(i)
        if(ja(i) < ia(i)) then !this control avoids malfunctions in the cases 
            ! where the matrix is declared symmetric but all its elements are
            ! explicitly stored see young1c.mtx from "Matrix Market".
            ! Nominally Matrix Market only stores lower triangle. 
          nel    = nel+1
          as(nnzero+nel)  = a(i)
          jas(nnzero+nel) = ia(i)
          ias(nnzero+nel) = ja(i)
        end if
      end do
        
      !     .... order with key ias ...
      nzr = nnzero + nel
      call mrgsrt(nzr,ias,aux,iret)
      if (iret == 0) call zreordvn(nzr,as,ias,jas,aux)
      !     .... order with key jas ...

      i    = 1
      j    = i
      do 
        if (i > nzr) exit
        do
          if (j > nzr) exit
          if (ias(j) /= ias(i)) exit
          j = j+1
        enddo
        nzl = j - i
        call mrgsrt(nzl,jas(i),aux,iret)
        if (iret.eq.0) call zreordvn(nzl,as(i),ias(i),jas(i),aux)
        i = j

      enddo
      
      return                                                                    
    end subroutine zdesym

end module mmio
