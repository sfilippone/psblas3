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
module psb_blacs_mod
  
  interface gebs2d
    module procedure igebs2ds, igebs2dv, igebs2dm,&
         &           dgebs2ds, dgebs2dv, dgebs2dm,&
         &           zgebs2ds, zgebs2dv, zgebs2dm
  end interface

  interface gebr2d
    module procedure igebr2ds, igebr2dv, igebr2dm,&
         &           dgebr2ds, dgebr2dv, dgebr2dm,&
         &           zgebr2ds, zgebr2dv, zgebr2dm    
  end interface


  interface gesd2d
    module procedure igesd2ds, igesd2dv, igesd2dm,&
         &           dgesd2ds, dgesd2dv, dgesd2dm,&
         &           zgesd2ds, zgesd2dv, zgesd2dm
  end interface

  interface gerv2d
    module procedure igerv2ds, igerv2dv, igerv2dm,&
         &           dgerv2ds, dgerv2dv, dgerv2dm,&
         &           zgerv2ds, zgerv2dv, zgerv2dm
  end interface

  interface gsum2d
    module procedure igsum2ds, igsum2dv, igsum2dm,&
         &           dgsum2ds, dgsum2dv, dgsum2dm,&
         &           zgsum2ds, zgsum2dv, zgsum2dm
  end interface

  interface gamx2d
    module procedure igamx2ds, igamx2dv, igamx2dm,&
         &           dgamx2ds, dgamx2dv, dgamx2dm,&
         &           zgamx2ds, zgamx2dv, zgamx2dm
  end interface


  interface gamn2d
    module procedure igamn2ds, igamn2dv, igamn2dm,&
         &           dgamn2ds, dgamn2dv, dgamn2dm,&
         &           zgamn2ds, zgamn2dv, zgamn2dm
  end interface

 
contains 
  
  subroutine igebs2ds(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt,dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    character             :: top_ 
    
    interface 
      subroutine igebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine igebs2d
    end interface

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call igebs2d(icontxt,scope,top_,1,1,dat,1)
    
  end subroutine igebs2ds

  subroutine igebs2dv(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt,dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    

    interface 
      subroutine igebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine igebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call igebs2d(icontxt,scope,top_,size(dat,1),1,dat,size(dat,1))
    
  end subroutine igebs2dv

  subroutine igebs2dm(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt,dat(:,:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine igebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine igebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call igebs2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine igebs2dm



  subroutine dgebs2ds(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(in)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(icontxt,scope,top_,1,1,dat,1)
    
  end subroutine dgebs2ds

  subroutine dgebs2dv(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(in)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(icontxt,scope,top_,size(dat),1,dat,size(dat))
    
  end subroutine dgebs2dv

  subroutine dgebs2dm(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(in)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine dgebs2dm



  subroutine zgebs2ds(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(in)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(icontxt,scope,top_,1,1,dat,1)
    
  end subroutine zgebs2ds

  subroutine zgebs2dv(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(in)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(icontxt,scope,top_,size(dat),1,dat,size(dat))
    
  end subroutine zgebs2dv

  subroutine zgebs2dm(icontxt,scope,dat,top)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(in)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(icontxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine zgebs2dm





  subroutine dgebr2ds(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine dgebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol
    
    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call dgebr2d(icontxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine dgebr2ds

  subroutine dgebr2dv(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call dgebr2d(icontxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine dgebr2dv

  subroutine dgebr2dm(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call dgebr2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine dgebr2dm




  subroutine zgebr2ds(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine zgebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol

    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call zgebr2d(icontxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine zgebr2ds

  subroutine zgebr2dv(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call zgebr2d(icontxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine zgebr2dv

  subroutine zgebr2dm(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call zgebr2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine zgebr2dm



  subroutine igebr2ds(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine igebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol

    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call igebr2d(icontxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine igebr2ds

  subroutine igebr2dv(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call igebr2d(icontxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine igebr2dv

  subroutine igebr2dm(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igebr2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = 0
    case('C','c')
      rrt_ = 0
      crt_ = mycol
    case('A','a')
      rrt_ = 0
      crt_ = 0
    case default
      rrt_ = 0
      crt_ = 0
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call igebr2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine igebr2dm



  subroutine dgesd2ds(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    call dgesd2d(icontxt,1,1,dat,1,rdst,cdst)
    
  end subroutine dgesd2ds


  subroutine dgesd2dv(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    call dgesd2d(icontxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine dgesd2dv

  subroutine dgesd2dm(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    call dgesd2d(icontxt,size(dat,1),size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine dgesd2dm


  subroutine igesd2ds(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    integer, intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    call igesd2d(icontxt,1,1,dat,1,rdst,cdst)
    
  end subroutine igesd2ds


  subroutine igesd2dv(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    integer, intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    call igesd2d(icontxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine igesd2dv

  subroutine igesd2dm(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    integer, intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    call igesd2d(icontxt,size(dat,1),size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine igesd2dm



  subroutine zgesd2ds(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    call zgesd2d(icontxt,1,1,dat,1,rdst,cdst)
    
  end subroutine zgesd2ds


  subroutine zgesd2dv(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    call zgesd2d(icontxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine zgesd2dv

  subroutine zgesd2dm(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgesd2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    call zgesd2d(icontxt,size(dat,1),size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine zgesd2dm



  subroutine dgerv2ds(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    call dgerv2d(icontxt,1,1,dat,1,rdst,cdst)
    
  end subroutine dgerv2ds


  subroutine dgerv2dv(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    call dgerv2d(icontxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine dgerv2dv

  subroutine dgerv2dm(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    call dgerv2d(icontxt,size(dat,1),size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine dgerv2dm


  subroutine igerv2ds(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface

    call igerv2d(icontxt,1,1,dat,1,rdst,cdst)
    
  end subroutine igerv2ds


  subroutine igerv2dv(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface

    call igerv2d(icontxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine igerv2dv

  subroutine igerv2dm(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface

    call igerv2d(icontxt,size(dat,1),size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine igerv2dm



  subroutine zgerv2ds(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface

    call zgerv2d(icontxt,1,1,dat,1,rdst,cdst)
    
  end subroutine zgerv2ds


  subroutine zgerv2dv(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface

    call zgerv2d(icontxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine zgerv2dv

  subroutine zgerv2dm(icontxt,dat,rdst,cdst)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgerv2d(icontxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface

    call zgerv2d(icontxt,size(dat,1),size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine zgerv2dm



  subroutine dgsum2ds(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine dgsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call dgsum2d(icontxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine dgsum2ds

  subroutine dgsum2dv(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call dgsum2d(icontxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine dgsum2dv

  subroutine dgsum2dm(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    real(kind(1.d0)), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call dgsum2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine dgsum2dm



  subroutine igsum2ds(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine igsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call igsum2d(icontxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine igsum2ds

  subroutine igsum2dv(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call igsum2d(icontxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine igsum2dv

  subroutine igsum2dm(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    integer, intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call igsum2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine igsum2dm



  subroutine zgsum2ds(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine zgsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    call zgsum2d(icontxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine zgsum2ds

  subroutine zgsum2dv(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call zgsum2d(icontxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine zgsum2dv

  subroutine zgsum2dm(icontxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: icontxt
    complex(kind(1.d0)), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgsum2d(icontxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    call zgsum2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine zgsum2dm




  subroutine dgamx2ds(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    real(kind(1.d0)), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine dgamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call dgamx2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call dgamx2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine dgamx2ds


  subroutine dgamx2dv(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).and.present(cia)) then 
      call dgamx2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamx2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine dgamx2dv

  subroutine dgamx2dm(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld,ldia
        real(kind(1.d0)), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine dgamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call dgamx2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamx2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine dgamx2dm



  subroutine igamx2ds(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    integer, intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine igamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        integer, intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine igamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call igamx2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call igamx2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine igamx2ds


  subroutine igamx2dv(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    integer, intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        integer, intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine igamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call igamx2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call igamx2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine igamx2dv

  subroutine igamx2dm(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    integer, intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld,ldia
        integer, intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine igamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call igamx2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call igamx2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine igamx2dm

  

  subroutine zgamx2ds(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    complex(kind(1.d0)), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine zgamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call zgamx2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call zgamx2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine zgamx2ds


  subroutine zgamx2dv(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    complex(kind(1.d0)), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call zgamx2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamx2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine zgamx2dv

  subroutine zgamx2dm(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    complex(kind(1.d0)), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamx2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld,ldia
        complex(kind(1.d0)), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine zgamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call zgamx2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamx2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine zgamx2dm
  

  subroutine dgamn2ds(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    real(kind(1.d0)), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine dgamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call dgamn2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call dgamn2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine dgamn2ds


  subroutine dgamn2dv(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    real(kind(1.d0)), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        real(kind(1.d0)), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call dgamn2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamn2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine dgamn2dv

  subroutine dgamn2dm(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    real(kind(1.d0)), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld,ldia
        real(kind(1.d0)), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine dgamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call dgamn2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamn2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine dgamn2dm



  subroutine igamn2ds(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    integer, intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine igamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        integer, intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine igamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call igamn2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call igamn2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine igamn2ds


  subroutine igamn2dv(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    integer, intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        integer, intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine igamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call igamn2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call igamn2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine igamn2dv

  subroutine igamn2dm(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    integer, intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld,ldia
        integer, intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine igamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call igamn2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call igamn2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine igamn2dm

  

  subroutine zgamn2ds(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    complex(kind(1.d0)), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine zgamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif
    
    if (present(ria).or.present(cia)) then 
      call zgamn2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call zgamn2d(icontxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine zgamn2ds


  subroutine zgamn2dv(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    complex(kind(1.d0)), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld
        complex(kind(1.d0)), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call zgamn2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamn2d(icontxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine zgamn2dv

  subroutine zgamn2dm(icontxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: icontxt
    complex(kind(1.d0)), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamn2d(icontxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: icontxt,m,n,ld,ldia
        complex(kind(1.d0)), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine zgamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(icontxt,nrows,ncols,myrow,mycol)
    select case(scope)
    case('R','r')
      rrt_ = myrow
      crt_ = -1
    case('C','c')
      rrt_ = -1
      crt_ = mycol
    case('A','a')
      rrt_ = -1
      crt_ = -1
    case default
      rrt_ = -1
      crt_ = -1
    end select


    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif
    if (present(rrt)) then
      rrt_ = rrt
    endif
    if (present(crt)) then 
      crt_ = crt
    endif

    if (present(ria).and.present(cia)) then 
      call zgamn2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamn2d(icontxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine zgamn2dm
  

end module psb_blacs_mod
