module psb_blacs_mod
  use psb_const_mod
  
  interface gebs2d
    module procedure igebs2ds, igebs2dv, igebs2dm,&
         &           dgebs2ds, dgebs2dv, dgebs2dm,&
         &           zgebs2ds, zgebs2dv, zgebs2dm,&
         &           sgebs2ds, sgebs2dv, sgebs2dm,&
         &           cgebs2ds, cgebs2dv, cgebs2dm
  end interface

  interface gebr2d
    module procedure igebr2ds, igebr2dv, igebr2dm,&
         &           dgebr2ds, dgebr2dv, dgebr2dm,&
         &           zgebr2ds, zgebr2dv, zgebr2dm,&
         &           sgebr2ds, sgebr2dv, sgebr2dm,&
         &           cgebr2ds, cgebr2dv, cgebr2dm    

  end interface

  interface gesd2d
    module procedure igesd2ds, igesd2dv, igesd2dm,&
         &           dgesd2ds, dgesd2dv, dgesd2dm,&
         &           zgesd2ds, zgesd2dv, zgesd2dm,&
         &           sgesd2ds, sgesd2dv, sgesd2dm,&
         &           cgesd2ds, cgesd2dv, cgesd2dm
  end interface

  interface gerv2d
    module procedure igerv2ds, igerv2dv, igerv2dm,&
         &           dgerv2ds, dgerv2dv, dgerv2dm,&
         &           zgerv2ds, zgerv2dv, zgerv2dm,&
         &           sgerv2ds, sgerv2dv, sgerv2dm,&
         &           cgerv2ds, cgerv2dv, cgerv2dm
  end interface

  interface gsum2d
    module procedure igsum2ds, igsum2dv, igsum2dm,&
         &           dgsum2ds, dgsum2dv, dgsum2dm,&
         &           zgsum2ds, zgsum2dv, zgsum2dm,&
         &           sgsum2ds, sgsum2dv, sgsum2dm,&
         &           cgsum2ds, cgsum2dv, cgsum2dm
  end interface

  interface gamx2d
    module procedure igamx2ds, igamx2dv, igamx2dm,&
         &           dgamx2ds, dgamx2dv, dgamx2dm,&
         &           zgamx2ds, zgamx2dv, zgamx2dm,&
         &           sgamx2ds, sgamx2dv, sgamx2dm,&
         &           cgamx2ds, cgamx2dv, cgamx2dm
  end interface


  interface gamn2d
    module procedure igamn2ds, igamn2dv, igamn2dm,&
         &           dgamn2ds, dgamn2dv, dgamn2dm,&
         &           zgamn2ds, zgamn2dv, zgamn2dm,&
         &           sgamn2ds, sgamn2dv, sgamn2dm,&
         &           cgamn2ds, cgamn2dv, cgamn2dm
  end interface
contains
  

  
  subroutine igebs2ds(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt,dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    character             :: top_ 
    
    interface 
      subroutine igebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine igebs2d
    end interface

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call igebs2d(ictxt,scope,top_,1,1,dat,1)
    
  end subroutine igebs2ds

  subroutine igebs2dv(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt,dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    

    interface 
      subroutine igebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
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

    call igebs2d(ictxt,scope,top_,size(dat,1),1,dat,size(dat,1))
    
  end subroutine igebs2dv

  subroutine igebs2dm(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt,dat(:,:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine igebs2d(ictxt,scope,top,m,n,v,ld)
        integer, intent(in)   :: ictxt,m,n,ld
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

    call igebs2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine igebs2dm



  subroutine dgebs2ds(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(in)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(ictxt,scope,top_,1,1,dat,1)
    
  end subroutine dgebs2ds

  subroutine dgebs2dv(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(in)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(ictxt,scope,top_,size(dat),1,dat,size(dat))
    
  end subroutine dgebs2dv

  subroutine dgebs2dm(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(in)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine dgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine dgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call dgebs2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine dgebs2dm


  subroutine sgebs2ds(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(in)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine sgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine sgebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call sgebs2d(ictxt,scope,top_,1,1,dat,1)
    
  end subroutine sgebs2ds

  subroutine sgebs2dv(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(in)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine sgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine sgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call sgebs2d(ictxt,scope,top_,size(dat),1,dat,size(dat))
    
  end subroutine sgebs2dv

  subroutine sgebs2dm(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(in)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine sgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine sgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call sgebs2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine sgebs2dm


  subroutine zgebs2ds(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(in)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(ictxt,scope,top_,1,1,dat,1)
    
  end subroutine zgebs2ds

  subroutine zgebs2dv(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(in)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(ictxt,scope,top_,size(dat),1,dat,size(dat))
    
  end subroutine zgebs2dv

  subroutine zgebs2dm(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(in)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine zgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine zgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call zgebs2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine zgebs2dm

  subroutine cgebs2ds(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(in)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine cgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(in)   :: v
        character, intent(in) :: scope, top
      end subroutine cgebs2d
    end interface
    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call cgebs2d(ictxt,scope,top_,1,1,dat,1)
    
  end subroutine cgebs2ds

  subroutine cgebs2dv(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(in)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine cgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(in)   :: v(*)
        character, intent(in) :: scope, top
      end subroutine cgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call cgebs2d(ictxt,scope,top_,size(dat),1,dat,size(dat))
    
  end subroutine cgebs2dv

  subroutine cgebs2dm(ictxt,scope,dat,top)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(in)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    
    interface 
      subroutine cgebs2d(ictxt,scope,top,m,n,v,ld)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(in)   :: v(ld,*)
        character, intent(in) :: scope, top
      end subroutine cgebs2d
    end interface

    character :: top_ 

    if (present(top)) then 
      top_ = top
    else
      top_ = ' '
    endif

    call cgebs2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1))
    
  end subroutine cgebs2dm





  subroutine dgebr2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine dgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol
    
    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call dgebr2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine dgebr2ds

  subroutine dgebr2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call dgebr2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine dgebr2dv

  subroutine dgebr2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call dgebr2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine dgebr2dm


  subroutine sgebr2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine sgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine sgebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol
    
    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call sgebr2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine sgebr2ds

  subroutine sgebr2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine sgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine sgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call sgebr2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine sgebr2dv

  subroutine sgebr2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine sgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine sgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call sgebr2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine sgebr2dm


  subroutine zgebr2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine zgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call zgebr2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine zgebr2ds

  subroutine zgebr2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call zgebr2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine zgebr2dv

  subroutine zgebr2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call zgebr2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine zgebr2dm

  subroutine cgebr2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine cgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine cgebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call cgebr2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine cgebr2ds

  subroutine cgebr2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine cgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine cgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call cgebr2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine cgebr2dv

  subroutine cgebr2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine cgebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine cgebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call cgebr2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine cgebr2dm



  subroutine igebr2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine igebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call igebr2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine igebr2ds

  subroutine igebr2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call igebr2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine igebr2dv

  subroutine igebr2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igebr2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igebr2d
    end interface

    character :: top_
    integer   :: nrows,ncols,myrow,mycol
    integer   :: rrt_, crt_

    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call igebr2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine igebr2dm



  subroutine sgesd2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine sgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine sgesd2d
    end interface

    call sgesd2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine sgesd2ds


  subroutine sgesd2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine sgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine sgesd2d
    end interface

    call sgesd2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine sgesd2dv

  subroutine sgesd2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_
    interface 
      subroutine sgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine sgesd2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call sgesd2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine sgesd2dm


  subroutine dgesd2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    call dgesd2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine dgesd2ds


  subroutine dgesd2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    call dgesd2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine dgesd2dv

  subroutine dgesd2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_
    interface 
      subroutine dgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine dgesd2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call dgesd2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine dgesd2dm


  subroutine igesd2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    integer, intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    call igesd2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine igesd2ds


  subroutine igesd2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    integer, intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    call igesd2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine igesd2dv

  subroutine igesd2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    integer, intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine igesd2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine igesd2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call igesd2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine igesd2dm



  subroutine cgesd2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine cgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine cgesd2d
    end interface

    call cgesd2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine cgesd2ds


  subroutine cgesd2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine cgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine cgesd2d
    end interface

    call cgesd2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine cgesd2dv

  subroutine cgesd2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine cgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine cgesd2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call cgesd2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine cgesd2dm


  subroutine zgesd2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(in)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(in)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    call zgesd2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine zgesd2ds


  subroutine zgesd2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(in)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(in)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    call zgesd2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine zgesd2dv

  subroutine zgesd2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(in)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine zgesd2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(in)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine zgesd2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call zgesd2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine zgesd2dm



  subroutine sgerv2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine sgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine sgerv2d
    end interface

    call sgerv2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine sgerv2ds


  subroutine sgerv2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine sgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine sgerv2d
    end interface

    call sgerv2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine sgerv2dv

  subroutine sgerv2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine sgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine sgerv2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call sgerv2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine sgerv2dm

  subroutine dgerv2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    call dgerv2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine dgerv2ds


  subroutine dgerv2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine dgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    call dgerv2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine dgerv2dv

  subroutine dgerv2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine dgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine dgerv2d
    end interface

    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call dgerv2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine dgerv2dm


  subroutine igerv2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface

    call igerv2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine igerv2ds


  subroutine igerv2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine igerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface

    call igerv2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine igerv2dv

  subroutine igerv2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_

    interface 
      subroutine igerv2d(ictxt,m,n,v,ld,rd,cd)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine igerv2d
    end interface


    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call igerv2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine igerv2dm



  subroutine cgerv2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine cgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine cgerv2d
    end interface

    call cgerv2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine cgerv2ds


  subroutine cgerv2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine cgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine cgerv2d
    end interface

    call cgerv2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine cgerv2dv

  subroutine cgerv2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_
    
    interface 
      subroutine cgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine cgerv2d
    end interface


    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call cgerv2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine cgerv2dm

  subroutine zgerv2ds(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface

    call zgerv2d(ictxt,1,1,dat,1,rdst,cdst)
    
  end subroutine zgerv2ds


  subroutine zgerv2dv(ictxt,dat,rdst,cdst)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat(:)
    integer, intent(in)  :: rdst,cdst
    
    interface 
      subroutine zgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v(*)
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface

    call zgerv2d(ictxt,size(dat),1,dat,size(dat),rdst,cdst)
    
  end subroutine zgerv2dv

  subroutine zgerv2dm(ictxt,dat,rdst,cdst,m)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat(:,:)
    integer, intent(in)  :: rdst,cdst
    integer, intent(in), optional :: m
    
    integer :: m_
    
    interface 
      subroutine zgerv2d(ictxt,m,n,v,ld,rd,cd)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v(ld,*)
        integer, intent(in)   :: rd,cd
      end subroutine zgerv2d
    end interface


    if (present(m)) then 
      m_ = m
    else
      m_ = size(dat,1)
    endif

    call zgerv2d(ictxt,m_,size(dat,2),dat,size(dat,1),rdst,cdst)
    
  end subroutine zgerv2dm



  subroutine sgsum2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine sgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine sgsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call sgsum2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine sgsum2ds

  subroutine sgsum2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine sgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine sgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call sgsum2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine sgsum2dv

  subroutine sgsum2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_spk_), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine sgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_spk_), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine sgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call sgsum2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine sgsum2dm

  subroutine dgsum2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine dgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call dgsum2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine dgsum2ds

  subroutine dgsum2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call dgsum2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine dgsum2dv

  subroutine dgsum2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    real(psb_dpk_), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine dgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine dgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call dgsum2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine dgsum2dm



  subroutine igsum2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine igsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call igsum2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine igsum2ds

  subroutine igsum2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call igsum2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine igsum2dv

  subroutine igsum2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    integer, intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine igsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        integer, intent(in)   :: ictxt,m,n,ld
        integer, intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine igsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call igsum2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine igsum2dm


  subroutine cgsum2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine cgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine cgsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call cgsum2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine cgsum2ds

  subroutine cgsum2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine cgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine cgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call cgsum2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine cgsum2dv

  subroutine cgsum2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_spk_), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine cgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine cgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call cgsum2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine cgsum2dm


  subroutine zgsum2ds(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt
    
    interface 
      subroutine zgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
    
    call zgsum2d(ictxt,scope,top_,1,1,dat,1,rrt_,crt_)
    
  end subroutine zgsum2ds

  subroutine zgsum2dv(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat(:)
    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v(*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call zgsum2d(ictxt,scope,top_,size(dat),1,dat,size(dat),rrt_,crt_)
    
  end subroutine zgsum2dv

  subroutine zgsum2dm(ictxt,scope,dat,top,rrt,crt)
    integer, intent(in)   :: ictxt
    complex(psb_dpk_), intent(inout)   :: dat(:,:)

    character, intent(in) :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional  :: rrt,crt

    interface 
      subroutine zgsum2d(ictxt,scope,top,m,n,v,ld,rrt,crt)
        use psb_const_mod
        integer, intent(in)   :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout)   :: v(ld,*)
        character, intent(in) :: scope, top
        integer, intent(in)   :: rrt,crt
      end subroutine zgsum2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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

    call zgsum2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),rrt_,crt_)
    
  end subroutine zgsum2dm


  subroutine sgamx2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine sgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        real(psb_spk_), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine sgamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call sgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call sgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine sgamx2ds


  subroutine sgamx2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine sgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        real(psb_spk_), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine sgamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call sgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call sgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine sgamx2dv

  subroutine sgamx2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine sgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        real(psb_spk_), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine sgamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call sgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call sgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine sgamx2dm


  subroutine dgamx2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine dgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call dgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call dgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine dgamx2ds


  subroutine dgamx2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call dgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine dgamx2dv

  subroutine dgamx2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        real(psb_dpk_), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine dgamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call dgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine dgamx2dm



  subroutine igamx2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine igamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
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


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call igamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call igamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine igamx2ds


  subroutine igamx2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
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


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call igamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call igamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine igamx2dv

  subroutine igamx2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld,ldia
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


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call igamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call igamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine igamx2dm

 
  subroutine cgamx2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine cgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine cgamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call cgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call cgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine cgamx2ds


  subroutine cgamx2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine cgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine cgamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call cgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call cgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine cgamx2dv

  subroutine cgamx2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine cgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        complex(psb_spk_), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine cgamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call cgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call cgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine cgamx2dm


  subroutine zgamx2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine zgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamx2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call zgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call zgamx2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine zgamx2ds


  subroutine zgamx2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamx2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call zgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamx2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine zgamx2dv

  subroutine zgamx2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamx2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        complex(psb_dpk_), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine zgamx2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call zgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamx2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine zgamx2dm
  

  subroutine sgamn2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine sgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        real(psb_spk_), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine sgamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call sgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call sgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine sgamn2ds


  subroutine sgamn2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine sgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        real(psb_spk_), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine sgamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call sgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call sgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine sgamn2dv

  subroutine sgamn2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_spk_), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine sgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        real(psb_spk_), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine sgamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call sgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call sgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine sgamn2dm

  subroutine dgamn2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine dgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call dgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call dgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine dgamn2ds


  subroutine dgamn2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        real(psb_dpk_), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine dgamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call dgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine dgamn2dv

  subroutine dgamn2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    real(psb_dpk_), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine dgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        real(psb_dpk_), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine dgamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call dgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call dgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine dgamn2dm


  subroutine igamn2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine igamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
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


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call igamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call igamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine igamn2ds


  subroutine igamn2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld
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


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call igamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call igamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine igamn2dv

  subroutine igamn2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    integer, intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine igamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        integer, intent(in)             :: ictxt,m,n,ld,ldia
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


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call igamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call igamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine igamn2dm


  subroutine cgamn2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine cgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine cgamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call cgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call cgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine cgamn2ds


  subroutine cgamn2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine cgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        complex(psb_spk_), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine cgamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call cgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call cgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine cgamn2dv

  subroutine cgamn2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_spk_), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine cgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        complex(psb_spk_), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine cgamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call cgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call cgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine cgamn2dm

  subroutine zgamn2ds(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt
    integer, intent(inout), optional :: ria,cia
    
    interface 
      subroutine zgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout) :: v
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamn2d
    end interface
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1),cia_(1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call zgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,1,rrt_,crt_)
      if (present(ria)) ria=ria_(1)
      if (present(cia)) cia=cia_(1)
    else
      call zgamn2d(ictxt,scope,top_,1,1,dat,1,ria_,cia_,-1,rrt_,crt_)
    endif
  end subroutine zgamn2ds


  subroutine zgamn2dv(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:)
    character, intent(in)            :: scope
    character, intent(in), optional  :: top
    integer, intent(inout), optional :: ria(:),cia(:)
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld
        complex(psb_dpk_), intent(inout) :: v(*)
        character, intent(in)           :: scope, top
        integer, intent(inout)          :: ria(*),cia(*)
        integer, intent(in)             :: rrt,crt,ldia
      end subroutine zgamn2d
    end interface

    integer   :: ldia_,ria_(1),cia_(1)
    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call zgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria,cia,min(size(ria),size(cia)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamn2d(ictxt,scope,top_,size(dat),1,dat,size(dat),&
           &  ria_,cia_,ldia_,rrt_,crt_)
    end if
    
  end subroutine zgamn2dv

  subroutine zgamn2dm(ictxt,scope,dat,top,ria,cia,rrt,crt)
    integer, intent(in)              :: ictxt
    complex(psb_dpk_), intent(inout)  :: dat(:,:)
    character, intent(in)            :: scope
    integer, intent(inout), optional :: ria(:,:),cia(:,:)
    character, intent(in), optional  :: top
    integer, intent(in), optional    :: rrt,crt

    interface 
      subroutine zgamn2d(ictxt,scope,top,m,n,v,ld,ria,cia,ldia,rrt,crt)
        use psb_const_mod
        integer, intent(in)             :: ictxt,m,n,ld,ldia
        complex(psb_dpk_), intent(inout) :: v(ld,*)
        integer, intent(inout)          :: ria(ldia,*),cia(ldia,*)
        character, intent(in)           :: scope, top
        integer, intent(in)             :: rrt,crt
      end subroutine zgamn2d
    end interface

    character :: top_ 
    integer   :: rrt_, crt_
    integer   :: ldia_,ria_(1,1),cia_(1,1)
    integer   :: nrows,ncols,myrow,mycol


    call blacs_gridinfo(ictxt,nrows,ncols,myrow,mycol)
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
      call zgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria,cia,min(size(ria,1),size(cia,1)),rrt_,crt_)
    else
      ldia_ = -1
      call zgamn2d(ictxt,scope,top_,size(dat,1),size(dat,2),dat,size(dat,1),&
           & ria_,cia_,ldia_,rrt_,crt_)
    end if
      
  end subroutine zgamn2dm

end module psb_blacs_mod
