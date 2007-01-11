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
module psb_spsb_mod

  interface
      real(kind(1.d0)) function ddot(n,dx,incx,dy,incy)
      real(kind(1.d0)) :: dx(*),dy(*),dtemp
      integer          :: i,incx,incy,ix,iy,m,mp1,n
    end function ddot
 end interface

  interface
      subroutine  daxpby(m, n, alpha, X, lldx, beta, Y, lldy, info)
      real(kind(1.d0)) :: X(lldx,*), Y(lldy,*)
      real(kind(1.d0)) :: alpha, beta
      integer          :: n, m, lldx, lldy, info
    end subroutine daxpby
 end interface

 interface
      subroutine dcsmm(trans,m,n,k,alpha,pl,fida,descra,a,ia1,ia2,&
     &   infoa,pr,b,ldb,beta,c,ldc,work,lwork,ierror)
      integer          ::  m,n,k,ldb,ldc,lwork, ierror
      character        ::  trans
      real(kind(1.d0)) ::  alpha,beta
      integer          ::  ia1(*),ia2(*),infoa(*),pl(*),pr(*)
      character        ::  descra*11, fida*5
      real(kind(1.d0)) ::  a(*),b(ldb,*),c(ldc,*),work(*)
    end subroutine dcsmm
 end interface

 interface
    real(kind(1.d0)) function dcsnmi(trans,m,n,fida,descra,a,ia1,ia2,&
         &                 infoa,ierror)
      integer          ::  m,n, ierror
      character        ::  trans
      integer          ::  ia1(*),ia2(*),infoa(*)
      character        ::  descra*11, fida*5
      real(kind(1.d0)) ::  a(*)
    end function dcsnmi
 end interface

 interface
      subroutine dcsrp(trans,m,n,fida,descra,ia1,ia2,infoa,&
     &  p,work,lwork,ierror)
      implicit none                                                      
      integer          :: lwork, m, n, ierror
      character        :: trans
      real(kind(1.d0)) :: work(lwork)
      integer          :: ia1(*), ia2(*), infoa(*), p(*), int_val(5)
      character        :: descra*11, fida*5
    end subroutine dcsrp
 end interface

 interface
      subroutine dcssm(trans,m,n,alpha,unitd,d,&
     &   pl,fidt,descrt,t,it1,it2,infot,pr,&
     &   b,ldb,beta,c,ldc,work,lwork,ierror)
      real(kind(1.d0)) :: alpha, beta
      integer          :: n, ldb, ldc, m, lwork, ierror
      character        :: unitd, trans
      real(kind(1.d0)) :: t(*), b(ldb,*), c(ldc,*), d(*), work(*)
      integer          :: it1(*), it2(*), infot(*), pl(*), pr(*)
      character        :: descrt*11, fidt*5
    end subroutine dcssm
 end interface

 interface
    subroutine dcsupd(m, n, fida, descra, a, ia1,&
         &  ia2, infoa, ia, ja, fidh, descrh, h, ih1, ih2,&
         &  infoh, ih, jh, flag, glob_to_loc,&
         &  iwork, liwork, ierror)
      integer          :: ia, ja, ih, jh, m, n,&
           &  ierror, flag, liwork 
      integer          :: ia1(*),ia2(*),ih1(*),ih2(*),&
           &  infoa(*),infoh(*),iwork(*),&
           &  glob_to_loc(*)
      character        :: descra*11,descrh*11, fida*5, fidh*5
      real(kind(1.d0)) :: a(*),h(*)    
    end subroutine dcsupd
 end interface
 
 interface
    subroutine dgelp(trans,m,n,p,b,ldb,work,lwork,ierror)
      integer          :: ldb, m, n, lwork, ierror
      character        :: trans
      real(kind(1.d0)) :: b(ldb,*), work(*)
      integer          :: p(*)
    end subroutine dgelp
 end interface
 
 interface
    subroutine dlpupd(m,n,perm,b,ldb,beta,c,ldc)
      integer          :: m, n, ldb, ldc
      real(kind(1.d0)) :: beta
      integer          :: perm(*)
      real(kind(1.d0)) :: b(ldb,*), c(ldc,*)
    end subroutine dlpupd
 end interface
 
 interface
      subroutine dswmm(trans,m,n,k,alpha,fida,descra,a,ia1,ia2,&
           & infoa,b,ldb,beta,c,ldc,work,lwork,ierror)
      integer          :: m,n,k,ldb,ldc,lwork,ierror
      character        :: trans
      real(kind(1.d0)) :: alpha,beta
      integer          :: ia1(*),ia2(*),infoa(*), int_val(5)
      character        :: descra*11, fida*5
      real(kind(1.d0)) :: a(*),b(ldb,*),c(ldc,*),work(*)
    end subroutine dswmm
 end interface


 interface
    subroutine dswprt(m,n,fida,descra,a,ia1,ia2,infoa,title,&
         & iout,ierror)
      integer       :: m,n,iout,ierror
      integer       :: ia1(*),ia2(*),infoa(*)
      character     :: descra*11, fida*5, title*(*)
    end subroutine dswprt
 end interface
 
 interface
      subroutine dswsm(trans,m,n,alpha,unitd,d,fidt,descrt,t,it1,it2,& 
           & infot,b,ldb,beta,c,ldc,work,lwork,ierror)
      integer          :: m, n, ldb, ldc, lwork, ierror
      character        :: unitd, trans
      real(kind(1.d0)) :: alpha, beta
      integer          :: it1(*), it2(*), infot(*)
      character        :: descrt*11, fidt*5
      real(kind(1.d0)) :: t(*), b(ldb,*), c(ldc,*), d(*), work(*)
    end subroutine dswsm
 end interface

 interface
      subroutine symbmm (n, m, l, ia, ja, diaga, ib,&
           & jb, diagb, ic, jc, diagc, index)
      integer          :: ia(*), ja(*), diaga, ib(*),&
           & jb(*), diagb, diagc, index(*)
      integer, pointer :: ic(:),jc(:)
      integer          :: nze, info
    end subroutine symbmm
 end interface

end module psb_spsb_mod
