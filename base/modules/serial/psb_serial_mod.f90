!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
module psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_realloc_mod
  use psb_string_mod
  use psb_sort_mod

  use psi_serial_mod, &
       & psb_gth => psi_gth,&
       & psb_sct => psi_sct

  use psb_s_serial_mod
  use psb_d_serial_mod
  use psb_c_serial_mod
  use psb_z_serial_mod

  interface psb_nrm1
    module procedure psb_snrm1, psb_dnrm1, psb_cnrm1, psb_znrm1
  end interface psb_nrm1

  interface psb_minreal
    module procedure psb_sminreal, psb_dminreal, psb_cminreal, psb_zminreal
  end interface psb_minreal


  interface psb_nspaxpby
    subroutine psb_d_nspaxpby(nz,iz,z,alpha, nx, ix, x, beta, ny,iy,y,info)
      use psb_const_mod, only : psb_ipk_, psb_spk_, psb_dpk_
      integer(psb_ipk_), intent(out)              :: nz
      integer(psb_ipk_), intent(out)              :: iz(:)
      real(psb_dpk_), intent (out)      :: z(:)
      integer(psb_ipk_), intent(in)               :: nx, ny
      integer(psb_ipk_), intent(in)               :: ix(:), iy(:)
      real(psb_dpk_), intent (in)       :: x(:), y(:)
      real(psb_dpk_), intent (in)       :: alpha, beta
      integer(psb_ipk_), intent(out)              :: info
    end subroutine psb_d_nspaxpby
  end interface psb_nspaxpby


contains

  elemental function psb_snrm1(x) result(res)
    real(psb_spk_), intent(in)  :: x
    real(psb_spk_)              :: res
    res = abs( x )
  end function psb_snrm1

  elemental function psb_dnrm1(x) result(res)
    real(psb_dpk_), intent(in) :: x
    real(psb_dpk_)             :: res
    res = abs(  x )
  end function psb_dnrm1

  elemental function psb_cnrm1(x) result(res)
    complex(psb_spk_), intent(in)  :: x
    real(psb_spk_)                 :: res
    res = abs( real( x ) ) + abs( aimag( x ) )  
  end function psb_cnrm1

  elemental function psb_znrm1(x) result(res)
    complex(psb_dpk_), intent(in)  :: x
    real(psb_dpk_)                 :: res
    res = abs( real( x ) ) + abs( aimag( x ) )  
  end function psb_znrm1

  elemental function psb_sminreal(x) result(res)
    real(psb_spk_), intent(in)  :: x
    real(psb_spk_)              :: res
    res = ( x )
  end function psb_sminreal

  elemental function psb_dminreal(x) result(res)
    real(psb_dpk_), intent(in) :: x
    real(psb_dpk_)             :: res
    res = (  x )
  end function psb_dminreal

  elemental function psb_cminreal(x) result(res)
    complex(psb_spk_), intent(in)  :: x
    real(psb_spk_)                 :: res
    res = min( real( x ) , aimag( x ) )  
  end function psb_cminreal

  elemental function psb_zminreal(x) result(res)
    complex(psb_dpk_), intent(in)  :: x
    real(psb_dpk_)                 :: res
    res = min( real( x ) , aimag( x ) )  
  end function psb_zminreal


  subroutine crot( n, cx, incx, cy, incy, c, s )
    !
    !  -- lapack auxiliary routine (version 3.0) --
    !     univ. of tennessee, univ. of california berkeley, nag ltd.,
    !     courant institute, argonne national lab, and rice university
    !     october 31, 1992
    !
    !     .. scalar arguments ..
    integer(psb_mpik_) :: incx, incy, n
    real(psb_spk_)    c
    complex(psb_spk_)   s
    !     ..
    !     .. array arguments ..
    complex(psb_spk_) cx( * ), cy( * )
    !     ..
    !
    !  purpose
    !  == = ====
    !
    !  zrot   applies a plane rotation, where the cos (c) is real and the
    !  sin (s) is complex, and the vectors cx and cy are complex.
    !
    !  arguments
    !  == = ======
    !
    !  n       (input) integer
    !          the number of elements in the vectors cx and cy.
    !
    !  cx      (input/output) complex*16 array, dimension (n)
    !          on input, the vector x.
    !          on output, cx is overwritten with c*x + s*y.
    !
    !  incx    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  cy      (input/output) complex*16 array, dimension (n)
    !          on input, the vector y.
    !          on output, cy is overwritten with -conjg(s)*x + c*y.
    !
    !  incy    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  c       (input) double precision
    !  s       (input) complex*16
    !          c and s define a rotation
    !             [  c          s  ]
    !             [ -conjg(s)   c  ]
    !          where c*c + s*conjg(s) = 1.0.
    !
    ! == = ==================================================================
    !
    !     .. local scalars ..
    integer(psb_mpik_) :: i, ix, iy
    complex(psb_spk_)         stemp
    !     ..
    !     .. intrinsic functions ..
    !     ..
    !     .. executable statements ..
    !
    if( n <= 0 ) return
    if( incx == 1 .and. incy == 1 ) then 
      !
      !     code for both increments equal to 1
      !
      do  i = 1, n
        stemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - conjg(s)*cx(i)
        cx(i) = stemp
      end do
    else
      !
      !     code for unequal increments or equal increments not equal to 1
      !
      ix = 1
      iy = 1
      if( incx < 0 )ix = ( -n+1 )*incx + 1
      if( incy < 0 )iy = ( -n+1 )*incy + 1
      do  i = 1, n
        stemp  = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - conjg(s)*cx(ix)
        cx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    return
  end subroutine crot
  !
  !
  subroutine crotg(ca,cb,c,s)
    complex(psb_spk_) ca,cb,s
    real(psb_spk_) c
    real(psb_spk_) norm,scale
    complex(psb_spk_) alpha
    !
    if (cabs(ca) == 0.0) then 
      !
      c = 0.0d0
      s = (1.0,0.0)
      ca = cb
      return
    end if
    !

    scale = cabs(ca) + cabs(cb)
    norm = scale*sqrt((cabs(ca/cmplx(scale,0.0)))**2 +&
         &   (cabs(cb/cmplx(scale,0.0)))**2)
    alpha = ca /cabs(ca)
    c = cabs(ca) / norm
    s = alpha * conjg(cb) / norm
    ca = alpha * norm
    !

    return
  end subroutine crotg



  subroutine zrot( n, cx, incx, cy, incy, c, s )
    !
    !  -- lapack auxiliary routine (version 3.0) --
    !     univ. of tennessee, univ. of california berkeley, nag ltd.,
    !     courant institute, argonne national lab, and rice university
    !     october 31, 1992
    !
    !     .. scalar arguments ..
    integer(psb_mpik_) :: incx, incy, n
    real(psb_dpk_)    c
    complex(psb_dpk_)   s
    !     ..
    !     .. array arguments ..
    complex(psb_dpk_) cx( * ), cy( * )
    !     ..
    !
    !  purpose
    !  == = ====
    !
    !  zrot   applies a plane rotation, where the cos (c) is real and the
    !  sin (s) is complex, and the vectors cx and cy are complex.
    !
    !  arguments
    !  == = ======
    !
    !  n       (input) integer
    !          the number of elements in the vectors cx and cy.
    !
    !  cx      (input/output) complex*16 array, dimension (n)
    !          on input, the vector x.
    !          on output, cx is overwritten with c*x + s*y.
    !
    !  incx    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  cy      (input/output) complex*16 array, dimension (n)
    !          on input, the vector y.
    !          on output, cy is overwritten with -conjg(s)*x + c*y.
    !
    !  incy    (input) integer
    !          the increment between successive values of cy.  incx <> 0.
    !
    !  c       (input) double precision
    !  s       (input) complex*16
    !          c and s define a rotation
    !             [  c          s  ]
    !             [ -conjg(s)   c  ]
    !          where c*c + s*conjg(s) = 1.0.
    !
    ! == = ==================================================================
    !
    !     .. local scalars ..
    integer(psb_mpik_) :: i, ix, iy
    complex(psb_dpk_)         stemp
    !     ..
    !     .. intrinsic functions ..
    intrinsic          dconjg
    !     ..
    !     .. executable statements ..
    !
    if( n <= 0 ) return
    if( incx == 1 .and. incy == 1 ) then 
      !
      !     code for both increments equal to 1
      !
      do  i = 1, n
        stemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - dconjg(s)*cx(i)
        cx(i) = stemp
      end do
    else
      !
      !     code for unequal increments or equal increments not equal to 1
      !
      ix = 1
      iy = 1
      if( incx < 0 )ix = ( -n+1 )*incx + 1
      if( incy < 0 )iy = ( -n+1 )*incy + 1
      do  i = 1, n
        stemp  = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - dconjg(s)*cx(ix)
        cx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
    end if
    return
  end subroutine zrot
  !
  !
  subroutine zrotg(ca,cb,c,s)
    complex(psb_dpk_) ca,cb,s
    real(psb_dpk_) c
    real(psb_dpk_) norm,scale
    complex(psb_dpk_) alpha
    !
    if (cdabs(ca) == 0.0d0) then 
      !
      c = 0.0d0
      s = (1.0d0,0.0d0)
      ca = cb
      return
    end if
    !

    scale = cdabs(ca) + cdabs(cb)
    norm = scale*dsqrt((cdabs(ca/cmplx(scale,0.0d0,kind=psb_dpk_)))**2 +&
         &   (cdabs(cb/cmplx(scale,0.0d0,kind=psb_dpk_)))**2)
    alpha = ca /cdabs(ca)
    c = cdabs(ca) / norm
    s = alpha * conjg(cb) / norm
    ca = alpha * norm
    !

    return
  end subroutine zrotg


end module psb_serial_mod

