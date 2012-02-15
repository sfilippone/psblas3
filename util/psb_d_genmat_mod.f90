module psb_d_genmat_mod
  use psb_base_mod
  interface 
    function d_func_3d(x,y,z) result(val)
      import :: psb_dpk_
      real(psb_dpk_), intent(in) :: x,y,z
      real(psb_dpk_) :: val
    end function d_func_3d
  end interface 

  interface 
    subroutine gen_prob3d(ictxt,idim,a,bv,xv,desc_a,afmt,a1,a2,a3,b1,b2,b3,c,g,info)
      !
      !   Discretizes the partial differential equation
      ! 
      !   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)  
      ! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = 0
      !      dxdx     dydy       dzdz        dx       dy         dz   
      !
      ! with Dirichlet boundary conditions
      !   u = g 
      !
      !  on the unit cube  0<=x,y,z<=1.
      !
      !
      ! Note that if a1=a2=a3=c=0., the PDE is the  Laplace equation.
      !
      import  :: psb_ipk_, psb_desc_type, psb_dspmat_type, psb_d_vect_type, d_func_3d
      implicit none
      procedure(d_func_3d)  :: a1,a2,a3,c,b1,b2,b3,g
      integer(psb_ipk_)     :: idim
      type(psb_dspmat_type) :: a
      type(psb_d_vect_type) :: xv,bv
      type(psb_desc_type)   :: desc_a
      integer(psb_ipk_)     :: ictxt, info
      character             :: afmt*5
    end subroutine gen_prob3d
  end interface


end module psb_d_genmat_mod
