module psb_s_genpde_mod

  use psb_base_mod, only : psb_spk_, psb_ipk_, psb_desc_type,&
       &  psb_sspmat_type, psb_s_vect_type, szero,&
       &  psb_s_base_sparse_mat, psb_s_base_vect_type

  interface 
    function s_func_3d(x,y,z) result(val)
      import :: psb_spk_
      real(psb_spk_), intent(in) :: x,y,z
      real(psb_spk_) :: val
    end function s_func_3d
  end interface 

  interface psb_gen_pde3d
    subroutine psb_s_gen_pde3d(ictxt,idim,a,bv,xv,desc_a,afmt,&
         & a1,a2,a3,b1,b2,b3,c,g,info,f,amold,vmold)
      !
      !   Discretizes the partial differential equation
      ! 
      !   a1 dd(u)  a2 dd(u)    a3 dd(u)    b1 d(u)   b2 d(u)  b3 d(u)  
      ! -   ------ -  ------ -  ------ +  -----  +  ------  +  ------ + c u = f
      !      dxdx     dydy       dzdz        dx       dy         dz   
      !
      ! with Dirichlet boundary conditions
      !   u = g 
      !
      !  on the unit cube  0<=x,y,z<=1.
      !
      !
      ! Note that if b1=b2=b3=c=0., the PDE is the  Laplace equation.
      !
      import  :: psb_ipk_, psb_desc_type, psb_sspmat_type, psb_s_vect_type, &
           & s_func_3d, psb_s_base_sparse_mat, psb_s_base_vect_type
      implicit none
      procedure(s_func_3d)  :: a1,a2,a3,c,b1,b2,b3,g
      integer(psb_ipk_)     :: idim
      type(psb_sspmat_type) :: a
      type(psb_s_vect_type) :: xv,bv
      type(psb_desc_type)   :: desc_a
      integer(psb_ipk_)     :: ictxt, info
      character(len=*)      :: afmt
      procedure(s_func_3d), optional :: f
      class(psb_s_base_sparse_mat), optional :: amold
      class(psb_s_base_vect_type), optional :: vmold
    end subroutine psb_s_gen_pde3d
  end interface


  interface 
    function s_func_2d(x,y) result(val)
      import :: psb_spk_
      real(psb_spk_), intent(in) :: x,y
      real(psb_spk_) :: val
    end function s_func_2d
  end interface 

  interface psb_gen_pde2d 
    subroutine psb_s_gen_pde2d(ictxt,idim,a,bv,xv,desc_a,afmt,&
         & a1,a2,b1,b2,c,g,info,f,amold,vmold)
      !
      !   Discretizes the partial differential equation
      ! 
      !   a1 dd(u)  a2 dd(u)    b1 d(u)   b2 d(u)   
      ! -   ------ -  ------  +  -----  +  ------  +  c u = f
      !      dxdx     dydy         dx       dy         
      !
      ! with Dirichlet boundary conditions
      !   u = g 
      !
      !  on the unit square  0<=x,y<=1.
      !
      !
      ! Note that if b1=b2=c=0., the PDE is the  Laplace equation.
      !
      import  :: psb_ipk_, psb_desc_type, psb_sspmat_type, psb_s_vect_type,&
           &  s_func_2d, psb_s_base_sparse_mat, psb_s_base_vect_type
      implicit none
      procedure(s_func_2d)  :: a1,a2,c,b1,b2,g
      integer(psb_ipk_)     :: idim
      type(psb_sspmat_type) :: a
      type(psb_s_vect_type) :: xv,bv
      type(psb_desc_type)   :: desc_a
      integer(psb_ipk_)     :: ictxt, info
      character(len=*)      :: afmt
      procedure(s_func_2d), optional :: f
      class(psb_s_base_sparse_mat), optional :: amold
      class(psb_s_base_vect_type), optional :: vmold
    end subroutine psb_s_gen_pde2d
  end interface

contains

  function s_null_func_3d(x,y,z) result(val)

    real(psb_spk_), intent(in) :: x,y,z
    real(psb_spk_) :: val
    
    val = szero

  end function s_null_func_3d

  function s_null_func_2d(x,y) result(val)

    real(psb_spk_), intent(in) :: x,y
    real(psb_spk_) :: val
    
    val = szero

  end function s_null_func_2d


end module psb_s_genpde_mod
