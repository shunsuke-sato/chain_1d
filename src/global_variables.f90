!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
module global_variables

! Parameters
  real(8),parameter :: pi=3.14159265358979323846d0
  complex(8),parameter :: zI=(0d0,1d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0

! Grid
  real(8) :: lattice_a, H, dKx
  integer :: Nx,Nk,NB,NBocc
  real(8),allocatable :: Lx(:),Kx(:)

! Electrons
  complex(8),allocatable :: zpsi(:,:,:),zpsi_GS(:,:,:),ztpsi(:),zhtpsi(:)
  real(8),allocatable :: Veff(:)
  real(8) :: L_coef(-2:2),G_coef(-2:2)
  integer :: Nelec
  real(8),allocatable :: spe(:,:)
  complex(8),allocatable :: zH_mat(:,:),zH0_mat(:,:),zF_mat(:,:,:)
  
! Ions
  integer,parameter :: Nion = 2
  real(8) :: Zion(Nion), Rion(Nion), Mion(Nion)
  real(8),allocatable :: Phi_FC(:,:,:,:),Fion(:,:)

! Energy
  real(8) :: E_ii

! Laser fields
  real(8) :: Acx

! I/O parameter
  character(50) :: filename

contains

  function int_pot(x)
    real(8) :: int_pot,x
    real(8),parameter :: sg = 1.0d0, v0 = 0.3d0
    
    int_pot = v0/(sqrt(2d0*pi)*sg)*exp(-0.5d0*(x/sg)**2)
   
    return
  end function int_pot

  function int_pot_drv1(x)
    real(8) :: int_pot_drv1,x
    real(8),parameter :: sg = 1.0d0, v0 = 0.3d0
    
    int_pot_drv1 = v0/(sqrt(2d0*pi)*sg)*(-x/sg**2)*exp(-0.5d0*(x/sg)**2)
   
    return
  end function int_pot_drv1

end module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
