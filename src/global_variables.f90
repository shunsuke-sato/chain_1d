!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
module global_variables

! Parameters
!  real(8),parameter :: pi=3.14159265358979323846d0
  real(8),parameter :: pi=4d0*atan(1d0)
  complex(8),parameter :: zI=(0d0,1d0)
  real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0

! Grid
  real(8) :: lattice_a, H, dKx
  integer :: Nx,Nk,NB,NBocc
  real(8),allocatable :: Lx(:),Kx(:)

! time propagation
  integer :: Nt
  real(8) :: dt

! Electrons
  complex(8),allocatable :: zpsi(:,:,:),zpsi_GS(:,:,:),ztpsi(:),zhtpsi(:)
  real(8),allocatable :: Veff(:)
  real(8) :: L_coef(-2:2),G_coef(-2:2)
  integer :: Nelec
  real(8),allocatable :: spe(:,:)
  complex(8),allocatable :: zpsi_Ct(:,:),zpsi_Ct_t(:,:),zhpsi_Ct_t(:,:)
  complex(8),allocatable :: zH_mat(:,:),zH0_mat(:,:),zW_mat(:,:),zP_mat(:,:,:)
  complex(8),allocatable :: zF_mat(:,:,:),zG_mat(:,:,:)
  complex(8),allocatable :: zF_mat_full(:,:,:,:),zG_mat_full(:,:,:,:)
  complex(8),allocatable :: zG3_mat(:,:,:)
  complex(8),allocatable :: zG3_mat_full(:,:,:,:)
! Lanczos method
  integer,parameter :: Nlanczos = 4
  complex(8),allocatable :: zpsi_Ct_t_Lan(:,:,:)
  
! Ions
  integer,parameter :: Nion = 2
  real(8) :: Zion(Nion), Rion(Nion), Mion(Nion)
  real(8),allocatable :: Phi_FC(:,:,:,:),Fion(:,:), Uion(:,:),Dph_mat(:,:)
  real(8),allocatable :: Uion_new(:,:),Uion_old(:,:),Vion(:,:)
  real(8),allocatable :: w2_ph(:,:),w_ph(:,:)
  real(8) :: beta_kB

! Energy
  real(8) :: E_tot, E_elec, E_ii,Ekin_ion
  real(8) :: E_tot0, E_elec0, E_ii0,Ekin_ion0 ! for t=0

! Laser fields
  real(8) :: Acx

! I/O parameter
  character(50) :: filename


! potential parameters
  real(8),parameter :: sg = 2.0d0, v0 = 1.00   !sg = 2.0d0, v0 = 1.0d0  
contains

  function int_pot(x)
    real(8) :: int_pot,x
    
    int_pot = v0/(sqrt(2d0*pi)*sg)*exp(-0.5d0*(x/sg)**2)
   
    return
  end function int_pot

  function int_pot_drv1(x)
    real(8) :: int_pot_drv1,x
    
    int_pot_drv1 = v0/(sqrt(2d0*pi)*sg)*(-x/sg**2)*exp(-0.5d0*(x/sg)**2)
   
    return
  end function int_pot_drv1


  function int_pot_drv2(x)
    real(8) :: int_pot_drv2,x
    
    int_pot_drv2 = v0/(sqrt(2d0*pi)*sg)*(-1d0/sg**2+x**2/sg**4)*exp(-0.5d0*(x/sg)**2)
!    int_pot_drv2 = (-1d0/sg**2+x**2/sg**4)*int_pot(x)
   
    return
  end function int_pot_drv2

  function int_pot_drv3(x)
    real(8) :: int_pot_drv3,x
    
    int_pot_drv3 = v0/(sqrt(2d0*pi)*sg)*(2d0*(x/sg**2)/sg**2-x/sg**2) &
      *exp(-0.5d0*(x/sg)**2)
   
    return
  end function int_pot_drv3

end module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
