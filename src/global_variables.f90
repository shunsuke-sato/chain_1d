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
  integer,parameter :: Nlanczos = 9
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
  real(8),parameter :: sg_ii = 2.5d0, v0_ii = 1d2   !sg = 2.0d0, v0 = 1.0d0  
  real(8),parameter :: sg_ii2 = 2d0*sg_ii, v0_ii2 = -0.75d0*v0_ii   !sg = 2.0d0, v0 = 1.0d0  
  real(8),parameter :: sg_ei = 2.5d0, v0_ei = 1.00   !sg = 2.0d0, v0 = 1.0d0  
contains

! Ion-ion interaction
  function int_pot_ii(x)
    real(8) :: int_pot_ii,x
    
    int_pot_ii = v0_ii/(sqrt(2d0*pi)*sg_ii)*exp(-0.5d0*(x/sg_ii)**2)
    int_pot_ii = int_pot_ii &
      +v0_ii2/(sqrt(2d0*pi)*sg_ii2)*exp(-0.5d0*(x/sg_ii2)**2)
   
    return
  end function int_pot_ii

  function int_pot_drv1_ii(x)
    real(8) :: int_pot_drv1_ii,x
    
    int_pot_drv1_ii = v0_ii/(sqrt(2d0*pi)*sg_ii)*(-x/sg_ii**2) &
      *exp(-0.5d0*(x/sg_ii)**2)

    int_pot_drv1_ii = int_pot_drv1_ii &
      +v0_ii2/(sqrt(2d0*pi)*sg_ii2)*(-x/sg_ii2**2) &
      *exp(-0.5d0*(x/sg_ii2)**2)
   
    return
  end function int_pot_drv1_ii


  function int_pot_drv2_ii(x)
    real(8) :: int_pot_drv2_ii,x
    
    int_pot_drv2_ii = v0_ii/(sqrt(2d0*pi)*sg_ii)*(-1d0/sg_ii**2+x**2/sg_ii**4) &
      *exp(-0.5d0*(x/sg_ii)**2)
    int_pot_drv2_ii = int_pot_drv2_ii &
      +v0_ii2/(sqrt(2d0*pi)*sg_ii2)*(-1d0/sg_ii2**2+x**2/sg_ii2**4) &
      *exp(-0.5d0*(x/sg_ii2)**2)
!    int_pot_drv2 = (-1d0/sg**2+x**2/sg**4)*int_pot(x)
   
    return
  end function int_pot_drv2_ii

  function int_pot_drv3_ii(x)
    real(8) :: int_pot_drv3_ii,x
    
    int_pot_drv3_ii = v0_ii/(sqrt(2d0*pi)*sg_ii)*(2d0*(x/sg_ii**2)/sg_ii**2-x/sg_ii**2) &
      *exp(-0.5d0*(x/sg_ii)**2)
   
    return
  end function int_pot_drv3_ii

! electron-ion interaction
  function int_pot_ei(x)
    real(8) :: int_pot_ei,x
    
    int_pot_ei = v0_ei/(sqrt(2d0*pi)*sg_ei)*exp(-0.5d0*(x/sg_ei)**2)
   
    return
  end function int_pot_ei

  function int_pot_drv1_ei(x)
    real(8) :: int_pot_drv1_ei,x
    
    int_pot_drv1_ei = v0_ei/(sqrt(2d0*pi)*sg_ei)*(-x/sg_ei**2) &
      *exp(-0.5d0*(x/sg_ei)**2)
   
    return
  end function int_pot_drv1_ei


  function int_pot_drv2_ei(x)
    real(8) :: int_pot_drv2_ei,x
    
    int_pot_drv2_ei = v0_ei/(sqrt(2d0*pi)*sg_ei)*(-1d0/sg_ei**2+x**2/sg_ei**4) &
      *exp(-0.5d0*(x/sg_ei)**2)
!    int_pot_drv2 = (-1d0/sg**2+x**2/sg**4)*int_pot(x)
   
    return
  end function int_pot_drv2_ei

  function int_pot_drv3_ei(x)
    real(8) :: int_pot_drv3_ei,x
    
    int_pot_drv3_ei = v0_ei/(sqrt(2d0*pi)*sg_ei)*(2d0*(x/sg_ei**2)/sg_ei**2-x/sg_ei**2) &
      *exp(-0.5d0*(x/sg_ei)**2)
   
    return
  end function int_pot_drv3_ei

end module global_variables
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
