!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine set_parameters
  use global_variables

! Set parameters
  NK = 64
  Nx = 256
  NB=4
  NBocc=1
  Nelec=1
  lattice_a = 10.30d0 ! 6.0d0

! ion parameters
  Zion(1) = 2d0*2d0/3d0; Mion(1) = 70d0*2000d0
  Zion(2) = 2d0*1d0/3d0; Mion(2) = 70d0*2000d0
  Rion(1) = 0.0d0
!  Rion(2) = 0.5d0*(0.4906511171875000d0 + 0.4906511175937500d0)*lattice_a 
  Rion(2) = 0.5d0*lattice_a 



! time-propagation
  dt = 1.0d0
  Nt = 200000

  filename = "tst"
  return
end subroutine set_parameters
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
