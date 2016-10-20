!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine set_parameters
  use global_variables

! Set parameters
  NK = 16
  Nx = 128
  NB=6
  NBocc=1
  Nelec=1
  lattice_a = 6.0d0

! ion parameters
  Zion(1) = 1.5d0; Mion(1) = 2000d0*3d0
  Zion(2) = 0.5d0; Mion(2) = 1000d0*1d0
  Rion(1) = 0.0d0
  Rion(2) = 0.5d0*lattice_a

  filename = "tst"
  return
end subroutine set_parameters
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
