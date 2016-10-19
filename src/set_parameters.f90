!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine set_parameters
  use global_variables

! Set parameters
  NK=32
  Nx = 64
  NB=6
  NBocc=1
  Nelec=1
  lattice_a = 10.0d0

  filename = "tst"
  return
end subroutine set_parameters
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
