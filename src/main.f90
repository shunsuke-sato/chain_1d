!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
program main
  use global_variables
  implicit none
  integer :: iter
  integer :: ik,ib,ix1,ix2,ib2
  integer :: N_ex_elec
  real(8) :: s
  character(50) :: cnum
  real(8),allocatable :: nex_k(:,:,:)

  call set_parameters
  call preparation
  call ground_state
  call Electron_matrix_element
  call phonon_band

  stop
end program main
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
