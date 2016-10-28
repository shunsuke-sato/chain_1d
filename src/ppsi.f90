!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine ppsi(ik)
  use global_variables
  implicit none
  integer :: ix,ik

  zhtpsi(1)=-zI*( &
    &+G_coef(1)*(ztpsi(2)-ztpsi(Nx)) &
    &+G_coef(2)*(ztpsi(3)-ztpsi(Nx-1)))
  zhtpsi(2)=-zI*( &
    &+G_coef(1)*(ztpsi(3)-ztpsi(1)) &
    &+G_coef(2)*(ztpsi(4)-ztpsi(Nx)))
  
  do ix=3,Nx-2
    zhtpsi(ix)=-zI*( &
      &+G_coef(1)*(ztpsi(ix+1)-ztpsi(ix-1)) &
      &+G_coef(2)*(ztpsi(ix+2)-ztpsi(ix-2)))
  end do

  zhtpsi(Nx-1)=-zI*( &
    &+G_coef(1)*(ztpsi(Nx)-ztpsi(Nx-2)) &
    &+G_coef(2)*(ztpsi(1)-ztpsi(Nx-3)))
  zhtpsi(Nx)=-zI*( &
    &+G_coef(1)*(ztpsi(1)-ztpsi(Nx-1)) &
    &+G_coef(2)*(ztpsi(2)-ztpsi(Nx-2)))

  return
end subroutine ppsi
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
