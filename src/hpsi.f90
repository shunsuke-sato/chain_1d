!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine hpsi(ik)
  use global_variables
  implicit none
  integer :: ix,ik

  zhtpsi(1)=-0.5d0*( &
    &+L_coef(0)*ztpsi(1) &
    &+L_coef(1)*(ztpsi(2)+ztpsi(Nx)) &
    &+L_coef(2)*(ztpsi(3)+ztpsi(Nx-1))) &
    &-zI*(Kx(ik)+Acx)*( &
    &+G_coef(1)*(ztpsi(2)-ztpsi(Nx)) &
    &+G_coef(2)*(ztpsi(3)-ztpsi(Nx-1)))
  zhtpsi(2)=-0.5d0*( &
    &+L_coef(0)*ztpsi(2) &
    &+L_coef(1)*(ztpsi(3)+ztpsi(1)) &
    &+L_coef(2)*(ztpsi(4)+ztpsi(Nx))) &
    &-zI*(Kx(ik)+Acx)*( &
    &+G_coef(1)*(ztpsi(3)-ztpsi(1)) &
    &+G_coef(2)*(ztpsi(4)-ztpsi(Nx)))
  
  do ix=3,Nx-2
    zhtpsi(ix)=-0.5d0*( &
      &+L_coef(0)*ztpsi(ix) &
      &+L_coef(1)*(ztpsi(ix+1)+ztpsi(ix-1)) &
      &+L_coef(2)*(ztpsi(ix+2)+ztpsi(ix-2))) &
      &-zI*(Kx(ik)+Acx)*( &
      &+G_coef(1)*(ztpsi(ix+1)-ztpsi(ix-1)) &
      &+G_coef(2)*(ztpsi(ix+2)-ztpsi(ix-2)))
  end do

  zhtpsi(Nx-1)=-0.5d0*( &
    &+L_coef(0)*ztpsi(Nx-1) &
    &+L_coef(1)*(ztpsi(Nx)+ztpsi(Nx-2)) &
    &+L_coef(2)*(ztpsi(1)+ztpsi(Nx-3))) &
    &-zI*(Kx(ik)+Acx)*( &
    &+G_coef(1)*(ztpsi(Nx)-ztpsi(Nx-2)) &
    &+G_coef(2)*(ztpsi(1)-ztpsi(Nx-3)))
  zhtpsi(Nx)=-0.5d0*( &
    &+L_coef(0)*ztpsi(Nx) &
    &+L_coef(1)*(ztpsi(1)+ztpsi(Nx-1)) &
    &+L_coef(2)*(ztpsi(2)+ztpsi(Nx-2))) &
    &-zI*(Kx(ik)+Acx)*( &
    &+G_coef(1)*(ztpsi(1)-ztpsi(Nx-1)) &
    &+G_coef(2)*(ztpsi(2)-ztpsi(Nx-2)))
  
  
  zhtpsi(:)=zhtpsi(:)+(0.5d0*(kx(ik)+Acx)**2+Veff(:))*ztpsi(:)
  
  return
end subroutine hpsi
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
