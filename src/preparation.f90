!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine preparation
  use global_variables
  implicit none
  integer :: ix,ik,aion,bion, icopy
  real(8) :: x
  real(8) :: lambda
  real(8) :: sigma_t,xp1,xp2,xp3,xp4,xm1,xm2,xm3,xm4
  real(8) :: Za,Zb,xshift


  H=lattice_a/dble(Nx)
  dKx=2d0*pi/(dble(Nk)*lattice_a)

! Laplacian coefficients
  L_coef(0)=-30d0/12d0/(H**2)
  L_coef(1)=16d0/12d0/(H**2)
  L_coef(2)=-1d0/12d0/(H**2)
  L_coef(-1)=16d0/12d0/(H**2)
  L_coef(-2)=-1d0/12d0/(H**2)
! Gradient coefficients
  G_coef(0)=0d0
  G_coef(1)=8d0/12d0/H
  G_coef(2)=-1d0/12d0/H
  G_coef(-1)=-8d0/12d0/H
  G_coef(-2)=1d0/12d0/H


  allocate(Lx(Nx),Kx(NK))


  allocate(ztpsi(Nx),zhtpsi(Nx))
  allocate(Veff(Nx),spe(NB,NK))
  allocate(zpsi_GS(Nx,NB,NK))    
  allocate(zpsi_Ct(NB*NK,NK),zpsi_Ct_t(NB*NK,NK),zhpsi_Ct_t(NB*NK,NK))
  allocate( zpsi_Ct_t_Lan(NB*NK,NK,0:Nlanczos)); zpsi_Ct_t_Lan = 0d0

  do ix=1,Nx
    Lx(ix)=dble(ix-1)*H
  end do

  do ik=1,NK
    kx(ik)=dKx*dble(ik-1-NK/2)
!    kx(ik)=dKx*dble(ik-1)
  end do


  
  Veff = 0d0
  do icopy = -10,10
    do aion = 1,Nion
      do ix=1,Nx
        x = Lx(ix) - Rion(aion) + dble(icopy)*lattice_a
        Veff(ix) = Veff(ix) -Zion(aion)*int_pot_ei(x)
      end do
    end do
  end do

! ion
  allocate(Phi_FC(NK,Nion,NK,Nion), Fion(NK,Nion), Uion(NK,Nion))
  allocate(Dph_mat(NK*Nion,NK*Nion))
  allocate(w2_ph(Nion,NK),w_ph(Nion,NK))
  
  E_ii = 0d0
  do aion = 1,Nion
    do bion = 1,Nion
      do icopy = -10,10
        if(aion == bion .and. icopy == 0)cycle
        E_ii = E_ii &
          + Zion(aion)*Zion(bion)*int_pot_ii(Rion(aion)-Rion(bion)+lattice_a*icopy)
      end do
    end do
  end do
  E_ii = 0.5d0*E_ii

  return
end subroutine preparation
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
