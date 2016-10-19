!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine electron_matrix_element
  use global_variables
  implicit none
  real(8) :: x
  integer :: ik1,ik2,ib1,ib2,is1,is2,aion,icell,ix
  real(8) :: la_full
  integer :: Nmat
  integer :: iflag_cell(NK,Nion)
  real(8) :: tmp
  complex(8) :: zs

  allocate(zH_mat(NK*NB,NK*NB),zH0_mat(NK*NB,NK*NB), zF_mat(NK*NB,NK*NB,Nion))


! H0
  zH0_mat = 0d0
  is1=0
  do ik1=1,NK; do ib1=1,NB
    is1=is1+1
    zH0_mat(is1,is1)=spe(ib1,ik1)
  end do; end do

! Force matrix

! Start: cell serach
  iflag_cell=0
  do aion=1,Nion
    do icell=1,NK
      do ix=1,Nx
        x = Lx(ix) + lattice_a*icell
        tmp= int_pot(x-Rion(aion)) &
              +int_pot(x-Rion(aion)+la_full) &
              +int_pot(x-Rion(aion)-la_full)
        tmp = Zion(aion)*tmp
        if(abs(tmp)>1d-16)iflag_cell(icell,aion)=1
             
      end do
    end do
  end do
! End: cell serach

  do aion=1,Nion

    is1=0
    do ik1=1,NK; do ib1=1,NB
      is1=is1+1
      write(*,*)aion,is1
      is2=0
      do ik2=1,NK; do ib2=1,NB
        is2=is2+1
        if(is2<is1)cycle

        zs = 0d0
        do icell=1,NK
          if(iflag_cell(icell,aion)==0)cycle
          do ix=1,Nx
            x = Lx(ix) + lattice_a*icell
            zs = zs - Zion(aion)*exp(zI*(kx(ik2)-kx(ik1))*x) &
              *conjg(zpsi_GS(ix,ib1,ik1))*zpsi_GS(ix,ib2,ik2) &
              *( &
              int_pot_drv1(x-Rion(aion)) &
              +int_pot_drv1(x-Rion(aion)+la_full) &
              +int_pot_drv1(x-Rion(aion)-la_full) &
              )
          end do
        end do
        zF_mat(is1,is2,aion)=-zs*H/dble(NK)
        zF_mat(is2,is1,aion)=conjg(zF_mat(is1,is2,aion))
        
      end do; end do
    end do; end do

  end do

  return
end subroutine electron_matrix_element
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
