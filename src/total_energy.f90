subroutine total_energy
  use global_variables
  implicit none
  complex(8) :: zCt(NK*NB)
  integer :: icell, icell2, aion, bion, ik
  real(8) :: la_full,x

  la_full = lattice_a*dble(NK)  

  E_elec = 0d0
  do ik = 1,NK
    zCt = matmul(zH_mat(:,:),zpsi_Ct(:,ik))
    E_elec = E_elec + 2d0*sum(conjg(zpsi_Ct(:,ik))*zCt(:))
  end do

  E_ii = 0d0
  do icell = 1,NK; do aion=1,Nion
    do icell2 = 1,NK; do bion=1,Nion
      if(icell2 == icell .and. aion == bion)cycle
      x = Rion(aion) + dble(icell-1)*lattice_a &
        - (Rion(bion) + dble(icell2-1)*lattice_a) &
        + Uion(icell,aion) - Uion(icell2,bion)
          
      E_ii = E_ii + Zion(aion)*Zion(bion)*( &
        int_pot(x) + int_pot(x + la_full) + int_pot(x - la_full) )
    end do; end do
  end do; end do

  E_tot = E_elec + E_ii
  
  return
end subroutine total_energy
