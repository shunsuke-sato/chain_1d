subroutine force_ion
  use global_variables
  implicit none
  complex(8) :: zCt(NK*NB),zAt(NK*NB,NK*NB)
  integer :: icell, icell2, aion, bion,ik
  real(8) :: la_full,x
!LAPACK ==
  complex(8),parameter :: zONE = (1d0,0d0),zZero = 0d0
  integer :: Nmax
  integer :: info
  Nmax = NB*NK
!LAPACK ==

  la_full = lattice_a*dble(NK)
  
! Force
  Fion = 0d0
  do aion = 1,Nion
    do icell=1,NK

! electron-ion interaction ==
      zAt = zF_mat_full(:,:,icell,aion) +2d0*zG_mat_full(:,:,icell,aion)*Uion(icell,aion)

      call zhemm('L','U',Nmax,NK,zONE,zAt,Nmax,zpsi_Ct,Nmax,zZero,zhpsi_Ct_t,Nmax)
      Fion(icell,aion) = Fion(icell,aion) + 2d0*sum(conjg(zpsi_Ct(:,:))*zhpsi_Ct_t(:,:))
!      do ik = 1,NK
!        zCt = matmul(zF_mat_full(:,:,icell,aion),zpsi_Ct(:,ik))
!        zCt = matmul(zAt(:,:),zpsi_Ct(:,ik))
!        Fion(icell,aion) = Fion(icell,aion) + 2d0*sum(conjg(zpsi_Ct(:,ik))*zCt(:))
!      end do

! ion-ion interaction ==
      do icell2 = 1,NK
        do bion=1,Nion
          if(icell2 == icell .and. bion == aion)cycle
          x = Rion(aion) + dble(icell-1)*lattice_a &
            - (Rion(bion) + dble(icell2-1)*lattice_a) &
            + Uion(icell,aion) - Uion(icell2,bion)
          
          Fion(icell,aion)=Fion(icell,aion) - Zion(aion)*Zion(bion)*( &
        int_pot_drv1_ii(x) + int_pot_drv1_ii(x + la_full) + int_pot_drv1_ii(x - la_full))

!            write(*,"(3I6,e26.16e3)")aion,bion,cion,- 2d0*Zion(bion)*Zion(cion)*( &
!              int_pot_drv1(x) + int_pot_drv1(x + la_full) + int_pot_drv1(x - la_full) &
!              + int_pot_drv1(x + 2d0*la_full) + int_pot_drv1(x - 2d0*la_full) )

        end do
      end do
    end do
  end do

!! zero-force
!  x = sum(Fion)/dble(NK*NB)
!  Fion = Fion -x
  
  return
end subroutine force_ion
