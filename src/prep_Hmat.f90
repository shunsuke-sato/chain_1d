subroutine prep_Hmat
  use global_variables
  implicit none
  integer :: aion,icell

  zH_mat = zH0_mat
  do aion = 1,Nion
    do icell = 1,NK

      zH_mat(:,:) = zH_mat(:,:) - zF_mat_full(:,:,icell,aion)*Uion(icell,aion) &
        -zG_mat_full(:,:,icell,aion)*Uion(icell,aion)**2 ! &
!        - zG3_mat_full(:,:,icell,aion)*Uion(icell,aion)**3


    end do
  end do
  
  return
end subroutine prep_Hmat
