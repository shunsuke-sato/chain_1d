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

  allocate(zH_mat(NK*NB,NK*NB),zH0_mat(NK*NB,NK*NB),zW_mat(NK*NB,NK*NB))
  allocate(zP_mat(NB,NB,NK))
  allocate(zF_mat(NK*NB,NK*NB,Nion),zG_mat(NK*NB,NK*NB,Nion))
  allocate(zF_mat_full(NK*NB,NK*NB,NK,Nion),zG_mat_full(NK*NB,NK*NB,NK,Nion))
  allocate(zG3_mat(NK*NB,NK*NB,Nion),zG3_mat_full(NK*NB,NK*NB,NK,Nion))



! H0
  zH0_mat = 0d0
  is1=0
  do ik1=1,NK; do ib1=1,NB
    is1=is1+1
    zH0_mat(is1,is1)=spe(ib1,ik1)
  end do; end do

! Pmat
  zP_mat = 0d0
  do ik1 = 1,NK
    do ib1 =1,NB
      ztpsi(:) = zpsi_GS(:,ib1,ik1)
      call ppsi(ik1)
      do ib2=ib1,NB
        zs = sum(conjg(zpsi_GS(:,ib2,ik1))*zhtpsi(:))*H
        zP_mat(ib2,ib1,ik1)=conjg(zs)
        zP_mat(ib1,ib2,ik1)=zs
      end do
    end do
  end do

  open(21,file="pmat.out")
  do ik1 = 1,NK
    write(21,"(999e26.16e3)")kx(ik1),abs(zP_mat(1,2,ik1))**2
  end do
  close(21)

! Force matrix

! Start: cell serach
  la_full = lattice_a*NK
  iflag_cell=0
  do aion=1,Nion
    do icell=1,NK
      do ix=1,Nx
        x = Lx(ix) + lattice_a*(icell-1)
        tmp= int_pot(x-Rion(aion)) &
              +int_pot(x-Rion(aion)+la_full) &
              +int_pot(x-Rion(aion)-la_full) &
              +int_pot(x-Rion(aion)+2d0*la_full) &
              +int_pot(x-Rion(aion)-2d0*la_full) 
        tmp = Zion(aion)*tmp
        if(abs(tmp)>1d-20)iflag_cell(icell,aion)=1
             
      end do
    end do
  end do
! End: cell serach

! Fmat
  do aion=1,Nion
    is1=0
    do ik1=1,NK; do ib1=1,NB
      is1=is1+1
      is2=0
      do ik2=1,NK; do ib2=1,NB
        is2=is2+1
        if(is2<is1)cycle

        zs = 0d0
        do icell=1,NK
          if(iflag_cell(icell,aion)==0)cycle
          do ix=1,Nx
            x = Lx(ix) + lattice_a*(icell-1)
            zs = zs + Zion(aion)*exp(zI*(kx(ik2)-kx(ik1))*x) &
              *conjg(zpsi_GS(ix,ib1,ik1))*zpsi_GS(ix,ib2,ik2) &
              *( &
              int_pot_drv1(x-Rion(aion)) &
              +int_pot_drv1(x-Rion(aion)+la_full) &
              +int_pot_drv1(x-Rion(aion)-la_full) &
              +int_pot_drv1(x-Rion(aion)+2d0*la_full) &
              +int_pot_drv1(x-Rion(aion)-2d0*la_full) &
              )
          end do
        end do
        zF_mat(is1,is2,aion)=zs*H/dble(NK)
        zF_mat(is2,is1,aion)=conjg(zF_mat(is1,is2,aion))
        
      end do; end do
    end do; end do
  end do


  do icell = 1,NK
    is1=0
    do ik1=1,NK; do ib1=1,NB
      is1=is1+1
      is2=0
      do ik2=1,NK; do ib2=1,NB
        is2=is2+1
        if(is2>is1)cycle

        zF_mat_full(is1,is2,icell,:) = &
          -exp(zI*(kx(ik2)-kx(ik1))*lattice_a*(icell-1))*zF_mat(is1,is2,:)
        zF_mat_full(is2,is1,icell,:) = conjg(zF_mat_full(is1,is2,icell,:) )

      end do; end do
    end do; end do
  end do

! Gmat
  do aion=1,Nion
    is1=0
    do ik1=1,NK; do ib1=1,NB
      is1=is1+1
!      write(*,*)aion,is1
      is2=0
      do ik2=1,NK; do ib2=1,NB
        is2=is2+1
        if(is2<is1)cycle

        zs = 0d0
        do icell=1,NK
          if(iflag_cell(icell,aion)==0)cycle
          do ix=1,Nx
            x = Lx(ix) + lattice_a*(icell-1)
            zs = zs + Zion(aion)*exp(zI*(kx(ik2)-kx(ik1))*x) &
              *conjg(zpsi_GS(ix,ib1,ik1))*zpsi_GS(ix,ib2,ik2) &
              *( &
              int_pot_drv2(x-Rion(aion)) &
              +int_pot_drv2(x-Rion(aion)+la_full) &
              +int_pot_drv2(x-Rion(aion)-la_full) &
              +int_pot_drv2(x-Rion(aion)+2d0*la_full) &
              +int_pot_drv2(x-Rion(aion)-2d0*la_full) &
              )
          end do
        end do
        zG_mat(is1,is2,aion)=zs*H/dble(NK)
        zG_mat(is2,is1,aion)=conjg(zG_mat(is1,is2,aion))
        
      end do; end do
    end do; end do
  end do

  do icell = 1,NK
    is1=0
    do ik1=1,NK; do ib1=1,NB
      is1=is1+1
      is2=0
      do ik2=1,NK; do ib2=1,NB
        is2=is2+1
        if(is2<is1)cycle

        zG_mat_full(is1,is2,icell,:) = &
          0.5d0*exp(zI*(kx(ik2)-kx(ik1))*lattice_a*(icell-1))*zG_mat(is1,is2,:)
        zG_mat_full(is2,is1,icell,:) = conjg(zG_mat_full(is1,is2,icell,:) )

      end do; end do
    end do; end do
  end do


  return
! G3mat
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
            x = Lx(ix) + lattice_a*(icell-1)
            zs = zs + Zion(aion)*exp(zI*(kx(ik2)-kx(ik1))*x) &
              *conjg(zpsi_GS(ix,ib1,ik1))*zpsi_GS(ix,ib2,ik2) &
              *( &
              int_pot_drv3(x-Rion(aion)) &
              +int_pot_drv3(x-Rion(aion)+la_full) &
              +int_pot_drv3(x-Rion(aion)-la_full) &
              +int_pot_drv3(x-Rion(aion)+2d0*la_full) &
              +int_pot_drv3(x-Rion(aion)-2d0*la_full) &
              )
          end do
        end do
        zG3_mat(is1,is2,aion)=zs*H/dble(NK)
        zG3_mat(is2,is1,aion)=conjg(zG3_mat(is1,is2,aion))
        
      end do; end do
    end do; end do
  end do

  do icell = 1,NK
    is1=0
    do ik1=1,NK; do ib1=1,NB
      is1=is1+1
      is2=0
      do ik2=1,NK; do ib2=1,NB
        is2=is2+1
        if(is2>is1)cycle

        zG3_mat_full(is1,is2,icell,:) = &
          -1d0/(2d0*3d0)*exp(zI*(kx(ik2)-kx(ik1))*lattice_a*(icell-1))*zG3_mat(is1,is2,:)
        zG3_mat_full(is2,is1,icell,:) = conjg(zG3_mat_full(is1,is2,icell,:) )

      end do; end do
    end do; end do
  end do

  return
end subroutine electron_matrix_element
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
