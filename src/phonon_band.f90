!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine phonon_band
  use global_variables
  implicit none
  real(8),parameter :: Udist = 0.10d0
  complex(8),allocatable :: zWp_mat(:,:), zF_mat_t(:,:), zCt(:)
  real(8) :: x
  integer :: ik1,ik2,ib1,ib2,is1,is2
  integer :: aion,bion,cion,icell,icell2,ix
  complex(8) :: zs
  real(8) :: la_full
  integer :: Nmat
  integer :: iflag_cell(NK,Nion)
  real(8) :: tmp,ss
  real(8) :: Eii_tmp
!LAPACK ==
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
  lwork=6*(NK*NB)**2
  allocate(work_lp(lwork),rwork(3*NK*NB-2),w(NK*NB))
!LAPACK ==

  write(*,"(A)")"!Start Phonon band calculation"

  Nmat = NK*NB
  allocate(Phi_FC(NK,Nion,NK,Nion), Fion(NK,Nion))
  allocate(zWp_mat(NK*NB,NK*NB),zCt(NK*NB), zF_mat_t(NK*NB,NK*NB))
  la_full = lattice_a*NK

! Start: cell serach
  iflag_cell=0
  do aion=1,Nion
    do icell=1,NK
      do ix=1,Nx
        x = Lx(ix) + lattice_a*(icell-1)
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

      is2=0
      do ik2=1,NK; do ib2=1,NB
        is2=is2+1
        if(is2<is1)cycle

        zs = 0d0
        do icell=1,NK
          if(iflag_cell(icell,aion)==0)cycle
          do ix=1,Nx
            x = Lx(ix) + lattice_a*(icell-1)
            zs = zs - exp(zI*(kx(ik2)-kx(ik1))*x) &
              *conjg(zpsi_GS(ix,ib1,ik1))*zpsi_GS(ix,ib2,ik2) &
              *Zion(aion)*( &
              int_pot(x-(Rion(aion)+Udist))-int_pot(x-Rion(aion)) &
              +int_pot(x-(Rion(aion)+Udist)+la_full)-int_pot(x-Rion(aion)+la_full) &
              +int_pot(x-(Rion(aion)+Udist)-la_full)-int_pot(x-Rion(aion)-la_full) &
              )
          end do
        end do
        zWp_mat(is1,is2)=zs*H/dble(NK)
        zWp_mat(is2,is1)=conjg(zWp_mat(is1,is2))
        
      end do; end do
    end do; end do


    zH_mat(:,:)=zH0_mat(:,:) + zWp_mat(:,:)
    Call zheev('V', 'U', Nmat, zH_mat, Nmat, w, work_lp, lwork, rwork, info)

! Force
    Fion = 0d0
    do bion = 1,Nion
      do icell=1,NK
        if(aion == bion .and. icell == 1)cycle

! electron-ion interaction ==
        is1=0
        do ik1=1,NK; do ib1=1,NB
          is1=is1+1
          is2=0
          do ik2=1,NK; do ib2=1,NB
            is2=is2+1
            if(is2>is1)cycle

            zF_mat_t(is1,is2) = exp(zI*(kx(ik2)-kx(ik1))*lattice_a*(icell-1))*zF_mat(is1,is2,bion)
            zF_mat_t(is2,is1) = conjg(zF_mat_t(is1,is2) )

          end do; end do
        end do; end do


        zs = 0d0
        do is1=1,NBocc*NK
          zCt= matmul(zF_mat_t,zH_mat(:,is1))
          ss = sum(abs(zH_mat(:,is1))**2) 
          zs = zs + sum(conjg(zH_mat(:,is1))*zCt(:))/ss
        end do
        Fion(icell,bion)=-2d0*real(zs)
! ion-ion interaction ==
        do icell2 = 1,NK
          do cion=1,Nion
            if(icell2 == icell .and. bion == cion)cycle
            x = Rion(bion) + (icell-1)*lattice_a - (Rion(cion) + (icell2-1)*lattice_a)
            if(icell2 == 1 .and. cion == aion) x = x - Udist
            Fion(icell,bion)=Fion(icell,bion) - 2d0*Zion(bion)*Zion(cion)*( &
              int_pot_drv1(x) + int_pot_drv1(x + la_full) + int_pot_drv1(x - la_full) )

          end do
        end do

      end do
    end do
    ss = sum(Fion)
    Fion(1,aion) = - ss
    Phi_FC(:,:,1,aion) = Fion(:,:) /Udist

! Eii calc
    Eii_tmp=0d0
    do icell = 1,NK; do bion=1,Nion
    do icell2 = 1,NK; do cion=1,Nion
      if(icell2 == icell .and. bion == cion)cycle
      x = Rion(bion) + (icell-1)*lattice_a - (Rion(cion) + (icell2-1)*lattice_a)
      if(icell2 == 1 .and. cion == aion) x = x - Udist
      if(icell == 1 .and. bion == aion) x = x + Udist
      Eii_tmp = Eii_tmp + Zion(bion)*Zion(cion)*( &
        int_pot(x) + int_pot(x + la_full) + int_pot(x - la_full) )
    end do; end do
    end do; end do

    write(*,"(A,2x,999e26.16e3)")"Etot,Eii",sum(w(1:NK))*2d0+Eii_tmp,Eii_tmp
    write(*,"(A,2x,999e26.16e3)")"Force,Udist",Fion(1,aion),Udist
    write(*,"(A,2x,999e26.16e3)")"res.",Udist,sum(w(1:NK))*2d0+Eii_tmp,Fion(1,aion)
!    write(*,"(A,2x,999e26.16e3)")"res2.",Udist,Eii_tmp,Fion(1,aion)
  end do


  write(*,*)sum(abs(zF_mat(:,:,:)))
  open(20,file="force.dat")
  do icell  = 1,NK
    write(20,"(I7,2x,999e26.16e3)")icell,Phi_FC(icell,:,1,1),Phi_FC(icell,:,1,2)
  end do
  close(20)
  write(*,"(A)")"!End Phonon band calculation"
  
  return
end subroutine phonon_band
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
