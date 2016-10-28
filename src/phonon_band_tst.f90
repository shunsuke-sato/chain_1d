!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine phonon_band_tst
  use global_variables
  implicit none
  real(8),parameter :: Udist = 1d-2
  integer,parameter :: Num=10
  complex(8),allocatable :: zWp_mat(:,:), zF_mat_t(:,:), zCt(:)
  complex(8),allocatable :: zDw(:,:,:)
  real(8) :: x
  integer :: ik1,ik2,ib1,ib2,is1,is2
  integer :: aion,bion,cion,icell,icell2,icell3,ix
  complex(8) :: zs
  real(8) :: la_full
  integer :: Nmat
  integer :: iflag_cell(NK,Nion)
  real(8) :: tmp,ss,a,b2,c
  real(8) :: Eii_tmp,E_gs
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
  allocate(zCt(NK*NB), zDw(Nion,Nion,NK))
  la_full = lattice_a*dble(NK)

! GS
  zpsi_Ct = 0d0
  Uion = 0d0

  call prep_Hmat
  Call zheev('V', 'U', Nmat, zH_mat, Nmat, w, work_lp, lwork, rwork, info)
  zpsi_Ct(1:Nmat,1:NK) = zH_mat(1:Nmat,1:NK)
  call prep_Hmat
  call total_energy
  call force_ion

  E_gs = E_tot

  open(21,file="ef_tst.out")
  do is1 = 0,Num
    Uion(1,1) = dble(is1)*Udist/dble(Num)

    call prep_Hmat
    Call zheev('V', 'U', Nmat, zH_mat, Nmat, w, work_lp, lwork, rwork, info)
    zpsi_Ct(1:Nmat,1:NK) = zH_mat(1:Nmat,1:NK)
    call prep_Hmat
    call total_energy
    call force_ion

    write(21,"(999e26.16e3)")Uion(1,1),E_tot-E_gs,Fion(1,1),E_tot,E_elec,E_ii
  end do
  close(21)

  write(*,*)"f-sum",sum(Fion(:,:)),sum(Fion(:,1)),sum(Fion(:,2))
  do is1 = 1,NK
    write(*,"(I6,2x,999e16.6e3)")is1,Fion(is1,1),Fion(is1,2)
  end do


! phonon calculation
  Phi_FC(:,:,:,:) = 0d0
  do aion = 1,Nion
    Uion(:,:) = 0d0
    Uion(1,aion) = Udist
    call prep_Hmat
    Call zheev('V', 'U', Nmat, zH_mat, Nmat, w, work_lp, lwork, rwork, info)
    zpsi_Ct(1:Nmat,1:NK) = zH_mat(1:Nmat,1:NK)
    call prep_Hmat
    call total_energy
    call force_ion

    Phi_FC(:,:,1,aion) = - Fion(:,:) /Udist
  end do


  do icell = 2,NK
    do icell2 = 1,NK
      icell3 = icell2 + (icell-1)
      if(icell3 > NK)icell3 = icell3-NK
      Phi_FC(icell3,:,icell,:) = Phi_FC(icell2,:,1,:) 
    end do
  end do

  do ik1 = 1,NK
    do aion = 1,Nion
      do bion = 1,Nion
        zs = 0d0
        do icell = 1,NK
          zs = zs + Phi_FC(1,aion,icell,bion) &
          *exp(-zI*kx(ik1)*(Rion(aion) - Rion(bion) - lattice_a*(icell-1)))
        end do
        zDw(aion,bion,ik1)=zs/sqrt(Mion(aion)*Mion(bion))
      end do
    end do
  end do
! phonon band; assuming Nion = 2

  do ik1 = 1,NK
    a = real(zDw(1,1,ik1)); c = real(zDw(2,2,ik1)); b2 = abs(zDw(1,2,ik1))**2
    w2_ph(1,ik1) = 0.5d0*(a+c-sqrt((a-c)**2+4d0*b2))
    w2_ph(2,ik1) = 0.5d0*(a+c+sqrt((a-c)**2+4d0*b2))
    if(w2_ph(1,ik1) < 0d0) w2_ph(1,ik1) = 0d0
    if(w2_ph(2,ik1) < 0d0) w2_ph(2,ik1) = 0d0
  end do
  w_ph = sqrt(w2_ph)

  open(10,file=trim(filename)//'_phonon_band_tst.out')
  do ik1=1,NK
    write(10,'(100e26.16E3)')kx(ik1),w_ph(:,ik1),w2_ph(:,ik1)
  end do
    write(10,'(100e26.16E3)')-kx(1),w_ph(:,1),w2_ph(:,1)
  close(10)



  return
end subroutine phonon_band_tst
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
