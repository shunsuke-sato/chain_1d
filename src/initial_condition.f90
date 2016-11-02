!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine initial_condition
  use global_variables
  implicit none
  integer :: aion,ik,ib,is,ikt,ibt,ist
  complex(8),allocatable :: zCt(:)
  integer :: Nmat
!LAPACK ==
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
  lwork=6*(NK*NB)**2
  allocate(work_lp(lwork),rwork(3*NK*NB-2),w(NK*NB))
!LAPACK ==
  Nmat = NK*NB
  allocate(zCt(NK*NB))


  allocate(Uion_new(NK,Nion),Uion_old(NK,Nion),Vion(NK,Nion))


! initial conditions for ions ==
  beta_KB = 25.7d-3/27.2d0; beta_KB = 1d0/beta_KB
  call phonon_dist
!  Vion = 0d0
!  Uion = 0d0
!  Uion(1,1)=1d-2
!  Vion(1,1)=1d-1/Mion(1)
!  Vion(1,2)=-1d-1/Mion(2)



! initial condtions for electrons
  call prep_Hmat
!  zH_mat = zH0_mat
!  ist = 0
!  do ikt=1,NK; do ibt=1,NB
!    ist = ist+1 
!    if(ikt == NK/4+1 .and. ibt==1)then
!      is=ist
!      ik=ikt
!      ib=ibt
!    end if
!  end do; end do
!  write(*,*)"in",is,ik,ib
!  zH_mat(is,is+1) =   zH_mat(is,is+1) + sqrt(0.2d0)*(zH0_mat(is+1,is+1)-zH0_mat(is,is))
!  zH_mat(is+1,is) =  conjg(zH_mat(is,is+1))
  Call zheev('V', 'U', Nmat, zH_mat, Nmat, w, work_lp, lwork, rwork, info)
  zpsi_Ct(1:Nmat,1:NK) = zH_mat(1:Nmat,1:NK)
  call prep_Hmat
  call distortion
!  call total_energy
  call force_ion


  do aion = 1,Nion
    Uion_new(:,aion) = Uion(:,aion) + dt*Vion(:,aion) &
      +0.5d0*Fion(:,aion)/Mion(aion)*dt**2
    Uion_old(:,aion) = Uion(:,aion) - dt*Vion(:,aion) &
      +0.5d0*Fion(:,aion)/Mion(aion)*dt**2
  end do

  return
end subroutine initial_condition
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine distortion
  use global_variables
  implicit none
  real(8),parameter :: dist = 0.6d0,sig=0.05d0
  integer :: ik1,ib1,is1,ik2,ib2,is2,iexp
  real(8) :: ss,k0
  complex(8) :: zCt_tmp(NK*NB,NK),zhCt_tmp(NK*NB,NK),zfact
  zW_mat = 0d0

  is1 = 0
  do ik1=1,NK; do ib1=1,NB
    is1 = is1 + 1
    is2 = 0
    do ik2=1,NK; do ib2=1,NB
      is2 = is2 + 1

      if(ik1 /= ik2)cycle
      if((ib1 == 1.and.ib2 ==2))then
        k0 = kx(ik1)+0.5d0*pi/lattice_a
        ss = exp(-0.5d0*(k0/sig)**2)
        k0 = kx(ik1)+0.5d0*pi/lattice_a +2d0*pi/lattice_a 
        ss = ss + exp(-0.5d0*(k0/sig)**2)
        k0 = kx(ik1)+0.5d0*pi/lattice_a -2d0*pi/lattice_a 
        ss = ss + exp(-0.5d0*(k0/sig)**2)
        zW_mat(is1,is2) = zP_mat(ib1,ib2,ik1)*ss
        zW_mat(is2,is1) = conjg(zW_mat(is1,is2))
        
      end if

    end do;end do
  end do; end do

  zCt_tmp(:,:) = zpsi_Ct
  zfact = 1d0
  do iexp =1,16
    zfact = zfact*(-zI*dist)/dble(iexp)
    zhCt_tmp = matmul(zW_mat,zCt_tmp)
    zpsi_Ct = zpsi_Ct + zfact*zhCt_tmp

    zCt_tmp = zhCt_tmp
  end do

  return
end subroutine distortion
