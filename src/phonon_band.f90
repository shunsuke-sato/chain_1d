!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine phonon_band
  use global_variables
  implicit none
  real(8),parameter :: Udist = 1d-4
  complex(8),allocatable :: zWp_mat(:,:)
  real(8) :: x
  integer :: ik1,ik2,ib1,ib2,is1,is2
  integer :: aion,bion,icell,ix
  complex(8) :: zs
  real(8) :: la_full
  integer :: Nmat
  integer :: iflag_cell(NK,Nion)
  real(8) :: tmp
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
  allocate(Phi_FC(NK,Nion,NK,Nion))
  allocate(zWp_mat(NK*NB,NK*NB))
  la_full = lattice_a*NK

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
            zs = zs - exp(zI*(kx(ik2)-kx(ik1))*x) &
              *conjg(zpsi_GS(ix,ib1,ik1))*zpsi_GS(ix,ib2,ik2) &
              *Zion(aion)*( &
              int_pot(x-Rion(aion)+Udist)-int_pot(x-Rion(aion)) &
              +int_pot(x-Rion(aion)+Udist+la_full)-int_pot(x-Rion(aion)+la_full) &
              +int_pot(x-Rion(aion)+Udist-la_full)-int_pot(x-Rion(aion)-la_full) &
              )
          end do
        end do
        zWp_mat(is1,is2)=zs*H/dble(NK)
        zWp_mat(is2,is1)=conjg(zWp_mat(is1,is2))
        
      end do; end do
    end do; end do


    zH_mat(:,:)=zH0_mat(:,:) + zWp_mat(:,:)
    write(*,*)"diag"
    Call zheev('V', 'U', Nmat, zH_mat, Nmat, w, work_lp, lwork, rwork, info)
    write(*,*)"diag-end"
  end do
  
  write(*,"(A)")"!End Phonon band calculation"
  
  return
end subroutine phonon_band
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
