!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine ground_state
  use global_variables
  implicit none
  complex(8) :: za(Nx,Nx),zv(Nx,Nx)
  real(8) :: ar(Nx,Nx),ai(Nx,Nx),vr(Nx,Nx),vi(Nx,Nx)
  real(8) :: e(Nx)
  integer :: ix,jx,jx2,ik,ib
  real(8) :: ce,ve,s
!LAPACK ==
  integer :: lwork
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
  lwork=6*Nx**2
  allocate(work_lp(lwork),rwork(3*Nx-2),w(Nx))
!LAPACK ==

  write(*,"(A)")"!Start Electronic ground state calculation"
  Acx = 0d0

  do ik=1,NK
    za=0d0
    do ix=1,Nx
      do jx=ix-2,ix+2
        jx2=jx
        if(jx<=0)jx2=Nx+jx
        if(jx>Nx)jx2=jx-Nx
        za(ix,jx2)=-0.5d0*(L_coef(jx-ix)+zI*2d0*Kx(ik)*G_coef(jx-ix))
      end do
      za(ix,ix)=za(ix,ix)+0.5d0*kx(ik)**2+Veff(ix)
    end do

    Call zheev('V', 'U', Nx, za, Nx, w, work_lp, lwork, rwork, info)
    do ib=1,NB
      s=sum(abs(za(:,ib))**2)*H
      zpsi_GS(:,ib,ik)=za(:,ib)/sqrt(s)
    end do

    do ib=1,NB
      spe(ib,ik)=w(ib)
    end do
    
  end do

  open(10,file=trim(filename)//'_gswf.out')
  do ix=1,Nx
    write(10,'(100e26.16E3)')Lx(ix),abs(zpsi_GS(ix,1,1))**2,abs(zpsi_GS(ix,1,Nk/2))**2 &
      ,abs(zpsi_GS(ix,2,1))**2,abs(zpsi_GS(ix,2,Nk/2))**2
  end do
  close(10)
  
  
  open(10,file=trim(filename)//'_band_map.out')
  do ik=1,NK
    write(10,'(100e26.16E3)')kx(ik),(spe(ib,ik),ib=1,NB)
  end do
    write(10,'(100e26.16E3)')-kx(1),(spe(ib,1),ib=1,NB)
  close(10)
  
  open(10,file=trim(filename)//'_pot.out')
  do ix=1,Nx
    write(10,'(100e26.16E3)')Lx(ix),Veff(ix)
  end do
  close(10)
  
  write(*,*)'Band gap       =',spe(2,NK)-spe(1,NK),(spe(2,NK)-spe(1,NK))*Ry*2d0
  
  write(*,"(A,2x,e26.16e3)")'Total energy =',2d0*sum(spe(1,:))/dble(NK)+E_ii
  write(*,"(A,2x,e26.16e3)")'Electronic energy =',2d0*sum(spe(1,:))/dble(NK)
  write(*,"(A,2x,e26.16e3)")'Ion-Ionc energy =',E_ii
  write(*,"(A,2x,e26.16e3)")'Band-gap (eV)=',2d0*Ry*(spe(2,1) - spe(1,1))

  write(*,"(A)")"!End Electronic ground state calculation"  
!  zpsi(:,1:NBocc,:)=zpsi_GS(:,1:NBocc,:)
!  stop  
  return
end subroutine ground_state
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
