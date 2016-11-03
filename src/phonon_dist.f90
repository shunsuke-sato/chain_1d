subroutine phonon_dist
  use global_variables
  implicit none
  integer :: is1,is2, aion,bion, icell, icell2, irank
  real(8) :: ss,x1,x2,tmp
!LAPACK ==
  integer :: lwork,Nmat
  real(8),allocatable :: work_lp(:),Amat(:,:)
  real(8),allocatable :: rwork(:),w(:),w2(:),w3(:)
  integer :: info
  lwork=6*(NK*Nion)**2
  Nmat =NK*Nion
  allocate(work_lp(lwork),rwork(3*Nmat-2),w(Nmat),w2(Nmat),Amat(Nmat,Nmat))
!LAPACK ==


!  Phi_FC(:,:,1,aion) = - Fion(:,:) /Udist
  is1 =0
  do aion=1,Nion; do icell =1,NK
    is1 = is1 + 1
    is2 = 0
    do bion=1,Nion; do icell2 =1,NK
      is2 = is2 + 1
      Dph_mat(is1,is2)=Phi_FC(icell2,bion,icell,aion)/sqrt(Mion(aion)*Mion(bion))
    end do; end do
  end do; end do

  Amat = Dph_mat
!  Call zheev('V', 'U', Nmat, zAmat, Nmat, w, work_lp, lwork, rwork, info)
  Call dsyev('V', 'U', Nmat, Amat, Nmat, w, work_lp, lwork, info)

  
  ss = 0d0
  is1 =0
  do aion=1,Nion; do icell =1,NK
    is1 = is1 + 1
    w2(is1) = w2_ph(aion,icell)
  end do; end do

  do is1 = 1,Nmat
    do is2 = 1,Nmat-1
      if(w2(is2) > w2(is2+1))then
        ss = w2(is2)
        w2(is2) = w2(is2+1)
        w2(is2+1) = ss
      end if
    end do
  end do

  Uion =0d0
  Vion = 0d0
  do irank = 0,Nprocs-1
    do is1 =2,NK*Nion

      call normal_random_number(x1,x2)
      call random_number(x2)
      if(irank /= myrank)cycle
!    ss = x1/sqrt(beta_KB*w(is1)) !Boltzmann distribution
      tmp = 2d0*tanh(beta_KB*sqrt(w(is1))/2d0)*sqrt(w(is1)) !Quantum distribution
      ss = x1/sqrt(tmp) !Quantum distribution

      is2 = 0
      do aion = 1,Nion
        do icell = 1,NK
          is2 = is2 + 1
          Uion(icell,aion) = Uion(icell,aion) &
            + ss*sin(2d0*pi*x2)/sqrt(Mion(aion))*Amat(is2,is1)
          Vion(icell,aion) = Vion(icell,aion) &
            + ss*sqrt(w(is1)/Mion(aion))*cos(2d0*pi*x2)*Amat(is2,is1)
        end do
      end do
    end do
  end do
  
  write(*,*)"rank,U",myrank,sqrt(sum(Uion**2)/dble(NK*Nion))
  write(*,*)"rank,V",myrank,sqrt(sum(Vion**2)/dble(NK*Nion))
  write(*,*)"rank,U-max",myrank,maxval(abs(Uion))
  write(*,*)"rank,V-max",myrank,maxval(abs(Vion))

!  stop
  return
end subroutine phonon_dist

subroutine normal_random_number(x1,x2)
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  real(8) :: r1,r2,tmp

  call random_number(r1)
  call random_number(r2)

  if(r1 == 0d0)then
    x1 = 0d0
    x2 = 0d0
  else 
    tmp = sqrt(-2d0*log(r1))
    x1 = tmp*cos(2d0*pi*r2)
    x2 = tmp*sin(2d0*pi*r2)
  end if
  return
end subroutine normal_random_number
