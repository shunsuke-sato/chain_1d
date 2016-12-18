!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine dt_evolve
  use global_variables
  implicit none
  complex(8),parameter :: zONE = (1d0,0d0),zZero = 0d0
  integer,parameter :: NTaylor = 16
  integer :: Nmax,iTaylor,aion,ist
  complex(8) :: zfact
  Nmax = NB*NK

  zW_mat = zH_mat - zH0_mat
!  zW_mat = zH_mat
!  call dt_half_evolve_Taylor
  call dt_half_evolve_Lanczos

  do ist = 1,Nmax
    zpsi_Ct(ist,:) = exp(-zI*zH0_mat(ist,ist)*dt)*zpsi_Ct(ist,:)
  end do

  Uion_old = Uion
  Uion = Uion_new
  call prep_Hmat

  zW_mat = zH_mat - zH0_mat
!  zW_mat = zH_mat
!  call dt_half_evolve_Taylor
  call dt_half_evolve_Lanczos

  call force_ion
  do aion = 1,Nion
    Uion_new(:,aion) = 2d0*Uion(:,aion) -Uion_old(:,aion)+Fion(:,aion)/Mion(aion)*dt**2
  end do

  return

  contains
    subroutine dt_half_evolve_Taylor
      implicit none

      zpsi_Ct_t = zpsi_Ct
      zfact = 1d0
      do itaylor = 1,NTaylor
        zfact=zfact*(-zI*0.5d0*dt)/dble(itaylor)
        call zhemm('L','U',Nmax,NK,zONE,zW_mat,Nmax,zpsi_Ct_t,Nmax,zZero,zhpsi_Ct_t,Nmax)

        zpsi_Ct = zpsi_Ct + zfact*zhpsi_Ct_t 
        zpsi_Ct_t = zhpsi_Ct_t
      end do

      return
    end subroutine dt_half_evolve_Taylor

    subroutine dt_half_evolve_Lanczos
      implicit none
      integer :: iLan, ik
      real(8) :: alpha(NLanczos,NK),beta(NLanczos,NK),Hmat_Lan(NLanczos,NLanczos)
      complex(8) :: zvec(Nlanczos)
      real(8) :: ss
!LAPACK ==
      integer :: lwork,Nmat
      real(8),allocatable :: work_lp(:),Amat(:,:)
      real(8),allocatable :: rwork(:),w(:)
      integer :: info
      Nmat = Nlanczos
      lwork=6*(Nmat)**2
      allocate(work_lp(lwork),rwork(3*Nmat-2),w(Nmat),Amat(Nmat,Nmat))
!LAPACK ==

      do ik = 1,NK
        ss = sum(abs(zpsi_Ct(:,ik))**2)
        zpsi_Ct_t_Lan(:,ik,1) = zpsi_Ct(:,ik)/sqrt(ss)
      end do
      beta(1,:) = 0d0

      do iLan = 1,NLanczos-1
        zpsi_Ct_t(:,:) = zpsi_Ct_t_Lan(:,:,iLan)
        call zhemm('L','U',Nmax,NK,zONE,zW_mat,Nmax,zpsi_Ct_t,Nmax,zZero,zhpsi_Ct_t,Nmax)
        do ik = 1,NK
          alpha(iLan,ik) = sum(conjg(zpsi_Ct_t_Lan(:,ik,iLan))*zhpsi_Ct_t(:,ik))
        end do

        do ik = 1,NK
          zpsi_Ct_t_Lan(:,ik,iLan+1) = zhpsi_Ct_t(:,ik) & 
            - alpha(iLan,ik)*zpsi_Ct_t_Lan(:,ik,iLan) &
            - beta(iLan,ik)*zpsi_Ct_t_Lan(:,ik,iLan-1) 
        end do


        do ik = 1,NK
          beta(iLan+1,ik) = sqrt(sum(abs( zpsi_Ct_t_Lan(:,ik,iLan+1) )**2))
          zpsi_Ct_t_Lan(:,ik,iLan+1) = zpsi_Ct_t_Lan(:,ik,iLan+1)/beta(iLan+1,ik)
        end do
      end do
      zpsi_Ct_t(:,:) = zpsi_Ct_t_Lan(:,:,Nlanczos)
      call zhemm('L','U',Nmax,NK,zONE,zW_mat,Nmax,zpsi_Ct_t,Nmax,zZero,zhpsi_Ct_t,Nmax)
      do ik = 1,NK
        alpha(Nlanczos,ik) = sum(conjg(zpsi_Ct_t_Lan(:,ik,Nlanczos))*zhpsi_Ct_t(:,ik))
      end do

      do ik = 1,NK
        Amat = 0d0
        do iLan = 1,Nlanczos
          Amat(iLan,iLan) = alpha(iLan,ik)
        end do
        do iLan = 2,Nlanczos
          Amat(iLan-1,iLan) = beta(iLan,ik)
          Amat(iLan,iLan-1) = beta(iLan,ik)
        end do
        Call dsyev('V', 'U', Nmat, Amat, Nmat, w, work_lp, lwork, info)
        zvec = 0d0; zvec(1) = 1d0
        zvec = matmul(transpose(Amat),zvec)
        zvec(:) = exp(-zI*0.5d0*dt*w(:))*zvec(:)
        zvec = matmul(Amat,zvec)
        zpsi_Ct(:,ik) = matmul(zpsi_Ct_t_Lan(:,ik,1:Nlanczos),zvec(1:Nlanczos))
      end do
      return
    end subroutine dt_half_evolve_Lanczos
end subroutine dt_evolve

!-------10--------20--------30--------40--------50--------60--------70--------80--------90
