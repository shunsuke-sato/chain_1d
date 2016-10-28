!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
  subroutine write_neh_dist(it)
    use global_variables
    implicit none
    integer :: it,ik,ib,ist
    character(50) :: cit,filename_t
    real(8) :: neh(NB,NK),n_tot

    write(cit,'(I7.7)')it
    filename_t = "neh_dist_"//trim(cit)//".out"

    ist = 0
    do ik=1,NK
      do ib=1,NB
        ist = ist + 1
        neh(ib,ik) = sum(abs(zpsi_Ct(ist,:))**2)
      end do
    end do

    n_tot = sum(neh)
    neh(1,:) = 1d0 - neh(1,:)
    open(101,file=filename_t)
    write(101,"(A,2x,999e26.16e3)")"n_tot/2,hole [%],ele [%],",sum(neh(1,:))/n_tot*1d2 &
      ,sum(neh(2:NB,:))/n_tot*1d2
    do ik = 1,NK
      write(101,"(9999e26.16e3)")kx(ik),neh(:,ik)
    end do
    close(101)

    return
  end subroutine write_neh_dist
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
