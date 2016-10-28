!---------------------------------------------------!
! Copyright (c) 2016 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine time_dependent_run
  use global_variables
  implicit none
  integer :: it,ist

  call initial_condition
  call total_energy_TD
  E_tot0=E_tot; E_elec0=E_elec; E_ii0=E_ii; Ekin_ion0=Ekin_ion
  open(41,file="energy_td.out")
  write(41,"(999e26.16e3)")0d0,E_tot,E_elec,E_ii,Ekin_ion &
    ,E_tot-E_tot0,E_elec-E_elec0,E_ii-E_ii0,Ekin_ion-Ekin_ion0 &
    ,sum(abs(zpsi_Ct)**2)/dble(NK*NB)
  open(42,file="pex_dist_td.out")
  write(42,"(999e26.16e3)")0d0,(sum(abs(zpsi_Ct(ist,:))**2),ist=1,NB*NK)
  it = 0
  call write_neh_dist(it)

!  do it = 1,Nt
  it = 0
  do 
    it = it + 1
    write(*,*)"it=",it,Nt,sum(abs(zpsi_Ct)**2)/dble(NK)
    call dt_evolve
    if(mod(it,10) == 0)then
      call total_energy_TD
      write(41,"(999e26.16e3)")dt*it,E_tot,E_elec,E_ii,Ekin_ion &
        ,E_tot-E_tot0,E_elec-E_elec0,E_ii-E_ii0,Ekin_ion-Ekin_ion0 &
        ,sum(abs(zpsi_Ct)**2)/dble(NK*NB)
      write(42,"(999e26.16e3)")dt*it,(sum(abs(zpsi_Ct(ist,:))**2),ist=1,NB*NK)
    end if

    if(mod(it,10) == 0) call write_neh_dist(it)

    if(dble(it)*dt*0.02418d0 > 50d0)exit
  end do
  close(41)

end subroutine time_dependent_run
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
