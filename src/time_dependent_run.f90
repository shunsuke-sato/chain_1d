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
  call energy_average_within_multi_trajectory
  E_tot0=E_tot; E_elec0=E_elec; E_ii0=E_ii; Ekin_ion0=Ekin_ion


  if(myrank == 0)then
    open(41,file="energy_td.out")
    write(41,"(999e26.16e3)")0d0,E_tot,E_elec,E_ii,Ekin_ion &
      ,E_tot-E_tot0,E_elec-E_elec0,E_ii-E_ii0,Ekin_ion-Ekin_ion0 &
      ,sum(abs(zpsi_Ct)**2)/dble(NK*NB)
    open(42,file="energy_td_MT.out")
    write(42,"(999e26.16e3)")0d0,E_tot_MT,E_elec_MT,E_ii_MT,Ekin_ion_MT
  end if

  it = 0
  call write_neh_dist(it)

!  do it = 1,Nt
  it = 0
  do 
    it = it + 1
    if(myrank == 0 .and. mod(it,10)==0)write(*,*)"it=",it,Nt,sum(abs(zpsi_Ct)**2)/dble(NK)
    call dt_evolve
    if(mod(it,20) == 0)then
      call total_energy_TD
      call energy_average_within_multi_trajectory
      if(myrank == 0)then
        write(41,"(999e26.16e3)")dt*it,E_tot,E_elec,E_ii,Ekin_ion &
          ,E_tot-E_tot0,E_elec-E_elec0,E_ii-E_ii0,Ekin_ion-Ekin_ion0 &
          ,sum(abs(zpsi_Ct)**2)/dble(NK*NB)
        write(42,"(999e26.16e3)")dt*it,E_tot_MT,E_elec_MT,E_ii_MT,Ekin_ion_MT
      end if
    end if

    if(mod(it,100) == 0) call write_neh_dist(it)

    if(dble(it)*dt*0.02418d0 > 500d0)exit
  end do
  close(41)

  contains
    subroutine energy_average_within_multi_trajectory
      call MPI_ALLREDUCE(E_tot,E_tot_MT,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(E_elec,E_elec_MT,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Ekin_ion,Ekin_ion_MT,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(E_ii,E_ii_MT,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      E_tot_MT = E_tot_MT/dble(Nprocs)
      E_elec_MT = E_elec_MT/dble(Nprocs)
      E_ii_MT = E_ii_MT/dble(Nprocs)
      Ekin_ion_MT = Ekin_ion_MT/dble(Nprocs)
    end subroutine energy_average_within_multi_trajectory

end subroutine time_dependent_run
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
