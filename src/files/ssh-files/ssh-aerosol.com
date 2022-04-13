!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

PROGRAM SSHaerosol

  use aInitialization
  use jAdaptstep
  use bCoefficientRepartition
  use mEmissions
  use netcdf
  use lDiscretization
  use Resultoutput
  use gCoagulation
  use mod_photolysis
  use mod_meteo
  
  implicit none

  integer :: t, j, s,jesp,day,stat,hour
  logical :: file_exists
  character (len=400) :: namelist_ssh  ! Configuration file
  double precision, dimension(:), allocatable :: timer,timer2 !zhizhao
  double precision, dimension(:), allocatable :: cinorg_use !zhizhao

  double precision :: t_since_update_photolysis, t0,t1, delta_t2

  ! Initial time (seconds)
  call cpu_time(t0)

  ! Initialisation: discretization and distribution
  if (iargc() == 0) then
     write(*,*) "usage: ssh-aerosol namelist.input"
     stop
  else
     call getarg(1, namelist_ssh)
  end if

  ! Read the number of gas-phase species and chemical reactions
  call ssh_dimensions(N_gas, n_reaction, n_photolysis)

  call ssh_read_namelist(namelist_ssh)

  call ssh_read_inputs()                                

  call ssh_read_meteo()
  
  ! Read the meteorological data.
  temperature = temperature_array(1)
  pressure = pressure_array(1)
  humidity = humidity_array(1)
  relative_humidity = relative_humidity_array(1) 

  call ssh_init_parameters()
  
  call ssh_init_coag()  

  call ssh_init_distributions()  

  if (output_type .ne. 0) call ssh_init_output_conc() 

  call ssh_save_report()

  if (output_type .ne. 0) call ssh_save_concentration() 

  if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
    ! Allocate arrays for photolysis
    call ssh_allocate_photolysis()    
    ! Read photolysis rate for current day
    current_time = initial_time
    call ssh_init_photolysis() 
    call ssh_interpol_photolysis()
  endif
  ! **** simulation starts 
  t_since_update_photolysis = 0.d0

  ! Initialization is finished
  allocate(timer(nt+3))
  allocate(timer2(nt*2)) !zhizhao
  !allocate(cinorg_use(size(cstindex))) !zhizhao
  allocate(cinorg_use(ncst_chem))!zhizhao

  timer(1) = t0
  call cpu_time(t0)
  timer(2) = t0

  do t = 1, nt

     current_time = initial_time + (t - 1) * delta_t
     t_since_update_photolysis = t_since_update_photolysis +  delta_t

     if (ssh_standalone) write(*,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))
     if (ssh_logger) write(logfile,*) "Performing iteration #" // trim(str(t)) // "/" // trim(str(nt))

     ! Read the photolysis rates.
     if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
       if (t_since_update_photolysis >= time_update_photolysis) then
        call ssh_interpol_photolysis()
        t_since_update_photolysis = 0.d0
       endif
     endif

     ! Emissions
     if (tag_emis .ne. 0) call ssh_emission(delta_t)

     ! Read the meteorological data.
     temperature = temperature_array(t)
     pressure = pressure_array(t)
     humidity = humidity_array(t)
     relative_humidity = relative_humidity_array(t) 

     ! Gas-phase chemistry

     ! initial and final physical parameters are set same
     ! no volumetric emission
     ! 0 : vertical gas volumetric emission    1 : with number 
     ! 0 : not take into account cloud    0.d0 : air water content fracion sets to 0

     call cpu_time(t0) !zhizhao
     if (tag_chem .ne. 0) then
       ! get hourly radical conc.
       ! hour
       hour = int(mod(current_time/3600.,24.)) !int((current_time/3600.)%24)
       !print*,'time',current_time,'current hour', hour, cinorg(:,hour+1)

       if (ncst_chem .gt. 0) then
          cinorg_use = 0.0
          do s =1, ncst_chem !size(cstindex)
             cinorg_use(s) = cinorg(s,hour+1)
             !if (s.eq.1) print*,s,hour,cinorg_use(s)
          enddo
       endif

       if (t==1) then
          !delta_t=1=0.001*delta_t
          call ssh_chem(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
               ns_source, source_index, conversionfactor, conversionfactorjacobian,&
               0, lwc_cloud_threshold, molecular_weight, &
               current_time, attenuation, &
               humidity, temperature,&
               pressure, source, &
               photolysis_rate, 0.001*delta_t, attenuation,&
               humidity, temperature,&
               pressure, source, &
               photolysis_rate, longitude,&
               latitude, concentration_gas_all,&
               0, with_heterogeneous, n_aerosol, n_size, n_fracmax,&
               0.d0,&
               diam_bound, fixed_density, &
               wet_diameter, &
               heterogeneous_reaction_index, &
               concentration_mass,&
               with_adaptive, adaptive_time_step_tolerance,&
               min_adaptive_time_step, option_photolysis, ind_jbiper, ind_kbiper,&
               1, not(with_fixed_density), concentration_number, &
               mass_density, &
               ncst_chem, cinorg_use, cstindex, & !zhizhao use constant gas conc.
               tag_RO2, nRO2_chem, iRO2, iRO2_cst, RO2index, &   !zhizhao treatment of RO2
               iSumM, 0, aerosol_species_interact(:) )

          ! re-calculate total_mass(N_aerosol) because mass change due to gas-phase chemistry
          total_aero_mass = 0.d0
          total_mass = 0.d0
          do s = 1, N_aerosol_layers
             jesp = List_species(s)
             do j=1,N_size
                total_aero_mass(jesp) = total_aero_mass(jesp) + concentration_mass(j,s)
             enddo
          enddo
          ! update mass conc. of aerosol precursors
          ! concentration_gas_all(precursor_index) -> concentration_gas(n_aerosol)
          do s = 1, N_aerosol
             if (aerosol_species_interact(s) .gt. 0) then
                concentration_gas(s) = concentration_gas_all(aerosol_species_interact(s))
             end if
             total_mass(s) = total_mass(s) + concentration_gas(s) + total_aero_mass(s)
          end do

          ! Aerosol dynamic
          CALL SSH_AERODYN(current_time,0.001*delta_t)

          ! update mass conc. of aerosol precursors
          ! concentration_gas(n_aerosol) -> concentration_gas_all(precursor_index)
          do s = 1, N_aerosol
             if (aerosol_species_interact(s) .gt. 0) then
                concentration_gas_all(aerosol_species_interact(s)) = concentration_gas(s)
             end if
          end do
          delta_t2=0.999*delta_t
       else
          delta_t2=delta_t
       endif

       call ssh_chem(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
          ns_source, source_index, conversionfactor, conversionfactorjacobian,&
          0, lwc_cloud_threshold, molecular_weight, &
          current_time, attenuation, &
          humidity, temperature,&
          pressure, source, &
          photolysis_rate, delta_t2, attenuation,&
          humidity, temperature,&
          pressure, source, &
          photolysis_rate, longitude,&
          latitude, concentration_gas_all,&
          0, with_heterogeneous, n_aerosol, n_size, n_fracmax,&
          0.d0,&
          diam_bound, fixed_density, &
          wet_diameter, &
          heterogeneous_reaction_index, &
          concentration_mass,&
          with_adaptive, adaptive_time_step_tolerance,&
          min_adaptive_time_step, option_photolysis, ind_jbiper, ind_kbiper,&
          1, not(with_fixed_density), concentration_number, &
          mass_density, &
          ncst_chem, cinorg_use, cstindex, & !zhizhao use constant gas conc.
          tag_RO2, nRO2_chem, iRO2, iRO2_cst, RO2index, &   !zhizhao treatment of RO2
          iSumM, 0, aerosol_species_interact(:) )!tag_SumM)!zhizhao add cst sumM

          !ncst_chem, nRO2_chem,iRO2,tag_inorg, &	!zhizhao
          !tag_RO2,cstindex,RO2index,cinorg_use,iHO2,iSumM_cst)!zhizhao
      end if
    !zhizhao
    call cpu_time(t1)
    timer2(t*2-1)=t1-t0

    ! zhizhao put cst_aero, sizebin and layer are assume to be 1
    if (ncst_aero .gt. 0) then
      do s = 1, ncst_aero
        jesp = cst_aero_index(s)
        do j=1,N_size
          concentration_mass(j,jesp) = cst_aero(s,hour+1)
        enddo
      enddo
    endif

    ! re-calculate total_mass(N_aerosol) because mass change due to gas-phase chemistry
    total_aero_mass = 0.d0
    total_mass = 0.d0
    do s = 1, N_aerosol_layers
       jesp = List_species(s)
       do j=1,N_size
         total_aero_mass(jesp) = total_aero_mass(jesp) + concentration_mass(j,s)
       enddo
    enddo
    ! update mass conc. of aerosol precursors
    ! concentration_gas_all(precursor_index) -> concentration_gas(n_aerosol)
    do s = 1, N_aerosol
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas(s) = concentration_gas_all(aerosol_species_interact(s))
       end if
          total_mass(s) = total_mass(s) + concentration_gas(s) + total_aero_mass(s)
    end do

    ! Aerosol dynamic
    CALL SSH_AERODYN(current_time,delta_t2)

    ! update mass conc. of aerosol precursors
    ! concentration_gas(n_aerosol) -> concentration_gas_all(precursor_index)
    do s = 1, N_aerosol
       if (aerosol_species_interact(s) .gt. 0) then
          concentration_gas_all(aerosol_species_interact(s)) = concentration_gas(s)
       end if
    end do
    !zhizhao
    call cpu_time(t0)
    timer2(t*2)=t0-t1
    !

    if (output_type .ne. 0) call ssh_save_concentration()         ! Text or Binary format outout

    ! Time step is finished
    call cpu_time(t0)
    timer(t+2) = t0

  end do ! finsh simulation


  !if (output_type .ne. 0)  call ssh_delete_empty_file() ! delete empty output files

  call ssh_free_allocated_memory()
  IF (with_coag.EQ.1) call ssh_DeallocateCoefficientRepartition()
 
    ! Desallocate arrays for photolysis
  if ((tag_chem .ne. 0).AND.(option_photolysis.eq.2)) then
    call ssh_deallocate_photolysis()    
  endif

  if (ssh_standalone) write(*,*) "============================================"
  if (ssh_standalone) write(*,*) "==== SSH-aerosol simulation completed  ====="
  if (ssh_standalone) write(*,*) "============================================"
  if (ssh_logger) write(logfile,*) "============================================"
  if (ssh_logger) write(logfile,*) "==== SSH-aerosol simulation completed  ====="
  if (ssh_logger) write(logfile,*) "============================================"

  ! Simulation is finished
  call cpu_time(t0)
  timer(nt+3) = t0

  ! Print various times
  !if (ssh_standalone) then
    write(*,*) ""
    write(*,*) "Total simulation time in seconds : ", timer(nt+3) - timer(1)
    write(*,*) "Initialization time in seconds : ", timer(2)-timer(1)
    write(*,*) "Average time per time step in seconds : ", sum(timer(3:nt+2) - timer(2:nt+1))/dble(nt)
    write(*,*) "Maximal time per time step in seconds : ", maxval(timer(3:nt+2) - timer(2:nt+1))
    write(*,*) "Minimal time per time step in seconds : ", minval(timer(3:nt+2) - timer(2:nt+1))
  !endif
  if (ssh_logger) then
    write(logfile,*) ""
    write(logfile,*) "Total simulation time in seconds : ", timer(nt+3) - timer(1)
    write(logfile,*) "Initialization time in seconds : ", timer(2)-timer(1)
    write(logfile,*) "Average time per time step in seconds : ", sum(timer(3:nt+2) - timer(2:nt+1))/dble(nt)
    write(logfile,*) "Maximal time per time step in seconds : ", maxval(timer(3:nt+2) - timer(2:nt+1))
    write(logfile,*) "Minimal time per time step in seconds : ", minval(timer(3:nt+2) - timer(2:nt+1))
  endif

!zhizhao
    inquire (file = trim(output_directory) // "/" // "timer.txt", exist = file_exists)
    if (file_exists) then
            open(unit=10, file = trim(output_directory) // "/" // "timer.txt", status='old', iostat=stat)
            if (stat == 0) close(10, status='delete')
    endif
    ! write the new report
    open(unit=10,file=trim(output_directory) // "/" // "timer.txt", status="new")
      write(unit=10,FMT=*),'time_step\ttime_chem\ttime_dynamic'
      do s=1,nt
         write(unit=10,FMT=*),timer(s+2)-timer(s+1),timer2(s*2-1),timer2(s*2)
      enddo
    close(10)
!
  ! Free memory
  deallocate(timer)
  deallocate(timer2) !zhizhao
  deallocate(cinorg_use) !zhizhao

end PROGRAM SSHaerosol
