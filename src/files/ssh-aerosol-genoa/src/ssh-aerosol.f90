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
  !t!double precision, dimension(:), allocatable :: timer
  double precision, dimension(:), allocatable :: cinorg_use !genoa
  double precision:: rerr, perr !genoa storage error
  double precision:: rerr_num, perr_num, rerr_deno, perr_deno

  double precision :: t_since_update_photolysis, t0,t1, delta_t2

  ! Initial time (seconds)
  !t!call cpu_time(t0)

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

  if (output_type .ne. 0) call ssh_init_output_conc_sim() 

  !call ssh_save_report()

  if (output_type .ne. 0) call ssh_save_concentration_sim(1) 

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
  !t!allocate(timer(nt+3))
  !allocate(cinorg_use(size(cstindex))) !genoa
  allocate(cinorg_use(ncst_chem))!genoa

  ! init error analysis
  rerr = 0.d0
  perr = 0.d0
  rerr_num = 0.d0
  perr_num = 0.d0
  rerr_deno = 0.d0
  perr_deno = 0.d0

  !t!timer(1) = t0
  !t!call cpu_time(t0)
  !t!timer(2) = t0

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

     !t!call cpu_time(t0) !genoa
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

       if (tag_RO2 .eq. 4) then
          if (iRO2_cst .ne. 0) then
              cinorg_use(iRO2_cst) = RO2_pool(t)
          else
            print*, 'tag_RO2 = 4 but iRO2_cst = 0. need to input cst gas file with RO2.'
            stop
          endif

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
               ncst_chem, cinorg_use, cstindex, & !genoa use constant gas conc.
               tag_RO2, nRO2_chem, iRO2, iRO2_cst, RO2index, &   !genoa treatment of RO2
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
          ncst_chem, cinorg_use, cstindex, & !genoa use constant gas conc.
          tag_RO2, nRO2_chem, iRO2, iRO2_cst, RO2index, &   !genoa treatment of RO2
          iSumM, 0, aerosol_species_interact(:) )!tag_SumM)!genoa add cst sumM

          !ncst_chem, nRO2_chem,iRO2,tag_inorg, &	!genoa
          !tag_RO2,cstindex,RO2index,cinorg_use,iHO2,iSumM_cst)!genoa
      end if

    ! genoa put cst_aero, sizebin and layer are assume to be 1
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

    if (output_type .ne. 0) call ssh_save_concentration_sim(t)         ! Text or Binary format outout

    ! Time step is finished
    !t!call cpu_time(t0)
    !t!timer(t+2) = t0

    ! error analysis
    if (ierr_ref) then
       rerr_num =  rerr_num + dabs(total_soa(t) - ref_soa(t))
       rerr_deno = rerr_deno + total_soa(t) + ref_soa(t)
       if (delta_t * t .eq. 864.d2.or.delta_t * t .eq. 432.d3) then ! set time
          rerr = max(rerr, 2 * rerr_num / (rerr_deno + TINYM))
          rerr_num = 0.d0
          rerr_deno = 0.d0
      endif
    endif
    if (ierr_pre) then
       perr_num =  perr_num + dabs(total_soa(t) - pre_soa(t))
       perr_deno = perr_deno + total_soa(t) + pre_soa(t)
      if (delta_t * t .eq. 864.d2.or.delta_t * t .eq. 432.d3) then
          perr = max(perr, 2 * perr_num / (perr_deno + TINYM))
          perr_num = 0.d0
          perr_deno = 0.d0
      endif
    endif

  end do ! finsh simulation

  ! save error
  if (ierr_ref) print*,'err_ref: ',rerr
  if (ierr_pre) print*,'err_pre: ',perr
  
  if (output_type .ne. 0) then
    call ssh_close_file_sim() !genoa close all files
    !call ssh_delete_empty_file() ! delete empty output files
  endif

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
  !t!call cpu_time(t0)
  !t!timer(nt+3) = t0

  ! Print various times
  !t!if (ssh_standalone) then
    !t!write(*,*) ""
    !t!write(*,*) "Total simulation time in seconds : ", timer(nt+3) - timer(1)
    !t!write(*,*) "Initialization time in seconds : ", timer(2)-timer(1)
    !t!write(*,*) "Average time per time step in seconds : ", sum(timer(3:nt+2) - timer(2:nt+1))/dble(nt)
    !t!write(*,*) "Maximal time per time step in seconds : ", maxval(timer(3:nt+2) - timer(2:nt+1))
    !t!write(*,*) "Minimal time per time step in seconds : ", minval(timer(3:nt+2) - timer(2:nt+1))
  !t!endif
  !t!if (ssh_logger) then
    !t!write(logfile,*) ""
    !t!write(logfile,*) "Total simulation time in seconds : ", timer(nt+3) - timer(1)
    !t!write(logfile,*) "Initialization time in seconds : ", timer(2)-timer(1)
    !t!write(logfile,*) "Average time per time step in seconds : ", sum(timer(3:nt+2) - timer(2:nt+1))/dble(nt)
    !t!write(logfile,*) "Maximal time per time step in seconds : ", maxval(timer(3:nt+2) - timer(2:nt+1))
    !t!write(logfile,*) "Minimal time per time step in seconds : ", minval(timer(3:nt+2) - timer(2:nt+1))
  !t!endif

  ! Free memory
  !t!deallocate(timer)
  deallocate(cinorg_use) !genoa

end PROGRAM SSHaerosol
