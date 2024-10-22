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

  integer :: t, j, s,jesp,day,hour
  character (len=10) :: ivoc0
  character (len=400) :: namelist_ssh  ! Configuration file
  double precision, dimension(:), allocatable :: timer
  ! need if use constant input concentrations
  double precision, dimension(:), allocatable :: cst_gas_use
  ! error calculation
  double precision, dimension(:), allocatable :: rerr, perr ! storage error
  double precision, dimension(:), allocatable :: rerr_num, perr_num, rerr_deno, perr_deno

  double precision :: t_since_update_photolysis, t0, delta_t2

  ! Initial time (seconds)
  !t!call cpu_time(t0)
  ivoc = 0 ! init

  ! Initialisation: discretization and distribution
  if (iargc() == 0) then
     write(*,*) "usage: ssh-aerosol namelist.input"
     stop
  else
     call getarg(1, namelist_ssh)
     if (iargc() == 2) then
         call getarg(2, ivoc0)
         read(ivoc0,'(I1)') ivoc ! read as integer, only allow one digit for now
         print*,'read ivoc: ', ivoc
     endif
  end if

  ! Read the number of gas-phase species and chemical reactions
  call ssh_dimensions(N_gas, n_reaction, n_photolysis)

  call ssh_read_namelist(namelist_ssh)

  call ssh_read_inputs()                                

  call ssh_read_meteo()
  
  ! Read the meteorological data before output
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
  !t!timer(1) = t0
  !t!call cpu_time(t0)
  !t!timer(2) = t0

  ! for constant concentrations
  if (ncst_gas.gt.0) then
     allocate(cst_gas_use(ncst_gas))
  else
     allocate(cst_gas_use(0))
  endif

  ! init error analysis
  if (ierr_ref) then
     allocate(rerr(nout_total))
     allocate(rerr_num(nout_total))
     allocate(rerr_deno(nout_total))
     rerr = 0.d0
     rerr_num = 0.d0
     rerr_deno = 0.d0
  endif

  if (ierr_pre) then
     allocate(perr(nout_total))
     allocate(perr_num(nout_total))
     allocate(perr_deno(nout_total))
     perr = 0.d0
     perr_num = 0.d0
     perr_deno = 0.d0
  endif
  !print*,'init',rerr,perr

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

     ! get hourly radical conc.
     ! hour
     hour = int(mod(current_time/3600.,24.)) !int((current_time/3600.)%24)
     !print*,'time',current_time,'current hour', hour

     if (tag_chem .ne. 0) then
       ! update constant concentrations if exist
       if (ncst_gas .gt. 0) then
          cst_gas_use = 0.0
          do s =1, ncst_gas !size(cst_gas_index)
             cst_gas_use(s) = cst_gas(s,hour+1)
          enddo
       endif

       if (t==1) then
          !delta_t=1=0.001*delta_t
          call ssh_chem_twostep(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
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
               ncst_gas, cst_gas_use, cst_gas_index, & !genoa use constant gas conc.
               tag_RO2, nRO2_chem, iRO2, iRO2_cst, RO2index, &
               aerosol_species_interact(:))

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

           ! solve chemistry with the two-step time numerical solver if tag_twostep .eq. 1
           call ssh_chem_twostep(n_gas, n_reaction, n_photolysis, photolysis_reaction_index,&
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
              ncst_gas, cst_gas_use, cst_gas_index, & !genoa use constant gas conc.
              tag_RO2, nRO2_chem, iRO2, iRO2_cst, RO2index, &
              aerosol_species_interact(:))
      end if

    ! check ngas
    do s=1, N_gas
       if (concentration_gas_all(s).gt.1d3) then
          print*,'sshError 1: gas conc > 1d3',s,concentration_gas_all(s)
          stop
       endif
    enddo

    ! set cst_aero(n_species,n_size,n_step) if need. N_sizebin is assumed to be N_size. Only for internal mixing.
    if (ncst_aero .gt. 0) then
      do s = 1, ncst_aero
        jesp = cst_aero_index(s)
        do j=1,N_size
          concentration_mass(j,jesp) = cst_aero(s,j,hour+1)
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
      do s = 1, nout_total
        !check conc. if > 1E6: something is wrong
        if (total_soa(s,t) .gt. 1d3) then
            print*,"sshError 2: conc > 1d3! ",total_soa(s,t),s,t
            stop
        endif
        !print*,s,t,total_soa(s,t),ref_soa(s,t)
        rerr_num(s) =  rerr_num(s) + dabs(total_soa(s,t) - ref_soa(s,t))
        rerr_deno(s) = rerr_deno(s) + total_soa(s,t) + ref_soa(s,t)
        if (delta_t * t .eq. 864.d2.or.delta_t * t .eq. 432.d3) then ! set time
           rerr(s) = max(rerr(s), 2d0 * rerr_num(s) / (rerr_deno(s) + TINYM))
           !print*,'rerr',s,rerr(s),rerr_num(s),rerr_deno(s),TINYM
           rerr_num(s) = 0.d0
           rerr_deno(s) = 0.d0
        endif
      enddo
    endif
    if (ierr_pre) then
      do s = 1, nout_total
        !print*,s,t,total_soa(s,t),pre_soa(s,t),perr_num(s),dabs(total_soa(s,t) - pre_soa(s,t))
        perr_num(s) =  perr_num(s) + dabs(total_soa(s,t) - pre_soa(s,t))
        perr_deno(s) = perr_deno(s) + total_soa(s,t) + pre_soa(s,t)
        if (delta_t * t .eq. 864.d2.or.delta_t * t .eq. 432.d3) then
           perr(s) = max(perr(s), 2d0 * perr_num(s) / (perr_deno(s) + TINYM))
           !print*,'perr',s,perr(s),perr_num(s),perr_deno(s),TINYM
           perr_num(s) = 0.d0
           perr_deno(s) = 0.d0
        endif
      enddo
    endif

  end do			! finsh simulation


  ! save error
  if (ierr_ref) then
    print*, rerr
    write(*,101) maxval(rerr) !rerr
  endif

  if (ierr_pre) then
    print*, perr
    write(*,102) maxval(perr) !perr
  endif
 
101 format('err_ref: ',f15.4)
102 format('err_pre: ',f15.4)

  if (output_type .ne. 0) then
    call ssh_close_file_sim() !zhizhao close all files
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
  !t!if (allocated(timer)) deallocate(timer)
  if (allocated(cst_gas_use)) deallocate(cst_gas_use) !genoa
  if (allocated(rerr))  deallocate(rerr)
  if (allocated(perr))  deallocate(perr)
  if (allocated(rerr_num))  deallocate(rerr_num)
  if (allocated(perr_num))  deallocate(perr_num)
  if (allocated(rerr_deno))  deallocate(rerr_deno)
  if (allocated(perr_deno))  deallocate(perr_deno)

end PROGRAM SSHaerosol
