!!-----------------------------------------------------------------------
!!     Copyright (C) 2019 CEREA (ENPC) - INERIS
!!     SSH-aerosol is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

Module Resultoutput
  use aInitialization
  use dPhysicalbalance

  implicit none

! out_aero : array of file names; outpout time variation of organic, inorganic, PM2.5, PM10 results
  character(20), save :: out_aero(4) 
  character(4), save :: out_type(2) = (/".txt",".bin"/) ! 1: text, 2: binary

contains

   subroutine ssh_save_report()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine write report file "report.txt", which records most
!     settings, physical conditions and options adopted in the simulation.
!
!     File "report.txt" is saved in the directory : output_directory/
!     provided by user (namelist.ssh).
!
!------------------------------------------------------------------------
!
!     -- OUTPUT 
!     "report.txt"
!
!------------------------------------------------------------------------

   implicit none
	integer :: j,stat
        logical :: file_exists

    ! delete the old report file if it exists under the current saving directory
    inquire (file = trim(output_directory) // "/" // "report.txt", exist = file_exists)
    if (file_exists) then
            open(unit=10, file = trim(output_directory) // "/" // "report.txt", status='old', iostat=stat)
            if (stat == 0) close(10, status='delete')
    endif
    ! write the new report
    open(unit=10,file=trim(output_directory) // "/" // "report.txt", status="new")

	write(unit=10,FMT=*) '<<<< Meteorological setup >>>>'
        write(unit=10,FMT=*) 'location', latitude, 'N','	', longitude,'E','	','Temperature', Temperature, 'K'
        write(unit=10,FMT=*) 'Pressure', Pressure, 'Pa', '	','Specific Humidity', Humidity, '	',&
				'Cloud attenuation field', attenuation, '	relative humidity', relative_humidity
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Simulation time setup >>>>'
        write(unit=10,FMT=*) 'Begining time (from Jan. 1st)', initial_time, 's','	',&
				'Simulation Time', final_time,'s','	','Initial Time Step', delta_t,'s',&
         			'	','Number of iterations:', nt
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Inition condition >>>>'
	if (tag_init == 0) write(unit=10,FMT=*) &
			'Internally mixed aerosol species are provided for initial condition.'
	write(unit=10,FMT=*)  'Gas-phase conc. input file :', init_gas_conc_file
	write(unit=10,FMT=*)  'Particle conc. input file :', init_aero_conc_mass_file
	write(unit=10,FMT=*)  'N_sizebin', N_sizebin
        if (with_init_num .eq. 1) write(unit=10,FMT=*)  'Aerosol number conc. is read from file :',&
							 init_aero_conc_num_file
	if (with_init_num .eq. 0) then
		write(unit=10,FMT=*) ' Aerosol number conc. is estimated from mass and diameter.' 
		write(unit=10,FMT=*)  '====== concentration_number : ======'
		write(unit=10,FMT=*)  concentration_number
	end if
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Mixing state >>>>'
	if (tag_external == 1) write(unit=10,FMT=*)  'simulation is mixing-state resolved.', &
					'N_frac', N_frac ,'	','N_groups', N_groups
	if (tag_external == 0) write(unit=10,FMT=*)  'simulation is internally mixed.' ,&
					'N_frac', N_frac ,'	','N_groups', N_groups
	write(unit=10,FMT=*) 'frac_bound', frac_bound
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Species lists >>>>'
	write(unit=10,FMT=*)  'gas phase species file :', species_list_file
	write(unit=10,FMT=*)  'particle species file :', aerosol_species_list_file
	write(unit=10,FMT=*)
	if (tag_chem == 0) write(unit=10,FMT=*) '<<<< Without Gas-phase chemistry >>>>'
	if (tag_chem == 1) then
		write(unit=10,FMT=*) '<<<< Gas-phase chemistry >>>>'
		write(unit=10,FMT=*)  'with_heterogeneous', with_heterogeneous,'	','with_adaptive', with_adaptive,&
					'	','adaptive time step tolerance', adaptive_time_step_tolerance,&
					'	','min adaptive time step', min_adaptive_time_step
	end if
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Emissions >>>>'
	if (tag_emis == 1) then
		write(unit=10,FMT=*)  'With internally-mixed emissions.'
		write(unit=10,FMT=*)  'Gas-phase conc. emission file :', emis_gas_file
		write(unit=10,FMT=*)  'Particle conc. emission file :', emis_aero_mass_file
		if (with_emis_num == 1) write(unit=10,FMT=*) &
			'Emitted aerosol number conc. is read from file :', emis_aero_num_file
		if (with_emis_num == 0) write(unit=10,FMT=*) &
			'Emitted aerosol number conc. is estimated from mass and diameter.'
	end if
        if (tag_emis == 0) write(unit=10,FMT=*) 'Without emission.'
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) '<<<< Particle Dynamic >>>>'
	if (with_cond == 1) then 
		write(unit=10,FMT=*)  'With condensation', '	','Cut_dim', Cut_dim,'	', 'ISOAPDYN', ISOAPDYN
	else
		write(unit=10,FMT=*)  'Without condensation'
	end if
	if (with_coag == 1) then 
		write(unit=10,FMT=*)  'With coagulation', '	','coefficient file :', Coefficient_file
	else
		write(unit=10,FMT=*)  'Without coagulation'
	end if
	if (with_nucl == 1) then 
		write(unit=10,FMT=*)  'With nucleation', '	','nucl_model_binary', nucl_model_binary
		write(unit=10,FMT=*)  'With nucleation', '	','nucl_model_ternary', nucl_model_ternary
		write(unit=10,FMT=*)  'With nucleation', '	','nucl_model_hetero', nucl_model_hetero
	else
		write(unit=10,FMT=*)  'Without nucleation'
	end if
	write(unit=10,FMT=*)  'DTAEROMIN', DTAEROMIN, '	','redistribution_method', redistribution_method
	write(unit=10,FMT=*)  'Method',dynamic_solver, '	 with_oligomerization', with_oligomerization
	write(unit=10,FMT=*)  'with_fixed_density', with_fixed_density,'	', 'fixed_density', fixed_density, 'kg/m^3'
	write(unit=10,FMT=*)
	write(unit=10,FMT=*) 'output directory :', trim(output_directory),'/'

   CLOSE(10)

   end subroutine ssh_save_report


  subroutine ssh_init_output_conc_sim()

    integer :: stat, s, b
    logical :: file_exists
    character (len=200) output_filename
    character (len=200) :: cmd

    ! Create directory if it does not exist.
    cmd = trim('mkdir -p '// trim(output_directory) //"/aero/")
    call system(cmd)

    ! organics
    output_filename = trim(output_directory) // "/aero/Organics_1" // trim(out_type(output_type))
    ! Remove if output files exist
    inquire (file = output_filename, exist = file_exists)
    if (file_exists) then
        open(unit=100, file = output_filename, status='old', iostat=stat)
        if (stat == 0) close(100, status='delete')
    endif
    ! creative new empty file 
    open(unit=11,file=output_filename, status="new")
    !close(100)

    ! aerosols
    do b = 1, N_size
        do s = 1, n_output_aero
            output_filename = trim(output_directory) // "/aero/"&
            //trim(output_aero(s)) &
            //"_"//trim(str(b)) // trim(out_type(output_type))
            ! Remove if output files exist
            inquire (file = output_filename, exist = file_exists)
            if (file_exists) then
                open(unit=100, file = output_filename, status='old', iostat=stat)
                if (stat == 0) close(100, status='delete')
            endif
            ! creative new empty file ! totoal: Nb*Ns
            open(unit=s+(b-1)*n_output_aero+11 ,file=output_filename, status="new")
            !close(100)
        enddo
    end do

    ! add output for gas
    if (n_output_gas.gt.0.) then
        cmd = trim('mkdir -p '// trim(output_directory) //"/gas/")
        call system(cmd)
    endif

    do s = 1, n_output_gas
        output_filename = trim(output_directory) // "/gas/" &
           //trim(output_gas(s))// trim(out_type(output_type))
        ! Remove if output files exist
        inquire (file = output_filename, exist = file_exists)
        if (file_exists) then
            open(unit=100, file = output_filename, status='old', iostat=stat)
            if (stat == 0) close(100, status='delete')
        endif
        ! creative new empty file 
        open(unit=s+N_size*n_output_aero+11,file=output_filename, status="new")
        !close(100)
    enddo
    ! gas

  end subroutine ssh_init_output_conc_sim

  subroutine ssh_save_concentration_sim(t)

    integer :: s, b, i, j, jesp, t
    !double precision :: tmp
    !character (len=100) output_filename

    concentration_mass_tmp = 0.d0

    do b = 1, N_size
       do s = 1, N_aerosol_layers
          jesp = List_species(s)
          concentration_mass_tmp(b ,jesp) = concentration_mass_tmp(b ,jesp) + concentration_mass(b ,s)
       enddo
    enddo

    ! **** output_directory/aero/
    if (nout_soa.eq.1) then ! only output total soa
        do b = 1, N_size
          i = 1
          ! compute total SOA
          total_soa(i,t) = 0d0
          do s = 1, N_aerosol ! remove water and no-organics
            !if (aerosol_species_name(s).eq.'PBiMT') cycle
            !if (aerosol_species_name(s).eq.'PSOAlP') cycle
            if (aerosol_type(s).ne.4.or.Index_groups(s).lt.0) cycle ! only keep organics from primary vocs
            
            total_soa(i,t) = total_soa(i,t) + concentration_mass_tmp(b, s)
          end do

          ! save organics
          !write(11,*) total_soa(i,t)

          ! save outputs
          do i = 1, n_output_aero
             s = output_aero_index(i)
             write(i+(b-1)*n_output_aero+11,*) concentration_mass_tmp(b, s)
          enddo
        end do
    else ! only for MT
        do b = 1, N_size
            do i = 1, nout_soa
              ! compute total SOA
              total_soa(i,t) = 0d0
              do s = 1, N_aerosol ! remove water and no-organics
                !if (aerosol_species_name(s).eq.'PBiMT') cycle
                !if (aerosol_species_name(s).eq.'PSOAlP') cycle
                if (aerosol_type(s).ne.4.or.Index_groups(s).lt.0) cycle ! only keep organics from primary vocs
                !if (i.ne.nout_soa.and.Index_groups(s).ne.i) cycle ! individual soa

                ! only for MTs
                ! nout_soa = 1: apinene, Index_groups(s) = 1,2,5,7; != 3,4,6
                !if (i.eq.2.and.any((/3,4,6/)==Index_groups(s))) cycle
                ! nout_soa = 2: bpinene, Index_groups(s) = 1,3,5,6; != 2,4,7
                !if (i.eq.3.and.any((/2,4,7/)==Index_groups(s))) cycle
                ! nout_soa = 3: limonene,Index_groups(s) = 1,4,6,7; != 2,3,5
                !if (i.eq.4.and.any((/2,3,5/)==Index_groups(s))) cycle
                if (i.ne.1 .and. Index_groups(s).ne.i-1) cycle
                total_soa(i,t) = total_soa(i,t) + concentration_mass_tmp(b, s)
              end do
            enddo
            ! save organics
            !output_filename = trim(output_directory) // "/aero/Organics_"// trim(str(b)) // trim(out_type(output_type))
            !open(unit=b+10,file=output_filename, status="old", position = "append")
            !write(11,*) (total_soa(i,t),i = 1, nout_total)
            !total_soa(t) = tmp ! save total soa for error analysis
            !close(100)

            ! save outputs
            do i = 1, n_output_aero
              s = output_aero_index(i)
              !output_filename = trim(output_directory) // "/aero/"&
              !      //trim(aerosol_species_name(s))//"_"&
              !      // trim(str(b)) // trim(out_type(output_type))
              !open(unit=100,file=output_filename, status="old", position = "append")
                write(i+(b-1)*n_output_aero+11,*) concentration_mass_tmp(b, s)
              !close(100)
            enddo
        end do
    !else
    !    print*,'nout_soa not in 1,4. please check.',nout_soa
    !    stop
    endif

    ! get RO2 groups concs
    if (nRO2_group.gt.0) then
        do b = 1, nRO2_group
            total_soa(b+nout_soa,t)=0d0
            do i = 1, nRO2_chem
              s = RO2groups(i) ! index of RO2 group
              j = RO2index(i)  ! index of RO2 species

              ! only for MTs
              if (b.ne.nRO2_group.and.s.ne.b) cycle
              ! nout_soa = 1: apinene, Index_groups(s) = 1,2,5,7; != 3,4,6
              !if (b.eq.1.and.any((/3,4,6/) == s)) cycle
              ! nout_soa = 2: bpinene, Index_groups(s) = 1,3,5,6; != 2,4,7
              !if (b.eq.2.and.any((/2,4,7/) == s)) cycle
              ! nout_soa = 3: limonene,Index_groups(s) = 1,4,6,7; != 2,3,5
              !if (b.eq.3.and.any((/2,3,5/) == s)) cycle
              !if (t.eq.1.and.concentration_gas_all(j).gt.0d0) print*,j,trim(species_name(j)),concentration_gas_all(j)
              total_soa(b+nout_soa,t)=total_soa(b+nout_soa,t)+ concentration_gas_all(j)

            enddo
        enddo
        ! total RO2
        !if (iRO2.ne.0) then
        !  total_soa(nout_total,t) = concentration_gas_all(iRO2)
        !else
        !  do i = 1, nRO2_chem
        !      j = RO2index(i)  ! index of RO2 species
        !      total_soa(nout_total,t)=total_soa(nout_total,t)+ concentration_gas_all(j)
        !  enddo
        !endif
    endif

    ! round values
    do i = 1, nout_total
        j = anint(total_soa(i,t)*1d6)
        total_soa(i,t) = j/1d6
    enddo
    ! output total SOA
    write(11,'(999E15.6)') (total_soa(i,t),i = 1, nout_total)
    !write(11,*) (total_soa(i,t),i = 1, nout_total)

    ! gas
    do i = 1, n_output_gas
      s = output_gas_index(i)
      !output_filename = trim(output_directory) // "/gas/" &
      !   //trim(species_name(s))// trim(out_type(output_type))
      !open(unit=100,file=output_filename, status='old', position = "append")
       write(i+N_size*n_output_aero+11,*) concentration_gas_all(s) 
      !close(100)
    enddo
    ! gas

  end subroutine ssh_save_concentration_sim

  subroutine ssh_close_file_sim()
    integer :: s, b, i, j, jesp
    double precision :: conc_save

    ! **** output_directory/aero/
    close(11)
    ! aero
    do b = 1, N_size
        do i = 1, n_output_aero
           close(i+(b-1)*n_output_aero+11)
        enddo
    end do
    ! gas
    do i = 1, n_output_gas
      close(i+N_size*n_output_aero+11)
    enddo
  end subroutine ssh_close_file_sim


  subroutine ssh_save_concentration()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine records simulation results over each time step. 
!
!------------------------------------------------------------------------
!
!     -- OUTPUTS 
!     Mass concentrations of each gas-phase species:
!     >>> output_directory/gas/species name/".txt"(".bin")
!
!     Mass concentrations of each aerosol species in each grid cell:
!     >>> output_directory/aero/aerosol species name_i/".txt"(".bin"), i = 1,2,3...N_size
!
!     Mass concentrations of organic, inorganic, Black_Carbon, Dust, PM2.5, PM10:
!     >>> output_directory/aero/name/".txt"(".bin")
!
!     Total mass concentrations of each aerosol species:
!     >>> output_directory/TM/aerosol species name_TM/".txt"(".bin")
!
!     Total mass concentrations of each aerosol species + precursor:
!     >>> output_directory/TM/aerosol species name_precursor name_TM/".txt"(".bin")
!
!     Number concentrations of each grid cell :
!     >>> output_directory/number/NUMBER_i/".txt"(".bin"), i = 1,2,3...N_size
!
!     Total number concentrations :
!     >>> output_directory/number/TNUM/".txt"(".bin"), i = 1,2,3...N_size
!
!     Average diameter of each grid cell :
!     >>> output_directory/diameter/DIAMETER_i/".txt"(".bin"), i = 1,2,3...N_size
!
!     -- unitS : 
!     mass concentration [ug/m3] 
!     number concentration [#/m3]
!     diameter [um]
!
!------------------------------------------------------------------------

    logical :: iPM25, iPM10
    integer :: s, b, i, j, jesp
    double precision :: conc_save, out_conc(4) ! for out_aero
    character (len=200) output_filename

      output_filename = trim(output_directory) // "/meteo.dat"
      open(unit=100,file=output_filename, status='old', position = "append")
       write(100,*) temperature, pressure, humidity, relative_humidity
      close(100)

    ! **** output_directory/gas/
    ! save gas concentration results over each time step
    do s = 1, n_gas
          output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // trim(out_type(output_type))
          open(unit=100,file=output_filename, status='old', position = "append")
	       write(100,*) concentration_gas_all(s) 
          close(100)
    enddo

    do s = 1, N_aerosol
       do b = 1, N_size
          concentration_mass_tmp(b ,s) = 0.d0
       enddo
    enddo
    
    do b = 1, N_size
       do s = 1, N_aerosol_layers
          jesp = List_species(s)
          concentration_mass_tmp(b ,jesp) = concentration_mass_tmp(b ,jesp) + concentration_mass(b ,s)
       enddo
    enddo

    ! **** output_directory/aero/
    ! save aerosol concentration results over each time step
    do s = 1, N_aerosol
       do b = 1, N_size
          output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &  
                  // "_" // trim(str(b)) // trim(out_type(output_type))
          open(unit=100,file=output_filename, status="old", position = "append")
               write(100,*) concentration_mass_tmp(b, s)
          close(100)
       end do
    end do


    !** save organic, inorganic, PM2.5, PM10 per each time step
    out_conc = 0.d0 ! init

    do b = 1, N_size
        ! check diameter for PM2.5 and PM10
        s = concentration_index(b,1)! get index of size bins
        if (diam_bound(s) .gt. 1d1) then ! d > 10
            iPM10 = .false.
            iPM25 = .false.
        else if (diam_bound(s) .gt. 2.5d0) then ! 2.5 < d <= 10
            iPM10 = .true.
            iPM25 = .false.
        else ! d <= 2.5
            iPM10 = .true.
            iPM25 = .true.
        end if
        ! compute concs
        do s = 1, N_aerosol ! remove water and no-organics
            ! OM and IM
            if (aerosol_type(s).eq.3) then
                out_conc(1)=out_conc(1)+concentration_mass_tmp(b,s)! add inorganics
            else if (aerosol_type(s).eq.4) then
                out_conc(2)=out_conc(2)+concentration_mass_tmp(b,s)! add organics
            end if
            ! PM2.5 and PM10
            if (aerosol_type(s).ne.9) then
                if (iPM25) out_conc(3)=out_conc(3)+concentration_mass_tmp(b,s)! add PM2.5
                if (iPM10) out_conc(4)=out_conc(4)+concentration_mass_tmp(b,s)! add PM10
            end if
        end do
    end do
    ! save concs
    do s = 1, 4
        output_filename=trim(output_directory)//"/aero/"//trim(out_aero(s))//trim(out_type(output_type))
        open(unit=100,file=output_filename,status="old",position="append")
            write(100,*) out_conc(s)
        close(100)
    end do

    ! save organics
    ! **** output_directory/aero/
    do b = 1, N_size
        ! compute total SOA
        conc_save = 0.d0
        do s = 1, N_aerosol ! remove water and no-organics
            if (aerosol_type(s).ne. 4) cycle
            if (aerosol_species_name(s).eq.'PSOAlP') cycle
            if (aerosol_species_name(s).eq.'PBiMT') cycle
            conc_save = conc_save + concentration_mass_tmp(b, s)
        end do
        output_filename = trim(output_directory) // "/aero/Organics_"// trim(str(b)) // trim(out_type(output_type))
        open(unit=100,file=output_filename, status="old", position = "append")
            write(100,*) conc_save
        close(100)
    end do

  end subroutine ssh_save_concentration



  subroutine ssh_init_output_conc()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine initiailize output files, which should be called 
!     before save_concentration()
!
!------------------------------------------------------------------------

    integer :: stat, s, b
    logical :: file_exists
    character (len=200) output_filename
    character (len=200) :: cmd
    character (len=10) :: out_dir(5) 
    out_dir(1) = "/aero/"
    out_dir(2) = "/gas/"	
    out_dir(3) = "/number/"
    out_dir(4) = "/TM/"
    out_dir(5) = "/diameter/"
    ! init ! do not change order
    out_aero(1) = 'Inorganic'
    out_aero(2) = 'Organic'
    out_aero(3) = 'PM2.5'
    out_aero(4) = 'PM10'

    ! Create directory if it does not exist.
    do s = 1, 2
       	cmd = trim('mkdir -p '// trim(output_directory) // out_dir(s))
       	call system(cmd)
    end do

      output_filename = trim(output_directory) // "/meteo.dat"
      ! Remove if output file exist
      inquire (file = output_filename, exist = file_exists)
      if (file_exists) then
         open(unit=100, file = output_filename, status='old', iostat=stat)
         if (stat == 0) close(100, status='delete')
      endif
      ! creative new empty file 
      open(unit=100,file=output_filename, status="new")
      close(100)

    ! gas phase
    do s = 1, n_gas
          output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // trim(out_type(output_type))
          ! Remove if output file exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif
          ! creative new empty file 
          open(unit=100,file=output_filename, status="new")
          close(100)
    enddo


    ! organic, inorganic, PM2.5, PM10
    do s = 1, 4
          output_filename = trim(output_directory) // "/aero/" // trim(out_aero(s))//trim(out_type(output_type))
          ! Remove if output file exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif
          ! creative new empty file 
          open(unit=100,file=output_filename, status="new")
          close(100)
    enddo

    ! aerosols
    do s = 1, N_aerosol
       do b = 1, N_size
	  output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &
                  // "_" // trim(str(b)) // trim(out_type(output_type))
          ! Remove if output files exist
          inquire (file = output_filename, exist = file_exists)
          if (file_exists) then
             open(unit=100, file = output_filename, status='old', iostat=stat)
             if (stat == 0) close(100, status='delete')
          endif
          ! creative new empty file 
          open(unit=100,file=output_filename, status="new")
          close(100)
       end do
    end do

    ! total organics
    do b = 1, N_size
        output_filename = trim(output_directory) // "/aero/Organics_"// trim(str(b)) // trim(out_type(output_type))
        ! Remove if output files exist
        inquire (file = output_filename, exist = file_exists)
        if (file_exists) then
            open(unit=100, file = output_filename, status='old', iostat=stat)
            if (stat == 0) close(100, status='delete')
        endif
        ! creative new empty file 
        open(unit=100,file=output_filename, status="new")
        close(100)
    end do

  end subroutine ssh_init_output_conc


  subroutine ssh_delete_empty_file()
!------------------------------------------------------------------------
!
!     -- DESCRIPTION
!
!     This subroutine delete output file in where all the values are zero.
!
!------------------------------------------------------------------------
    integer :: stat, s, b
    real :: conc_value
    logical :: file_exists
    character (len=200) output_filename
    character (len=200) :: cmd
    character (len=10) :: out_dir(5) 
   ! gas phase
    do s = 1, n_gas
       output_filename = trim(output_directory) // "/gas/" // trim(species_name(s)) // trim(out_type(output_type))
       open(unit=100, file = output_filename, status='old', iostat=stat)
           conc_value = 0.0
           do while(stat .eq. 0)
              read(100, *,iostat=stat) conc_value
              if (conc_value .gt. 0.0) exit
           end do
       if (conc_value .eq. 0.0) close(100, status='delete') ! if all values are zero, delete file
       if (conc_value .ne. 0.0) close(100)
    enddo

    do s = 1, N_aerosol
    ! aerosols mass conc. of each cell
       do b = 1, N_size
	  output_filename = trim(output_directory) // "/aero/" // trim(aerosol_species_name(s)) &
                  // "_" // trim(str(b)) // trim(out_type(output_type))
          open(unit=100, file = output_filename, status='old', iostat=stat)
              conc_value = 0.0
              do while(stat .eq. 0)
                 read(100, *,iostat=stat) conc_value
                 if (conc_value .gt. 0.0) exit
              end do
          if (conc_value .eq. 0.0) close(100, status='delete')
          if (conc_value .ne. 0.0) close(100)
       end do
    enddo

  end subroutine ssh_delete_empty_file

  character(len=20) function str(k)
    !   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str
  
  end module Resultoutput
