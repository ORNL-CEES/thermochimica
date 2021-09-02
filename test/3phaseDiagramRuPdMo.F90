program thermo

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo

    implicit none

    integer :: i, j, nx, nSim
    character(16) :: intStr

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'
    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Specify values:
    dPressure        = 1D0
    dTemperature     = 900D0

    ! Specify output mode:
    iPrintResultsMode = 0

    nx = 30

    open(1, file= DATA_DIRECTORY // '../thermoout.json', &
        status='REPLACE', action='write')
    write(1,*) '{'
    close (1)

    nSim = 1
    do i = 0, nx
      do j =  0, nx
        dElementMass(42) = REAL(i)/REAL(nx) !+ REAL(MODULO(j,10)) / (10D0 * REAL(nx))
        dElementMass(44) = REAL(j)/REAL(nx) * MAX(1D0 - (dElementMass(42)),0D0)
        dElementMass(46) = MAX(1D0 - (dElementMass(42) + dElementMass(44)),0D0)

        ! Call Thermochimica:
        if (INFOThermo == 0)        call Thermochimica

        if (INFOThermo == 0) then
          open(1, file= DATA_DIRECTORY // '../thermoout.json', &
              status='OLD', position='append', action='write')
          write(intStr,*) nSim
          write(1,*) '"', TRIM(ADJUSTL(intStr)) ,'":'
          close (1)
          call WriteJSON(.TRUE.)
          open(1, file= DATA_DIRECTORY // '../thermoout.json', &
              status='OLD', position='append', action='write')
          if ((i < nx) .OR. (j < nx)) write(1,*) ','
          close (1)
          nSim = nSim + 1
        else
          call ThermoDebug
        end if

        ! Perform post-processing of results:
        if (iPrintResultsMode > 0)  call PrintResults

        ! Destruct everything:
        call ResetThermo
      end do
    end do

    open(1, file= DATA_DIRECTORY // '../thermoout.json', &
        status='OLD', position='append', action='write')
    write(1,*) '}'
    close (1)

end program thermo
