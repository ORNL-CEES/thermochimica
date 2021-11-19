program thermo

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo

    implicit none

    integer :: i, j, nt, nx, nSim
    real(8) :: tlo, thi, xlo, xhi, dTbase
    logical :: lJson
    character(16) :: intStr

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'Kaye_NobleMetals.dat'
    iPrintResultsMode     = 2

    call ParseCSDataFile(cThermoFileName)

    ! Specify values:
    dPressure              = 1D0
    tlo                    = 900
    thi                    = 2700
    xlo                    = 0
    xhi                    = 1
    nt                     = 100
    nx                     = 100

    open(1, file= DATA_DIRECTORY // '../thermoout.json', &
        status='REPLACE', action='write')
    write(1,*) '{'
    close (1)

    lJson = .FALSE.
    dTemperature = (thi + tlo) / 2
    dElementMass(44) = 0.5        ! Ru
    dElementMass(46) = 0.5        ! Pd
    nSim = 1
    call Thermochimica
    open(1, file= DATA_DIRECTORY // '../thermoout.json', &
        status='OLD', position='append', action='write')
    write(intStr,*) nSim
    write(1,*) '"', TRIM(ADJUSTL(intStr)) ,'":'
    close (1)
    call WriteJSON(.TRUE.)
    open(1, file= DATA_DIRECTORY // '../thermoout.json', &
        status='OLD', position='append', action='write')
    write(1,*) ','
    close (1)
    ! call PrintResults
    call ResetThermo
    ! call exit
    lJson = .TRUE.
    do i = 0, nt
      dTbase = tlo + (REAL(i)/REAL(nt))*(thi-tlo)
      do j = 0, nx
        dTemperature = dTbase + (REAL(MODULO(j,10))/10D0)*(1D0/REAL(nt))*(thi-tlo)
        nSim = nSim + 1
        dElementMass(44) = xlo + REAL(j)/REAL(nx)*(xhi-xlo)
        dElementMass(46) = 1D0-dElementMass(44)
        call Thermochimica
        open(1, file= DATA_DIRECTORY // '../thermoout.json', &
            status='OLD', position='append', action='write')
        write(intStr,*) nSim
        write(1,*) '"', TRIM(ADJUSTL(intStr)) ,'":'
        close (1)
        call WriteJSON(.TRUE.)
        open(1, file= DATA_DIRECTORY // '../thermoout.json', &
            status='OLD', position='append', action='write')
        if ((i < nt) .OR. (j < nx)) write(1,*) ','
        close (1)
        ! call PrintResults
        call ResetThermo
        lJson = .TRUE.
        ! call ThermoDebug
      end do
    end do

    open(1, file= DATA_DIRECTORY // '../thermoout.json', &
        status='OLD', position='append', action='write')
    write(1,*) '}'
    close (1)

end program thermo
