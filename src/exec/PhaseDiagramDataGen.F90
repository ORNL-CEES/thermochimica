program PhaseDiagramDataGen

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo

    implicit none
    character(1024) :: cInputFile
    integer :: i, j, nt, nx, nSim, iEl1, iEl2
    real(8) :: tlo, thi, xlo, xhi, dTbase, dDeltaT, dDeltaX, dPress
    character(16) :: intStr

    ! Read input argument to get filename
    call get_command_argument(1, cInputFile)
    if (len_trim(cInputFile) == 0) then
      print *, 'No input file specified'
      call EXIT(1)
    endif

    call ParseInputPhaseDiagram(cInputFile,tlo,thi,dDeltaT,xlo,xhi,dDeltaX,iEl1,iEl2,dPress)

    call ParseCSDataFile(cThermoFileName)

    ! Specify values:
    dPressure              = dPress
    if ((thi == tlo) .OR. dDeltaT == 0D0) then
      nt = 0
    else
      nt = CEILING((thi - tlo) / dDeltaT)
    end if
    if ((xhi == xlo) .OR. dDeltaX == 0D0) then
      nx = 0
    else
      nx = CEILING((xhi - xlo) / dDeltaX)
    end if

    open(1, file= DATA_DIRECTORY // '../outputs/thermoout.json', &
        status='REPLACE', action='write')
    write(1,*) '{'
    close (1)

    nSim = 1
    do i = 0, nt
      dTbase = tlo
      if (nt > 0) dTbase = dTbase + (REAL(i)/REAL(nt))*(thi-tlo)
      do j = 0, nx
        dTemperature = dTbase
        if ((i > 0) .AND. (i < nt)) dTemperature = dTemperature + ((REAL(MODULO(j,10))-5D0)/10D0)*(1D0/REAL(nt))*(thi-tlo)
        dTemperature = MAX(dTemperature,tlo)
        dTemperature = MIN(dTemperature,thi)
        dElementMass(iEl2) = xlo
        if (nx > 0) dElementMass(iEl2) = dElementMass(iEl2) + REAL(j)/REAL(nx)*(xhi-xlo)
        dElementMass(iEl1) = 1D0-dElementMass(iEl2)
        call Thermochimica
        if (INFOThermo == 0) then
          open(1, file= DATA_DIRECTORY // '../outputs/thermoout.json', &
              status='OLD', position='append', action='write')
          if ((i > 0) .OR. (j > 0)) write(1,*) ','
          write(intStr,*) nSim
          write(1,*) '"', TRIM(ADJUSTL(intStr)) ,'":'
          close (1)
          call WriteJSON(.TRUE.)
          nSim = nSim + 1
        else
          call ResetThermoAll
          INFOThermo = 0
          call ParseCSDataFile(cThermoFileName)
        end if
        call ResetThermo
      end do
    end do

    open(1, file= DATA_DIRECTORY // '../outputs/thermoout.json', &
        status='OLD', position='append', action='write')
    write(1,*) '}'
    close (1)

end program PhaseDiagramDataGen
