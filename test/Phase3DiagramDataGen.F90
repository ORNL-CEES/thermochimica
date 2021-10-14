program Phase3DiagramDataGen

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo

    implicit none
    character(1024) :: cInputFile
    integer :: i, j, nx1, nx2, nSim, iEl1, iEl2, iEl3
    real(8) :: t, x1lo, x1hi, dDeltaX1, x2lo, x2hi, dDeltaX2, dPress
    character(16) :: intStr

    ! Read input argument to get filename
    call get_command_argument(1, cInputFile)
    if (len_trim(cInputFile) == 0) then
      print *, 'No input file specified'
      call EXIT(1)
    endif

    call ParseInput3PhaseDiagram(cInputFile,t,x1lo,x1hi,dDeltaX1,x2lo,x2hi,dDeltaX2,iEl1,iEl2,iEl3,dPress)

    call ParseCSDataFile(cThermoFileName)

    ! Specify values:
    dPressure = dPress
    dTemperature = t
    if ((x1hi == x1lo) .OR. dDeltaX1 == 0D0) then
      nx1 = 0
    else
      nx1 = CEILING((x1hi - x1lo) / dDeltaX1)
    end if
    if ((x2hi == x2lo) .OR. dDeltaX2 == 0D0) then
      nx2 = 0
    else
      nx2 = CEILING((x2hi - x2lo) / dDeltaX2)
    end if

    open(1, file= DATA_DIRECTORY // '../thermoout.json', &
        status='REPLACE', action='write')
    write(1,*) '{'
    close (1)

    nSim = 1
    do i = 0, nx1
      do j =  0, nx2
        dElementMass(iEl1) = x1lo
        if (nx1 > 0) dElementMass(iEl1) = x1lo + (x1hi - x1lo) * REAL(i)/REAL(nx1)
        dElementMass(iEl2) = x2lo
        if (nx2 > 0) dElementMass(iEl2) = x2lo + (x2hi - x2lo) * REAL(j)/REAL(nx2)
        if ((x1hi == 1D0) .AND. (x1lo == 0) .AND. (x2hi == 1D0) .AND. (x2lo == 0)) then
            dElementMass(iEl2) = REAL(j)/REAL(nx2) * MAX(1D0 - (dElementMass(iEl1)),0D0)
        end if
        dElementMass(iEl3) = MAX(1D0 - (dElementMass(iEl1) + dElementMass(iEl2)),0D0)

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
          if ((i < nx1) .OR. (j < nx2)) write(1,*) ','
          close (1)
          nSim = nSim + 1
        else
          INFOThermo = 0
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

end program Phase3DiagramDataGen
