program RunCalculationList

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none
    character(1024) :: cInputFile
    integer :: i, j, nElIn, nCalc
    integer, dimension(:), allocatable :: iEls
    real(8), dimension(:), allocatable :: dEls
    character(:), allocatable :: cLine, cErrMsg, cTag, cValue, cElementNumber
    integer :: iDelimiterPosition, iOpenPosition, iClosePosition, iElementNumber
    character(1024) :: cLineInit
    logical :: lEnd, lPressureUnit, lTemperatureUnit, lMassUnit, lData, lEl, lNel

    character(16) :: intStr

    ! Initialize INFO
    INFO = 0

    ! Read input argument to get filename
    call get_command_argument(1, cInputFile)
    if (len_trim(cInputFile) == 0) then
      print *, 'No input file specified'
      call EXIT(1)
    endif
    ! Open input file
    open (UNIT = 3, FILE = cInputFile, STATUS = 'old', ACTION = 'read', IOSTAT = INFO)
    ! Check for error on attempt to open
    if (INFO /= 0) then
      INFOThermo = 40
      print *, 'Cannot open input file ', cInputFile
      return
    endif

    ! Initialize for read loop
    lEnd = .FALSE.
    iCounter = 0
    ! Read all line of input file
    LOOP_ReadFile: do while (INFO == 0)
      ! Keep track of line number
      iCounter = iCounter + 1
      ! Read a line
      READ(3,'(A)',IOSTAT = INFO) cLineInit
      ! If there was an error on read, give line number and return
      if (INFO > 0) then
        INFOThermo = 41
        WRITE(cErrMsg, '(A35,I10)') 'Reading input file failed on line: ', iCounter
        print *,  trim(cErrMsg)
        return
      ! If file end reached, break loop
      elseif (INFO < 0) then
        exit LOOP_ReadFile
      endif
      ! Remove leading then trailing spaces on line
      cLine = trim(adjustl(cLineInit))
      ! Check for comment line (going to be liberal with choices of comment indicators)
      if (scan(cLine,'!@#$%&*/\?|') == 1) then
        cycle LOOP_ReadFile
      endif
      ! Also look for lines without "="
      if (scan(cLine,'=') == 0) then
        cycle LOOP_ReadFile
      endif
      ! If we get to here, should be a data line, and therefore contain '=' delimiter separating tag and value terms
      iDelimiterPosition = scan(cLine,'=')
      ! Tag is on LHS of delimiter, extract this and trim whitespace
      cTag = trim(adjustl(cLine(1 : (iDelimiterPosition - 1))))
      ! Value if on RHS of delimiter, do same as above
      cValue = trim(adjustl(cLine((iDelimiterPosition + 1) : len(cLine))))
      ! Check if line contains a mass, need to treat these separately
      ! Masses will be the only lines to contain '()', so look for these
      iOpenPosition = scan(cLine,'(')
      if (iOpenPosition > 0) then
        iClosePosition = scan(cLine,')')
        ! Check for no close ')'
        if (iClosePosition == 0) then
          INFOThermo = 42
          write (cErrMsg, '(A31,I10)') 'Open ( but no close ) on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
        cElementNumber = trim(adjustl(cTag((iOpenPosition + 1) : (iClosePosition - 1))))
        read(cElementNumber,*,IOSTAT = INFO) iElementNumber
        if (INFO /= 0) then
          INFOThermo = 43
          write (cErrMsg, '(A36,I10)') 'Cannot read element number on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
        cTag = trim(adjustl(cTag(1 : (iOpenPosition - 1))))
      endif

      ! Now look through possible tags to assign variables
      select case (cTag)
        case ('nel','nelement','nEl','nElement','nelements','nElements')
          read(cValue,*,IOSTAT = INFO) nElIn
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A26,I10)') 'Cannot read number of elements on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          end if
          if (allocated(iEls)) deallocate(iEls,dEls)
          allocate(iEls(nElIn),dEls(nElIn))
          lNel = .TRUE.
        case ('iel','ielement','iEl','iElement')
          if (.NOT. lNel) then
            INFOThermo = 44
            write (cErrMsg, '(A50,I10)') 'Need number of elements before element indices at ', iCounter
            print *,  trim(cErrMsg)
            return
          end if
          read(cValue,*,IOSTAT = INFO) iEls
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A26,I10)') 'Cannot read elements on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          lEl = .TRUE.
        case ('ncalc','nCalc')
          read(cValue,*,IOSTAT = INFO) nCalc
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A44,I10)') 'Cannot read number of calculations on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          exit LOOP_ReadFile
        case ('p_unit','pressure_unit','Pressure_unit','Pressure_Unit','P_unit','pressure unit',&
          'Pressure Unit','Pressure unit','press_unit','press unit','Press unit','Press Unit',&
          'p unit','P unit','P Unit',&
          'p_units','pressure_units','Pressure_units','Pressure_Units','P_units','P_Units','pressure units',&
            'Pressure Units','Pressure units','press_units','press units','Press units','Press Units'&
            'p units','P units','P Units')
          read(cValue,*,IOSTAT = INFO) cInputUnitPressure
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A35,I10)') 'Cannot read pressure unit on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          lPressureUnit = .TRUE.
        case ('t_unit','temperature_unit','Temperature_unit','Temperature_Unit','T_unit','T_Unit',&
          'temperature unit','Temperature Unit','Temperature unit','temp_unit','temp unit',&
          'Temp unit','Temp Unit','t unit','T unit','T Unit',&
          't_units','temperature_units','Temperature_units','Temperature_Units','T_units','T_Units',&
            'temperature units','Temperature Units','Temperature units','temp_units','temp units',&
            'Temp units','Temp Units','t units','T units','T Units')
          read(cValue,*,IOSTAT = INFO) cInputUnitTemperature
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A38,I10)') 'Cannot read temperature unit on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          lTemperatureUnit = .TRUE.
        case ('m_unit','mass_unit','Mass_unit','Mass_Unit','m unit','mass unit','Mass unit','Mass Unit',&
          'm_units','mass_units','Mass_units','Mass_Units','m units','mass units','Mass units','Mass Units')
          read(cValue,*,IOSTAT = INFO) cInputUnitMass
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A31,I10)') 'Cannot read mass unit on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          lMassUnit = .TRUE.
        case ('data','Data','data_file','Data_file','data file','Data file','Data File',&
          'dat','Dat','dat_file','Dat_file','dat file','Dat file','Dat File')
          read(cValue,'(A)',IOSTAT = INFO) cThermoFileName
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A35,I10)') 'Cannot read data filename on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          lData = .TRUE.
        case ('print_mode','Print_mode','Print_Mode',&
          'print mode','Print mode','Print Mode')
          read(cValue,*,IOSTAT = INFO) iPrintResultsMode
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A32,I10)') 'Cannot read print mode on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
        case ('debug_mode','Debug_mode','Debug_Mode',&
          'debug mode','Debug mode','Debug Mode')
          read(cValue,*,IOSTAT = INFO) lDebugMode
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A32,I10)') 'Cannot read debug mode on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
        case ('reinit','Reinit','reinitialization','Reinitialization',&
          'reinit_mode','Reinit_mode','Reinit_Mode','reinitialization_mode','Reinitialization_mode','Reinitialization_Mode',&
          'reinit mode','Reinit mode','Reinit Mode','reinitialization mode','Reinitialization mode','Reinitialization Mode')
          read(cValue,*,IOSTAT = INFO) lReinitRequested
          if (INFO /= 0) then
            INFOThermo = 44
            write (cErrMsg, '(A43,I10)') 'Cannot read reinitialization mode on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
        case default
          write (cErrMsg, '(A34,I10)') 'Input tag not recognized on line: ', iCounter
          print *,  trim(cErrMsg)
      endselect
    end do LOOP_ReadFile

    ! Now check that all required variables have been set
    if (.NOT. lNel) then
      INFOThermo = 45
      print *, 'Number of elements not set'
      return
    endif
    if (.NOT. lEl) then
      INFOThermo = 45
      print *, 'No elements set'
      return
    endif
    if (.NOT. lPressureUnit) then
      INFOThermo = 45
      cInputUnitPressure = 'atm'
      return
    endif
    if (.NOT. lTemperatureUnit) then
      INFOThermo = 45
      cInputUnitTemperature = 'K'
      return
    endif
    if (.NOT. lMassUnit) then
      INFOThermo = 45
      cInputUnitMass = 'moles'
      return
    endif

    call ParseCSDataFile(cThermoFileName)

    ! Specify values:

    OPEN(2, file= DATA_DIRECTORY // '../thermoout.json', &
        status='REPLACE', action='write')
    WRITE(2,*) '{'
    CLOSE(2)

    do i = 1, nCalc
      READ(3,*,IOSTAT = INFO) dTemperature, dPressure, dEls
      dElementMass = 0D0
      do j = 1, nElIn
        dElementMass(iEls(j)) = dEls(j)
      end do
      call Thermochimica
      call PrintResults
      open(2, file= DATA_DIRECTORY // '../thermoout.json', &
          status='OLD', position='append', action='write')
      write(intStr,*) i
      write(2,*) '"', TRIM(ADJUSTL(intStr)) ,'":'
      close (2)
      call WriteJSON(.TRUE.)
      open(2, file= DATA_DIRECTORY // '../thermoout.json', &
          status='OLD', position='append', action='write')
      if (i < nCalc) write(2,*) ','
      close (2)
      INFOThermo = 0
      call ResetThermo
    end do
    CLOSE(3)

    open(2, file= DATA_DIRECTORY // '../thermoout.json', &
        status='OLD', position='append', action='write')
    write(2,*) '}'
    close (2)

end program RunCalculationList
