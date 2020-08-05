
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseInput.f90
    !> \brief   Parse the contents of an input script file.
    !> \author  M. Poschmann
    !> \date    Nov. 28, 2018
    !> \sa      CheckThermoInput.f90
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified and is completely
    ! independent of ChemApp and related products, including Solgas, Solgasmix, Fact, FactSage
    ! and ChemSage.
    !
    !
    ! Revisions:
    ! ==========
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   28/11/2018      M. Poschmann         Original code
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this subroutine is to parse an input file for Thermochimica
    !
    !
    ! Pertinent variables:
    ! ====================
    !> \param[in]   cInputFileName     The name of the input file to be read.
    !
    ! dTemperature              Temperature (converted to K)
    ! dPresssure                Absolute hydrostatic pressure (converted to atm)
    ! dElementMass              Total mass of an element, where the coefficient corresponds to the
    !                            atomic number (e.g., dMolesElement(92) refers to uranium).
    ! cThermoInputUnits         A character vector containing the units for temperature, pressure and
    !                            mass.
    ! INFOThermo                A scalar integer that indicates a successful exit or identifies an error.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine ParseInput(cInputFileName,dTempLow,dTempHigh,dDeltaT)

  USE ModuleThermoIO
  USE ModuleGEMSolver
  USE ModuleParseCS

  implicit none

  character(*)                :: cInputFileName
  integer                     :: iDelimiterPosition, iOpenPosition, iClosePosition, iElementNumber, iColon1, iColon2
  logical                     :: lEnd, lPressure, lTemperature, lMass, lPressureUnit, lTemperatureUnit, lMassUnit, lData
  character(:), allocatable   :: cLine, cErrMsg, cTag, cValue, cElementNumber
  character(1024)             :: cLineInit
  real(8), intent(out)        :: dTempLow, dTempHigh, dDeltaT

  ! Initialize INFO
  INFO = 0
  ! Open input file
  open (UNIT = 1, FILE = cInputFileName, STATUS = 'old', ACTION = 'read', IOSTAT = INFO)
  ! Check for error on attempt to open
  if (INFO /= 0) then
    INFOThermo = 40
    print *, 'Cannot open input file ' // cInputFileName
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
    read (1,'(A)',IOSTAT = INFO) cLineInit
    ! If there was an error on read, give line number and return
    if (INFO > 0) then
      INFOThermo = 41
      write (cErrMsg, '(A35,I10)') 'Reading input file failed on line: ', iCounter
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
      case ('p','P','pressure','press','Pressure','Press')
        read(cValue,*,IOSTAT = INFO) dPressure
        if (INFO /= 0) then
          INFOThermo = 44
          write (cErrMsg, '(A30,I10)') 'Cannot read pressure on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
        lPressure = .TRUE.
      case ('t','temp','temperature','T','Temp','Temperature')
        iColon1 = 0
        iColon2 = 0
        iColon1 = INDEX(cValue,':')
        if (iColon1 > 0) then
          iColon2 = INDEX(cValue(iColon1+1:),':') + iColon1
          read(cValue(1:iColon1-1),*,IOSTAT = INFO) dTempLow
          if (iColon2 > 0) then
            read(cValue(iColon1+1:iColon2-1),*,IOSTAT = INFO) dTempHigh
            read(cValue(iColon2+1:),*,IOSTAT = INFO) dDeltaT
          else
            read(cValue(iColon1+1:),*,IOSTAT = INFO) dTempHigh
            if (dTempHigh >  dTempLow) dDeltaT =  10
            if (dTempHigh <  dTempLow) dDeltaT = -10
            if (dTempHigh == dTempLow) dDeltaT =  0
          end if
        else
          read(cValue,*,IOSTAT = INFO) dTempLow
          dTempHigh = dTempLow
          dDeltaT = 0
        end if
        if (INFO /= 0) then
          INFOThermo = 44
          print *, 'Cannot read temperature on line: ', iCounter
          return
        endif
        lTemperature = .TRUE.
      case ('m','mass','M','Mass')
        read(cValue,*,IOSTAT = INFO) dElementMass(iElementNumber)
        if (INFO /= 0) then
          INFOThermo = 44
          write (cErrMsg, '(A26,I10)') 'Cannot read mass on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
        lMass = .TRUE.
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
      case default
        write (cErrMsg, '(A34,I10)') 'Input tag not recognized on line: ', iCounter
        print *,  trim(cErrMsg)
    endselect
  enddo LOOP_ReadFile

  ! Now check that all required variables have been set
  if (.NOT. lPressure) then
    INFOThermo = 45
    print *,  'Pressure not set'
    return
  endif
  if (.NOT. lTemperature) then
    INFOThermo = 45
    print *,  'Temperature not set'
    return
  endif
  if (.NOT. lMass) then
    INFOThermo = 45
    print *,  'No masses set'
    return
  endif
  if (.NOT. lPressureUnit) then
    INFOThermo = 45
    print *,  'Pressure unit not set'
    return
  endif
  if (.NOT. lTemperatureUnit) then
    INFOThermo = 45
    print *,  'Temperature unit not set'
    return
  endif
  if (.NOT. lMassUnit) then
    INFOThermo = 45
    print *,  'Mass unit not set'
    return
  endif

end subroutine ParseInput
