
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseInput3PhaseDiagram.f90
    !> \brief   Parse the contents of an input script file.
    !> \author  M. Poschmann
    !> \date    Sep. 2, 2021
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
    !   02/09/2021      M. Poschmann         Original code
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this subroutine is to parse an input file for Thermochimica generation of
    !   phase diagram data.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine ParseInput3PhaseDiagram(cInputFileName,dTemp,dX1lo,dX1hi,dDeltaX1,dX2lo,dX2hi,dDeltaX2,iEl1,iEl2,iEl3,dPress)

  USE ModuleThermoIO
  USE ModuleGEMSolver
  USE ModuleParseCS

  implicit none

  character(*)                :: cInputFileName
  integer                     :: iDelimiterPosition, iOpenPosition, iClosePosition, iElementNumber, iColon1, iColon2
  logical                     :: lEnd, lPressure, lTemperature, lPressureUnit, lTemperatureUnit, lData, lEl, lMassUnit
  character(:), allocatable   :: cLine, cErrMsg, cTag, cValue, cElementNumber
  character(1024)             :: cLineInit, cThermoFileNameTemp
  real(8), intent(out)        :: dTemp, dPress, dX1lo, dX1hi, dDeltaX1, dX2lo, dX2hi, dDeltaX2
  integer, intent(out)        :: iEl1, iEl2, iEl3

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

  dX1lo = 0D0
  dX1hi = 1D0
  dDeltaX1 = 1D-1
  dX2lo = 0D0
  dX2hi = 1D0
  dDeltaX2 = 1D-1

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
        iColon1 = 0
        iColon2 = 0
        iColon1 = INDEX(cValue,':')
        if (iColon1 > 0) then
        else
          read(cValue,*,IOSTAT = INFO) dPress
        end if
        if (INFO /= 0) then
          INFOThermo = 44
          write (cErrMsg, '(A30,I10)') 'Cannot read pressure on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
        lPressure = .TRUE.
      case ('t','temp','temperature','T','Temp','Temperature')
        read(cValue,*,IOSTAT = INFO) dTemp
        if (INFO /= 0) then
          INFOThermo = 44
          print *, 'Cannot read temperature on line: ', iCounter
          return
        endif
        lTemperature = .TRUE.
      case ('x1','X1')
        iColon1 = 0
        iColon2 = 0
        iColon1 = INDEX(cValue,':')
        if (iColon1 > 0) then
          iColon2 = INDEX(cValue(iColon1+1:),':') + iColon1
          read(cValue(1:iColon1-1),*,IOSTAT = INFO) dX1lo
          if (iColon2 > iColon1) then
            read(cValue(iColon1+1:iColon2-1),*,IOSTAT = INFO) dX1hi
            read(cValue(iColon2+1:),*,IOSTAT = INFO) dDeltaX1
          else
            read(cValue(iColon1+1:),*,IOSTAT = INFO) dX1hi
            if (dX1hi >  dX1lo) dDeltaX1 =  0.1
            if (dX1hi <  dX1lo) dDeltaX1 = -0.1
            if (dX1hi == dX1lo) dDeltaX1 =  0
          end if
        else
          read(cValue,*,IOSTAT = INFO) dX1lo
          dX1hi = dX1lo
          dDeltaX1 = 0
        end if
        if (INFO /= 0) then
          INFOThermo = 44
          print *, 'Cannot read temperature on line: ', iCounter
          return
        endif
      case ('x2','X2')
        iColon1 = 0
        iColon2 = 0
        iColon1 = INDEX(cValue,':')
        if (iColon1 > 0) then
          iColon2 = INDEX(cValue(iColon1+1:),':') + iColon1
          read(cValue(1:iColon1-1),*,IOSTAT = INFO) dX2lo
          if (iColon2 > iColon1) then
            read(cValue(iColon1+1:iColon2-1),*,IOSTAT = INFO) dX2hi
            read(cValue(iColon2+1:),*,IOSTAT = INFO) dDeltaX2
          else
            read(cValue(iColon1+1:),*,IOSTAT = INFO) dX2hi
            if (dX2hi >  dX2lo) dDeltaX2 =  0.1
            if (dX2hi <  dX2lo) dDeltaX2 = -0.1
            if (dX2hi == dX2lo) dDeltaX2 =  0
          end if
        else
          read(cValue,*,IOSTAT = INFO) dX2lo
          dX2hi = dX2lo
          dDeltaX2 = 0
        end if
        if (INFO /= 0) then
          INFOThermo = 44
          print *, 'Cannot read temperature on line: ', iCounter
          return
        endif
      case ('iel','ielement','iEl','iElement')
        read(cValue,*,IOSTAT = INFO) iEl1, iEl2, iEl3
        if (INFO /= 0) then
          INFOThermo = 44
          write (cErrMsg, '(A26,I10)') 'Cannot read elements on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
        lEl = .TRUE.
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
        read(cValue,'(A)',IOSTAT = INFO) cThermoFileNameTemp
        if (INFO /= 0) then
          INFOThermo = 44
          write (cErrMsg, '(A35,I10)') 'Cannot read data filename on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
        cThermoFileName = cThermoFileNameTemp
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
      case ('stepTogether','steptogether','StepTogether','step_Together',&
        'step_together','Step_Together','step Together','step together','Step Together')
        read(cValue,*,IOSTAT = INFO) lStepTogether
        if (INFO /= 0) then
          INFOThermo = 44
          write (cErrMsg, '(A40,I10)') 'Cannot read step-together mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
      case ('writeJSON','writejson','WriteJSON','write_JSON',&
        'write_json','Write_JSON','write JSON','write json','Write JSON')
        read(cValue,*,IOSTAT = INFO) lWriteJSON
        if (INFO /= 0) then
          INFOThermo = 44
          write (cErrMsg, '(A37,I10)') 'Cannot read write JSON mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
      case ('fuzzy','fuzzy stoichiometry','fuzzystoichiometry','fuzzy_stoichiometry',&
        'Fuzzy','Fuzzy Stoichiometry','FuzzyStoichiometry','Fuzzy_Stoichiometry')
        read(cValue,*,IOSTAT = INFO) lFuzzyStoich
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A47,I10)') 'Cannot read fuzzy stoichiometry mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('fuzzy mag','fuzzy magnitude','fuzzymagnitude','fuzzy_magnitude',&
        'Fuzzy Mag','Fuzzy Magnitude','FuzzyMagnitude','Fuzzy_Magnitude')
        read(cValue,*,IOSTAT = INFO) dFuzzMag
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A52,I10)') 'Cannot read fuzzy stoichiometry magnitude on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('gibbs', 'gibbs min', 'gibbs minimum', 'gibbs_minimum', 'gibbsmin',&
        'Gibbs', 'Gibbs Min', 'Gibbs Minimum', 'Gibbs_Minimum', 'GibbsMin')
        read(cValue,*,IOSTAT = INFO) lGibbsMinCheck
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A46,I10)') 'Cannot read Gibbs energy check mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case default
        write (cErrMsg, '(A34,I10)') 'Input tag not recognized on line: ', iCounter
        print *,  trim(cErrMsg)
    endselect
  enddo LOOP_ReadFile

  ! Now check that all required variables have been set
  if (.NOT. lPressure) then
    INFOThermo = 45
    print *, 'Pressure not set'
    return
  endif
  if (.NOT. lTemperature) then
    INFOThermo = 45
    print *, 'Temperature not set'
    return
  endif
  if (.NOT. lEl) then
    INFOThermo = 45
    print *, 'No elements set'
    return
  endif
  if (.NOT. lPressureUnit) then
    INFOThermo = 45
    ! print *,  'Pressure unit not set'
    cInputUnitPressure = 'atm'
    return
  endif
  if (.NOT. lTemperatureUnit) then
    INFOThermo = 45
    ! print *,  'Temperature unit not set'
    cInputUnitTemperature = 'K'
    return
  endif
  if (.NOT. lMassUnit) then
    INFOThermo = 45
    ! print *,  'Mass unit not set'
    cInputUnitMass = 'moles'
    return
  endif

  select case (cInputUnitTemperature)
      case ('K')
          ! Do nothing.
      case ('C')
          dTemp  = dTemp + 273.15D0
          cInputUnitTemperature = 'K'
      case ('F')
          dTemp  = (dTemp  + 459.67D0) * (5D0/9D0)
          cInputUnitTemperature = 'K'
      case ('R')
          dTemp  = dTemp  * (5D0/9D0)
          cInputUnitTemperature = 'K'
      case default
          ! Temperature unit not recognized.
          INFOThermo = 4
          return
  end select

  select case (cInputUnitPressure)
      case ('atm')
          ! Do nothing.
      case ('psi')
          dPress  = dPress  * 0.068045957064D0
          cInputUnitPressure = 'atm'
      case ('bar')
          dPress  = dPress  * 0.98692316931D0
          cInputUnitPressure = 'atm'
      case ('Pa')
          dPress  = dPress  * 0.009869231693D0 * 1D-3
          cInputUnitPressure = 'atm'
      case ('kPa')
          dPress  = dPress  * 0.009869231693D0
          cInputUnitPressure = 'atm'
      case default
          ! Pressure unit not recognized.
          INFOThermo = 4
          return
  end select

end subroutine ParseInput3PhaseDiagram
