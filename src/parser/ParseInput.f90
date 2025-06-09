
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

subroutine ParseInput(cInputFileName,dTempLow,dTempHigh,dDeltaT,dPressLow,dPressHigh,dDeltaP)

  USE ModuleThermoIO
  USE ModuleGEMSolver
  USE ModuleParseCS

  implicit none

  character(*)                :: cInputFileName
  integer                     :: iDelimiterPosition, iOpenPosition, iClosePosition, iElementNumber, iEqualPosition
  integer                     :: iColon1, iColon2
  logical                     :: lEnd, lPressure, lTemperature, lMass, lPressureUnit, lTemperatureUnit, lMassUnit, lData
  character(:), allocatable   :: cLine, cErrMsg, cTag, cValue, cElementNumber
  character(1024)             :: cLineInit, cThermoFileNameTemp, cOutputFilePathTemp
  real(8), intent(out)        :: dTempLow, dTempHigh, dDeltaT, dPressLow, dPressHigh, dDeltaP
  !character(1024)             :: cOutputFileName Nothing to see here
  ! Initialize INFO
  INFO = 0

  ! lWriteJSON true by default
  lWriteJSON = .TRUE.
  cOutputFilePathTemp = '../outputs/thermoout.json'
  ! Open input file
  open (UNIT = 1, FILE = cInputFileName, STATUS = 'old', ACTION = 'read', IOSTAT = INFO)
  ! Check for error on attempt to open
  if (INFO /= 0) then
    INFOThermo = 50
    print *, 'Cannot open input file ' // cInputFileName
    call EXIT(1)
  end if

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
      INFOThermo = 51
      write (cErrMsg, '(A35,I10)') 'Reading input file failed on line: ', iCounter
      print *,  trim(cErrMsg)
      return
    ! If file end reached, break loop
    elseif (INFO < 0) then
      exit LOOP_ReadFile
    end if
    ! Remove leading then trailing spaces on line
    cLine = trim(adjustl(cLineInit))
    ! Check for comment line (going to be liberal with choices of comment indicators)
    if (scan(cLine,'!@#$%&*/\?|') == 1) then
      cycle LOOP_ReadFile
    end if
    ! Also look for lines without "="
    if (scan(cLine,'=') == 0) then
      cycle LOOP_ReadFile
    end if
    ! If we get to here, should be a data line, and therefore contain '=' delimiter separating tag and value terms
    iDelimiterPosition = scan(cLine,'=')
    ! Tag is on LHS of delimiter, extract this and trim whitespace
    cTag = trim(adjustl(cLine(1 : (iDelimiterPosition - 1))))
    ! Value if on RHS of delimiter, do same as above
    cValue = trim(adjustl(cLine((iDelimiterPosition + 1) : len(cLine))))
    ! Check if line contains a mass, need to treat these separately
    ! Masses will be the only lines to contain '()' on the LHS, so look for these
    iOpenPosition = scan(cLine,'(')
    iEqualPosition = scan(cLine,'=')
    if ((iOpenPosition > 0) .AND. (iOpenPosition < iEqualPosition)) then
      iClosePosition = scan(cLine,')')
      ! Check for no close ')'
      if (iClosePosition == 0) then
        INFOThermo = 52
        write (cErrMsg, '(A31,I10)') 'Open ( but no close ) on line: ', iCounter
        print *,  trim(cErrMsg)
        return
      end if
      cElementNumber = trim(adjustl(cTag((iOpenPosition + 1) : (iClosePosition - 1))))
      read(cElementNumber,*,IOSTAT = INFO) iElementNumber
      if (INFO /= 0) then
        INFOThermo = 53
        write (cErrMsg, '(A36,I10)') 'Cannot read element number on line: ', iCounter
        print *,  trim(cErrMsg)
        return
      end if
      cTag = trim(adjustl(cTag(1 : (iOpenPosition - 1))))
    end if

    ! Now look through possible tags to assign variables
    select case (cTag)
      case ('p','P','pressure','press','Pressure','Press')
        iColon1 = 0
        iColon2 = 0
        iColon1 = INDEX(cValue,':')
        if (iColon1 > 0) then
          iColon2 = INDEX(cValue(iColon1+1:),':') + iColon1
          read(cValue(1:iColon1-1),*,IOSTAT = INFO) dPressLow
          if (iColon2 > iColon1) then
            read(cValue(iColon1+1:iColon2-1),*,IOSTAT = INFO) dPressHigh
            read(cValue(iColon2+1:),*,IOSTAT = INFO) dDeltaP
          else
            read(cValue(iColon1+1:),*,IOSTAT = INFO) dPressHigh
            if (dPressHigh >  dPressLow) dDeltaP =  1
            if (dPressHigh <  dPressLow) dDeltaP = -1
            if (dPressHigh == dPressLow) dDeltaP =  0
          end if
        else
          read(cValue,*,IOSTAT = INFO) dPressLow
          dPressHigh = dPressLow
          dDeltaP = 0
        end if
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A30,I10)') 'Cannot read pressure on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
        lPressure = .TRUE.
      case ('t','temp','temperature','T','Temp','Temperature')
        iColon1 = 0
        iColon2 = 0
        iColon1 = INDEX(cValue,':')
        if (iColon1 > 0) then
          iColon2 = INDEX(cValue(iColon1+1:),':') + iColon1
          read(cValue(1:iColon1-1),*,IOSTAT = INFO) dTempLow
          if (iColon2 > iColon1) then
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
          INFOThermo = 54
          print *, 'Cannot read temperature on line: ', iCounter
          return
        end if
        lTemperature = .TRUE.
      case ('m','mass','M','Mass')
        read(cValue,*,IOSTAT = INFO) dElementMass(iElementNumber)
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A26,I10)') 'Cannot read mass on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
        lMass = .TRUE.
      case ('p_unit','pressure_unit','Pressure_unit','Pressure_Unit','P_unit','pressure unit',&
        'Pressure Unit','Pressure unit','press_unit','press unit','Press unit','Press Unit',&
        'p unit','P unit','P Unit',&
        'p_units','pressure_units','Pressure_units','Pressure_Units','P_units','P_Units','pressure units',&
          'Pressure Units','Pressure units','press_units','press units','Press units','Press Units'&
          'p units','P units','P Units')
        read(cValue,*,IOSTAT = INFO) cInputUnitPressure
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A35,I10)') 'Cannot read pressure unit on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
        lPressureUnit = .TRUE.
      case ('t_unit','temperature_unit','Temperature_unit','Temperature_Unit','T_unit','T_Unit',&
        'temperature unit','Temperature Unit','Temperature unit','temp_unit','temp unit',&
        'Temp unit','Temp Unit','t unit','T unit','T Unit',&
        't_units','temperature_units','Temperature_units','Temperature_Units','T_units','T_Units',&
          'temperature units','Temperature Units','Temperature units','temp_units','temp units',&
          'Temp units','Temp Units','t units','T units','T Units')
        read(cValue,*,IOSTAT = INFO) cInputUnitTemperature
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A38,I10)') 'Cannot read temperature unit on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
        lTemperatureUnit = .TRUE.
      case ('m_unit','mass_unit','Mass_unit','Mass_Unit','m unit','mass unit','Mass unit','Mass Unit',&
        'm_units','mass_units','Mass_units','Mass_Units','m units','mass units','Mass units','Mass Units')
        read(cValue,*,IOSTAT = INFO) cInputUnitMass
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A31,I10)') 'Cannot read mass unit on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
        lMassUnit = .TRUE.
      case ('output_file','output file','Output File','Output file')
        read(cValue,'(A)',IOSTAT=INFO) cOutputFilePathTemp
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A35,I10)') 'That is some crap buddy', iCounter !Need to create an error here, temp is the same as data
          print *, trim(cErrMsg)
          return
        end if
        cOutputFilePath = cOutputFilePathTemp
      case ('data','Data','data_file','Data_file','data file','Data file','Data File',&
        'dat','Dat','dat_file','Dat_file','dat file','Dat file','Dat File')
        read(cValue,'(A)',IOSTAT = INFO) cThermoFileNameTemp
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A35,I10)') 'Cannot read data filename on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
        cThermoFileName = cThermoFileNameTemp
        lData = .TRUE.
      case ('print_mode','Print_mode','Print_Mode',&
        'print mode','Print mode','Print Mode')
        read(cValue,*,IOSTAT = INFO) iPrintResultsMode
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A32,I10)') 'Cannot read print mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('debug_mode','Debug_mode','Debug_Mode',&
        'debug mode','Debug mode','Debug Mode')
        read(cValue,*,IOSTAT = INFO) lDebugMode
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A32,I10)') 'Cannot read debug mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('reinit','Reinit','reinitialization','Reinitialization',&
        'reinit_mode','Reinit_mode','Reinit_Mode','reinitialization_mode','Reinitialization_mode','Reinitialization_Mode',&
        'reinit mode','Reinit mode','Reinit Mode','reinitialization mode','Reinitialization mode','Reinitialization Mode')
        read(cValue,*,IOSTAT = INFO) lReinitRequested
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A43,I10)') 'Cannot read reinitialization mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('heat capacity','entropy','enthalpy','Heat Capacity','Entropy','Enthalpy',&
        'heatCapacityEntropyEnthalpy','HeatCapacityEntropyEnthalpy')
        read(cValue,*,IOSTAT = INFO) lHeatCapacityEntropyEnthalpy
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A43,I10)') 'Cannot read heat capacity / entropy / enthalpy mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('stepTogether','steptogether','StepTogether','step_Together',&
        'step_together','Step_Together','step Together','step together','Step Together')
        read(cValue,*,IOSTAT = INFO) lStepTogether
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A40,I10)') 'Cannot read step-together mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('writeJSON','writejson','WriteJSON','write_JSON',&
        'write_json','Write_JSON','write JSON','write json','Write JSON')
        read(cValue,*,IOSTAT = INFO) lWriteJSON
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A37,I10)') 'Cannot read write JSON mode on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('nMinSpeciesPerPhase','species per phase','min species','minimum species per phase')
        read(cValue,*,IOSTAT = INFO) nMinSpeciesPerPhase
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A47,I10)') 'Cannot read minimum species per phase on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('nPhasesExcluded','nphasesexcluded','number of phases excluded','number excluded')
        read(cValue,*,IOSTAT = INFO) nPhasesExcluded
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A47,I10)') 'Cannot read number of phases excluded on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('phasesExcluded','phasesexcluded','phases excluded')
        if (nPhasesExcluded == 0) then
          INFOThermo = 54
          write (cErrMsg, '(A63,I10)') 'Need (nonzero) number of phases excluded before phase list at ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
        read(cValue,*,IOSTAT = INFO) cPhasesExcluded(1:nPhasesExcluded)
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A41,I10)') 'Cannot parse phases exclusions on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        endif
      case ('nPhasesExcludedExcept','nphasesexcludedexcept','number of phases excluded except','number excluded except')
        read(cValue,*,IOSTAT = INFO) nPhasesExcludedExcept
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A47,I10)') 'Cannot read number of phases excluded on line: ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
      case ('phasesExcludedExcept','phasesexcludedexcept','phases excluded except')
        if (nPhasesExcludedExcept == 0) then
          INFOThermo = 54
          write (cErrMsg, '(A75,I10)') 'Need (nonzero) number of phases exclusion exceptions before phase list at ', iCounter
          print *,  trim(cErrMsg)
          return
        end if
        read(cValue,*,IOSTAT = INFO) cPhasesExcludedExcept(1:nPhasesExcludedExcept)
        if (INFO /= 0) then
          INFOThermo = 54
          write (cErrMsg, '(A50,I10)') 'Cannot parse phase exclusion exceptions on line: ', iCounter
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
    INFOThermo = 55
    print *,  'Pressure not set'
    return
  end if
  if (.NOT. lTemperature) then
    INFOThermo = 55
    print *,  'Temperature not set'
    return
  end if
  if (.NOT. lMass) then
    INFOThermo = 55
    print *,  'No masses set'
    return
  end if
  if (.NOT. lPressureUnit) then
    INFOThermo = 55
    ! print *,  'Pressure unit not set'
    cInputUnitPressure = 'atm'
    return
  end if
  if (.NOT. lTemperatureUnit) then
    INFOThermo = 55
    ! print *,  'Temperature unit not set'
    cInputUnitTemperature = 'K'
    return
  end if
  if (.NOT. lMassUnit) then
    INFOThermo = 55
    ! print *,  'Mass unit not set'
    cInputUnitMass = 'moles'
    return
  end if

  select case (cInputUnitTemperature)
      case ('K')
          ! Do nothing.
      case ('C')
          dTempLow  = dTempLow  + 273.15D0
          dTempHigh = dTempHigh + 273.15D0
          cInputUnitTemperature = 'K'
      case ('F')
          dTempLow  = (dTempLow  + 459.67D0) * (5D0/9D0)
          dTempHigh = (dTempHigh + 459.67D0) * (5D0/9D0)
          dDeltaT   = dDeltaT * (5D0/9D0)
          cInputUnitTemperature = 'K'
      case ('R')
          dTempLow  = dTempLow  * (5D0/9D0)
          dTempHigh = dTempHigh * (5D0/9D0)
          dDeltaT   = dDeltaT * (5D0/9D0)
          cInputUnitTemperature = 'K'
      case default
          ! Temperature unit not recognized.
          INFOThermo = 5
          return
  end select

  select case (cInputUnitPressure)
      case ('atm')
          ! Do nothing.
      case ('psi')
          dPressLow  = dPressLow  * 0.068045957064D0
          dPressHigh = dPressHigh * 0.068045957064D0
          dDeltaP    = dDeltaP * 0.068045957064D0
          cInputUnitPressure = 'atm'
      case ('bar')
          dPressLow  = dPressLow  * 0.98692316931D0
          dPressHigh = dPressHigh * 0.98692316931D0
          dDeltaP    = dDeltaP * 0.98692316931D0
          cInputUnitPressure = 'atm'
      case ('Pa')
          dPressLow  = dPressLow  * 0.009869231693D0 * 1D-3
          dPressHigh = dPressHigh * 0.009869231693D0 * 1D-3
          dDeltaP    = dDeltaP * 0.009869231693D0 * 1D-3
          cInputUnitPressure = 'atm'
      case ('kPa')
          dPressLow  = dPressLow  * 0.009869231693D0
          dPressHigh = dPressHigh * 0.009869231693D0
          dDeltaP    = dDeltaP * 0.009869231693D0
          cInputUnitPressure = 'atm'
      case default
          ! Pressure unit not recognized.
          INFOThermo = 5
          return
  end select

end subroutine ParseInput
