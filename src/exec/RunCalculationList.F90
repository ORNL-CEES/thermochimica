program RunCalculationList

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS
    implicit none
#ifdef USE_MPI
  include 'mpif.h'
#endif

    character(1024) :: cInputFile
    integer :: i, j, nElIn, nCalc
    integer, dimension(:), allocatable :: iEls
    real(8), dimension(:), allocatable :: dEls
    character(:), allocatable :: cLine, cErrMsg, cTag, cValue, cElementNumber
    integer :: iDelimiterPosition, iOpenPosition, iClosePosition, iElementNumber, iEqualPosition
    character(1024) :: cLineInit, cThermoFileNameTemp, cOutputFilePathTemp
    logical :: lEnd, lPressureUnit, lTemperatureUnit, lMassUnit, lData, lEl, lNel
    character(15) :: cRunUnitTemperature, cRunUnitPressure, cRunUnitMass
    integer :: MPI_rank,MPI_size
    character(16) :: intStr
    integer :: fileCheck
    character(1024) :: fileOut
    character(3) :: integerString
#ifdef USE_MPI
    integer :: ierr
#endif
    ! Initialize INFO
    INFO = 0
    MPI_rank = 1
    MPI_size = 1
    ! lWriteJSON true by default
    lWriteJSON = .TRUE.
    ! Read input argument to get filename
    call get_command_argument(1, cInputFile)
    if (len_trim(cInputFile) == 0) then
      print *, 'No input file specified'
      call EXIT(1)
    endif
    ! Open input file
    open (UNIT = 2561, FILE = cInputFile, STATUS = 'old', ACTION = 'read', IOSTAT = INFO)
    ! Check for error on attempt to open
    if (INFO /= 0) then
      INFOThermo = 50
      print *, 'Cannot open input file ', cInputFile
      return
    endif
    ! Set default output file path
    fileOut = 'thermoout'
    
    ! Initialize for read loop
    lEnd = .FALSE.
    iCounter = 0
    ! Read all line of input file
    LOOP_ReadFile: do while (INFO == 0)
      ! Keep track of line number
      iCounter = iCounter + 1
      ! Read a line
      READ(2561,'(A)',IOSTAT = INFO) cLineInit
      ! If there was an error on read, give line number and return
      if (INFO > 0) then
        INFOThermo = 51
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
        endif
        cElementNumber = trim(adjustl(cTag((iOpenPosition + 1) : (iClosePosition - 1))))
        read(cElementNumber,*,IOSTAT = INFO) iElementNumber
        if (INFO /= 0) then
          INFOThermo = 53
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
            INFOThermo = 54
            write (cErrMsg, '(A26,I10)') 'Cannot read number of elements on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          end if
          if (allocated(iEls)) deallocate(iEls,dEls)
          allocate(iEls(nElIn),dEls(nElIn))
          lNel = .TRUE.
        case ('iel','ielement','iEl','iElement')
          if (.NOT. lNel) then
            INFOThermo = 54
            write (cErrMsg, '(A50,I10)') 'Need number of elements before element indices at ', iCounter
            print *,  trim(cErrMsg)
            return
          end if
          read(cValue,*,IOSTAT = INFO) iEls
          if (INFO /= 0) then
            INFOThermo = 54
            write (cErrMsg, '(A26,I10)') 'Cannot read elements on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          lEl = .TRUE.
        case ('ncalc','nCalc')
          read(cValue,*,IOSTAT = INFO) nCalc
          if (INFO /= 0) then
            INFOThermo = 54
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
          read(cValue,*,IOSTAT = INFO) cRunUnitPressure
          if (INFO /= 0) then
            INFOThermo = 54
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
          read(cValue,*,IOSTAT = INFO) cRunUnitTemperature
          if (INFO /= 0) then
            INFOThermo = 54
            write (cErrMsg, '(A38,I10)') 'Cannot read temperature unit on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          lTemperatureUnit = .TRUE.
        case ('m_unit','mass_unit','Mass_unit','Mass_Unit','m unit','mass unit','Mass unit','Mass Unit',&
          'm_units','mass_units','Mass_units','Mass_Units','m units','mass units','Mass units','Mass Units')
          read(cValue,*,IOSTAT = INFO) cRunUnitMass
          if (INFO /= 0) then
            INFOThermo = 54
            write (cErrMsg, '(A31,I10)') 'Cannot read mass unit on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
          lMassUnit = .TRUE.
        case('output','output file','output_file','Output File','Output file','output File', 'path','output path',&
          'Output Path','Output path','outputpath','outputfile','filepath','Output_File','out','json','JSON','JSON File',&
          'jsonout','JSON out','JSON output file')
          read(cValue,'(A)',IOSTAT = INFO) cOutputFilePathTemp
          if(INFO /= 0) then
            INFOThermo = 54
            write(cErrMsg,'(A35,I10)') 'Cannot read output file on line', iCounter
            return
          endif
          fileOut = cOutputFilePathTemp
        case ('data','Data','data_file','Data_file','data file','Data file','Data File',&
          'dat','Dat','dat_file','Dat_file','dat file','Dat file','Dat File')
          read(cValue,'(A)',IOSTAT = INFO) cThermoFileNameTemp
          if (INFO /= 0) then
            INFOThermo = 54
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
            INFOThermo = 54
            write (cErrMsg, '(A32,I10)') 'Cannot read print mode on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
        case ('debug_mode','Debug_mode','Debug_Mode',&
          'debug mode','Debug mode','Debug Mode')
          read(cValue,*,IOSTAT = INFO) lDebugMode
          if (INFO /= 0) then
            INFOThermo = 54
            write (cErrMsg, '(A32,I10)') 'Cannot read debug mode on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
        case ('reinit','Reinit','reinitialization','Reinitialization',&
          'reinit_mode','Reinit_mode','Reinit_Mode','reinitialization_mode','Reinitialization_mode','Reinitialization_Mode',&
          'reinit mode','Reinit mode','Reinit Mode','reinitialization mode','Reinitialization mode','Reinitialization Mode')
          read(cValue,*,IOSTAT = INFO) lReinitRequested
          if (INFO /= 0) then
            INFOThermo = 54
            write (cErrMsg, '(A43,I10)') 'Cannot read reinitialization mode on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          endif
        case ('heat capacity','entropy','enthalpy','Heat Capacity','Entropy','Enthalpy',&
          'heatCapacityEntropyEnthalpy','HeatCapacityEntropyEnthalpy')
          read(cValue,*,IOSTAT = INFO) lHeatCapacityEntropyEnthalpy
          if (INFO /= 0) then
            INFOThermo = 54
            write (cErrMsg, '(A43,I10)') 'Cannot read heat capacity / entropy / enthalpy mode on line: ', iCounter
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
            write (cErrMsg, '(A65,I10)') 'Need (nonzero) number of phases excluded before phase list at ', iCounter
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
        case ('writeJSON','writejson','WriteJSON','write_JSON',&
          'write_json','Write_JSON','write JSON','write json','Write JSON')
          read(cValue,*,IOSTAT = INFO) lWriteJSON
          if (INFO /= 0) then
            INFOThermo = 54
            write (cErrMsg, '(A38,I10)') 'Cannot read write JSON mode on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          end if
        case ('fuzzy','fuzzy stoichiometry','fuzzystoichiometry','fuzzy_stoichiometry',&
          'Fuzzy','Fuzzy Stoichiometry','FuzzyStoichiometry','Fuzzy_Stoichiometry',&
          'lFuzzyStoich','fuzz','Fuzz','fuzzy stoich')
          read(cValue,*,IOSTAT = INFO) lFuzzyStoich
          if (INFO /= 0) then
            INFOThermo = 54
            write (cErrMsg, '(A47,I10)') 'Cannot read fuzzy stoichiometry mode on line: ', iCounter
            print *,  trim(cErrMsg)
            return
          end if
        case ('fuzzy mag','fuzzy magnitude','fuzzymagnitude','fuzzy_magnitude',&
          'Fuzzy Mag','Fuzzy Magnitude','FuzzyMagnitude','Fuzzy_Magnitude',&
          'dFuzzMag','fuzz magnitude','Fuzz Magnitude')
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
          write (cErrMsg, '(A35,I10)') 'Input tag not recognized on line: ', iCounter
          print *,  trim(cErrMsg)
      endselect
    end do LOOP_ReadFile

    ! Now check that all required variables have been set
    if (.NOT. lNel) then
      INFOThermo = 55
      print *, 'Number of elements not set'
      return
    endif
    if (.NOT. lEl) then
      INFOThermo = 55
      print *, 'No elements set'
      return
    endif
    if (.NOT. lPressureUnit) then
      INFOThermo = 55
      cRunUnitPressure = 'atm'
      return
    endif
    if (.NOT. lTemperatureUnit) then
      INFOThermo = 55
      cRunUnitTemperature = 'K'
      return
    endif
    if (.NOT. lMassUnit) then
      INFOThermo = 55
      cRunUnitMass = 'moles'
      return
    endif
    ! print *, USE_MPI +
    call ParseCSDataFile(cThermoFileName)
    ! Specify values:
#ifdef USE_MPI
    if(USE_MPI > 0) then
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,MPI_rank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,MPI_size,ierr)
    ! print *, "MPI Size: ", MPI_size
    ! print *, "MPI Rank: ", MPI_rank
    endif
#endif
    ! print *, MPI_rank
    ! print *, integerString
    write(integerString,'(I0)') MPI_rank
    cOutputFilePath = trim(DATA_DIRECTORY) // '../outputs/' // trim(fileOut) // '_' // trim(adjustl(integerString)) // '.json'
    if (lWriteJSON) then
      do i = 0, MPI_size - 1
        OPEN(2 + i, file= cOutputFilePath, &
            status='REPLACE', action='write')
        WRITE(2+i,*) '{'
        CLOSE(2+i)
      end do
    end if
    
    do i = 0, nCalc - 1
      if (modulo(i,MPI_size) /= MPI_rank) CYCLE
      cInputUnitPressure = cRunUnitPressure
      cInputUnitTemperature = cRunUnitTemperature
      cInputUnitMass = cRunUnitMass
      READ(2561,*,IOSTAT = INFO) dTemperature, dPressure, dEls
      dElementMass = 0D0
      do j = 1, nElIn
        dElementMass(iEls(j)) = dEls(j)
      end do
      call Thermochimica
      call PrintResults
      if (iPrintResultsMode > 0) call ThermoDebug
      open(2+fileCheck, file= cOutputFilePath, &
          status='OLD', position='append', action='write')
      if (i > 1) write(2+fileCheck,*) ','
      write(intStr,*) i + 1
      write(2+fileCheck,*) '"', TRIM(ADJUSTL(intStr)) ,'":'
      close (2+fileCheck)
      if (lWriteJSON) then
          call WriteJSON(.TRUE.)
      end if
      ! Reset Thermochimica:
      if (INFOThermo == 0) then
          call ResetThermo
      else
          call ResetThermoAll
          INFOThermo = 0
          call ParseCSDataFile(cThermoFileName)
      end if
    end do
    CLOSE(2561)

    if (lWriteJSON) then
      do i = 0, MPI_size-1
        open(2+i, file= cOutputFilePath, &
            status='OLD', position='append', action='write')
        write(2+i,*) '}'
        close (2+i)
      end do
    end if
#ifdef USE_MPI    
    if(USE_MPI > 0) then
      call MPI_FINALIZE(ierr)
    endif
#endif
end program RunCalculationList