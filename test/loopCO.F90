program loopCO

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
  
    implicit none
    integer :: i
  
    ! Initialize variables:
    dTemperature            = 3.0D2
    dPressure               = 1D0
    dElementMass(6)         = 1.0d0 ! Carbon
    dElementMass(8)         = 1.0d0 ! Oxygen
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    CInputUnitMass          = 'moles'
    cThermoFileName         = DATA_DIRECTORY // 'C-O.dat'
  
    ! Specify output and debug modes:
    iPrintResultsMode     = 2
    lDebugMode            = .FALSE.
  
    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)
  
    ! Call Thermochimica:
    call Thermochimica
  
    ! Print first run output
    if (iPrintResultsMode > 0)  call PrintResults
  
    if (INFOThermo == 0) then
      LOOP_time: do i = 1,10
        call ResetThermo
        dTemperature            = 3.0D2+1.0D2*i
        call Thermochimica
        if (iPrintResultsMode > 0)  call PrintResults
      end do LOOP_time
    end if
  
    ! Reset Thermochimica:
    call ResetThermo
  
    ! Call the debugger:
    call ThermoDebug
  
end program loopCO
  