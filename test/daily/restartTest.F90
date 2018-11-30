
    !-------------------------------------------------------------------------------------------------------------
    !
    ! DISCLAIMER
    ! ==========
    !
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    30/11/2018    M. Poschmann        File creation.
    !
    ! Purpose:
    ! ========
    !
    ! Test implementation of restart data usage.
    !
    !-------------------------------------------------------------------------------------------------------------

program RestartTest

  USE ModuleThermoIO
  USE ModuleThermo
  USE ModuleGEMSolver

  implicit none
  integer :: i
  real    :: start, finish

  ! Initialize variables:
  dPressure             = 2D0
  dTemperature          = 900D0
  dElementMass(6)       = 1D0
  dElementMass(8)       = 1D0
  cInputUnitTemperature = 'K'
  cInputUnitPressure    = 'atm'
  cInputUnitMass        = 'moles'
  cThermoFileName       = DATA_DIRECTORY // 'C-O.dat'

  ! Specify output and debug modes:
  iPrintResultsMode     = 2
  lDebugMode            = .FALSE.

  ! Parse the ChemSage data-file:
  call ParseCSDataFile(cThermoFileName)

  ! Call Thermochimica:
  call Thermochimica
  ! Save restart data
  call SaveRestartData
  ! Load restart data
  if(lRestartAvailable) call LoadRestartData  
  ! Call Thermochimica a bunch more for timing
  call cpu_time(start)
  LOOP_time: do i = 1,1
    call Thermochimica
  end do LOOP_time
  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start

  ! Perform post-processing of results:
  if (iPrintResultsMode > 0)  call PrintResults

  ! Reset Thermochimica:
  call ResetThermo

  ! Call the debugger:
  call ThermoDebug

end program RestartTest
