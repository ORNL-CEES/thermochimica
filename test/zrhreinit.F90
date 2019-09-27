
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

program zrhrestart

  USE ModuleThermoIO
  USE ModuleThermo
  USE ModuleGEMSolver

  implicit none
  integer :: i, j
  real    :: start, finish

  ! Initialize variables:
  cThermoFileName       = DATA_DIRECTORY // 'ZRHD_MHP.dat'

  ! Specify values:
  cInputUnitTemperature = 'K'
  cInputUnitPressure    = 'atm'
  cInputUnitMass        = 'moles'
  dTemperature          = 1000D0
  dPressure             = 1.0D0
  dElementMass          = 0D0
  dElementMass(1)       = 1.0D0                              ! H
  dElementMass(40)      = 0.9D0                              ! Zr
  iPrintResultsMode       = 2
  lDebugMode              = .FALSE.

  ! Parse the ChemSage data-file:
  call ParseCSDataFile(cThermoFileName)

  ! Call Thermochimica:
  call Thermochimica

  ! Print first run output
  if (iPrintResultsMode > 0)  call PrintResults

  call ThermoDebug

  LOOP_restart: do j = 1,1
    ! Save restart data
    call SaveRestartData

    call ResetThermo

    ! Re-state input variables
    dTemperature          = 2000D0
    dPressure             = 1.0D0
    dElementMass          = 0D0
    dElementMass(1)       = 1.0D0                              ! H
    dElementMass(40)      = 0.9D0                              ! Zr
    iPrintResultsMode       = 2
    lDebugMode              = .FALSE.
    ! Load restart data
    lRestartRequested = .TRUE.

    ! Call Thermochimica a bunch more for timing
    call cpu_time(start)
    LOOP_time: do i = 1,1
      call Thermochimica
    end do LOOP_time
    call cpu_time(finish)
    ! print '("Time = ",f6.3," seconds.")',finish-start

    ! Print second run output
    if (iPrintResultsMode > 0)  call PrintResults

  end do LOOP_restart

  ! Call the debugger:
  call ThermoDebug

  ! Reset Thermochimica:
  call ResetThermoAll


end program zrhrestart
