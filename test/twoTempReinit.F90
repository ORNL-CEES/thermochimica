
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

program twotemprestart

  USE ModuleThermoIO
  USE ModuleThermo
  USE ModuleGEMSolver

  implicit none
  integer :: i, j
  real    :: start, finish

  ! Initialize variables:
  dTemperature            = 2073.15D0
  dPressure               = 1D0
  dElementMass(94)        = 0D0
  dElementMass(93)        = 0D0
  dElementMass(92)        = 0.5*0.56849D-2
  dElementMass(60)        = 0.18974d-6
  dElementMass(59)        = 0D0
  dElementMass(58)        = 0.75956d-7
  dElementMass(57)        = 3*0.3831d-7
  dElementMass(56)        = 0.99634d-7
  dElementMass(55)        = 0.12584d-8
  dElementMass(54)        = 0D0
  dElementMass(53)        = 0D0
  dElementMass(52)        = 0D0
  dElementMass(46)        = 0.10003d-4
  dElementMass(45)        = 2*0.73875d-6
  dElementMass(44)        = 0.42121d-5
  dElementMass(43)        = 0D0
  dElementMass(42)        = 0.14263d-4
  dElementMass(40)        = 0.11667d-6
  dElementMass(39)        = 0D0
  dElementMass(38)        = 0.86762d-7
  dElementMass(37)        = 0D0
  dElementMass(8)         = 2*1.14467207D-2
  dElementMass(1)         = 0D0
  cInputUnitTemperature   = 'K'
  cInputUnitPressure      = 'atm'
  cInputUnitMass          = 'moles'
  cThermoFileName         = DATA_DIRECTORY // 'Example7c.dat'
  iPrintResultsMode       = 2
  lDebugMode              = .FALSE.

  ! Parse the ChemSage data-file:
  call ParseCSDataFile(cThermoFileName)

  ! Call Thermochimica:
  call Thermochimica

  ! Print first run output
  if (iPrintResultsMode > 0)  call PrintResults

  LOOP_restart: do j = 1,1
    ! Save restart data
    call SaveRestartData

    call ResetThermo

    ! Re-state input variables
    dTemperature            = 1173.15D0
    dElementMass(94)        = 0D0
    dElementMass(93)        = 0D0
    dElementMass(92)        = 0.56849D-2
    dElementMass(60)        = 0.18974d-6
    dElementMass(59)        = 0D0
    dElementMass(58)        = 0.75956d-7
    dElementMass(57)        = 0.3831d-7
    dElementMass(56)        = 0.99634d-7
    dElementMass(55)        = 0.12584d-8
    dElementMass(54)        = 0D0
    dElementMass(53)        = 0D0
    dElementMass(52)        = 0D0
    dElementMass(46)        = 0.10003d-4
    dElementMass(45)        = 0.73875d-6
    dElementMass(44)        = 0.42121d-5
    dElementMass(43)        = 0D0
    dElementMass(42)        = 0.14263d-4
    dElementMass(40)        = 0.11667d-6
    dElementMass(39)        = 0D0
    dElementMass(38)        = 0.86762d-7
    dElementMass(37)        = 0D0
    dElementMass(8)         = 1.14467207D-2
    dElementMass(1)         = 0D0
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
  ! Reset Thermochimica:
  call ResetThermoAll

  ! Call the debugger:
  call ThermoDebug

end program twotemprestart
