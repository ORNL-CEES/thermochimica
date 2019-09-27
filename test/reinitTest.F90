
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
    ! Test implementation of reinit data usage.
    !
    !-------------------------------------------------------------------------------------------------------------

program ReinitTest

  USE ModuleThermoIO
  USE ModuleThermo
  USE ModuleGEMSolver

  implicit none
  integer :: i, j
  real(8)    :: start, finish, gibbs

  ! Initialize variables:
  dPressure             = 2D0
  dTemperature          = 1900D0
  dElementMass(74)      = 20D0        ! W
  dElementMass(79)      = 2D0         ! Au
  dElementMass(18)      = 7D0         ! Ar
  dElementMass(8)       = 5D0         ! O
  dElementMass(10)      = 1D0         ! Ne
  cInputUnitTemperature = 'K'
  cInputUnitPressure    = 'atm'
  cInputUnitMass        = 'moles'
  cThermoFileName       = DATA_DIRECTORY // 'W-Au-Ar-Ne-O_04.dat'

  ! Specify output and debug modes:
  iPrintResultsMode     = 2
  lDebugMode            = .FALSE.

  ! Parse the ChemSage data-file:
  call ParseCSDataFile(cThermoFileName)

  ! Call Thermochimica:
  call Thermochimica

  ! Print first run output
  if (iPrintResultsMode > 0)  call PrintResults

  LOOP_reinit: do j = 1,1
    ! Save reinit data
    call SaveReinitData

    ! Re-state input variables
    dPressure              = 2D0
    dTemperature           = 650D0
    dElementMass(74)       = 20D0        ! W
    dElementMass(79)       = 2D0         ! Au
    dElementMass(18)       = 7D0         ! Ar
    dElementMass(8)        = 5D0         ! O
    dElementMass(10)       = 1D0         ! Ne
    ! Load reinit data
    lReinitRequested = .FALSE.

    ! Call Thermochimica a bunch more for timing
    call cpu_time(start)
    LOOP_time: do i = 1,1
      call Thermochimica
    end do LOOP_time
    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start

    ! Print second run output
    if (iPrintResultsMode > 0)  call PrintResults

    call GibbsEnergyOfReinitData(gibbs)

  end do LOOP_reinit
  ! Reset Thermochimica:
  call ResetThermo

  ! Call the debugger:
  call ThermoDebug

end program ReinitTest
