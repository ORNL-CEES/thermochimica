
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo84.F90
    !> \brief   SUBM second sublattice ternary mixing.
    !> \author  M.H.A. Piro, M. Poschmann
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    05/14/2013    M.H.A. Piro         Original code
    !    11/04/2022    M. Poschmann        SUBM Test Case
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results a SUBM phase with a ternary excess mixing term on the second sublattice.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo84

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: gibbscheck
    logical :: s1pass


    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'subm-2sublatticeTest-mixing.dat'

    ! Specify values:
    dTemperature          = 1000
    dPressure             = 1.0D0
    dElementMass          = 0D0

    ! dElementMass(3)       = 1.2D0                              ! Li
    ! dElementMass(11)      = 0.4D0                              ! Na
    dElementMass(17)      = 0.5D0                              ! Cl
    dElementMass(9)       = 0.3D0                              ! F
    ! dElementMass(26)      = 0.3D0                              ! Fe
    dElementMass(8)       = 0.2D0                              ! O
    dElementMass(19)      = 1.2D0                              ! K

    gibbscheck = -1.76944D04

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    s1pass = .FALSE.

    ! Check results:
    if (INFOThermo == 0) then
        if ((DABS((dGibbsEnergySys - (gibbscheck))/(gibbscheck)) < 1D-3)) then
            s1pass = .TRUE.
        end if
    end if

    if (s1pass) then
        ! The test passed:
        print *, 'TestThermo84: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo84: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo84
