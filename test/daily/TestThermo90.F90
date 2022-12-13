
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo90.F90
    !> \brief   MQMQA reciprocal excess terms.
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
    !    12/13/2022    M. Poschmann        SUBM Test Case
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for a SUBQ phase with one center and one corner excess mixing term.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo90

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: gibbscheck
    logical :: s1pass


    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'reciprocal-test.dat'

    ! Specify values:
    dTemperature          = 1000
    dPressure             = 1.0D0
    dElementMass          = 0D0

    dTemperature           = 1000D0
    dElementMass(3)        = 0.5D0           ! Li
    dElementMass(11)       = 0.5D0           ! Na
    ! dElementMass(19)       = 0.2D0           ! K
    dElementMass(9)        = 0.5D0           ! F
    dElementMass(17)       = 0.5D0           ! Cl
    ! dElementMass(53)       = 0.8D0           ! I

    gibbscheck = -7.70822D03

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
        print *, 'TestThermo90: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo90: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo90
