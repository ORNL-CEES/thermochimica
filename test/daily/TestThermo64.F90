
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo89.F90
    !> \brief   SUBM worst case scenario.
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
    !    11/11/2022    M. Poschmann        SUBM Test Case
    !    04/17/2024    A.E.F Fitzsimmons   CsI data bug test
    !    05/06/2024    A.E.F Fitzsimmons   Naming convention update 
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results of CsI at low pressure. This test is consistent with the test in the tutorial manual.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo64

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    real(8) :: gibbscheck
    logical :: s1pass


    ! Specify units:
    cInputUnitTemperature = 'C'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'CsI-Pham.dat'

    ! Specify values:
    dTemperature          = 400
    dPressure             = 1D-5

    dElementMass(53)       = 1D0                              ! I
    dElementMass(55)      = 1D0                              ! cs


    gibbscheck = -4.41869D05

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
        print *, 'TestThermo64: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo64: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo64
