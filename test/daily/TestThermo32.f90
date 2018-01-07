

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
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    05/14/2013    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that Thermochimica computes the correct results for a
    ! ternary system representing U-Gd-O system developed by J. McMurray et al:
    !
    !   J.W. McMurray, D. Shin, B.W. Slone and T.M. Besmann, "Thermodynamic reassessment of U-Gd-O system,"
    !   Journal of Nuclear Materials.
    !
    ! This particular test considers a very narrow region of the UO2-GdO1.5 phase diagram where the fluorite
    ! phase is in equilibrium with a pure condensed phase.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo32

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = '../data/UGdO_ORNL.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2452D0
    dElementMass           = 0d0
    dElementMass(8)        = 9.5D0      ! O
    dElementMass(64)       = 5D0        ! Gd
    dElementMass(92)       = 1D0        ! U

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(2) - 5.0283D-11)/5.0283D-11 < 1D-3).AND. &
        (DABS(dMolFraction(29) - 0.14217D0)/0.14217D0 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-8.06205D6))/(-8.06205D6) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo32: PASS'
        else
            ! The test failed.
            print *, 'TestThermo32: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo32: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo32
