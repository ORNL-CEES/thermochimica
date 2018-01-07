

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
    ! binary system representing Am-O developed by T.M. Besmann:
    !
    !   T.M. Besmann, "Modeling the Thermochemical Behavior of Am$_{2-x}$," Journal of Nuclear Materials,
    !   402 (2010) 25-29.
    !
    ! This particular treatment only considers the americia solid solution phase, which is represented by the
    ! compound energy formalism.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo30

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = '../data/AmO_TMB.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1455D0
    dElementMass           = 0d0
    dElementMass(8)        = 1.95D0        ! O
    dElementMass(95)       = 1D0           ! Am

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The tetragonal & FCC phases should be stable at equilibrium.
        if ((DABS(dMolFraction(2) - 3.3348D-2)/3.3348D-2 < 1D-3).AND. &
        (DABS(dMolFraction(7) - 2.25D-2)/2.25D-2 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-1.1223D6))/(-1.1223D6) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo30: PASS'
        else
            ! The test failed.
            print *, 'TestThermo30: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo30: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo30
