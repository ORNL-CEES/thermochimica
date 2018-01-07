

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
    ! sample problem representing the Th-U-Pu ternary system.  This particular data-file was kindly provided by
    ! Ondrej Benes from ITU (Germany), which is described in the following paper:
    !
    !   O. Benes, D. Manara and R.J.M. Konings, "Thermodynamic Assessment of the Th-U-Pu System," Journal
    !   of Nuclear Materials, 449 (2014) 15-22.
    !
    ! This particular test is similar to test 28, except that the Th-U-Pu ternary is considered (as opposed to
    ! the Th-Pu binary).  A very narrow three phase region is considered.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo29

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = '../data/ThPuU_ITU.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 900D0
    dElementMass           = 0d0
    dElementMass(90)       = 0.555D0        ! Th
    dElementMass(92)       = 0.0889D0       ! U
    dElementMass(94)       = 0.345D0        ! Pu

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The  FCC phases should be stable at equilibrium.
        if ((DABS(dMolFraction(8) - 7.5514D-3)/7.5514D-3 < 1D-3).AND. &
        (DABS(dMolFraction(6) - 7.1517D-3)/7.1517D-3 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-6.45968D4))/(-6.45968D4) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo29: PASS'
        else
            ! The test failed.
            print *, 'TestThermo29: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo29: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo29
