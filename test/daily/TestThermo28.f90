

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
    ! sample problem representing the Th-Pu binary system.  This particular data-file was kindly provided by
    ! Ondrej Benes from ITU (Germany), which is described in the following paper:
    !
    !   O. Benes, D. Manara and R.J.M. Konings, "Thermodynamic Assessment of the Th-U-Pu System," Journal
    !   of Nuclear Materials, 449 (2014) 15-22.
    !
    ! This particular test considers a very narrow region of the phase diagram, whereby the FCC and tetragonal
    ! phases coexist.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo28

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
    dTemperature           = 763D0
    dElementMass           = 0d0
    dElementMass(90)       = 0.0115D0
    dElementMass(94)       = 1D0 - dElementMass(90)

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The tetragonal & FCC phases should be stable at equilibrium.
        if ((DABS(dMolFraction(8) - 1.1183D-2)/1.1183D-2 < 1D-3).AND. &
        (DABS(dMolFraction(6) - 1.1950D-2)/1.1950D-2 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-5.40590D4))/(-5.40590D4) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo28: PASS'
        else
            ! The test failed.
            print *, 'TestThermo28: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo28: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo28
