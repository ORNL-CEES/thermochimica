

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
    ! ternary system representing Pu-U-O system developed by C. Gueneau et al:
    !
    !   C. Gueneau, N. Dupin, B. Sundman, C. Martial, J.-C. Dumas, S. Gosse, S. Chatain, F. De Bruycker,
    !   D. Manara and R.J.M. Konings, "Thermodynamic Modelling of Advanced Oxide and Carbide Nuclear Fuels:
    !   Description of the U-Pu-O-C Systems," Journal of Nuclear Materials, 419 (2011) 145-167.
    !
    ! This particular test considers a particular combination of temperature, pressure and composition that can
    ! be compared with Figure 12b from the figure above.  This particular point was also chosen because there is
    ! experimental data that supports the model.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo34

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = '../data/Pu_U_O_CEA.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1073D0
    dElementMass           = 0d0
    dElementMass(8)        = 2.05D0     ! O
    dElementMass(92)       = 0.7D0      ! U
    dElementMass(94)       = 0.3D0      ! Pu

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(2) - 1.482D-10)/1.482D-10 < 1D-3).AND. &
        (DABS(dMolFraction(17) - 0.28497D0)/0.28497D0 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-1.22662D6))/(-1.22662D6) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo34: PASS'
        else
            ! The test failed.
            print *, 'TestThermo34: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo34: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo34
