

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
    ! This particular test considers a unique composition, temperature and pressure that can be directly
    ! compared to a figure in this paper, which includes a series of experimental data points.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo31

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
    dTemperature           = 1773D0
    dElementMass           = 0d0
    dElementMass(8)        = 2.05D0         ! O
    dElementMass(64)       = 0.169D0        ! Gd
    dElementMass(92)       = 0.831D0        ! U

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(2) - 6.7995D-4)/6.7995D-4 < 1D-3).AND. &
        (DABS(dMolFraction(15) - 0.53694D0)/0.53694D0 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-1.38222D6))/(-1.38222D6) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo31: PASS'
        else
            ! The test failed.
            print *, 'TestThermo31: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo31: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo31
