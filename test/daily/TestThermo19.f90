
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
    !    08/22/2012    M.H.A. Piro       Original code
    !    05/13/2014    M.H.A. Piro       I changed one of the conditions for a correct solution to the activity
    !                                     of Pd(g) instead of the mole fractions of stable species.  The
    !                                     reason for doing this is because a miscibility gap is predicted to
    !                                     be stable, which means that two computational phases of the same
    !                                     physical phase are simultaneously stable.  Since there are two comp-
    !                                     utational phases of the same phase, the indexes might be different,
    !                                     but the calculation is correct.
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that Thermochimica computes the correct results for a
    ! scenario involving a miscibility gap.  This particular test considers the Pd-Rh system from Kaye et al
    ! with a miscibility gap represented by a face centred cubic (FCC) phase.
    !
    !       M.H. Kaye, B.J. Lewis and W.T. Thompson, "Thermodynamic Treatment of Noble Metal Fission
    !       Products in Nuclear Fuel," Journal of Nuclear Materials, 366 (2007) 8-27.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo18

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0


    ! Initialize variables:
    dTemperature            = 400D0
    dPressure               = 1D0
    dElementMass(46)        = 0.4D0         ! Pd
    dElementMass(45)        = 0.6D0         ! Rh

    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/example6d.dat'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica

    if ((INFOThermo == 0).AND.(DABS(dMolFraction(2) - 1.024D-43)/1.024D-43 < 1D-3)) then
        ! The unit test passed: the correct error code was reported and exited gracefully.
        print *, 'TestThermo19: PASS'
    else
        ! The unit test failed.
        print *, 'TestThermo19: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo18
