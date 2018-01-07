

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
    !    06/10/2016    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that the subroutine GetOutputSiteFraction.f90 works
    ! correctly. The purpose of that subroutine is to get the site fraction of a particular constituent on
    ! a particular sublattice of a particular phase, which could be used by some other multi-physics
    ! code.  For example, simulating solid state diffusion in some multi-physics code may need the site
    ! fraction of interstitials of some component in a multi-subllatice phase.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo59

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0

    integer       :: INFO
    integer       :: iSublatticeOut, iConstituentOut
    character(25) :: cSolnOut
    real(8)       :: dSiteFractionOut, A, B

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = '../data/Pu_U_O_CEA.dat'

    ! Initialize variables:
    dTemperature            = 1234D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(8)         = 2.05D0
    dElementMass(92)        = 1D0

    ! Specify output variables to be fetched:
    cSolnOut          = 'O2ZRU_C'
    dSiteFractionOut  = 0D0

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Get the site fraction of vacanies on the normal anion sublattice:
    iSublatticeOut  = 2
    iConstituentOut = 2
    call GetOutputSiteFraction(cSolnOut, iSublatticeOut, iConstituentOut, dSiteFractionOut, INFO)
    A = dSiteFractionOut

    ! Get the site fraction of oxygen interstitials:
    iSublatticeOut  = 3
    iConstituentOut = 1
    call GetOutputSiteFraction(cSolnOut, iSublatticeOut, iConstituentOut, dSiteFractionOut, INFO)
    B = dSiteFractionOut

    !-------------------------------------------------------------------------------------------------------------
    !
    ! Commented out the following lines for the purposes of running application tests.
    ! These lines should indicate how one could use these variables.
    !
    !-------------------------------------------------------------------------------------------------------------

    !print *, ' The site fraction of vacanies on the normal anion sublattice is ', A, ', and the site fraction'
    !print *, ' of oxygen interstitials is ', B, '.'
    !print *

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(A - 1.939D-11)/1.939D-11 < 1D-3).AND. &
        (DABS(B - 5D-2)/5D-2 < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo59: PASS'
        else
            ! The test failed.
            print *, 'TestThermo59: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo59: FAIL <---'
        STATUS  = 1
    end if

    ! iPrintResultsMode = 2
    ! Perform post-processing of results:
    ! if (iPrintResultsMode > 0) call PrintResults

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo59
