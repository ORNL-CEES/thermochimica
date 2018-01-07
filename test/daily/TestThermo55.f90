

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
    !    09/19/2015    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this application test is to ensure that the subroutines GetOutputChemPot and
    ! GetOutputSolnSpecies function correctly.  The subroutine GetOutputChemPot takes as input
    ! a character representing a chemical element and returns the chemical potential of that
    ! element.  The subroutine GetOutputSolnSpecies takes as input the name of a solution species and
    ! the name of the solution phase that species belongs to, and returns the mole fraction and chemical
    ! potential of that species.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo55

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer :: STATUS = 0

    integer       :: INFO
    character(3)  :: cElementNameRequest
    character(25) :: cSpeciesOut, cSolnOut
    real(8)       :: dElementChemPot, dMolFractionOut, dChemPotSpecies


    ! Specify units:
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = '../data/UPMZ_MHP.dat'

    ! Initialize variables:
    dTemperature            = 456D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(94)        = 0.07D0
    dElementMass(92)        = 0.9D0
    dElementMass(42)        = 2D-2
    dElementMass(40)        = 1D-2

    ! Specify output variables to be fetched:
    cElementNameRequest = 'Pu'
    cSolnOut            = 'ZETA_PHASE'
    cSpeciesOut         = 'U'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Get the chemical potential as output:
    call GetOutputChemPot(cElementNameRequest,dElementChemPot,INFO)

    ! Get solution species information as output:
    call GetOutputSolnSpecies(cSolnOut, cSpeciesOut, dMolFractionOut, dChemPotSpecies, INFO)


    !-------------------------------------------------------------------------------------------------------------
    !
    ! Commented out the following lines for the purposes of running application tests.
    ! These lines should indicate how one could use these variables.
    !
    !-------------------------------------------------------------------------------------------------------------

    !print *, ' The chemical potential of element ', cElementNameRequest, ' is ', dElementChemPot, ' [J/g-at].'
    !print *
    !print *, 'Species ', cSpeciesOut, ' in solution phase ', cSolnOut, ' has a mole fraction of ', dMolFractionOut, &
    !      ' and a chemical potential of ', dChemPotSpecies, ' [J/mol].'
    !print *

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(67) - 1.66238D-3)/1.66238D-3 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-2.54833D4))/(-2.54833D4) < 1D-3).AND.  &
        (DABS(dMolFractionOut - 0.6501D0)/(0.6501D0) < 1D-3)) then
            ! The test passed:
            print *, 'TestThermo55: PASS'
        else
            ! The test failed.
            print *, 'TestThermo55: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo55: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo55
