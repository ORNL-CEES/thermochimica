

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
    ! GetOutputSolnSpecies return the correct error when given inappropriate input.  A species name
    ! is requested that is included in the database, but is not stable at equilibrium.
    !
    !-------------------------------------------------------------------------------------------------------------


program TestThermo56

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
    dTemperature            = 328D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(94)        = 360D0
    dElementMass(92)        = 635D0
    dElementMass(42)        = 395D0
    dElementMass(40)        = 587D0

    ! Specify output variables to be fetched:
    cElementNameRequest = 'Pu'
    cSolnOut            = 'HCP_A3'
    cSpeciesOut         = 'ZR'


    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0) call Thermochimica

    ! Get the chemical potential as output:
    call GetOutputChemPot(cElementNameRequest, dElementChemPot, INFO)

    ! Get solution species information as output:
    call GetOutputSolnSpecies(cSolnOut, cSpeciesOut, dMolFractionOut, dChemPotSpecies, INFO)

    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dMolFraction(69) - 0.40908D0)/0.40908D0 < 1D-3).AND. &
        (DABS(dGibbsEnergySys - (-3.92483D7))/(-3.92483D7) < 1D-3).AND.  &
        (INFO == 1)) then
            ! The test passed:
            print *, 'TestThermo56: PASS'
        else
            ! The test failed.
            print *, 'TestThermo56: FAIL <---'
            STATUS  = 1
        end if
    else
        ! The test failed.
        print *, 'TestThermo56: FAIL <---'
        STATUS  = 1
    end if

    ! Reset Thermochimica:
    call ResetThermo

    call EXIT(STATUS)
end program TestThermo56
