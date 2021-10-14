
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo64.F90
    !> \brief   Spot test - Nb-Zr-O-H 600 K.
    !> \author  M. Poschmann
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    08/31/2021    M. Poschmann         Original code
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that ternary mixing SUBL cases with only
    !!  one specified coefficient are handled correctly.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo64

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer :: i, j
    logical :: bccPass

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC-test64.dat'

    ! Specify values:
    dPressure              = 100D0
    dTemperature           = 600D0
    dElementMass(41)       = 1D0            ! Nb
    dElementMass(40)       = 1D0            ! Zr
    dElementMass(8)        = 1D0            ! O
    dElementMass(1)        = 0.1D0          ! H

    ! Specify output mode:
    iPrintResultsMode     = 2

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica

    bccPass = .FALSE.
    ! Check results:
    if (INFOThermo == 0) then
        ! The fluorite oxide phase should be the only one stable at equilibrium.
        if ((DABS(dGibbsEnergySys - (-5.24838E+05))/((-5.24838E+05))) < 1D-3) then
            do i = 1, nSolnPhases
                j = -iAssemblage(nElements + 1 - i)
                if (cSolnPhaseName(j) == 'BCC_A2') then
                    bccPass = .TRUE.
                end if
            end do
        end if
    end if

    if (bccPass) then
        ! The test passed:
        print *, 'TestThermo64: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo64: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo64
