
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo61.F90
    !> \brief   Tricky SUBL vacancy-vacancy
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
    !    24/06/2021    M. Poschmann         Original code
    !    04/17/2024    A.E.F. Fitzsimmons  Naming convention change
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica includes vacancy-vacancy
    !! species in the SUBL model if appropriate.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo40

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer :: i, j, k
    real(8) :: gibbsCheck, p1check, p2check, s1check, dHeatCapacityCheck
    logical :: subqPass, gasPass, cppass

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'FeTiVO.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1000D0
    dElementMass(8)        = 2D0              ! O
    dElementMass(22)       = 0.5D0            ! Ti

    gibbsCheck = -6.24557D05
    p1check    = 0.50000D0
    s1check    = 0.99772D0
    p2check    = 0.50114D0
    dHeatCapacityCheck = 54.8676

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica
    call HeatCapacity

    subqPass = .FALSE.
    gasPass  = .FALSE.
    cppass   = .FALSE.
    if (INFOThermo == 0) then
        if (DABS((dGibbsEnergySys - gibbsCheck)/gibbsCheck) < 1D-3) then
            do i = 1, nSolnPhases
                k = nElements + 1 - i
                j = -iAssemblage(k)
                if (cSolnPhaseName(j) == 'Rutilesoln') then
                    if (DABS((dMolesPhase(k)-p1check)/p1check) < 1D-3) then
                        if (DABS((dMolFraction(nSpeciesPhase(j-1)+1)-s1check)/s1check) < 1D-3) then
                            subqPass = .TRUE.
                        end if
                    end if
                else if (cSolnPhaseName(j) == 'gas_ideal') then
                    if (DABS((dMolesPhase(k)-p2check)/p2check) < 1D-3) then
                        gasPass = .TRUE.
                    end if
                end if
            end do
            if (ABS(dHeatCapacity - dHeatCapacityCheck)/dHeatCapacityCheck < 1D-3) cppass = .TRUE.
        end if
    end if

    if (subqPass .AND. gasPass .AND. cppass) then
        ! The test passed:
        print *, 'TestThermo40: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo40: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo40
