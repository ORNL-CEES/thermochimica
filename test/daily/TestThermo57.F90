
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo57.F90
    !> \brief   Spot test - Fe-Ti-V-O 2000 K.
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
    !    08/10/2020    M. Poschmann         Original code
    !    17/04/2024    A.E.F. Fitzsimmons   Naming convention change
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a system including SUBQ solution phase with compound constituent.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo57

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer :: i, j, k
    real(8) :: gibbsCheck, p1check, p2check, s1check, s2check, dHeatCapacityCheck
    logical :: subqPass, hemaPass, cppass

    ! Specify units:
    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'FeTiVO.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 2000D0
    dElementMass(8)        = 2D0              ! O
    dElementMass(22)       = 0.5D0            ! Ti
    dElementMass(23)       = 0.5D0            ! V
    dElementMass(26)       = 0.5D0            ! Fe

    gibbsCheck = -1.21336D06
    p1check    = 0.36803D0
    p2check    = 0.47842D0
    s1check    = 1.0885D-2
    s2check    = 5.6813D-2
    dHeatCapacityCheck = 186.933

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    if (INFOThermo == 0)        call Thermochimica
    call HeatCapacity

    subqPass = .FALSE.
    hemaPass = .FALSE.
    cppass   = .FALSE.
    if (INFOThermo == 0) then
        if (DABS((dGibbsEnergySys - gibbsCheck)/gibbsCheck) < 1D-3) then
            do i = 1, nSolnPhases
                k = nElements + 1 - i
                j = -iAssemblage(k)
                if (cSolnPhaseName(j) == 'SlagBsoln') then
                    if (DABS((dMolesPhase(k)-p1check)/p1check) < 1D-3) then
                        if (DABS((dMolFraction(nSpeciesPhase(j-1)+1)-s1check)/s1check) < 1D-3) then
                            subqPass = .TRUE.
                        end if
                    end if
                else if (cSolnPhaseName(j) == 'Hemasoln') then
                    if (DABS((dMolesPhase(k)-p2check)/p2check) < 1D-3) then
                        if (DABS((dMolFraction(nSpeciesPhase(j-1)+1)-s2check)/s2check) < 1D-3) then
                            hemaPass = .TRUE.
                        end if
                    end if
                end if
            end do
            if (ABS(dHeatCapacity - dHeatCapacityCheck)/dHeatCapacityCheck < 1D-3) cppass = .TRUE.
        end if
    end if

    if (subqPass .AND. hemaPass .AND. cppass) then
        ! The test passed:
        print *, 'TestThermo57: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo57: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo57
