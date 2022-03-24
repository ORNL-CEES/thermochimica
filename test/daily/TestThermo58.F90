
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo58.F90
    !> \brief   Spot test - Fe-Ti-O 2000 K.
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
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !! results for a system including SUBQ solution phase.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo58

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer :: i, j, k, l
    real(8) :: gibbsCheck, p1check, p2check, s1check, s2check, dHeatCapacityCheck
    logical :: subqPass, gasPass, cppass

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
    dElementMass(26)       = 0.5D0            ! Fe

    gibbsCheck = -1.00057D06
    p1check    = 1.1768D0
    p2check    = 0.14544D0
    s1check    = 5.1547D-2
    s2check    = 2.9406D-7
    dHeatCapacityCheck = 105.954

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
                if (cSolnPhaseName(j) == 'SlagBsoln') then
                    if (DABS((dMolesPhase(k)-p1check)/p1check) < 1D-3) then
                        if (DABS((dMolFraction(nSpeciesPhase(j-1)+1)-s1check)/s1check) < 1D-3) then
                            subqPass = .TRUE.
                        end if
                    end if
                else if (cSolnPhaseName(j) == 'gas_ideal') then
                    if (DABS((dMolesPhase(k)-p2check)/p2check) < 1D-3) then
                        do l = nSpeciesPhase(j-1) + 1, nSpeciesPhase(j)
                            if (DABS((dMolFraction(nSpeciesPhase(j-1)+5)-s2check)/s2check) < 1D-3) then
                                gasPass = .TRUE.
                            end if
                        end do
                    end if
                end if
            end do
            if (ABS(dHeatCapacity - dHeatCapacityCheck)/dHeatCapacityCheck < 1D-3) cppass = .TRUE.
        end if
    end if

    if (subqPass .AND. gasPass .AND. cppass) then
        ! The test passed:
        print *, 'TestThermo58: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo58: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

    ! Destruct everything:
    if (INFOThermo == 0)        call ResetThermoAll

    ! Call the debugger:
    call ThermoDebug

end program TestThermo58
