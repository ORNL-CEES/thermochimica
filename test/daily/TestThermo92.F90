    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo92.F90
    !> \brief   Unit test - phase fraction constraints (duplicate update + invalid fraction).
    !> \author  Thermochimica contributors
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo92

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModulePhaseConstraints

    implicit none

    integer :: i, kSoln, kCon, lSoln, lCon, nReal, infoConstraint
    real(8) :: totalElem, sumStoichSoln, sumStoichCon, fracSoln, fracCon
    real(8) :: targetConstraint, achievedConstraint, residualConstraint, lambdaConstraint
    logical :: pass

    pass = .TRUE.

    ! Case 1: Duplicate constraint should update target and remain valid.
    call ClearPhaseConstraints
    call SetupCOSystem

    call AddPhaseFractionConstraint('gas_ideal', 0.2D0)
    call AddPhaseFractionConstraint('gas_ideal', 0.7D0)
    call AddPhaseFractionConstraint('C_Graphite(s)', 0.3D0)

    call Thermochimica

    if (INFOThermo /= 0) pass = .FALSE.
    if (dGEMFunctionNorm > dTolerance(1)) pass = .FALSE.

    if (pass) then
        nReal = nElements - nChargedConstraints
        if (nReal < 1) nReal = nElements

        totalElem = 0D0
        do i = 1, nReal
            totalElem = totalElem + dMolesElement(i)
        end do

        ! Find solution phase index for gas_ideal
        kSoln = 0
        do i = 1, nSolnPhasesSys
            if (trim(adjustl(cSolnPhaseName(i))) == 'gas_ideal') then
                kSoln = i
                exit
            end if
        end do
        lSoln = 0
        if (kSoln > 0) then
            do i = nElements - nSolnPhases + 1, nElements
                if (iAssemblage(i) == -kSoln) then
                    lSoln = i
                    exit
                end if
            end do
        end if

        ! Find pure condensed phase index for graphite
        kCon = 0
        do i = 1, nSpecies
            if (iPhase(i) == 0) then
                if (trim(adjustl(cSpeciesName(i))) == 'C_Graphite(s)') then
                    kCon = i
                    exit
                end if
            end if
        end do
        lCon = 0
        if (kCon > 0) then
            do i = 1, nConPhases
                if (iAssemblage(i) == kCon) then
                    lCon = i
                    exit
                end if
            end do
        end if

        if ((kSoln == 0) .OR. (lSoln == 0) .OR. (kCon == 0) .OR. (lCon == 0)) pass = .FALSE.

        if (pass) then
            call CompStoichSolnPhase(kSoln)
            sumStoichSoln = 0D0
            do i = 1, nReal
                sumStoichSoln = sumStoichSoln + dEffStoichSolnPhase(kSoln,i)
            end do
            fracSoln = dMolesPhase(lSoln) * sumStoichSoln / totalElem

            sumStoichCon = 0D0
            do i = 1, nReal
                sumStoichCon = sumStoichCon + dStoichSpecies(kCon,i)
            end do
            fracCon = dMolesPhase(lCon) * sumStoichCon / totalElem

            if (DABS(fracSoln - 0.7D0) > dTolerance(1)) pass = .FALSE.
            if (DABS(fracCon  - 0.3D0) > dTolerance(1)) pass = .FALSE.
        end if
    end if

    ! Case 2: Invalid fraction values should be rejected.
    call ResetThermoAll
    INFOThermo = 0

    if (pass) then
        call SetupCOSystem
        call AddPhaseFractionConstraint('gas_ideal', -0.1D0)
        call AddPhaseFractionConstraint('C_Graphite(s)', 1.1D0)
        call Thermochimica

        if (INFOThermo /= 87) pass = .FALSE.
        call GetPhaseConstraintData(1, targetConstraint, achievedConstraint, residualConstraint, &
            lambdaConstraint, infoConstraint)
        if (infoConstraint /= 2) pass = .FALSE.
    end if

    ! Case 3: Error recovery must rebuild solver/database state without losing constraints.
    if (pass) then
        call ResetThermoAllPreservePhaseConstraints
        INFOThermo = 0
        call GetPhaseConstraintData(1, targetConstraint, achievedConstraint, residualConstraint, &
            lambdaConstraint, infoConstraint)
        if (infoConstraint /= 2) pass = .FALSE.
        call ParseCSDataFile(cThermoFileName)

        if (nPhaseConstraints /= 2) pass = .FALSE.
        call AddPhaseFractionConstraint('gas_ideal', 0.6D0)
        call AddPhaseFractionConstraint('C_Graphite(s)', 0.4D0)
        call Thermochimica

        if (INFOThermo /= 0) pass = .FALSE.
        if (dGEMFunctionNorm > dTolerance(1)) pass = .FALSE.
    end if

    ! Case 4: A single pure condensed phase can satisfy a phase constraint
    ! without a solution phase being present in the equilibrium assemblage.
    if (pass) then
        call ResetThermoAll
        INFOThermo = 0
        dTemperature            = 1000D0
        dPressure               = 1D0
        dElementMass            = 0D0
        dElementMass(6)         = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'
        call ParseCSDataFile(cThermoFileName)
        call AddPhaseFractionConstraint('C_Graphite(s)', 1D0)
        call Thermochimica

        if (INFOThermo /= 0) pass = .FALSE.
        if (nSolnPhases /= 0) pass = .FALSE.
        if (nConPhases /= 1) pass = .FALSE.
        call GetPhaseConstraintResult(1, achievedConstraint, residualConstraint, infoConstraint)
        if (infoConstraint /= 0) pass = .FALSE.
        if (DABS(achievedConstraint - 1D0) > dTolerance(1)) pass = .FALSE.
        if (DABS(residualConstraint) > dTolerance(1)) pass = .FALSE.
        if (dGEMFunctionNorm > dTolerance(1)) pass = .FALSE.
    end if

    if (pass) then
        print *, 'TestThermo92: PASS'
        call ResetThermoAll
        call EXIT(0)
    else
        print *, 'TestThermo92: FAIL <---'
        call ResetThermoAll
        call EXIT(1)
    end if

contains

    subroutine SetupCOSystem
        call ClearPhaseConstraints
        dTemperature            = 1000D0
        dPressure               = 1D0
        dElementMass            = 0D0
        dElementMass(6)         = 1D0
        dElementMass(8)         = 1D0
        cInputUnitTemperature   = 'K'
        cInputUnitPressure      = 'atm'
        cInputUnitMass          = 'moles'
        cThermoFileName         = DATA_DIRECTORY // 'CO.dat'
        call ParseCSDataFile(cThermoFileName)
    end subroutine SetupCOSystem

end program TestThermo92
