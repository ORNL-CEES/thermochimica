
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo91.F90
    !> \brief   Unit test - hard phase fraction constraints.
    !> \author  Thermochimica contributors
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo91

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModulePhaseConstraints

    implicit none

    integer :: i, kSoln, kCon, lSoln, lCon, nReal
    real(8) :: totalElem, sumStoichSoln, sumStoichCon, fracSoln, fracCon
    logical :: pass

    pass = .TRUE.

    call ClearPhaseConstraints

    ! Initialize variables:
    dTemperature            = 1000D0
    dPressure               = 1D0
    dElementMass            = 0D0
    dElementMass(6)         = 1D0
    dElementMass(8)         = 1D0
    cInputUnitTemperature   = 'K'
    cInputUnitPressure      = 'atm'
    cInputUnitMass          = 'moles'
    cThermoFileName         = DATA_DIRECTORY // 'CO.dat'

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Constrain phases to 50/50 element-moles split:
    call AddPhaseFractionConstraint('gas_ideal', 0.5D0)
    call AddPhaseFractionConstraint('C_Graphite(s)', 0.5D0)

    ! Call Thermochimica:
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

            if (DABS(fracSoln - 0.5D0) > dTolerance(1)) pass = .FALSE.
            if (DABS(fracCon  - 0.5D0) > dTolerance(1)) pass = .FALSE.
        end if
    end if

    if (pass) then
        print *, 'TestThermo91: PASS'
        call ResetThermoAll
        call EXIT(0)
    else
        print *, 'TestThermo91: FAIL <---'
        call ResetThermoAll
        call EXIT(1)
    end if

end program TestThermo91
