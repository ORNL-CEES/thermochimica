!> \file TestThermo91.F90
!> \brief Indexed thermodynamic result API checks.

program TestThermo91

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    interface
        subroutine GetSpeciesPotential(iPhase, iSpecies, dPotential, INFO) &
            bind(C, name="TCAPI_getSpeciesChemicalPotentialByIndex")
            import C_DOUBLE, C_INT
            integer(C_INT), intent(in)  :: iPhase, iSpecies
            real(C_DOUBLE), intent(out) :: dPotential
            integer(C_INT), intent(out) :: INFO
        end subroutine GetSpeciesPotential

        subroutine GetEndmemberPotential(iPhase, iEndmember, dPotential, INFO) &
            bind(C, name="TCAPI_getMqmqaEndmemberStoichiometricPotentialByIndex")
            import C_DOUBLE, C_INT
            integer(C_INT), intent(in)  :: iPhase, iEndmember
            real(C_DOUBLE), intent(out) :: dPotential
            integer(C_INT), intent(out) :: INFO
        end subroutine GetEndmemberPotential

        subroutine GetPhaseGibbsEnergy(iPhase, dTotal, dMolar, INFO) &
            bind(C, name="TCAPI_getPhaseGibbsEnergyBySystemIndex")
            import C_DOUBLE, C_INT
            integer(C_INT), intent(in)  :: iPhase
            real(C_DOUBLE), intent(out) :: dTotal, dMolar
            integer(C_INT), intent(out) :: INFO
        end subroutine GetPhaseGibbsEnergy

        subroutine GetPhaseDrivingForce(iPhase, dDrivingForce, INFO) &
            bind(C, name="TCAPI_getPhaseDrivingForceBySystemIndex")
            import C_DOUBLE, C_INT
            integer(C_INT), intent(in)  :: iPhase
            real(C_DOUBLE), intent(out) :: dDrivingForce
            integer(C_INT), intent(out) :: INFO
        end subroutine GetPhaseDrivingForce

        subroutine GetSystemGibbsEnergy(dGibbsEnergy, INFO) bind(C, name="TCAPI_getSystemGibbsEnergy")
            import C_DOUBLE, C_INT
            real(C_DOUBLE), intent(out) :: dGibbsEnergy
            integer(C_INT), intent(out) :: INFO
        end subroutine GetSystemGibbsEnergy
    end interface

    integer(C_INT) :: i, iFirst, iLast, iSpeciesLocal, iTargetPhase, iTargetSpecies, iSublattice, INFO
    real(C_DOUBLE) :: dActual, dAtomSum, dExpected, dFractionSum, dMolar, dPhaseMoles
    real(C_DOUBLE) :: dReferencePotential, dResidual, dTotal
    logical        :: lPass

    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'
    cThermoFileName       = DATA_DIRECTORY // 'FeTiVO.dat'

    dPressure       = 1D0
    dTemperature    = 2000D0
    dElementMass    = 0D0
    dElementMass(8) = 2D0
    dElementMass(22) = 0.5D0
    dElementMass(23) = 0.5D0
    dElementMass(26) = 0.5D0

    call ParseCSDataFile(cThermoFileName)
    call Thermochimica

    lPass = INFOThermo == 0
    iTargetPhase = 0
    do i = 1, nSolnPhasesSys
        if (cSolnPhaseName(i) == 'SlagBsoln') iTargetPhase = i
    end do
    lPass = lPass .AND. iTargetPhase > 0

    call GetSystemGibbsEnergy(dActual, INFO)
    lPass = lPass .AND. INFO == 0 .AND. Close(dActual, dGibbsEnergySys)

    iTargetSpecies = 1
    call GetSpeciesPotential(iTargetPhase, iTargetSpecies, dActual, INFO)
    dExpected = dChemicalPotential(nSpeciesPhase(iTargetPhase - 1) + iTargetSpecies) * &
                dIdealConstant * dTemperature
    lPass = lPass .AND. INFO == 0 .AND. Close(dActual, dExpected)

    iSublattice = iPhaseSublattice(iTargetPhase)
    call GetEndmemberPotential(iTargetPhase, iTargetSpecies, dActual, INFO)
    dExpected = SUM(dElementPotential(1:nElements) * &
                    dStoichPairs(iSublattice, iTargetSpecies, 1:nElements))
    dExpected = dExpected * dIdealConstant * dTemperature
    lPass = lPass .AND. INFO == 0 .AND. Close(dActual, dExpected)

    iFirst = nSpeciesPhase(iTargetPhase - 1) + 1
    iLast = nSpeciesPhase(iTargetPhase)
    dFractionSum = SUM(dMolFraction(iFirst:iLast))
    dExpected = SUM(dChemicalPotential(iFirst:iLast) * dMolFraction(iFirst:iLast)) / dFractionSum
    dExpected = dExpected * dIdealConstant * dTemperature

    dPhaseMoles = 0D0
    do i = 1, nElements
        if (iAssemblage(i) == -iTargetPhase) dPhaseMoles = dMolesPhase(i)
    end do
    call GetPhaseGibbsEnergy(iTargetPhase, dTotal, dMolar, INFO)
    lPass = lPass .AND. INFO == 0 .AND. Close(dMolar, dExpected)
    lPass = lPass .AND. Close(dTotal, dExpected * dPhaseMoles)

    dAtomSum = 0D0
    dResidual = 0D0
    do iSpeciesLocal = iFirst, iLast
        if (dMolFraction(iSpeciesLocal) > 0D0) then
            dReferencePotential = SUM(dElementPotential(1:nElements) * &
                                      dStoichSpecies(iSpeciesLocal, 1:nElements))
            dResidual = dResidual + dMolFraction(iSpeciesLocal) * &
                        (dChemicalPotential(iSpeciesLocal) - dReferencePotential)
            dAtomSum = dAtomSum + dMolFraction(iSpeciesLocal) * dSpeciesTotalAtoms(iSpeciesLocal)
        end if
    end do
    dExpected = dResidual / dAtomSum * dIdealConstant * dTemperature
    call GetPhaseDrivingForce(iTargetPhase, dActual, INFO)
    lPass = lPass .AND. INFO == 0 .AND. Close(dActual, dExpected)

    i = 0
    call GetPhaseGibbsEnergy(i, dTotal, dMolar, INFO)
    lPass = lPass .AND. INFO == 2

    if (lPass) then
        print *, 'TestThermo91: PASS'
        call ResetThermoAll
        call EXIT(0)
    else
        print *, 'TestThermo91: FAIL <---'
        call ResetThermoAll
        call EXIT(1)
    end if

contains

    logical function Close(dLeft, dRight)
        real(C_DOUBLE), intent(in) :: dLeft, dRight
        real(C_DOUBLE)             :: dScale

        dScale = MAX(1D0, ABS(dLeft), ABS(dRight))
        Close = ABS(dLeft - dRight) <= 1D-12 * dScale
    end function Close

end program TestThermo91
