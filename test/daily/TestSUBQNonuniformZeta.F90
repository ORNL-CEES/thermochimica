!-------------------------------------------------------------------------------------------------------------
!> \file    TestSUBQNonuniformZeta.F90
!> \brief   Regression test for zeta-weighted SUBQ pair fractions.
!>
!> \details Replace the uniform zeta values in the FeTiVO SUBQ phase with
!!          deterministic nonuniform values and check the resulting equilibrium
!!          calculation.
!-------------------------------------------------------------------------------------------------------------

program TestSUBQNonuniformZeta

    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer :: i, iPhaseIndex, iSPI, j, k
    logical :: lPass, lSlagFound
    real(8), parameter :: dGibbsExpected = -1.213336160774011D6
    real(8), parameter :: dSlagAmountExpected = 3.660813286750651D-1
    real(8), parameter :: dGibbsRelativeTolerance = 1D-6
    real(8), parameter :: dSlagAmountTolerance = 1D-5
    real(8) :: dSlagAmount

    ! Define the FeTiVO state used by TestThermo57.
    cInputUnitTemperature = 'K'
    cInputUnitPressure = 'atm'
    cInputUnitMass = 'moles'
    cThermoFileName = DATA_DIRECTORY // 'FeTiVO.dat'

    dPressure = 1D0
    dTemperature = 2000D0
    dElementMass(8) = 2D0
    dElementMass(22) = 0.5D0
    dElementMass(23) = 0.5D0
    dElementMass(26) = 0.5D0

    lPass = .TRUE.
    iPhaseIndex = 0
    iSPI = 0
    call ParseCSDataFile(cThermoFileName)
    lPass = lPass .AND. (INFOThermo == 0)

    ! Locate the SUBQ phase and install discriminating nonuniform zeta values.
    if (INFOThermo == 0) then
        do i = 1, nSolnPhasesSysCS
            if (TRIM(cSolnPhaseNameCS(i)) == 'SlagBsoln') then
                iPhaseIndex = i
                iSPI = iPhaseSublatticeCS(i)
                exit
            end if
        end do
        lPass = lPass .AND. (iPhaseIndex > 0)
        lPass = lPass .AND. (iSPI > 0)
    end if

    if (iSPI > 0) then
        lPass = lPass .AND. (TRIM(cSolnPhaseTypeCS(iPhaseIndex)) == 'SUBQ')
        lPass = lPass .AND. (nPairsSROCS(iSPI,1) == 5)
        if (nPairsSROCS(iSPI,1) == 5) then
            lPass = lPass .AND. ALL(DABS(dZetaSpeciesCS(iSPI,1:5)-2.4D0) < 1D-12)
            do i = 1, 5
                dZetaSpeciesCS(iSPI,i) = 1.8D0 + 0.15D0*DFLOAT(i)
            end do
        end if
    end if

    if (lPass) then
        call Thermochimica
        lPass = lPass .AND. (INFOThermo == 0)
    end if

    ! Extract the active Slag phase amount.
    dSlagAmount = 0D0
    lSlagFound = .FALSE.
    if (lPass) then
        do i = 1, nSolnPhases
            k = nElements + 1 - i
            j = -iAssemblage(k)
            if (TRIM(cSolnPhaseName(j)) == 'SlagBsoln') then
                lSlagFound = .TRUE.
                dSlagAmount = dMolesPhase(k)
                exit
            end if
        end do
        lPass = lPass .AND. lSlagFound
    end if

    ! Check the two outputs that distinguish the corrected calculation.
    if (lPass) then
        lPass = lPass .AND. &
            (DABS((dGibbsEnergySys-dGibbsExpected)/dGibbsExpected) < dGibbsRelativeTolerance)
        lPass = lPass .AND. (DABS(dSlagAmount-dSlagAmountExpected) < dSlagAmountTolerance)
    end if

    if (lPass) then
        print *, 'TestSUBQNonuniformZeta: PASS'
        call ResetThermoAll
        call EXIT(0)
    else
        if (INFOThermo == 0) then
            write(*,'(A,2ES24.15)') 'Gibbs actual/expected = ',dGibbsEnergySys,dGibbsExpected
            write(*,'(A,2ES24.15)') 'Slag amount actual/expected = ',dSlagAmount,dSlagAmountExpected
        end if
        print *, 'TestSUBQNonuniformZeta: FAIL <---'
        call ResetThermoAll
        call EXIT(1)
    end if

end program TestSUBQNonuniformZeta
