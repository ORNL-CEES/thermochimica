
subroutine RetryCalculationFirstPhase
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleReinit

    implicit none

    logical :: oldReinit

    oldReinit = lReinitRequested
    if ((INFOThermo == 12) .AND. (nSolnPhasesSys > 0)) then
        if (nSpeciesPhase(1) > 0) then
            INFOThermo = 0
            lReinitRequested = .TRUE.
            call SaveReinitData
            if (lReinitAvailable) then
                iAssemblage_Old = 0
                iAssemblage_Old(nElements) = -1
                dMolFraction_Old = 1D0 / nSpeciesPhase(1)
                dMolesPhase_Old = 0D0
                dMolesPhase_Old(nElements) = 1D0
                dChemicalPotential_Old = 0D0
                call ResetThermo
                lRetryAttempted = .TRUE.
                call Thermochimica
            end if
            lReinitRequested = oldReinit
        end if
    end if

end subroutine RetryCalculationFirstPhase
