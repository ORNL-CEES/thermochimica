    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ModuleTesting.F90
    !> \brief   Module for testing.
    !> \author  A.E.F. Fitzsimmons
    !
    !> \todo    Fix duplicates
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    04/10/2024    A.E.F. Fitzsimmons   Original Code
    !
    ! Purpose:
    ! ========
    !> \details Using a bisection search, this subroutine will find and record temperatures of which a phase has
    !! made a transition. Returning the list of transition temperatures.
    !
    !-------------------------------------------------------------------------------------------------------------
recursive subroutine PhaseTransitionrecur(dTempMin, dTempMax, dTempTolerance, dPhaseTransitionTemp, iTransitions, lFirst)
    USE ModuleThermoIO
    USE ModuleThermo
    implicit none

    real(8), intent(in) :: dTempMin, dTempMax, dTempTolerance
    real(8), dimension(10), intent(out) :: dPhaseTransitionTemp
    real(8) :: dTempUpper, dTempLower
    integer, intent(out) :: iTransitions
    integer :: i, iMaxIteration
    integer, dimension(nElements) :: iAssemblageUpper, iAssemblageLower
    integer, dimension(10,nElements) :: iPhaseTransitionAssemblage
    logical :: lBisectionConvergence, lCompareAssemblage, lCompareAssemblage2
    
    !Init
    iMaxIteration = 50
    dTempUpper = dTempMax
    dTempLower = dTempMin
    lBisectionConvergence = .FALSE.
    lCompareAssemblage = .FALSE.
    lCompareAssemblage2 = .FALSE.
    
    if(dTempMin>=dTempMax) then
        INFOThermo = 56
        return
    end if

    if(nElements < 2) then
        INFOThermo = 59
        return
    end if
    
    !Init Upper and Lower assemblage
    dTemperature = dTempUpper
    call Thermochimica 
    iAssemblageUpper = iAssemblage
    call SortInt(iAssemblageUpper, nElements)

    dTemperature = dTempLower
    call Thermochimica
    iAssemblageLower = iAssemblage
    call SortInt(iAssemblageLower, nElements)

    lCompareAssemblage = ALL(iAssemblageLower == iAssemblageUpper)

    !Check if any transition is made
    if(lCompareAssemblage) then 
        INFOThermo = 58
        return
    end if 
    
    !Bisection
    LOOP_TEMP: do i = 1, iMaxIteration
        !reinit
        lCompareAssemblage = .FALSE.
        lCompareAssemblage2 = .FALSE.
        dTemperature = (dTempUpper + dTempLower) / 2D0

        !Check if reached convergence
        if (dTempUpper-dTempLower <= dTempTolerance) then
            !Save convergence Temperature and Assemblage
            
            iTransitions = iTransitions + 1
            dPhaseTransitionTemp(iTransitions) = dTemperature
            if (iTransitions > 1) then
                if (abs(dPhaseTransitionTemp(iTransitions) - dPhaseTransitionTemp(iTransitions-1)) < dTempTolerance) then
                    iTransitions = iTransitions - 1
                end if
            end if
            
            call SortInt(iAssemblage, nElements)
            iPhaseTransitionAssemblage(iTransitions, :) = iAssemblage 
            lBisectionConvergence = .TRUE.
        !    print *, "AssemblageUpper: ", iAssemblageUpper, NEW_LINE(" "), " AssemblageLower: ", iAssemblageLower, NEW_LINE(" "), " AssemblageTransition: ", iAssemblage, NEW_LINE(" ")
            
            EXIT LOOP_TEMP
        end if

        !Midpoint calculation
        call Thermochimica
        call SortInt(iAssemblage, nElements)

        print *, "Iteration: ", i, " dTempLower: ", dTempLower, " dTempUpper: ", dTempUpper, " dTempGuess: ", dTemperature, NEW_LINE(" ")
        lCompareAssemblage = ALL(iAssemblage == iAssemblageUpper)
        lCompareAssemblage2 = ALL(iAssemblage == iAssemblageLower)

        !Bisection step
        if (lCompareAssemblage) then 
            dTempUpper = dTemperature
            iAssemblageUpper = iAssemblage
            call SortInt(iAssemblageUpper, nElements)
        elseif (lCompareAssemblage2) then
            dTempLower = dTemperature
            iAssemblageLower = iAssemblage
            call SortInt(iAssemblageLower, nElements)
        else
            !else if upper-lower<
            print *, "CALL NEW UPPER"
            call PhaseTransition(dTempLower, dTemperature - dTempTolerance, dTempTolerance, dPhaseTransitionTemp, iTransitions)
            print *, "CALL NEW LOWER"
            call PhaseTransition(dTemperature + dTempTolerance, dTempUpper, dTempTolerance, dPhaseTransitionTemp, iTransitions)
            EXIT LOOP_TEMP
        end if
            
    end do LOOP_TEMP

end subroutine PhaseTransition



