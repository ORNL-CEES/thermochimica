    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    PhaseTransition.F90
    !> \brief   Find Phase Transitions
    !> \author  A.E.F. Fitzsimmons
    !
    !> \todo    Input file setup, Return assemblages
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    04/10/2024    A.E.F. Fitzsimmons   Original Code
    !
    ! Purpose:
    ! ========
    !> \details Using a bisection search, this subroutine will find and record phase transition temperatures.
    !
    !
    !-------------------------------------------------------------------------------------------------------------
subroutine PhaseTransition(dTempMin, dTempMax, dTempTolerance, dPhaseTransitionTemp, iTransitions)
    USE ModuleThermoIO
    USE ModuleThermo
    implicit none

    real(8), intent(in) :: dTempMin, dTempMax, dTempTolerance
    real(8), dimension(10), intent(out) :: dPhaseTransitionTemp
    integer, intent(out) :: iTransitions

    real(8) :: dTempUpper, dTempLower
    integer :: i, iMaxIteration
    integer, dimension(nElements) :: iAssemblageUpper, iAssemblageLower
    logical :: lCompareAssemblage, lCompareAssemblage2

    ! Init
    i = 1
    iMaxIteration = 50
    iTransitions = 0
    dPhaseTransitionTemp = 0D0
    dTempLower = dTempMin
    dTempUpper = dTempMax
    lCompareAssemblage = .FALSE.
    lCompareAssemblage2 = .FALSE.

    if(dTempMin>=dTempMax) then
        INFOThermo = 56
        return
    end if

    ! Initialize Lower assemblage
    dTemperature = dTempLower
    call Thermochimica
    iAssemblageLower = iAssemblage
    call SortInt(iAssemblageLower, nElements)

    ! Loop until we reach the max temperature
    do while (dTempLower < dTempMax)

        ! Initialize Upper assemblage
        dTempUpper = dTempMax
        dTemperature = dTempUpper
        call Thermochimica
        iAssemblageUpper = iAssemblage
        call SortInt(iAssemblageUpper, nElements)

        ! Check if there's a transition between lower and upper
        if (ALL(iAssemblageLower == iAssemblageUpper)) then
            !No transitions in this range
            exit
        end if

        ! Bisection loop to find transition within [dTempLower, dTempUpper]
        do i = 1, iMaxIteration
            dTemperature = (dTempUpper + dTempLower) / 2D0

            ! Check convergence
            if ((dTempUpper - dTempLower) <= dTempTolerance) then
                ! Found transition
                iTransitions = iTransitions + 1
                dPhaseTransitionTemp(iTransitions) = dTemperature
                call SortInt(iAssemblage, nElements)

                ! New Lower
                dTempLower = dTemperature + dTempTolerance
                iAssemblageLower = iAssemblageUpper
                exit  ! Found a transition, break bisection loop to restart from new lower
            end if

            ! Midpoint calculation
            call Thermochimica
            call SortInt(iAssemblage, nElements)

            lCompareAssemblage = ALL(iAssemblage == iAssemblageUpper)
            lCompareAssemblage2 = ALL(iAssemblage == iAssemblageLower)

            ! Bisection check
            if (lCompareAssemblage) then
                dTempUpper = dTemperature
                iAssemblageUpper = iAssemblage
            elseif (lCompareAssemblage2) then
                dTempLower = dTemperature
                iAssemblageLower = iAssemblage
            else
                ! Phase not equal to lower or higher assemblage. Treat midpoint as new upper bound to keep bisecting
                dTempUpper = dTemperature
                iAssemblageUpper = iAssemblage
            end if
        end do

        ! If no transition found within tolerance after bisection, stop
        if ((dTempUpper - dTempLower) > dTempTolerance) then
            INFOThermo = 57
            exit
        end if

    end do

end subroutine PhaseTransition
