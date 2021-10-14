subroutine HeatCapacity(dHeatCapacity)

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer               :: k, nMaxHeatCapAttempt
    real(8), dimension(2) :: dGibbsEnergies
    real(8), intent(out)  :: dHeatCapacity
    real(8)               :: dTStepSize, dSecondDer
    real(8)               :: dtemp0, dGibbs0, dGibbsDiff, dTargetDiff, dIncrease, dDecrease

    dtemp0         = dTemperature
    dGibbs0        = dGibbsEnergySys
    dTStepSize     = 2.5e-1
    dGibbsEnergies = 0D0
    dTargetDiff    = 2.5e-4
    dGibbsDiff     = dTargetDiff
    dIncrease      = 1.5D0
    dDecrease      = 0.6D0
    nMaxHeatCapAttempt = 5

    lReinitRequested = .FALSE.
    if(lReinitRequested) call SaveReinitData

    HeatCapTrial: do k = 1, nMaxHeatCapAttempt
        if (dGibbsDiff > 2D0*dTargetDiff) then
            dTStepSize = dTStepSize * dDecrease
        else if (dGibbsDiff < 1D-1*dTargetDiff) then
            dTStepSize = dTStepSize * dIncrease
        else if (k > 1) then
            exit HeatCapTrial
        end if
        ! Lower T
        dTemperature = dtemp0 - dTStepSize
        call ResetThermo
        call Thermochimica
        if (INFOThermo /= 0) then
            INFOThermo = 0
            cycle HeatCapTrial
        end if
        dGibbsEnergies(1) = dGibbsEnergySys
        ! Higher T
        dTemperature = dtemp0 + dTStepSize
        call ResetThermo
        call Thermochimica
        if (INFOThermo /= 0) then
            INFOThermo = 0
            cycle HeatCapTrial
        end if
        dGibbsEnergies(2)     = dGibbsEnergySys
        dGibbsDiff = ABS(dGibbsEnergies(2) - dGibbsEnergies(1)) / ABS(dGibbs0)
    end do HeatCapTrial

    ! Second derivative from  3-point stencil
    dSecondDer = (dGibbsEnergies(1)+dGibbsEnergies(2)-2D0*dGibbs0)/dTStepSize**2
    ! Second derivative from  5-point stencil
    ! dSecondDer = (-dGibbsEnergies(1)-dGibbsEnergies(5)-30D0*dGibbsEnergies(3)+16*dGibbsEnergies(2)+16*dGibbsEnergies(4)) &
    !               / (12*dTStepSize**2)

    ! Heat capacity from 2nd derivative of G
    dHeatCapacity = -dtemp0*dSecondDer

    call ResetThermo
    dTemperature = dtemp0
    call Thermochimica

end subroutine HeatCapacity
