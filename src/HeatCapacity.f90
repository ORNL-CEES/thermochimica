subroutine HeatCapacity(dHeatCapacity)

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer   :: j
    real(8), dimension(3) :: dGibbsEnergies
    real(8), intent(out)  :: dHeatCapacity
    real(8)               :: dTStepSize, dSecondDer
    real(8)               :: dtemp0

    dtemp0 = dTemperature
    dTStepSize = 5e-1
    dGibbsEnergies = 0D0

    lReinitRequested = .FALSE.
    if(lReinitRequested) call SaveReinitData

    do j=-1,1
        dTemperature = dtemp0 + j*dTStepSize
        call ResetThermo
        call Thermochimica
        dGibbsEnergies(j+2) = dGibbsEnergySys
    end do

    ! Second derivative from  3-point stencil
    dSecondDer = (dGibbsEnergies(1)+dGibbsEnergies(3)-2D0*dGibbsEnergies(2))/dTStepSize**2
    ! Second derivative from  5-point stencil
    ! dSecondDer = (-dGibbsEnergies(1)-dGibbsEnergies(5)-30D0*dGibbsEnergies(3)+16*dGibbsEnergies(2)+16*dGibbsEnergies(4)) &
    !               / (12*dTStepSize**2)

    ! Heat capacity from 2nd derivative of G
    dHeatCapacity = -dtemp0*dSecondDer

end subroutine HeatCapacity
