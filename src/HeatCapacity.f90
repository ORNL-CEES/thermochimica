subroutine HeatCapacity(dHeatCapacity)

    USE ModuleThermoIO
    USE ModuleThermo

    implicit none

    integer   :: j
    real(8), dimension(3) :: dGibbsEnergies
    real(8), intent(out)  :: dHeatCapacity
    real(8)               :: dTStepSize
    real(8)               :: dtemp0

    dtemp0 = dTemperature
    dTStepSize = 1e-2
    dGibbsEnergies = 0D0

    lReinitRequested = .TRUE.
    if(lReinitRequested) call SaveReinitData

    do j=-1,1
        dTemperature = dtemp0 + j*dTStepSize
        call ResetThermo
        call Thermochimica
        dGibbsEnergies(j+2) = dGibbsEnergySys
    end do

    ! Heat capacity from 3-point stencil 2nd derivative of G
    dHeatCapacity = -dtemp0*(dGibbsEnergies(1)+dGibbsEnergies(3)-2D0*dGibbsEnergies(2))/dTStepSize**2

end subroutine HeatCapacity
