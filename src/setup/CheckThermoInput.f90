
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckThermoInput.f90
    !> \brief   Check the input quantities and character units.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      Thermochimica.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    08/30/2011    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to apply a unit conversion to the input variables (if necessary)
    !! and to ensure that the input variables are appropriate.  The working units for temperature,
    !! absolute hydrostaic pressure and mass are Kelvin [K], atmospheres [atm] and moles [mol], respectively.
    !! These variables are tested to ensure that they are within an acceptable range and that they are
    !! real.  An integer scalar INFOThermo is returned in a similar style as LAPACK to identify an error.
    !!
    !! A description of each value of INFOThermo that could be returned from this subroutine is given below:
    !!
    !   INFOThermo Value:       Description:
    !   -----------------       -------------
    !> \details
    !! <ul>
    !! <li> 0 - Successful exit, </li>
    !! <li> 1 - Temperature is out of range or a NAN, </li>
    !! <li> 2 - Hydrostatic pressure is out of range or a NAN, </li>
    !! <li> 3 - The mass of any element is out of range or a NAN, and </li>
    !! <li> 4 - The character string representing the input units is unrecognizable. </li>
    !! </ul>
    !!
    !! The variable cThermoInputUnits can take on the following values:
    !!
    !> \details
    !! <table border="1" width="400">
    !! <tr>
    !!    <td> <b> cThermoInputUnits </td> <td> Description </b> </td>
    !! </tr>
    !! <tr>
    !!    <td> "K" </td>
    !!    <td> Temperature in Kelvin </td>
    !! </tr>
    !! <tr>
    !!    <td> "C" </td>
    !!    <td> Temperature in Celsius </td>
    !! </tr>
    !! <tr>
    !!    <td> "F" </td>
    !!    <td> Temperature in Fahrenheit </td>
    !! </tr>
    !! <tr>
    !!    <td> "R" </td>
    !!    <td> Temperature in Rankine </td>
    !! </tr>
    !! <tr>
    !!    <td> "atm" </td>
    !!    <td> Pressure in atmospheres </td>
    !! </tr>
    !! <tr>
    !!    <td> "psi" </td>
    !!    <td> Pressure in pounds per square inch </td>
    !! </tr>
    !! <tr>
    !!    <td> "bar" </td>
    !!    <td> Pressure in bars </td>
    !! </tr>
    !! <tr>
    !!    <td> "Pa" </td>
    !!    <td> Pressure in Pascals </td>
    !! </tr>
    !! <tr>
    !!    <td> "kPa" </td>
    !!    <td> Pressure in kiloPascals </td>
    !! </tr>
    !! </table>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dTemperature              Temperature (converted to K)
    ! dPresssure                Absolute hydrostatic pressure (converted to atm)
    ! dElementMass              Total mass of an element, where the coefficient corresponds to the
    !                            atomic number (e.g., dMolesElement(92) refers to uranium).
    ! cThermoInputUnits         A character vector containing the units for temperature, pressure and
    !                            mass.
    ! INFOThermo                A scalar integer that indicates a successful exit or identifies an error.
    !
    !
    ! Note: the mass conversion of each element is performed in the CheckSystem.f90 subroutine.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine CheckThermoInput

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer::   i, j


    ! Convert termperature to Kelvin:
    select case (cInputUnitTemperature)
    case ('K')
        ! Do nothing.
    case ('C')
        dTemperature = dTemperature + 273.15D0
    case ('F')
        dTemperature = (dTemperature + 459.67D0) * (5D0/9D0)
    case ('R')
        dTemperature = dTemperature * (5D0/9D0)
    case default
        ! Temperature unit not recognized.
        INFOThermo = 4
        return
    end select

    cInputUnitTemperature = 'K'

    ! Check that the absolute temperature [K] is within an acceptable range and real:
    if ((dTemperature < 295D0).OR.(dTemperature > 6000D0).OR.(dTemperature /= dTemperature)) then
        INFOThermo = 1
        return
    end if

    ! Convert absolute hydrostatic pressure to atm:
    select case (cInputUnitPressure)
    case ('atm')
        ! Do nothing.
    case ('psi')
        dPressure = dPressure * 0.068045957064D0
    case ('bar')
        dPressure = dPressure * 0.98692316931D0
    case ('Pa')
        dPressure = dPressure * 0.009869231693D0 * 1D-3
    case ('kPa')
        dPressure = dPressure * 0.009869231693D0
    case default
        ! Pressure unit not recognized.
        INFOThermo = 4
        return
    end select

    cInputUnitPressure = 'atm'

    ! Check that mass units are valid and set scaling (scale everything to grams)
    dMassScale = 1D0
    select case (cInputUnitMass)
    case ('mass fraction','mole fraction','atom fraction','gram-atoms','moles','mol','grams','g')
        ! Assume any fractional mass units are just grams,
        ! the user will have to set a mass unit if they care.
        dMassScale = 1D0
    case ('atoms')
        ! This seems silly, but here we are.
        dMassScale = 1D0 / 6.0221409D23
    case ('kilograms','kg')
        dMassScale = dMassScale * 1D3
    case ('pounds','lbs')
        dMassScale = dMassScale * 4.53592D2
    case default
        ! The character string representing input units is not recognized.
        INFOThermo = 4
        return
    end select

    ! Check that the absolute hydrostatic pressure [atm] is within an acceptable range and real:
    if ((dPressure < 1D-6).OR.(dPressure > 1D6).OR.(dPressure /= dPressure)) then
        INFOThermo = 2
        return
    end if

    ! Check that all input quantities are real and non-negative:
    do i = 0,nElementsPT
        if ((dElementMass(i) /= dElementMass(i)).OR.(dElementMass(i) < 0D0)) then
            INFOThermo = 3
            return
        end if
    end do

    ! Check compound masses and stoichiometries
    do i = 1,nCompounds
        if ((dCompoundMass(i) /= dCompoundMass(i)).OR.(dCompoundMass(i) < 0D0)) then
            INFOThermo = 3
            return
        end if
        do j = 0,nElementsPT
            if ((dCompoundStoich(i,j) /= dCompoundStoich(i,j)).OR.(dCompoundStoich(i,j) < 0D0)) then
                INFOThermo = 3
                return
            end if
        end do
    end do

    ! Note: the mass conversion of each element is performed in the CheckSystem.f90 subroutine.

    return

end subroutine CheckThermoInput
