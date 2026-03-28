#ifndef THERMOCHIMICA_DEFAULT_TOLERANCE_EPSILON
#define THERMOCHIMICA_DEFAULT_TOLERANCE_EPSILON 1D-14
#endif

module ModuleThermoTolerance

    USE ModuleThermo, ONLY: dTolerance, dNormalizeInput

    implicit none

    SAVE

    integer, parameter :: iToleranceOverrideNone = 0
    integer, parameter :: iToleranceOverrideMassBalance = 1
    integer, parameter :: iToleranceOverrideMinMoleFraction = 2

    real(8), parameter :: dStandaloneToleranceEpsilon = THERMOCHIMICA_DEFAULT_TOLERANCE_EPSILON
    real(8)            :: dToleranceEpsilon = dStandaloneToleranceEpsilon
    real(8)            :: dToleranceOverrideValue = 0D0

    integer :: iToleranceOverrideMode = iToleranceOverrideNone
    logical :: lToleranceInitialized = .FALSE.

contains

    subroutine SetStandaloneDefaultTolerances

        implicit none

        call SetDefaultTolerances(dStandaloneToleranceEpsilon)

        return
    end subroutine SetStandaloneDefaultTolerances

    subroutine SetMachineDefaultTolerances

        implicit none

        call SetDefaultTolerances(EPSILON(1D0))

        return
    end subroutine SetMachineDefaultTolerances

    subroutine SetDefaultTolerances(dEpsilon)

        implicit none

        real(8), intent(in) :: dEpsilon

        dToleranceEpsilon   = dEpsilon
        lToleranceInitialized = .TRUE.

        dTolerance     = 0D0

        ! Tolerance of relative errors of mass balance equations (dimensionless):
        dTolerance(1)  = 1D-5

        ! Tolerance for the sum of mole fractions in a solution phase at equilibrium (dimensionless):
        dTolerance(2)  = 1D-5

        ! Tolerance for mass balance equations in LevelingSolver:
        dTolerance(3)  = -dToleranceEpsilon * 1D3

        ! Tolerance for Gibbs' Criteria applied to pure separate phases (dimensionless):
        dTolerance(4)  = DLOG(1D0 - dTolerance(2))

        ! Tolerance for residual of chemical potentials:
        dTolerance(5)  = DABS(dTolerance(4))

        ! Tolerance for minimum number of moles of an element of the system:
        dTolerance(6)  = dNormalizeInput * dToleranceEpsilon / dTolerance(1)

        ! Tolerance for minimum number of moles of a phase:
        dTolerance(7)  = dNormalizeInput * dToleranceEpsilon

        ! The minimum number of moles of a solution phase constituent:
        dTolerance(8)  = 1D-200

        ! The minimum number of moles of a solution phase that is added to the system:
        dTolerance(9)  = dNormalizeInput * 1D-3

        ! The maximum number of moles of a solution phase that is added to the system:
        dTolerance(10) = dNormalizeInput * 1D-1

        ! If the functional norm is less than this value, then skip the iteration history check:
        dTolerance(11) = 1D-5

        ! If the functional norm is less than this value and the system has exceeded a certain number of iterations:
        dTolerance(12) = 1D-3

        ! Tolerance for the maximum functional norm to add a pure condensed phase or solution phase:
        dTolerance(13) = 1D3

        ! Extreme-case solution-phase removal threshold:
        dTolerance(14) = dTolerance(1) * dTolerance(7)

        ! Tolerance for the maximum functional norm to check for a miscibility gap:
        dTolerance(15) = dTolerance(1)

        return
    end subroutine SetDefaultTolerances

    subroutine SetToleranceOverrideValue(dMassBalanceTolerance)

        implicit none

        real(8), intent(in) :: dMassBalanceTolerance

        iToleranceOverrideMode = iToleranceOverrideMassBalance
        dToleranceOverrideValue = dMassBalanceTolerance

        if (lToleranceInitialized) call ApplyToleranceOverride

        return
    end subroutine SetToleranceOverrideValue

    subroutine SetMinMoleFractionOverrideValue(dMinMolFrac)

        implicit none

        real(8), intent(in) :: dMinMolFrac

        iToleranceOverrideMode = iToleranceOverrideMinMoleFraction
        dToleranceOverrideValue = dMinMolFrac

        if (lToleranceInitialized) call ApplyToleranceOverride

        return
    end subroutine SetMinMoleFractionOverrideValue

    subroutine ApplyToleranceOverride

        implicit none

        real(8) :: dMinMolFrac, dTolClamped

        select case (iToleranceOverrideMode)
        case (iToleranceOverrideNone)
            return
        case (iToleranceOverrideMassBalance)
            dTolClamped = MAX(dToleranceOverrideValue, dToleranceEpsilon)
            dTolerance(1) = dTolClamped
        case (iToleranceOverrideMinMoleFraction)
            dMinMolFrac = MAX(dToleranceOverrideValue, dToleranceEpsilon)
            dTolerance(1) = dToleranceEpsilon / dMinMolFrac
        end select

        dTolerance(6)  = dNormalizeInput * dToleranceEpsilon / dTolerance(1)
        dTolerance(14) = dTolerance(1) * dTolerance(7)
        dTolerance(15) = dTolerance(1)

        return
    end subroutine ApplyToleranceOverride

end module ModuleThermoTolerance
