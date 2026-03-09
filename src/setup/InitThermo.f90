
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    InitThermo.f90
    !> \brief   Initialize Thermochimica
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \return  dNormalizeInput
    !> \sa      Thermochimica.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   09/09/2011      M.H.A. Piro         Created the ability to input different types of units.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to initialize Thermochimica. Various physical constants and
    !! numerical constants are defined.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dNormalizeInput       A double real scalar used to normalize the input moles of each element to
    !                        minimize numerical errors.
    !
    ! dIdealConstant        The ideal gas constant.
    ! dTolerance            A double real vector representing numerical tolerances, which is used in other
    !                        areas of the software.
    ! dEPS                  A double real scalar representing machine precision (approximate).
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine InitThermo

    USE ModuleThermo
    USE ModuleParseCS, ONLY: iParamPassCS, nParamMax

    implicit none

    real(8)::   dEPS


    ! Initialize variables:
    iParamPassCS   = 0
    nElements      = 0
    nSpecies       = 0
    nSolnPhases    = 0
    nSolnPhasesSys = 0
    nConPhasesSys  = 0
    nParam         = 0
    nMagParam      = 0
    nMaxParam      = nParamMax

    ! The universal gas constant used below is calculated using Avogadro's # (6.0221415D23) to 8
    ! sig. fig.'s and Boltzmann's constant (1.3806503D-23) to 8 sig. fig.'s.  dIdealConstant is
    ! therefore computed to 7 sig. fig.'s.  Note: SOLGASMIX uses 8.31433.
    dIdealConstant = 8.314472D0

    ! Machine precision (actual double-precision machine epsilon, ~2.22D-16):
    dEPS = EPSILON(1D0)

    ! Normalizing constant (applied to mass of each element) to minimize numerical errors in mass balance
    ! equations:
    dNormalizeSum   = 1000D0
    dNormalizeInput = dNormalizeSum


    ! Numerical tolerances:
    ! ---------------------
    !
    ! The minimum supportable mole fraction for an element is determined by the formula:
    !     x_min = EPSILON(1D0) / dTolerance(1)
    !
    ! With the default dTolerance(1) = 1D-5, this yields x_min ~ 2.22D-11 (~22 ppt).
    ! To support smaller concentrations at the cost of looser mass balance accuracy, the
    ! calling code may reduce dTolerance(1) after InitThermo returns.  For example:
    !     dTolerance(1) = 1D-3  ->  x_min ~ 2.22D-13  (~0.2 ppt)
    !     dTolerance(1) = 1D-1  ->  x_min ~ 2.22D-15  (sub-femtomole fractions)
    ! Note that dTolerance(2..15) are derived from dTolerance(1) here and must be
    ! recomputed if dTolerance(1) is changed after this routine returns.

    ! Initialize variables:
    dTolerance     = 0D0

    ! Tolerance of relative errors of mass balance equations (dimensionless):
    dTolerance(1)  = 1D-5

    ! Tolerance for the sum of mole fractions in a solution phase at equilibrium (dimensionless):
    dTolerance(2)  = 1D-5

    ! Tolerance for mass balance equations in LevelingSolver:
    dTolerance(3)  = -dEPS * 1D3

    ! Tolerance for Gibbs' Criteria (see reference at top) applied to pure seperate phases (dimensionless).
    ! Note that this tolerance criterion is consistent with the second.
    dTolerance(4)  = DLOG(1D0 - dTolerance(2))

    ! Tolerance for residual of chemical potentials:
    dTolerance(5)  = DABS(DLOG(1D0 - dTolerance(2)))

    ! Tolerance for minimum number of moles of an element of the system:
    ! (= dNormalizeInput * EPSILON(1D0) / dTolerance(1); as a mole fraction: EPSILON(1D0)/dTolerance(1))
    dTolerance(6)  = dNormalizeInput * dEPS / dTolerance(1)

    ! Tolerance for minimum number of moles of a phase:
    dTolerance(7)  = dNormalizeInput * dEPS

    ! The minimum number of moles of a solution phase constituent:
    dTolerance(8)  = 1D-200

    ! The minimum number of moles of a solution phase that is added to the system (note that this is redefined
    ! in CheckSystem.f90 as the minimum of this value and ten times the minimum number of moles of any element:
    dTolerance(9) = dNormalizeInput * 1D-3

    ! The maximum number of moles of a solution phase that is added to the system:
    dTolerance(10) = dNormalizeInput * 1D-1

    ! If the functional norm is less than this value, then skip the iteration history check:
    dTolerance(11) = 1D-5

    ! If the functional norm is less than this value and the system has exceeded a certian number of iterations,
    ! call it a day:
    dTolerance(12) = 1D-3

    ! Tolerance for the maximum functional norm to add a pure condensed phase or solution phase to the system:
    dTolerance(13) = 1D3

    ! If the number of moles of a solution phase is sufficiently small and it cannot be removed directly, swapped
    ! for a pure condensed phase or a solution phase, then Thermochimica will try to remove a pure condensed phase.
    ! This is only tested when the number of moles of this phase is extremely small, as represented by:
    dTolerance(14) = dTolerance(1) * dTolerance(7)

    ! Tolerance for the maximum functional norm to check for a miscibility gap:
    dTolerance(15) = dTolerance(1)

    return

end subroutine InitThermo
