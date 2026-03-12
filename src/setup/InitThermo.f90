
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
    USE ModuleThermoTolerance, ONLY: SetStandaloneDefaultTolerances
    USE ModuleParseCS, ONLY: iParamPassCS, nParamMax

    implicit none


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

    ! Normalizing constant (applied to mass of each element) to minimize numerical errors in mass balance
    ! equations:
    dNormalizeSum   = 1000D0
    dNormalizeInput = dNormalizeSum


    ! Numerical tolerances:
    ! ---------------------
    !
    ! Standalone Thermochimica retains the legacy numerical epsilon by default.  This may be overridden at
    ! compile time with THERMOCHIMICA_DEFAULT_TOLERANCE_EPSILON.
    !
    ! The minimum supportable mole fraction for an element is determined by the formula:
    !     x_min = dToleranceEpsilon / dTolerance(1)
    !
    ! With the default dTolerance(1) = 1D-5 and the legacy dToleranceEpsilon = 1D-14,
    ! this yields x_min = 1D-9.
    !
    ! To support smaller concentrations at the cost of looser mass balance accuracy, the
    ! calling code may reduce dTolerance(1) after InitThermo returns.
    ! Note that dTolerance(2..15) are derived from dTolerance(1) and dToleranceEpsilon
    ! and must be recomputed consistently if either quantity is changed.

    call SetStandaloneDefaultTolerances

    return

end subroutine InitThermo
