
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file        ModuleThermoIO.f90
    !> \brief       Fortran module for input/output of Thermochimica
    !> \details     The purpose of this module is to provide the means to share information
    !!              between Thermochimica and any software that calls it.
    !> \author      M.H.A. Piro
    !> \sa          ThermoDebug.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   10/17/2011      M.H.A. Piro     Original code
    !   04/23/2012      M.H.A. Piro     Add dOxygen documentation
    !
    !
    !
    !> \param       INFOThermo      An integer scalar identifying whether the program exits successfully or if
    !!                              it encounters an error.  Details are provided in ThermoDebug.f90.
    !> \param       dTemperature    A double real scalar representing the absolute temperature.
    !> \param       dPressure       A double real scalar representing the absolute hydrostatic pressure.
    !> \param       dElementMass    A double real vector representing the mass of each chemical element for all
    !!                                  the elements on the periodic table.
    !> \param       cThermoFileName  Name of a ChemSage data-file (e.g., 'UO2fuelthermo.dat').
    !!                            NOTE: this has a maximum of 120 characters, which includes the path and the
    !!                            file extension.
    !> \param       cInputUnitTemperature:  A character scalar representing the temperature units
    !!                                       ['K', 'C', 'F', 'R'];
    !> \param       cInputUnitPressure:  A character scalar representing the pressure units
    !!                                    ['atm', 'psi', 'bar', 'Pa', 'kPa'];
    !> \param       cInputUnitMass:  A character scalar representing the mass units
    !!                                ['mass fraction', 'kilograms', 'grams', 'pounds', 'mole fraction',
    !!                                 'atom fraction', 'atoms', 'moles'].
    !> \param       iPostProcessMode An integer scalar indicating the level of detail for post-processing.
    !!                              - 0: no post-processing;
    !!                              - 1: basic output, similar to FactSage;
    !!                              - 2: advanced output, similar to FactSage with additional data.
    !
    !
    !-------------------------------------------------------------------------------------------------------------


module ModuleThermoIO

    implicit none

    SAVE

    ! INPUT VARIABLES:
    integer                                  :: iCounter, iPrintResultsMode, nMinSpeciesPerPhase = 2
    real(8)                                  :: dTemperature, dPressure, dFuzzMag = 1D-12
    real(8),       dimension(0:168)          :: dElementMass
    logical,       dimension(0:118)          :: lPreset = .FALSE.
    character(15)                            :: cInputUnitTemperature, cInputUnitPressure, cInputUnitMass
    character(:), allocatable                :: cThermoFileName
    logical                                  :: lReinitAvailable = .FALSE., lReinitLoaded = .FALSE., lReinitRequested = .FALSE.
    logical                                  :: lStepTogether = .FALSE., lWriteJSON = .FALSE.
    logical                                  :: lFuzzyStoich = .FALSE., lGibbsMinCheck = .FALSE.
    integer                                  :: nPhasesExcluded = 0, nPhasesExcludedExcept = 0
    character(25), dimension(1000)           :: cPhasesExcluded = '', cPhasesExcludedExcept = ''

    ! Compound variables:
    integer                                  :: nCompounds = 0
    real(8),       dimension(118)            :: dCompoundMass
    real(8),       dimension(118,0:118)      :: dCompoundStoich
    character(12), dimension(118)            :: cCompoundNames
    logical                                  :: lCompoundStoichCalculated = .FALSE., lRetryAttempted = .FALSE.
    logical                                  :: lHeatCapacityEntropyEnthalpy = .FALSE.

    ! OUTPUT VARIABLES:
    integer                                  :: INFOThermo, nSolnPhasesOut, nPureConPhaseOut, nSpeciesOut
    real(8)                                  :: dGibbsEnergySys
    real(8), dimension(:), allocatable       :: dSolnPhaseMolesOut, dPureConPhaseMolesOut, dSpeciesMoleFractionOut
    character(25), dimension(:), allocatable :: cSolnPhaseNameOut, cPureConPhaseNameOut, cSpeciesNameOut, cSpeciesPhaseOut
    logical, dimension(:), allocatable       :: lSpeciesStable
    real(8)                                  :: dHeatCapacity = 0D0, dEntropy = 0D0, dEnthalpy = 0D0

end module ModuleThermoIO
