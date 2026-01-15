
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
    character(:), allocatable                :: cOutputFilePath
    character(:), allocatable                :: cResolvedOutputFilePath
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

contains

    logical function PathIsAbsolute(path)

        character(*), intent(in) :: path
        character(len(path)) :: trimmed
        integer :: n

        print *, 'PathIsAbsolute'
        trimmed = adjustl(path)
        n = len_trim(trimmed)

        print *, "trimmed: ", trimmed, " n: ", n

        if (n <= 0) then
            print *, 'Cond1'
            PathIsAbsolute = .FALSE.
        else if ((trimmed(1:1) == '/') .OR. (trimmed(1:1) == '\\')) then
            print *, 'Cond2'
            print *, 'trimmed(1:1)', trimmed(1:1)
            PathIsAbsolute = .TRUE.
        else if (n >= 2) then
            print *, 'Cond3'
            PathIsAbsolute = (trimmed(1:2) == '..')
            print *, 'trimmed(2:2): ', trimmed(1:2)
            print *, PathIsAbsolute
        else
            print *, 'Cond4'
            PathIsAbsolute = .FALSE.
        end if

    end function PathIsAbsolute

    subroutine SetDefaultOutputFilePath()

        call UpdateOutputFilePath('../outputs/thermoout.json')

        print *, "SetDefaultOutputPath"

    end subroutine SetDefaultOutputFilePath

    subroutine UpdateOutputFilePath(rawPath)

        character(*), intent(in) :: rawPath
        character(len(rawPath)) :: cleaned
        integer :: n

        print *, "UpdateOutputFilePath"
        print *, rawPath

        cleaned = adjustl(rawPath)
        n = len_trim(cleaned)

        if (n <= 0) then
            cOutputFilePath = '../outputs/thermoout.json'
            cResolvedOutputFilePath = DATA_DIRECTORY // '../outputs/thermoout.json'
            return
        end if

        cOutputFilePath = cleaned(1:n)

        if (PathIsAbsolute(cleaned(1:n))) then
            cResolvedOutputFilePath = cleaned(1:n)
        else if ((INDEX(cleaned(1:n), '/') == 0) .AND. (INDEX(cleaned(1:n), '\\') == 0)) then
            cResolvedOutputFilePath = DATA_DIRECTORY // '../outputs/' // cleaned(1:n)
        else
            cResolvedOutputFilePath = DATA_DIRECTORY // cleaned(1:n)
        end if

        print *, "UpdateOutputFilePath: ", cResolvedOutputFilePath

    end subroutine UpdateOutputFilePath

    ! function GetResolvedOutputFilePath() result(path)

    !     character(len=4096) :: path   ! choose a size you’re comfortable with

    !     print *, "GetResolvedOutputFilePath"

    !     if (.not. allocated(cResolvedOutputFilePath)) call SetDefaultOutputFilePath()
    !     if (.not. allocated(cResolvedOutputFilePath)) error stop "Resolved path not set"

    !     path = cResolvedOutputFilePath  ! will pad with spaces

    !     print *, "ResolvedPath: ", path

    ! end function GetResolvedOutputFilePath

    FUNCTION GetResolvedOutputFilePath() RESULT(path)
        IMPLICIT NONE
        CHARACTER(len=:), ALLOCATABLE :: path

        ! In modern compilers, assignment automatically handles allocation
        if (.not. allocated(cResolvedOutputFilePath)) call SetDefaultOutputFilePath()
        if (.not. allocated(cResolvedOutputFilePath)) error stop "Resolved path not set"

        ! path = cResolvedOutputFilePath

        ! print *, "path: ", path
        ! print *, "POOP!"

        path = cResolvedOutputFilePath

    END FUNCTION GetResolvedOutputFilePath


end module ModuleThermoIO
