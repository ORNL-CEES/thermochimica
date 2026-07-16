subroutine SetThermoFileNameISO(cFileName, lcFileName) &
    bind(C, name="TCAPI_setThermoFilename")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    character(kind=c_char,len=1), target, intent(in) :: cFileName(*)
    integer(c_size_t), intent(in), value             :: lcFileName
    character(kind=c_char,len=lcFileName), pointer :: fFileName

    call c_f_pointer(cptr=c_loc(cFileName), fptr=fFileName)

    call SetThermoFileName(fFileName)

    return

end subroutine SetThermoFileNameISO

subroutine ClearPhaseSelectionISO() &
    bind(C, name="TCAPI_clearPhaseSelection")

    USE ModuleThermoIO, ONLY: nPhasesExcluded, nPhasesExcludedExcept, &
                              cPhasesExcluded, cPhasesExcludedExcept

    implicit none

    nPhasesExcluded = 0
    nPhasesExcludedExcept = 0
    cPhasesExcluded = ''
    cPhasesExcludedExcept = ''

end subroutine ClearPhaseSelectionISO

subroutine AddPhaseSelectionISO(cPhase, lcPhase, lInclude, iInfo) &
    bind(C, name="TCAPI_addPhaseSelection")

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermoIO, ONLY: nPhasesExcluded, nPhasesExcludedExcept, &
                              cPhasesExcluded, cPhasesExcludedExcept

    implicit none

    character(kind=c_char,len=1), target, intent(in) :: cPhase(*)
    integer(c_size_t), intent(in), value :: lcPhase
    logical(c_bool), intent(in) :: lInclude
    integer(c_int), intent(out) :: iInfo
    character(kind=c_char,len=lcPhase), pointer :: fPhase

    iInfo = 0
    if (lcPhase == 0 .OR. lcPhase > len(cPhasesExcluded(1))) then
        iInfo = 1
        return
    else if ((lInclude .AND. nPhasesExcludedExcept >= size(cPhasesExcludedExcept)) .OR. &
             (.NOT. lInclude .AND. nPhasesExcluded >= size(cPhasesExcluded))) then
        iInfo = 2
        return
    end if

    call c_f_pointer(cptr=c_loc(cPhase), fptr=fPhase)
    if (lInclude) then
        nPhasesExcludedExcept = nPhasesExcludedExcept + 1
        cPhasesExcludedExcept(nPhasesExcludedExcept) = fPhase
    else
        nPhasesExcluded = nPhasesExcluded + 1
        cPhasesExcluded(nPhasesExcluded) = fPhase
    end if

end subroutine AddPhaseSelectionISO

subroutine SetUnitTemperatureISO(cUnitTemperature, lcUnitTemperature) &
    bind(C, name="TCAPI_setUnitTemperature")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    character(kind=c_char,len=1), target, intent(in)      :: cUnitTemperature(*)
    integer(c_size_t), intent(in), value                  :: lcUnitTemperature
    character(kind=c_char,len=lcUnitTemperature), pointer :: fUnitTemperature

    call c_f_pointer(cptr=c_loc(cUnitTemperature), fptr=fUnitTemperature)
    call SetUnitTemperature(fUnitTemperature)

    return

end subroutine SetUnitTemperatureISO

subroutine SetUnitPressureISO(cUnitPressure, lcUnitPressure) &
    bind(C, name="TCAPI_setUnitPressure")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    character(kind=c_char,len=1), target, intent(in)   :: cUnitPressure(*)
    integer(c_size_t), intent(in), value               :: lcUnitPressure
    character(kind=c_char,len=lcUnitPressure), pointer :: fUnitPressure

    call c_f_pointer(cptr=c_loc(cUnitPressure), fptr=fUnitPressure)
    call SetUnitPressure(fUnitPressure)

    return

end subroutine SetUnitPressureISO

subroutine SetUnitMassISO(cUnitMass, lcUnitMass) &
    bind(C, name="TCAPI_setUnitMass")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    character(kind=c_char,len=1), target, intent(in)   :: cUnitMass(*)
    integer(c_size_t), intent(in), value               :: lcUnitMass
    character(kind=c_char,len=lcUnitMass), pointer :: fUnitMass

    call c_f_pointer(cptr=c_loc(cUnitMass), fptr=fUnitMass)
    call SetUnitMass(fUnitMass)

    return

end subroutine SetUnitMassISO

subroutine SetStandardUnitsISO() &
    bind(C, name="TCAPI_setStandardUnits")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call SetStandardUnits

    return

end subroutine SetStandardUnitsISO

subroutine SetModelicaUnitsISO() &
    bind(C, name="TCAPI_setModelicaUnits")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call SetModelicaUnits

    return

end subroutine SetModelicaUnitsISO

subroutine SetUnitsISO(cTemperature, lcTemperature, cPressure, lcPressure, cMass, lcMass) &
    bind(C, name="TCAPI_setUnits")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    character(kind=c_char,len=1), target, intent(in)  ::  cTemperature(*), cPressure(*), cMass(*)
    integer(c_size_t), intent(in), value              :: lcTemperature,   lcPressure,   lcMass
    character(kind=c_char,len=lcTemperature), pointer ::  fTemperature
    character(kind=c_char,len=lcPressure),    pointer ::  fPressure
    character(kind=c_char,len=lcMass),        pointer ::  fMass

    call c_f_pointer(cptr=c_loc(cTemperature), fptr=fTemperature)
    call c_f_pointer(cptr=c_loc(cPressure),    fptr=fPressure)
    call c_f_pointer(cptr=c_loc(cMass),        fptr=fMass)

    call SetUnits(fTemperature,fPressure,fMass)

    return

end subroutine SetUnitsISO

subroutine GetNumberElementsDatabaseISO(dElements) &
    bind(C, name="TCAPI_getNumberElementsDatabase")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleParseCS, ONLY: nElementsCS

    implicit none

    integer(C_INT), intent(out):: dElements

    dElements = nElementsCS

    return

end subroutine GetNumberElementsDatabaseISO

function GetElementAtIndexISO(index, len) &
    bind(C, name='TCAPI_getElementAtIndex')

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleParseCS, ONLY: cElementNameCS

    implicit none

    integer(c_int), intent(in) :: index
    integer(c_int), intent(out) :: len
    type(c_ptr) :: GetElementAtIndexISO

    GetElementAtIndexISO = c_loc(cElementNameCS(index))
    len = len_trim(cElementNameCS(index))

    return

end function GetElementAtIndexISO

subroutine SetTemperaturePressureISO(dTemp, dPress) &
    bind(C, name="TCAPI_setTemperaturePressure")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    real(C_DOUBLE), intent(in) :: dTemp, dPress

    call SetTemperaturePressure(dTemp, dPress)

    return

end subroutine SetTemperaturePressureISO

subroutine SetPrintResultsModeISO(Pinfo) &
    bind(C, name="TCAPI_setPrintResultsMode")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(out) :: Pinfo

    call SetPrintResultsMode(Pinfo)

    return

end subroutine SetPrintResultsModeISO

subroutine PresetElementMassISO(iAtom, dMass) &
    bind(C, name="TCAPI_presetElementMass")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(in) :: iAtom
    real(C_DOUBLE), intent(in) :: dMass

    call PresetElementMass(iAtom, dMass)

    return

end subroutine PresetElementMassISO

subroutine SetElementMassISO(iAtom, dMass) &
    bind(C, name="TCAPI_setElementMass")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none
    integer(C_INT), intent(in) :: iAtom
    real(C_DOUBLE), intent(in) :: dMass

    call SetElementMass(iAtom, dMass)

    return

end subroutine SetElementMassISO

subroutine GetElementMassISO(iAtom, dMass) &
    bind(C, name="TCAPI_getElementMass")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(in)  :: iAtom
    real(C_DOUBLE), intent(out) :: dMass

    call GetElementMass(iAtom, dMass)

    return

end subroutine GetElementMassISO

subroutine CheckINFOThermoISO(dbginfo) &
    bind(C, name="TCAPI_checkInfoThermo")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(out) :: dbginfo

    call CheckINFOThermo(dbginfo)

    return

end subroutine CheckINFOThermoISO

subroutine ResetINFOThermoISO() &
    bind(C, name="TCAPI_resetInfoThermo")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call ResetINFOThermo

    return

end subroutine ResetINFOThermoISO

subroutine ResetThermoISO() &
    bind(C, name="TCAPI_resetThermo")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call ResetThermo

    return

end subroutine ResetThermoISO

subroutine ResetThermoAllISO() &
    bind(C, name="TCAPI_resetThermoAll")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call ResetThermoAll

    return

end subroutine ResetThermoAllISO


subroutine GetNumberPhasesSystemISO(iSolnPhases, iConPhases) &
    bind(C, name="TCAPI_getNumberPhasesSystem")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(out):: iSolnPhases, iConPhases

    call GetNumberPhasesSystem(iSolnPhases, iConPhases)

    return

end subroutine GetNumberPhasesSystemISO

function GetPhaseNameAtIndexISO(phase_index, phase_name_len) &
    bind(C, name='TCAPI_getPhaseNameAtIndex')

    USE, INTRINSIC :: ISO_C_BINDING

    USE ModuleThermo

    implicit none

    integer(c_int), intent(in) :: phase_index
    integer(c_int), intent(out) :: phase_name_len
    type(c_ptr) :: GetPhaseNameAtIndexISO

    if (phase_index <=  nSolnPhasesSys) then
        GetPhaseNameAtIndexISO = c_loc(cSolnPhaseName(phase_index))
        phase_name_len = len_trim(cSolnPhaseName(phase_index))
    else
        GetPhaseNameAtIndexISO = c_loc(cSpeciesName(MAXVAL(nSpeciesPhase) + phase_index - nSolnPhasesSys))
        phase_name_len = len_trim(cSpeciesName(MAXVAL(nSpeciesPhase) + phase_index - nSolnPhasesSys))
    end if

    return

end function GetPhaseNameAtIndexISO

subroutine GetNumberSpeciesSystemISO(nSpeciesDB) &
    bind(C, name='TCAPI_getNumberSpeciesSystem')

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleParseCS, ONLY: nSolnPhasesSysCS

    implicit none

    integer(c_int), intent(out), dimension(nSolnPhasesSysCS) :: nSpeciesDB

    call GetNumberSpeciesSystem(nSpeciesDB)

    return

end subroutine GetNumberSpeciesSystemISO


function GetSpeciesAtIndexISO(index, len) &
    bind(C, name='TCAPI_getSpeciesAtIndex')

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo

    implicit none

    integer(c_int), intent(in) :: index
    integer(c_int), intent(out) :: len
    type(c_ptr) :: GetSpeciesAtIndexISO

    GetSpeciesAtIndexISO = c_loc(cSpeciesName(index))
    len = len_trim(cSpeciesName(index))


    return

end function GetSpeciesAtIndexISO

subroutine IsPhaseGasISO(phase_index, isGas) &
    bind(C, name='TCAPI_isPhaseGas')

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: cSolnPhaseType

    implicit none
    integer(c_int), intent(in) :: phase_index
    logical(c_bool), intent(out) :: isGas

    isGas = .FALSE.

    if ( cSolnPhaseType(phase_index) == 'IDMX    ' ) isGas = .TRUE.

    return

end subroutine IsPhaseGasISO

subroutine IsPhaseMQMISO(phase_index, isMQM) &
    bind(C, name='TCAPI_isPhaseMQM')

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: cSolnPhaseType

    implicit none
    integer(c_int), intent(in) :: phase_index
    logical(c_bool), intent(out) :: isMQM

    isMQM = .FALSE.

    if ( cSolnPhaseType(phase_index) == 'SUBG    ' .OR. cSolnPhaseType(phase_index) == 'SUBQ    ' ) isMQM = .TRUE.

    return

end subroutine IsPhaseMQMISO

function GetMqmqaPairNameISO(phase_index, pair_index, len) &
    bind(C, name='TCAPI_getMqmqaPairAtIndex')

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo

    implicit none

    integer(c_int), intent(in) :: phase_index, pair_index
    integer(c_int), intent(out) :: len
    type(c_ptr) :: GetMqmqaPairNameISO

    integer :: iSPI

    iSPI = iPhaseSublattice(phase_index)    ! Following the name convention used in other MQMQA functions

    GetMqmqaPairNameISO = c_loc(cPairName(iSPI, pair_index))
    len = len_trim(cPairName(iSPI, pair_index))

    return

end function GetMqmqaPairNameISO

subroutine ThermoDebugISO() &
    bind(C, name="TCAPI_thermoDebug")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call ThermoDebug

    return

end subroutine ThermoDebugISO

subroutine PrintResultsISO() &
    bind(C, name="TCAPI_printResults")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call PrintResults

    return

end subroutine PrintResultsISO

subroutine SaveReinitDataISO() &
    bind(C, name="TCAPI_saveReinitData")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call SaveReinitData

    return

end subroutine SaveReinitDataISO

subroutine ResetReinitISO() &
    bind(C, name="TCAPI_resetReinit")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call ResetReinit

    return

end subroutine ResetReinitISO

subroutine SolPhaseParseISO(iElem, dMolSum) &
    bind(C, name="TCAPI_solPhaseParse")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(in)  :: iElem
    real(C_DOUBLE), intent(out) :: dMolSum

    call SolPhaseParse(iElem, dMolSum)

    return
end subroutine SolPhaseParseISO

subroutine ParseCSDataFileISO() &
    bind(C, name="TCAPI_parseCSDataFile")

    USE,INTRINSIC :: ISO_C_BINDING

    USE :: ModuleThermoIO

    implicit none

    call ParseCSDataFile(cThermoFileName)
    ! call ParseCSDataFile

    return

end subroutine ParseCSDataFileISO

subroutine ThermochimicaISO() &
    bind(C, name="TCAPI_thermochimica")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call ThermochimicaAPI

    return

end subroutine ThermochimicaISO

subroutine ThermochimicaSetupISO() &
    bind(C, name="TCAPI_setup")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call ThermochimicaSetupAPI

    return

end subroutine ThermochimicaSetupISO

subroutine ThermochimicaSolverISO() &
    bind(C, name="TCAPI_solve")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call ThermochimicaSolver

    return

end subroutine ThermochimicaSolverISO

subroutine ThermochimicaCheckSystemISO() &
    bind(C, name="TCAPI_checkSystem")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call CheckSystem

    return

end subroutine ThermochimicaCheckSystemISO

subroutine ThermochimicaInitISO() &
    bind(C, name="TCAPI_init")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call InitThermoAPI

    return

end subroutine ThermochimicaInitISO

subroutine CompThermoDataISO() &
    bind(C, name="TCAPI_compThermoData")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call CompThermoData

    return

end subroutine CompThermoDataISO

subroutine getMolFractionISO(i, value, ierr) &
    bind(C, name="TCAPI_getMolFraction")

    USE,INTRINSIC :: ISO_C_BINDING

    USE ModuleThermo
    implicit none

    integer(C_INT), intent(in)  ::  i
    integer(C_INT), intent(out) :: ierr
    real(C_DOUBLE), intent(out) :: value

    call getMolFraction(i, value, ierr)

    return
end subroutine getMolFractionISO

subroutine getChemicalPotentialISO(i, value, ierr) &
    bind(C, name="TCAPI_getChemicalPotential")

    USE,INTRINSIC :: ISO_C_BINDING

    USE ModuleThermo
    implicit none

    integer(C_INT), intent(in)  ::  i
    integer(C_INT), intent(out) :: ierr
    real(C_DOUBLE), intent(out) :: value

    call getChemicalPotential(i, value, ierr)

    return
end subroutine getChemicalPotentialISO

subroutine getElementPotentialISO(i, value, ierr) &
    bind(C, name="TCAPI_getElementPotential")

    USE,INTRINSIC :: ISO_C_BINDING

    USE ModuleThermoIO
    USE ModuleThermo
    implicit none

    integer(C_INT), intent(in)  ::  i
    integer(C_INT), intent(out) :: ierr
    real(C_DOUBLE), intent(out) :: value

    call getElementPotential(i, value, ierr)

    return

end subroutine getElementPotentialISO

subroutine SetReinitRequestedISO(iRequested) &
    bind(C, name="TCAPI_setReinitRequested")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(in) :: iRequested

    call SetReinitRequested(iRequested)

    return

end subroutine SetReinitRequestedISO

subroutine getReinitDataSizesISO(mElements, mSpecies) &
    bind(C, name="TCAPI_getReinitDataSizes")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(out) :: mElements, mSpecies

    call getReinitDataSizes(mElements,mSpecies)

    return

end subroutine getReinitDataSizesISO

subroutine getReinitDataISO(mAssemblage,mMolesPhase,mElementPotential, &
    mChemicalPotential,mMolFraction,mElementsUsed,mReinitAvailable, mIterations) &
    bind(C, name="TCAPI_getReinitData")

    USE,INTRINSIC :: ISO_C_BINDING

    USE ModuleThermo, ONLY: nElements, nSpecies
    implicit none

    integer(C_INT), intent(out)                           :: mReinitAvailable, mIterations
    integer(C_INT), intent(out), dimension(nElements)     :: mAssemblage
    real(C_DOUBLE), intent(out), dimension(nElements)     :: mMolesPhase, mElementPotential
    real(C_DOUBLE), intent(out), dimension(nSpecies)      :: mChemicalPotential, mMolFraction
    integer(C_INT), intent(out), dimension(0:168)         :: mElementsUsed

    call getReinitData(mAssemblage,mMolesPhase,mElementPotential, &
            mChemicalPotential,mMolFraction,mElementsUsed,mReinitAvailable,mIterations)

    return

end subroutine getReinitDataISO

subroutine setReinitDataISO(mElements,mSpecies,mAssemblage,mMolesPhase, &
    mElementPotential,mChemicalPotential,mMolFraction,mElementsUsed) &
    bind(C, name="TCAPI_setReinitData")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(in)                            :: mElements, mSpecies
    integer(C_INT), intent(in), dimension(mElements)      :: mAssemblage
    real(C_DOUBLE), intent(in), dimension(mElements)      :: mMolesPhase, mElementPotential
    real(C_DOUBLE), intent(in), dimension(mSpecies)       :: mChemicalPotential, mMolFraction
    integer(C_INT), intent(in), dimension(0:168)          :: mElementsUsed

    call setReinitData(mElements,mSpecies,mAssemblage,mMolesPhase, &
              mElementPotential,mChemicalPotential,mMolFraction,mElementsUsed)

    return

end subroutine setReinitDataISO

subroutine GetMolesPhaseISO(mMolesPhase) &
    bind(C, name="TCAPI_getMolesPhase")

    USE,INTRINSIC :: ISO_C_BINDING

    USE ModuleThermo, ONLY: nElements

    implicit none

    real(C_DOUBLE), intent(out), dimension(nElements) :: mMolesPhase

    call GetMolesPhase(mMolesPhase)

    return

end subroutine GetMolesPhaseISO

subroutine GetAssemblageISO(mAssemblage) &
    bind(C, name="TCAPI_getAssemblage")

    USE,INTRINSIC :: ISO_C_BINDING

    USE ModuleThermo, ONLY: nElements
    implicit none

    integer(C_INT), intent(out), dimension(nElements) :: mAssemblage

    call GetAssemblage(mAssemblage)

    return

end subroutine GetAssemblageISO

subroutine GetAllElementPotentialISO(mElementPotential) &
    bind(C, name="TCAPI_getAllElementPotential")

    USE,INTRINSIC :: ISO_C_BINDING

    USE ModuleThermo, ONLY: nElements
    implicit none

    real(C_DOUBLE), intent(out), dimension(nElements) :: mElementPotential

    call GetAllElementPotential(mElementPotential)

    return

end subroutine GetAllElementPotentialISO

subroutine GetElementFractionISO(iAtom, dFrac) &
    bind(C, name="TCAPI_getElementFraction")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(in)  ::  iAtom
    real(C_DOUBLE), intent(out) ::  dFrac

    call GetElementFraction(iAtom, dFrac)

    return

end subroutine GetElementFractionISO

subroutine PrintStateISO() &
    bind(C, name="TCAPI_printState")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    call PrintState

end subroutine PrintStateISO

subroutine GetElementMoleFractionInPhaseISO(cElement, lcElement, cPhase, lcPhase, dMolesOut, INFO) &
    bind(C, name="TCAPI_getElementMoleFractionInPhase")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO
    real(C_DOUBLE),       intent(out)                :: dMolesOut
    character(kind=c_char,len=1), target, intent(in) :: cPhase(*), cElement(*)
    integer(c_size_t), intent(in), value             :: lcPhase, lcElement
    character(kind=c_char,len=lcPhase), pointer      :: fPhase
    character(kind=c_char,len=lcElement), pointer    :: fElement

    call c_f_pointer(cptr=c_loc(cPhase), fptr=fPhase)
    call c_f_pointer(cptr=c_loc(cElement), fptr=fElement)

    call GetElementMoleFractionInPhase(fElement, lcElement, fPhase, lcPhase, dMolesOut, INFO)

    return

end subroutine GetElementMoleFractionInPhaseISO

subroutine GetElementMolesInPhaseISO(cElement, lcElement, cPhase, lcPhase, dMolesOut, INFO) &
    bind(C, name="TCAPI_getElementMolesInPhase")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO
    real(C_DOUBLE),       intent(out)                :: dMolesOut
    character(kind=c_char,len=1), target, intent(in) :: cPhase(*), cElement(*)
    integer(c_size_t), intent(in), value             :: lcPhase, lcElement
    character(kind=c_char,len=lcPhase), pointer      :: fPhase
    character(kind=c_char,len=lcElement), pointer    :: fElement

    call c_f_pointer(cptr=c_loc(cPhase), fptr=fPhase)
    call c_f_pointer(cptr=c_loc(cElement), fptr=fElement)

    call GetElementMolesInPhase(fElement, lcElement, fPhase, lcPhase, dMolesOut, INFO)

    return

end subroutine GetElementMolesInPhaseISO

subroutine GetOutputChemPotISO(cElementNameRequest, lcElementNameRequest, dElementChemPot, INFO) &
    bind(C, name="TCAPI_getOutputChemPot")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),      intent(out)                         :: INFO
    real(C_DOUBLE),      intent(out)                         :: dElementChemPot
    character(kind=c_char,len=1), target, intent(in)         :: cElementNameRequest(*)
    integer(c_size_t), intent(in), value                     :: lcElementNameRequest
    character(kind=c_char,len=lcElementNameRequest), pointer :: fElementNameRequest

    call c_f_pointer(cptr=c_loc(cElementNameRequest), fptr=fElementNameRequest)

    call GetOutputChemPot(fElementNameRequest, dElementChemPot, INFO)

    return

end subroutine GetOutputChemPotISO

subroutine GetOutputMolSpeciesISO(cSpeciesOut, lcSpeciesOut, dMolFractionOut, dMolesOut, INFO) &
    bind(C, name="TCAPI_getOutputMolSpecies")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO
    real(C_DOUBLE),       intent(out)                :: dMolFractionOut, dMolesOut
    character(kind=c_char,len=1), target, intent(in) :: cSpeciesOut(*)
    integer(c_size_t), intent(in), value             :: lcSpeciesOut
    character(kind=c_char,len=lcSpeciesOut), pointer :: fSpeciesOut

    call c_f_pointer(cptr=c_loc(cSpeciesOut), fptr=fSpeciesOut)

    call GetOutputMolSpecies(fSpeciesOut, lcSpeciesOut, dMolFractionOut, dMolesOut, INFO)

    return

end subroutine GetOutputMolSpeciesISO

subroutine GetOutputMolSpeciesPhaseISO(cPhase, lcPhase, cSpecies, lcSpecies, dMolFractionOut, INFO) &
    bind(C, name="TCAPI_getOutputMolSpeciesPhase")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(out)                      :: INFO
    real(C_DOUBLE), intent(out)                      :: dMolFractionOut
    character(kind=c_char,len=1), target, intent(in) :: cPhase(*), cSpecies(*)
    integer(c_size_t), intent(in), value             :: lcPhase, lcSpecies
    character(kind=c_char,len=lcPhase), pointer      :: fPhase
    character(kind=c_char,len=lcSpecies), pointer    :: fSpecies

    call c_f_pointer(cptr=c_loc(cPhase), fptr=fPhase)
    call c_f_pointer(cptr=c_loc(cSpecies), fptr=fSpecies)

    call GetOutputMolSpeciesPhase(fPhase, lcPhase, fSpecies, lcSpecies, dMolFractionOut, INFO)

    return

end subroutine GetOutputMolSpeciesPhaseISO

subroutine GetOutputMolSpeciesPhaseByIndexISO(iPhaseSystem, iSpeciesPhase, dMolFractionOut, INFO) &
    bind(C, name="TCAPI_getOutputMolSpeciesPhaseByIndex")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: cSolnPhaseName, dMolFraction, iAssemblage, nElements, &
                            nSolnPhasesSys, nSpeciesPhase
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer(C_INT), intent(in)    :: iPhaseSystem, iSpeciesPhase
    real(C_DOUBLE), intent(out)   :: dMolFractionOut
    integer(C_INT), intent(out)   :: INFO
    integer(C_INT)                :: i, iActualPhase, iActualSpecies

    INFO = 0
    dMolFractionOut = 0D0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if

    if (iPhaseSystem <= 0 .OR. iPhaseSystem > nSolnPhasesSys .OR. iSpeciesPhase <= 0 .OR. &
        iSpeciesPhase > nSpeciesPhase(iPhaseSystem) - nSpeciesPhase(iPhaseSystem - 1)) then
        INFO = 2
        return
    end if

    do i = 1, nElements
        if (iAssemblage(i) < 0) then
            if (cSolnPhaseName(-iAssemblage(i)) == cSolnPhaseName(iPhaseSystem)) then
                iActualPhase = -iAssemblage(i)
                iActualSpecies = nSpeciesPhase(iActualPhase - 1) + iSpeciesPhase
                dMolFractionOut = dMolFraction(iActualSpecies)
                return
            end if
        end if
    end do

    INFO = 1

end subroutine GetOutputMolSpeciesPhaseByIndexISO

subroutine GetElementMolesInPhaseByIndexISO(iElementSystem, iPhaseSystem, dMolesOut, INFO) &
    bind(C, name="TCAPI_getElementMolesInPhaseByIndex")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: cSolnPhaseName, dMolesPhase, dMolesSpecies, dStoichSpecies, &
                            iAssemblage, nElements, nSolnPhasesSys, nSpeciesPhase
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer(C_INT), intent(in)  :: iElementSystem, iPhaseSystem
    real(C_DOUBLE), intent(out) :: dMolesOut
    integer(C_INT), intent(out) :: INFO
    integer(C_INT)              :: i, iActualPhase, iTarget, j

    INFO = 0
    dMolesOut = 0D0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if
    if (iElementSystem <= 0 .OR. iElementSystem > nElements) then
        INFO = 1
        return
    end if
    if (iPhaseSystem <= 0) then
        INFO = 2
        return
    end if

    if (iPhaseSystem <= nSolnPhasesSys) then
        do i = 1, nElements
            if (iAssemblage(i) < 0) then
                if (cSolnPhaseName(-iAssemblage(i)) == cSolnPhaseName(iPhaseSystem)) then
                    iActualPhase = -iAssemblage(i)
                    do j = nSpeciesPhase(iActualPhase - 1) + 1, nSpeciesPhase(iActualPhase)
                        dMolesOut = dMolesOut + dMolesSpecies(j) * dStoichSpecies(j, iElementSystem)
                    end do
                    return
                end if
            end if
        end do
    else
        iTarget = MAXVAL(nSpeciesPhase) + iPhaseSystem - nSolnPhasesSys
        do i = 1, nElements
            if (iAssemblage(i) == iTarget) then
                dMolesOut = dMolesPhase(i) * dStoichSpecies(iTarget, iElementSystem)
                return
            end if
        end do
    end if

    INFO = 2

end subroutine GetElementMolesInPhaseByIndexISO

subroutine GetOutputSiteFractionISO(cSolnOut, lcSolnOut, iSublatticeOut, iConstituentOut, dSiteFractionOut, INFO) &
    bind(C, name="TCAPI_getOutputSiteFraction")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)   :: INFO, iSublatticeOut, iConstituentOut
    real(C_DOUBLE),       intent(out)   :: dSiteFractionOut
    character(kind=c_char,len=1), target, intent(in) :: cSolnOut(*)
    integer(c_size_t), intent(in), value             :: lcSolnOut
    character(kind=c_char,len=lcSolnOut), pointer    :: fSolnOut

    call c_f_pointer(cptr=c_loc(cSolnOut), fptr=fSolnOut)

    call GetOutputSiteFraction(fSolnOut, lcSolnOut, iSublatticeOut, iConstituentOut, dSiteFractionOut, INFO)

    return

end subroutine GetOutputSiteFractionISO

subroutine GetOutputSolnSpeciesISO(cSolnOut, lcSolnOut, cSpeciesOut, lcSpeciesOut, dMolFractionOut, dChemPotSpecies, INFO) &
    bind(C, name="TCAPI_getOutputSolnSpecies")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO
    real(C_DOUBLE),       intent(out)                :: dMolFractionOut, dChemPotSpecies
    character(kind=c_char,len=1), target, intent(in) :: cSolnOut(*), cSpeciesOut(*)
    integer(c_size_t), intent(in), value             :: lcSolnOut, lcSpeciesOut
    character(kind=c_char,len=lcSolnOut), pointer    :: fSolnOut
    character(kind=c_char,len=lcSpeciesOut), pointer :: fSpeciesOut

    call c_f_pointer(cptr=c_loc(cSolnOut), fptr=fSolnOut)
    call c_f_pointer(cptr=c_loc(cSpeciesOut), fptr=fSpeciesOut)

    call GetOutputSolnSpecies(fSolnOut, lcSolnOut, fSpeciesOut, lcSpeciesOut, dMolFractionOut, dChemPotSpecies, INFO)

    return

end subroutine GetOutputSolnSpeciesISO

subroutine GetPhaseIndexISO(cPhaseName, lcPhaseName, iIndexOut, INFO) &
    bind(C, name="TCAPI_getPhaseIndex")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO, iIndexOut
    character(kind=c_char,len=1), target, intent(in) :: cPhaseName(*)
    integer(c_size_t), intent(in), value             :: lcPhaseName
    character(kind=c_char,len=lcPhaseName), pointer  :: fPhaseName

    call c_f_pointer(cptr=c_loc(cPhaseName), fptr=fPhaseName)

    call GetPhaseIndex(fPhaseName, lcPhaseName, iIndexOut, INFO)

    return

end subroutine GetPhaseIndexISO

subroutine GetPhaseIndexBySystemIndexISO(iPhaseSystem, iIndexOut, INFO) &
    bind(C, name="TCAPI_getPhaseIndexBySystemIndex")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: cSolnPhaseName, iAssemblage, nElements, nSolnPhasesSys, nSpeciesPhase
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer(C_INT), intent(in)  :: iPhaseSystem
    integer(C_INT), intent(out) :: iIndexOut, INFO
    integer(C_INT)              :: i, iTarget

    INFO = 0
    iIndexOut = 0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if

    if (iPhaseSystem <= 0) then
        INFO = 1
        return
    end if

    if (iPhaseSystem <= nSolnPhasesSys) then
        do i = 1, nElements
            if (iAssemblage(i) < 0) then
                if (cSolnPhaseName(-iAssemblage(i)) == cSolnPhaseName(iPhaseSystem)) then
                    iIndexOut = i
                    return
                end if
            end if
        end do
    else
        iTarget = MAXVAL(nSpeciesPhase) + iPhaseSystem - nSolnPhasesSys
        do i = 1, nElements
            if (iAssemblage(i) == iTarget) then
                iIndexOut = i
                return
            end if
        end do
    end if

end subroutine GetPhaseIndexBySystemIndexISO

subroutine GetPhaseMolesBySystemIndexISO(iPhaseSystem, dMolesOut, INFO) &
    bind(C, name="TCAPI_getPhaseMolesBySystemIndex")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: cSolnPhaseName, dMolesPhase, iAssemblage, nElements, &
                            nSolnPhasesSys, nSpeciesPhase
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer(C_INT), intent(in)  :: iPhaseSystem
    real(C_DOUBLE), intent(out) :: dMolesOut
    integer(C_INT), intent(out) :: INFO
    integer(C_INT)              :: i, iTarget

    INFO = 0
    dMolesOut = 0D0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if
    if (iPhaseSystem <= 0) then
        INFO = 1
        return
    end if

    if (iPhaseSystem <= nSolnPhasesSys) then
        do i = 1, nElements
            if (iAssemblage(i) < 0) then
                if (cSolnPhaseName(-iAssemblage(i)) == cSolnPhaseName(iPhaseSystem)) then
                    dMolesOut = dMolesPhase(i)
                    return
                end if
            end if
        end do
    else
        iTarget = MAXVAL(nSpeciesPhase) + iPhaseSystem - nSolnPhasesSys
        do i = 1, nElements
            if (iAssemblage(i) == iTarget) then
                dMolesOut = dMolesPhase(i)
                return
            end if
        end do
    end if

end subroutine GetPhaseMolesBySystemIndexISO

subroutine GetSpeciesChemicalPotentialByIndexISO(iPhaseSystem, iSpeciesPhase, dPotentialOut, INFO) &
    bind(C, name="TCAPI_getSpeciesChemicalPotentialByIndex")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleGEMSolver, ONLY: lSolnPhases
    USE ModuleThermo, ONLY: dChemicalPotential, dIdealConstant, &
                            nSolnPhasesSys, nSpeciesPhase
    USE ModuleThermoIO, ONLY: dTemperature, INFOThermo

    implicit none

    integer(C_INT), intent(in)  :: iPhaseSystem, iSpeciesPhase
    real(C_DOUBLE), intent(out) :: dPotentialOut
    integer(C_INT), intent(out) :: INFO
    integer(C_INT)              :: iSpecies

    INFO = 0
    dPotentialOut = 0D0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if
    if (iPhaseSystem <= 0 .OR. iPhaseSystem > nSolnPhasesSys .OR. iSpeciesPhase <= 0 .OR. &
        iSpeciesPhase > nSpeciesPhase(iPhaseSystem) - nSpeciesPhase(iPhaseSystem - 1)) then
        INFO = 2
        return
    end if
    if (.NOT. lSolnPhases(iPhaseSystem)) then
        INFO = 1
        return
    end if

    iSpecies = nSpeciesPhase(iPhaseSystem - 1) + iSpeciesPhase
    dPotentialOut = dChemicalPotential(iSpecies) * dIdealConstant * dTemperature

end subroutine GetSpeciesChemicalPotentialByIndexISO

subroutine GetMqmqaEndmemberStoichiometricPotentialByIndexISO(iPhaseSystem, iEndmember, dPotentialOut, INFO) &
    bind(C, name="TCAPI_getMqmqaEndmemberStoichiometricPotentialByIndex")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: cSolnPhaseType, dElementPotential, dIdealConstant, dStoichPairs, &
                            iPhaseSublattice, nElements, nPairsSRO, nSolnPhasesSys
    USE ModuleThermoIO, ONLY: dTemperature, INFOThermo

    implicit none

    integer(C_INT), intent(in)  :: iPhaseSystem, iEndmember
    real(C_DOUBLE), intent(out) :: dPotentialOut
    integer(C_INT), intent(out) :: INFO
    integer(C_INT)              :: i, iSublatticePhase

    INFO = 0
    dPotentialOut = 0D0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if
    if (iPhaseSystem <= 0 .OR. iPhaseSystem > nSolnPhasesSys) then
        INFO = 2
        return
    end if
    if (cSolnPhaseType(iPhaseSystem) /= 'SUBG' .AND. cSolnPhaseType(iPhaseSystem) /= 'SUBQ') then
        INFO = 2
        return
    end if

    iSublatticePhase = iPhaseSublattice(iPhaseSystem)
    if (iEndmember <= 0 .OR. iEndmember > nPairsSRO(iSublatticePhase, 1)) then
        INFO = 2
        return
    end if

    do i = 1, nElements
        dPotentialOut = dPotentialOut + dElementPotential(i) * dStoichPairs(iSublatticePhase, iEndmember, i)
    end do
    dPotentialOut = dPotentialOut * dIdealConstant * dTemperature

end subroutine GetMqmqaEndmemberStoichiometricPotentialByIndexISO

subroutine GetPhaseGibbsEnergyBySystemIndexISO(iPhaseSystem, dTotalOut, dMolarOut, INFO) &
    bind(C, name="TCAPI_getPhaseGibbsEnergyBySystemIndex")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: dChemicalPotential, dIdealConstant, dMolesPhase, dMolFraction, dStdGibbsEnergy, &
                            iAssemblage, nConPhasesSys, nElements, nSolnPhasesSys, nSpeciesPhase
    USE ModuleThermoIO, ONLY: dTemperature, INFOThermo

    implicit none

    integer(C_INT), intent(in)  :: iPhaseSystem
    real(C_DOUBLE), intent(out) :: dTotalOut, dMolarOut
    integer(C_INT), intent(out) :: INFO
    integer(C_INT)              :: i, iFirst, iLast, iSpecies
    real(C_DOUBLE)              :: dFractionSum

    INFO = 0
    dTotalOut = 0D0
    dMolarOut = 0D0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if
    if (iPhaseSystem <= 0 .OR. iPhaseSystem > nSolnPhasesSys + nConPhasesSys) then
        INFO = 2
        return
    end if

    if (iPhaseSystem <= nSolnPhasesSys) then
        do i = 1, nElements
            if (iAssemblage(i) == -iPhaseSystem) then
                iFirst = nSpeciesPhase(iPhaseSystem - 1) + 1
                iLast = nSpeciesPhase(iPhaseSystem)
                dFractionSum = SUM(dMolFraction(iFirst:iLast))
                if (dFractionSum <= 0D0) then
                    INFO = 1
                    return
                end if
                dMolarOut = SUM(dChemicalPotential(iFirst:iLast) * dMolFraction(iFirst:iLast)) / dFractionSum
                dMolarOut = dMolarOut * dIdealConstant * dTemperature
                dTotalOut = dMolarOut * dMolesPhase(i)
                return
            end if
        end do
    else
        iSpecies = nSpeciesPhase(nSolnPhasesSys) + iPhaseSystem - nSolnPhasesSys
        do i = 1, nElements
            if (iAssemblage(i) == iSpecies) then
                dMolarOut = dStdGibbsEnergy(iSpecies) * dIdealConstant * dTemperature
                dTotalOut = dMolarOut * dMolesPhase(i)
                return
            end if
        end do
    end if

    INFO = 1

end subroutine GetPhaseGibbsEnergyBySystemIndexISO

subroutine GetPhaseDrivingForceBySystemIndexISO(iPhaseSystem, dDrivingForceOut, INFO) &
    bind(C, name="TCAPI_getPhaseDrivingForceBySystemIndex")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: dChemicalPotential, dElementPotential, dIdealConstant, dMolFraction, &
                            dSpeciesTotalAtoms, dStdGibbsEnergy, dStoichSpecies, &
                            nConPhasesSys, nElements, nSolnPhasesSys, nSpeciesPhase
    USE ModuleThermoIO, ONLY: dTemperature, INFOThermo

    implicit none

    integer(C_INT), intent(in)  :: iPhaseSystem
    real(C_DOUBLE), intent(out) :: dDrivingForceOut
    integer(C_INT), intent(out) :: INFO
    integer(C_INT)              :: iFirst, iLast, iSpecies
    real(C_DOUBLE)              :: dAtomSum, dReferencePotential, dResidual

    INFO = 0
    dDrivingForceOut = 0D0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if
    if (iPhaseSystem <= 0 .OR. iPhaseSystem > nSolnPhasesSys + nConPhasesSys) then
        INFO = 2
        return
    end if

    if (iPhaseSystem <= nSolnPhasesSys) then
        iFirst = nSpeciesPhase(iPhaseSystem - 1) + 1
        iLast = nSpeciesPhase(iPhaseSystem)
        dAtomSum = 0D0
        dResidual = 0D0
        do iSpecies = iFirst, iLast
            if (dMolFraction(iSpecies) > 0D0) then
                dReferencePotential = SUM(dElementPotential(1:nElements) * dStoichSpecies(iSpecies, 1:nElements))
                dResidual = dResidual + dMolFraction(iSpecies) * &
                            (dChemicalPotential(iSpecies) - dReferencePotential)
                dAtomSum = dAtomSum + dMolFraction(iSpecies) * dSpeciesTotalAtoms(iSpecies)
            end if
        end do
        if (dAtomSum <= 0D0) then
            INFO = 1
            return
        end if
        dDrivingForceOut = dResidual / dAtomSum
    else
        iSpecies = nSpeciesPhase(nSolnPhasesSys) + iPhaseSystem - nSolnPhasesSys
        dReferencePotential = SUM(dElementPotential(1:nElements) * dStoichSpecies(iSpecies, 1:nElements))
        dDrivingForceOut = (dStdGibbsEnergy(iSpecies) - dReferencePotential) / dSpeciesTotalAtoms(iSpecies)
    end if

    dDrivingForceOut = dDrivingForceOut * dIdealConstant * dTemperature

end subroutine GetPhaseDrivingForceBySystemIndexISO

subroutine GetSystemGibbsEnergyISO(dGibbsEnergyOut, INFO) bind(C, name="TCAPI_getSystemGibbsEnergy")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermoIO, ONLY: dGibbsEnergySys, INFOThermo

    implicit none

    real(C_DOUBLE), intent(out) :: dGibbsEnergyOut
    integer(C_INT), intent(out) :: INFO

    INFO = 0
    dGibbsEnergyOut = 0D0
    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if
    dGibbsEnergyOut = dGibbsEnergySys

end subroutine GetSystemGibbsEnergyISO

subroutine GetElementIndexByAtomicNumberISO(iAtomicNumber, iIndexOut, INFO) &
    bind(C, name="TCAPI_getElementIndexByAtomicNumber")

    USE,INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: iElementSystem
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer(C_INT), intent(in)  :: iAtomicNumber
    integer(C_INT), intent(out) :: iIndexOut, INFO
    integer(C_INT)              :: i

    INFO = 0
    iIndexOut = 0

    if (INFOThermo /= 0) then
        INFO = -1
        return
    end if

    do i = 1, SIZE(iElementSystem)
        if (iElementSystem(i) /= 0) iIndexOut = iIndexOut + 1
        if (iElementSystem(i) == iAtomicNumber) return
    end do

    iIndexOut = 0
    INFO = 1

end subroutine GetElementIndexByAtomicNumberISO

subroutine GetPureConPhaseMolISO(cPureConOut, lcPureConOut, dPureConMolOut, INFO) &
    bind(C, name="TCAPI_getPureConPhaseMol")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO
    real(C_DOUBLE),       intent(out)                :: dPureConMolOut
    character(kind=c_char,len=1), target, intent(in) :: cPureConOut(*)
    integer(c_size_t), intent(in), value             :: lcPureConOut
    character(kind=c_char,len=lcPureConOut), pointer :: fPureConOut

    call c_f_pointer(cptr=c_loc(cPureConOut), fptr=fPureConOut)

    call GetPureConPhaseMol(fPureConOut, dPureConMolOut, INFO)

    return

end subroutine GetPureConPhaseMolISO

subroutine GetSUBLSiteMolISO(cSolnOut, lcSolnOut, iSublatticeOut, iConstituentOut, dSiteMolOut, INFO) &
    bind(C, name="TCAPI_getSublSiteMol")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO, iSublatticeOut, iConstituentOut
    real(C_DOUBLE),       intent(out)                :: dSiteMolOut
    character(kind=c_char,len=1), target, intent(in) :: cSolnOut(*)
    integer(c_size_t), intent(in), value             :: lcSolnOut
    character(kind=c_char,len=lcSolnOut), pointer    :: fSolnOut

    call c_f_pointer(cptr=c_loc(cSolnOut), fptr=fSolnOut)

    call GetSUBLSiteMol(fSolnOut, lcSolnOut, iSublatticeOut, iConstituentOut, dSiteMolOut, INFO)

    return

end subroutine GetSUBLSiteMolISO

subroutine GetSolnPhaseMolISO(cSolnOut, lcSolnOut, dSolnMolOut, INFO) &
    bind(C, name="TCAPI_getSolnPhaseMol")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO
    real(C_DOUBLE),       intent(out)                :: dSolnMolOut
    character(kind=c_char,len=1), target, intent(in) :: cSolnOut(*)
    integer(c_size_t), intent(in), value             :: lcSolnOut
    character(kind=c_char,len=lcSolnOut), pointer    :: fSolnOut

    call c_f_pointer(cptr=c_loc(cSolnOut), fptr=fSolnOut)

    call GetSolnPhaseMol(fSolnOut, dSolnMolOut, INFO)

    return

end subroutine GetSolnPhaseMolISO

subroutine GetMqmqaMolesPairsISO(cPhaseName, lcPhaseName, dMolesPairsOut, INFO) &
    bind(C, name="TCAPI_getMqmqaMolesPairs")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO
    real(C_DOUBLE),       intent(out)                :: dMolesPairsOut
    character(kind=c_char,len=1), target, intent(in) :: cPhaseName(*)
    integer(c_size_t), intent(in), value             :: lcPhaseName
    character(kind=c_char,len=lcPhaseName), pointer  :: fPhaseName

    call c_f_pointer(cptr=c_loc(cPhaseName), fptr=fPhaseName)

    call GetMqmqaMolesPairs(fPhaseName, dMolesPairsOut, INFO)

    return

end subroutine GetMqmqaMolesPairsISO

subroutine GetMqmqaPairMolFractionISO(cPhaseName, lcPhaseName, cPairIn, lcPairIn, dMolesPairOut, INFO) &
    bind(C, name="TCAPI_getMqmqaPairMolFraction")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO
    real(C_DOUBLE),       intent(out)                :: dMolesPairOut
    character(kind=c_char,len=1), target, intent(in) :: cPhaseName(*), cPairIn(*)
    integer(c_size_t), intent(in), value             :: lcPhaseName, lcPairIn
    character(kind=c_char,len=lcPhaseName), pointer  :: fPhaseName
    character(kind=c_char,len=lcPairIn), pointer     :: fPairIn


    call c_f_pointer(cptr=c_loc(cPhaseName), fptr=fPhaseName)
    call c_f_pointer(cptr=c_loc(cPairIn), fptr=fPairIn)

    call GetMqmqaPairMolFraction(fPhaseName, lcPhaseName, fPairIn, lcPairIn, dMolesPairOut, INFO)

    return

end subroutine GetMqmqaPairMolFractionISO

subroutine GetMqmqaNumberPairsQuadsISO(cPhaseName,lcPhaseName, nPairs, nQuads, INFO) &
    bind(C, name="TCAPI_getMqmqaNumberPairsQuads")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT),       intent(out)                :: INFO, nPairs, nQuads
    character(kind=c_char,len=1), target, intent(in) :: cPhaseName(*)
    integer(c_size_t), intent(in), value             :: lcPhaseName
    character(kind=c_char,len=lcPhaseName), pointer  :: fPhaseName

    call c_f_pointer(cptr=c_loc(cPhaseName), fptr=fPhaseName)

    call GetMqmqaNumberPairsQuads(fPhaseName, nPairs, nQuads, INFO)

    return

end subroutine GetMqmqaNumberPairsQuadsISO

subroutine GetMqmqaConstituentFractionISO(cPhase,lcPhase,iSublattice,cConstituent,lcConstituent,dConstituentFractionOut,INFO) &
    bind(C, name="TCAPI_getMqmqaConstituentFraction")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(out)                      :: INFO
    integer(C_INT), intent(in)                       :: iSublattice
    real(C_DOUBLE), intent(out)                      :: dConstituentFractionOut
    character(kind=c_char,len=1), target, intent(in) :: cPhase(*), cConstituent(*)
    integer(c_size_t), intent(in), value             :: lcPhase,   lcConstituent
    character(kind=c_char,len=lcPhase),       pointer  :: fPhase
    character(kind=c_char,len=lcConstituent), pointer  :: fConstituent

    call c_f_pointer(cptr=c_loc(cPhase),       fptr=fPhase)
    call c_f_pointer(cptr=c_loc(cConstituent), fptr=fConstituent)

    call GetMqmqaConstituentFraction(fPhase, iSublattice, fConstituent, dConstituentFractionOut, INFO)

    return

end subroutine GetMqmqaConstituentFractionISO

subroutine SetHeatCapacityEnthalpyEntropyRequestedISO(iRequested) &
    bind(C, name="TCAPI_setHeatCapacityEnthalpyEntropyRequested")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer(C_INT), intent(in)::  iRequested

    call SetHeatCapacityEnthalpyEntropyRequested(iRequested)

    return

end subroutine SetHeatCapacityEnthalpyEntropyRequestedISO

subroutine GetHeatCapacityEnthalpyEntropyISO(dHeatCapacityOut, dEnthalpyOut, dEntropyOut) &
    bind(C, name="TCAPI_getHeatCapacityEnthalpyEntropy")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none
    real(C_DOUBLE), intent(out):: dHeatCapacityOut, dEnthalpyOut, dEntropyOut

    call GetHeatCapacityEnthalpyEntropy(dHeatCapacityOut, dEnthalpyOut, dEntropyOut)

    return

end subroutine GetHeatCapacityEnthalpyEntropyISO


subroutine SetFuzzyStoichISO(lFuzzyStoichIn) &
    bind(C, name="TCAPI_setFuzzyStoich")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none
    logical(C_BOOL), intent(in) :: lFuzzyStoichIn

    call SetFuzzyStoich(lFuzzyStoichIn)

    return
  end subroutine SetFuzzyStoichISO

  subroutine SetFuzzyMagnitudeISO(dFuzzMagIn) &
    bind(C, name="TCAPI_setFuzzyMagnitude")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none
    real(C_DOUBLE), intent(in) :: dFuzzMagIn

    call SetFuzzyMagnitude(dFuzzMagIn)

    return
  end subroutine SetFuzzyMagnitudeISO

  subroutine SetGibbsMinCheckISO(lGibbsMinCheckIn) &
    bind(C, name="TCAPI_setGibbsMinCheck")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none
    logical(C_BOOL), intent(in) :: lGibbsMinCheckIn

    call SetGibbsMinCheck(lGibbsMinCheckIn)

    return
  end subroutine SetGibbsMinCheckISO

  subroutine SetMinMoleFractionISO(dMinMolFrac) &
    bind(C, name="TCAPI_setMinMoleFraction")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none
    real(C_DOUBLE), intent(in) :: dMinMolFrac

    call SetMinMoleFraction(dMinMolFrac)

    return
  end subroutine SetMinMoleFractionISO

  subroutine SetMassBalanceToleranceISO(dMassBalanceTolerance) &
    bind(C, name="TCAPI_setMassBalanceTolerance")

    USE,INTRINSIC :: ISO_C_BINDING

    implicit none
    real(C_DOUBLE), intent(in) :: dMassBalanceTolerance

    call SetMassBalanceTolerance(dMassBalanceTolerance)

    return
  end subroutine SetMassBalanceToleranceISO

subroutine GetNumberSublatticesISO(iPhase, iCount, iInfo) &
    bind(C, name="TCAPI_getNumberSublattices")

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: nSolnPhasesSys, iPhaseSublattice, nSublatticePhase

    implicit none

    integer(c_int), intent(in) :: iPhase
    integer(c_int), intent(out) :: iCount, iInfo
    integer :: iSublatticePhase

    iCount = 0
    iInfo = 0
    if (iPhase < 1 .OR. iPhase > nSolnPhasesSys) then
        iInfo = 1
        return
    end if

    iSublatticePhase = iPhaseSublattice(iPhase)
    if (iSublatticePhase < 1) then
        iInfo = 2
        return
    end if
    iCount = nSublatticePhase(iSublatticePhase)

end subroutine GetNumberSublatticesISO

subroutine GetNumberConstituentsISO(iPhase, iSublattice, iCount, iInfo) &
    bind(C, name="TCAPI_getNumberConstituents")

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: nSolnPhasesSys, iPhaseSublattice, nSublatticePhase, &
                            nConstituentSublattice

    implicit none

    integer(c_int), intent(in) :: iPhase, iSublattice
    integer(c_int), intent(out) :: iCount, iInfo
    integer :: iSublatticePhase

    iCount = 0
    iInfo = 0
    if (iPhase < 1 .OR. iPhase > nSolnPhasesSys) then
        iInfo = 1
        return
    end if
    iSublatticePhase = iPhaseSublattice(iPhase)
    if (iSublatticePhase < 1 .OR. iSublattice < 1 .OR. &
        iSublattice > nSublatticePhase(iSublatticePhase)) then
        iInfo = 2
        return
    end if
    iCount = nConstituentSublattice(iSublatticePhase, iSublattice)

end subroutine GetNumberConstituentsISO

function GetConstituentNameAtIndexISO(iPhase, iSublattice, iConstituent, iLength, iInfo) &
    bind(C, name="TCAPI_getConstituentNameAtIndex")

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: nSolnPhasesSys, iPhaseSublattice, nSublatticePhase, &
                            nConstituentSublattice, cConstituentNameSUB

    implicit none

    integer(c_int), intent(in) :: iPhase, iSublattice, iConstituent
    integer(c_int), intent(out) :: iLength, iInfo
    integer :: i, iSublatticePhase
    character(kind=c_char), dimension(8), save, target :: cName
    type(c_ptr) :: GetConstituentNameAtIndexISO

    iLength = 0
    iInfo = 0
    GetConstituentNameAtIndexISO = c_null_ptr
    if (iPhase < 1 .OR. iPhase > nSolnPhasesSys) then
        iInfo = 1
        return
    end if
    iSublatticePhase = iPhaseSublattice(iPhase)
    if (iSublatticePhase < 1 .OR. iSublattice < 1 .OR. &
        iSublattice > nSublatticePhase(iSublatticePhase)) then
        iInfo = 2
        return
    else if (iConstituent < 1 .OR. &
             iConstituent > nConstituentSublattice(iSublatticePhase, iSublattice)) then
        iInfo = 3
        return
    end if

    iLength = len_trim(cConstituentNameSUB(iSublatticePhase, iSublattice, iConstituent))
    cName = c_null_char
    do i = 1, iLength
        cName(i) = cConstituentNameSUB(iSublatticePhase, iSublattice, iConstituent)(i:i)
    end do
    GetConstituentNameAtIndexISO = c_loc(cName(1))

end function GetConstituentNameAtIndexISO

subroutine GetConstituentFractionByIndexISO(iPhase, iSublattice, iConstituent, dValue, iInfo) &
    bind(C, name="TCAPI_getConstituentFractionByIndex")

    USE, INTRINSIC :: ISO_C_BINDING
    USE ModuleThermo, ONLY: nSolnPhasesSys, cSolnPhaseName, cSolnPhaseType, iPhaseSublattice, &
                            nSublatticePhase, nConstituentSublattice, cConstituentNameSUB

    implicit none

    integer(c_int), intent(in) :: iPhase, iSublattice, iConstituent
    real(c_double), intent(out) :: dValue
    integer(c_int), intent(out) :: iInfo
    integer :: iSublatticePhase
    integer(c_size_t) :: iPhaseNameLength

    dValue = 0D0
    iInfo = 0
    if (iPhase < 1 .OR. iPhase > nSolnPhasesSys) then
        iInfo = 1
        return
    end if
    iSublatticePhase = iPhaseSublattice(iPhase)
    if (iSublatticePhase < 1 .OR. iSublattice < 1 .OR. &
        iSublattice > nSublatticePhase(iSublatticePhase)) then
        iInfo = 2
        return
    else if (iConstituent < 1 .OR. &
             iConstituent > nConstituentSublattice(iSublatticePhase, iSublattice)) then
        iInfo = 3
        return
    end if

    if (cSolnPhaseType(iPhase) == 'SUBL' .OR. cSolnPhaseType(iPhase) == 'SUBLM') then
        iPhaseNameLength = len_trim(cSolnPhaseName(iPhase))
        call GetOutputSiteFraction(cSolnPhaseName(iPhase), iPhaseNameLength, &
                                   iSublattice, iConstituent, dValue, iInfo)
    else if (cSolnPhaseType(iPhase) == 'SUBG' .OR. cSolnPhaseType(iPhase) == 'SUBQ') then
        call GetMqmqaConstituentFraction(cSolnPhaseName(iPhase), iSublattice, &
            cConstituentNameSUB(iSublatticePhase, iSublattice, iConstituent), dValue, iInfo)
    else
        iInfo = 4
    end if

end subroutine GetConstituentFractionByIndexISO
