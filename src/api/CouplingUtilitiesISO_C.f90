subroutine SetThermoFileNameISO(cFileName, lcFileName) &
    bind(C, name="TCAPI_setThermoFilename")

    USE ModuleThermoIO, ONLY: cThermoFileName
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    character(kind=c_char,len=1), target, intent(in) :: cFileName(*)
    integer(c_size_t), intent(in), value             :: lcFileName
    character(kind=c_char,len=lcFileName), pointer :: fFileName

    call c_f_pointer(cptr=c_loc(cFileName), fptr=fFileName)
    
    call SetThermoFileName(fFileName,lcFileName)

    return

end subroutine SetThermoFileNameISO

subroutine SetUnitTemperatureISO(cUnitTemperature, lcUnitTemperature) &
    bind(C, name="TCAPI_setUnitTemperature")

    USE ModuleThermoIO, ONLY: cInputUnitTemperature
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

    USE ModuleThermoIO, ONLY: cInputUnitPressure
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

    USE ModuleThermoIO, ONLY: cInputUnitMass
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    character(kind=c_char,len=1), target, intent(in)   :: cUnitMass(*)
    integer(c_size_t), intent(in), value               :: lcUnitMass
    character(kind=c_char,len=lcUnitMass), pointer :: fUnitMass

    call c_f_pointer(cptr=c_loc(cUnitMass), fptr=fUnitMass)
    cInputUnitMass = fUnitMass

    return

end subroutine SetUnitMassISO

subroutine SetStandardUnitsISO() &
    bind(C, name="TCAPI_setStandardUnits")

    USE ModuleThermoIO, ONLY: cInputUnitTemperature, cInputUnitPressure, cInputUnitMass

    implicit none

    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'

    return

end subroutine SetStandardUnitsISO

subroutine SetModelicaUnitsISO() &
    bind(C, name="TCAPI_setModelicaUnits")

    USE ModuleThermoIO, ONLY: cInputUnitTemperature, cInputUnitPressure, cInputUnitMass

    implicit none

    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'Pa'
    cInputUnitMass        = 'moles'

    return

end subroutine SetModelicaUnitsISO

subroutine SetUnitsISO(cTemperature, cPressure, cMass)

    USE ModuleThermoIO, ONLY: cInputUnitTemperature, cInputUnitPressure, cInputUnitMass

    implicit none

    character(*), intent(in)::  cTemperature
    character(*), intent(in)::  cPressure
    character(*), intent(in)::  cMass

    character(15) :: cTemperatureLen
    character(15) :: cPressureLen
    character(15) :: cMassLen

    cInputUnitTemperature = 'K'
    cInputUnitPressure    = 'atm'
    cInputUnitMass        = 'moles'


    if(len_trim(cTemperature) > 0)then
        cTemperatureLen = cTemperature(1:min(15,len(cTemperature)))
        cInputUnitTemperature = trim(cTemperatureLen)
    end if
    if(len_trim(cPressure) > 0)then
        cPressureLen = cPressure(1:min(15,len(cPressure)))
        cInputUnitPressure = trim(cPressureLen)
    end if
    if(len_trim(cMass) > 0)then
        cMassLen = cMass(1:min(15,len(cMass)))
        cInputUnitMass = trim(cMassLen)
    end if

    return

end subroutine SetUnitsISO

subroutine SetTemperaturePressureISO(dTemp, dPress) &
    bind(C, name="TCAPI_setTemperaturePressure")

    USE ModuleThermoIO, ONLY: dTemperature, dPressure

    implicit none

    real(8), intent(in)::  dTemp
    real(8), intent(in)::  dPress

    dTemperature = dTemp
    dPressure = dPress

    return

end subroutine SetTemperaturePressureISO

subroutine SetPrintResultsModeISO(Pinfo) &
    bind(C, name="TCAPI_setPrintResultsMode")

    USE ModuleThermoIO, ONLY: iPrintResultsMode

    implicit none

    integer Pinfo

    iPrintResultsMode = Pinfo

    return

end subroutine SetPrintResultsModeISO

subroutine PresetElementMassISO(iAtom, dMass) &
    bind(C, name="TCAPI_presetElementMass")

    USE ModuleThermoIO, ONLY: dElementMass, lPreset

    implicit none

    integer, intent(in)::  iAtom
    real(8), intent(in)::  dMass

    if( iAtom == 0 )then
        dElementMass = dMass
    else if( iAtom < 0 .or. iAtom > 118 )then
        write(*,*) 'Error in PresetElementMass ', iAtom, dMass
        stop
    else
        dElementMass(iAtom)      = dMass
        lPreset(iAtom)           = .TRUE.
    end if

    return

end subroutine PresetElementMassISO

subroutine SetElementMassISO(iAtom, dMass) &
    bind(C, name="TCAPI_setElementMass")


    USE ModuleThermoIO, ONLY: dElementMass, lPreset

    implicit none

    integer, intent(in)::  iAtom
    real(8), intent(in)::  dMass
    integer :: i

    if( iAtom == 0 ) then
        do i = 0, 118
            if (.NOT. lPreset(i)) dElementMass(i) = dMass
        end do
    else if( iAtom < 0 .or. iAtom > 118 )then
        write(*,*) 'Error in SetElementMass ', iAtom, dMass
        stop
    else
        if (.NOT. lPreset(iAtom)) dElementMass(iAtom) = dMass
    end if

    return

end subroutine SetElementMassISO

subroutine GetElementMassISO(iAtom, dMass)

    USE ModuleThermoIO, ONLY: dElementMass

    implicit none

    integer, intent(in)::  iAtom
    real(8), intent(out)::  dMass

    dMass = 0D0
    if (iAtom > 0 .AND. iAtom <= 118) then
        dMass =  dElementMass(iAtom)
    end if

    return

end subroutine GetElementMassISO

subroutine CheckINFOThermoISO(dbginfo) &
    bind(C, name="TCAPI_checkInfoThermo")

    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer, intent(out)::  dbginfo

    dbginfo = INFOThermo

    return

end subroutine CheckINFOThermoISO

subroutine ResetINFOThermoISO() &
    bind(C, name="TCAPI_resetInfoThermo")

    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    INFOThermo=0

    return

end subroutine ResetINFOThermoISO

subroutine ResetThermoISO() &
    bind(C, name="TCAPI_resetThermo")

    implicit none

    call ResetThermo

    return

end subroutine ResetThermoISO

subroutine ResetThermoAllISO() &
    bind(C, name="TCAPI_resetThermoAll")

    implicit none

    call ResetThermoAll

    return

end subroutine ResetThermoAllISO

subroutine ThermoDebugISO() &
    bind(C, name="TCAPI_thermoDebug")

    implicit none

    call ThermoDebug

    return

end subroutine ThermoDebugISO

subroutine PrintResultsISO() &
    bind(C, name="TCAPI_printResults")

    implicit none

    call PrintResults

    return

end subroutine PrintResultsISO

subroutine SaveReinitDataISO() &
    bind(C, name="TCAPI_saveReinitData")

    implicit none

    call SaveReinitData

    return

end subroutine SaveReinitDataISO

subroutine ResetReinitISO() &
    bind(C, name="TCAPI_resetReinit")

    implicit none

    call ResetReinit

    return

end subroutine ResetReinitISO


subroutine SolPhaseParseISO(iElem, dMolSum) &
    bind(C, name="TCAPI_solPhaseParse")

    ! quick hack for bison, ZrH, H in
    ! needs checking of input, intents, etc

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in):: iElem
    real(8), intent(out):: dMolSum
    real(8) :: dMolTemp
    integer                               :: i, j, k

    real(8),    dimension(:),   allocatable :: dTempVec

    ! Allocate arrays to sort solution phases:
    if (allocated(dTempVec)) deallocate(dTempVec)

    allocate(dTempVec(nSolnPhases))

    do i = 1, nSolnPhases
        j = nElements - i + 1
        dTempVec(i) = dMolesPhase(j)
    end do

    dMolSum = 0D0
    do j = 1, nSolnPhases

        ! Absolute solution phase index:
        k = -iAssemblage(nElements - j + 1)

        dMolTemp=dTempVec(j) * dEffStoichSolnPhase(k,iElem)
        dMolSum = dMolSum + dMolTemp
        !  write(*,"(A,A,e13.6)", ADVANCE="NO") ' -- ', trim(cSolnPhaseName(k)), dMolTemp
    end do
    ! write(*,*)

    return
end subroutine SolPhaseParseISO

subroutine SSParseCSDataFileISO() &
    bind(C, name="TCAPI_sSParseCSDataFile")

    USE ModuleThermoIO
    USE ModuleSS

    implicit none

    call ParseCSDataFile(cThermoFileName)

    return

end subroutine SSParseCSDataFileISO

subroutine ThermochimicaISO() &
    bind(C, name="TCAPI_thermochimica")

    implicit none

    call Thermochimica()

    return

end subroutine ThermochimicaISO

subroutine getMolFractionISO(i, value, ierr) &
    bind(C, name="TCAPI_getMolFraction")

    USE ModuleThermo
    implicit none

    integer, intent(in)::  i
    integer, intent(out):: ierr
    real(8), intent(out):: value

    ierr=0
    value=0D0
    if( i < 1 .OR. i > nSpecies )then
        ierr = 1
    else
        value=dMolFraction(i)
    endif

    return
end subroutine getMolFractionISO

subroutine getChemicalPotentialISO(i, value, ierr) &
    bind(C, name="TCAPI_getChemicalPotential")

    USE ModuleThermo
    implicit none

    integer, intent(in)::  i
    integer, intent(out):: ierr
    real(8), intent(out):: value

    ierr=0
    value=0D0
    if( i < 1 .OR. i > nSpecies )then
        ierr = 1
    else
        value=dChemicalPotential(i)
    endif

    return
end subroutine getChemicalPotentialISO

subroutine getElementPotentialISO(i, value, ierr) &
    bind(C, name="TCAPI_getElementPotential")

    USE ModuleThermoIO
    USE ModuleThermo
    implicit none

    integer, intent(in)::  i
    integer, intent(out):: ierr
    real(8), intent(out):: value

    integer k

    ierr=0
    value=0D0
    if( i < 1 .OR. i > nElements )then
        ierr = 1
        write(*,*) 'Element out of range ', i, nElements
        do k=1,nElements
            write(*,*) 'Element idx',k,' ',cElementName(k)
        enddo

    else
        value=dElementPotential(i)*dTemperature*dIdealConstant
    endif

    return

end subroutine getElementPotentialISO

subroutine SetReinitRequestedISO(iRequested) &
    bind(C, name="TCAPI_setReinitRequested")

    USE ModuleThermoIO, ONLY: lReinitRequested

    implicit none

    ! passing bool/logical was sketchy so just going with an int here
    integer, intent(in)::  iRequested
    if (iRequested == 0) then
        lReinitRequested = .FALSE.
    else
        lReinitRequested = .TRUE.
    end if

    return

end subroutine SetReinitRequestedISO

subroutine getReinitDataSizesISO(mElements, mSpecies) &
    bind(C, name="TCAPI_getReinitDataSizes")

    USE ModuleThermo, ONLY: nElements, nSpecies
    implicit none

    integer, intent(out)                           :: mElements, mSpecies

    mElements = nElements
    mSpecies = nSpecies

    return

end subroutine getReinitDataSizesISO

subroutine GetMolesPhaseISO(mMolesPhase) &
    bind(C, name="TCAPI_getMolesPhase")

    USE ModuleThermo, ONLY: nElements, dMolesPhase
    implicit none

    real(8), intent(out), dimension(nElements)     :: mMolesPhase

    mMolesPhase = dMolesPhase

    return

end subroutine GetMolesPhaseISO

subroutine GetAssemblageISO(mAssemblage) &
    bind(C, name="TCAPI_getAssemblage")

    USE ModuleThermo, ONLY: nElements, iAssemblage
    implicit none

    integer, intent(out), dimension(nElements)     :: mAssemblage

    mAssemblage = iAssemblage

    return

end subroutine GetAssemblageISO

subroutine GetAllElementPotentialISO(mElementPotential) &
    bind(C, name="TCAPI_getAllElementPotential")

    USE ModuleThermo, ONLY: nElements, dElementPotential
    implicit none

    real(8), intent(out), dimension(nElements)     :: mElementPotential

    mElementPotential = dElementPotential

    return

end subroutine GetAllElementPotentialISO

subroutine GetElementFractionISO(iAtom, dFrac) &
    bind(C, name="TCAPI_getElementFraction")

    USE ModuleThermoIO, ONLY: dElementMass
    USE ModuleThermo,   ONLY: nElementsPT

    implicit none

    integer, intent(in) ::  iAtom
    real(8), intent(out)::  dFrac
    real(8)             ::  dTotalElementMass
    integer :: i

    if( iAtom <= 0 .or. iAtom > 118 )then
        write(*,*) 'Error in GetElementFraction ', iAtom
        stop
    else
        dTotalElementMass = 0D0
        do i = 1, nElementsPT
            dTotalElementMass = dTotalElementMass + dElementMass(i)
        end do
        dFrac = dElementMass(iAtom) / dTotalElementMass
    end if

    return

end subroutine GetElementFractionISO
