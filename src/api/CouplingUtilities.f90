subroutine SetThermoFileName(cFileName)

  USE ModuleThermoIO, ONLY: cThermoFileName

  implicit none

  character(*), intent(in)::  cFileName

  cThermoFileName = cFileName

  return

end subroutine SetThermoFileName

subroutine SetUnitTemperature(cUnitTemperature)

  USE ModuleThermoIO, ONLY: cInputUnitTemperature

  implicit none

  character(*), intent(in)::  cUnitTemperature
  character(15) :: cUnitTemperatureLen

  cUnitTemperatureLen = cUnitTemperature(1:min(15,len(cUnitTemperature)))
  cInputUnitTemperature       = trim(cUnitTemperatureLen)

  return

end subroutine SetUnitTemperature

subroutine SetUnitPressure(cUnitPressure)

  USE ModuleThermoIO, ONLY: cInputUnitPressure

  implicit none

  character(*), intent(in)::  cUnitPressure
  character(15) :: cUnitPressureLen

  cUnitPressureLen = cUnitPressure(1:min(15,len(cUnitPressure)))
  cInputUnitPressure       = trim(cUnitPressureLen)

  return

end subroutine SetUnitPressure

subroutine SetUnitMass(cUnitMass)

  USE ModuleThermoIO, ONLY: cInputUnitMass

  implicit none

  character(*), intent(in)::  cUnitMass
  character(15) :: cUnitMassLen

  cUnitMassLen = cUnitMass(1:min(15,len(cUnitMass)))
  cInputUnitMass       = trim(cUnitMassLen)

  return

end subroutine SetUnitMass

subroutine SetStandardUnits

  USE ModuleThermoIO, ONLY: cInputUnitTemperature, cInputUnitPressure, cInputUnitMass

  implicit none

  cInputUnitTemperature = 'K'
  cInputUnitPressure    = 'atm'
  cInputUnitMass        = 'moles'

  return

end subroutine SetStandardUnits

subroutine SetModelicaUnits

  USE ModuleThermoIO, ONLY: cInputUnitTemperature, cInputUnitPressure, cInputUnitMass

  implicit none

  cInputUnitTemperature = 'K'
  cInputUnitPressure    = 'Pa'
  cInputUnitMass        = 'moles'

  return

end subroutine SetModelicaUnits

subroutine SetUnits(cTemperature, cPressure, cMass)

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

end subroutine SetUnits

subroutine SetTemperaturePressure(dTemp, dPress)

  USE ModuleThermoIO, ONLY: dTemperature, dPressure

  implicit none

  real(8), intent(in)::  dTemp
  real(8), intent(in)::  dPress

  dTemperature = dTemp
  dPressure = dPress

  return

end subroutine SetTemperaturePressure

subroutine SetPrintResultsMode(Pinfo)

  USE ModuleThermoIO, ONLY: iPrintResultsMode

  implicit none

  integer Pinfo

  iPrintResultsMode = Pinfo

  return

end subroutine SetPrintResultsMode

subroutine PresetElementMass(iAtom, dMass)

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

end subroutine PresetElementMass

subroutine SetElementMass(iAtom, dMass)

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
      ! write(*,*) 'SetElementMass ', iAtom, dElementMass(iAtom), dMass
  end if
  ! print *, "Set element mass in Thermochimica: Element=", iAtom, " Moles=", dElementMass(iAtom)

  return

end subroutine SetElementMass


subroutine GetNumberPhasesSystem(iNumSolnPhases, iNumConPhases)

  USE ModuleThermo, ONLY: nSolnPhasesSys, nConPhasesSys

  implicit none

  integer, intent(out):: iNumSolnPhases, iNumConPhases

  iNumConPhases = nConPhasesSys
  iNumSolnPhases = nSolnPhasesSys

  return

end subroutine GetNumberPhasesSystem

subroutine GetNumberSpeciesSystem(nSpeciesDB)
  USE ModuleThermo, ONLY: nSolnPhasesSys, nSpeciesPhase
  implicit none

  integer, intent(out), dimension(nSolnPhasesSys)     :: nSpeciesDB

  nSpeciesDB = nSpeciesPhase(1:size(nSpeciesDB))

  return

end subroutine GetNumberSpeciesSystem

subroutine GetElementMass(iAtom, dMass)

  USE ModuleThermoIO, ONLY: dElementMass

  implicit none

  integer, intent(in)::  iAtom
  real(8), intent(out)::  dMass

  dMass = 0D0
  if (iAtom > 0 .AND. iAtom <= 118) then
      dMass =  dElementMass(iAtom)
  end if

  return

end subroutine GetElementMass

subroutine CheckINFOThermo(dbginfo)

  USE ModuleThermoIO, ONLY: INFOThermo

  implicit none

  integer, intent(out)::  dbginfo

  dbginfo = INFOThermo

  return

end subroutine CheckINFOThermo

subroutine ResetINFOThermo

  USE ModuleThermoIO, ONLY: INFOThermo

  implicit none

  INFOThermo=0

  return

end subroutine ResetINFOThermo

subroutine SolPhaseParse(iElem, dMolSum)

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
end subroutine SolPhaseParse

! subroutine SSParseCSDataFile

!     USE ModuleThermoIO

!     implicit none

!     call ParseCSDataFile(cThermoFileName)

!     return

! end subroutine SSParseCSDataFile

subroutine APpmInBToMolInVol(dAppm, dAMassPerMol, dBMassPerMol, dBDens, dVol, iMolScale, dAMol, dBMol)

  ! Input
  ! dAppm          = element A, ppm, 0.000001 MU/MU
  ! dAMassPerMol   = element A, MU / mol
  ! dBMassPerMol   = element B, MU / mol
  ! dBDens         = element B, Density, Mass unit per unit volume  MU/LU^3
  ! dVol           = Volume, LU^3
  ! iMolScale      = Scale mol values wrt to 1, 2, or total
  ! Output
  ! dAMol          = Mol of element A in dVol
  ! dBMol          = Mol of element B in dVol

  implicit none

  integer, intent(in)  :: iMolScale
  real(8), intent(in)  :: dAppm, dAMassPerMol, dBMassPerMol, dBDens, dVol
  real(8), intent(out) :: dAMol, dBMol
  real(8)              :: dATotalMass, dBTotalMass

  ! ppm is 0.000001 MU/MU
  ! assume total density is equal to density of solvent
  dBTotalMass = dBDens * dVol
  dATotalMass = dBTotalMass * 0.000001 * dAppm

  dAMol = dATotalMass / dAMassPerMol
  dBMol = dBTotalMass / dBMassPerMol

  if( iMolScale == 1 )then
     dBMol = dBMol / dAMol
     dAmol = 1D0
  else if( iMolScale == 2 )then
     dAMol = dAMol / dBMol
     dBmol = 1D0
  else if( iMolScale == 3 )then
     dAMol = dAMol / (dAMol+dBMol)
     dBmol = 1D0 - dAMol
  end if

  return

end subroutine APpmInBToMolInVol

subroutine SSInitiateZRHD

  call SetThermoFileName('ZRHD_MHP.dat')
  call SetUnits('K','atm','moles')

  return

end subroutine SSInitiateZRHD

subroutine SSInitiateUO2PX

  call SetThermoFileName('DBV6_TMB_modified.dat')
  call SetUnits('K','atm','moles')

  return

end subroutine SSInitiateUO2PX

subroutine tokenize(str, delim, word, lword, n)

  implicit none

  character (len=*) :: str
  character(1)      :: delim
  integer           :: lword
  character(lword)  :: word(*)     ! need to fix the maximum n
  integer           :: n

  integer :: pos1, pos2, i

  n = 0
  pos1 = 1
  pos2 = 0

  DO
     pos2 = INDEX(str(pos1:), delim)
     IF (pos2 == 0) THEN
        n = n + 1
        word(n)=''
        word(n) = str(pos1:)
        EXIT
     END IF
     n = n + 1
     word(n)=''
     word(n) = str(pos1:pos1+pos2-2)
     pos1 = pos2+pos1
  END DO

  ! write(*,"(3A)") ' tokenize ', str,';'
  DO i = 1, n
      ! WRITE(*,"(2A)", ADVANCE="NO") word(i), "."
  END DO
  ! write(*,*)

end subroutine tokenize

subroutine chomp(str, len)
  ! remove \0 from c string
  ! must pass len that was result from strlen
  implicit none

  character (len=*) ::  str
  ! character(1)      :: str(*)
  integer           :: len,lchop

  ! write(*,*) 'chomp ',str,' len ',len

  lchop=len+1
  str(lchop:lchop)=""
  ! str=trim(str)

  return
end subroutine chomp

subroutine matchdict( word, dictionary, nwords, lenword, imatch )
  !
  !    Match word in a dictionary
  !    nwords - number of words in dictionary
  !    lenwords - array holding length of the words
  !    imatch - 0 = no match, 1 = match

  implicit none

  character (len=*) ::  word
  integer           :: nwords, lenword
  character (len=lenword) :: dictionary(nwords)
  character(25) :: cWord

  integer :: imatch

  integer :: i,lword

  imatch=0

  ! write(*,*) 'matchdict ',word

  lword=len(word)
  ! write(*,*) 'word ',word,'lword ', lword
  if(lword > 25)then
     write(*,"(A,i5)") "matchdict: word to match is too big ", lword
     stop
  end if

  cWord=""
  cWord(1:lword)=word(1:lword)

  do i=1,nwords
     if ( cWord == dictionary(i) ) then
        imatch = imatch + 1
     end if
  end do

  return
end subroutine matchdict

subroutine chopnull(str)

  implicit none

  !
  character (len=*) ::  str
  integer           ::  iloc

  iloc=scan(str,char(0))

  if(iloc > 0)then
     str(iloc:iloc)=""
  end if

  return
end subroutine chopnull

subroutine getMolFraction(i, value, ierr)
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
end subroutine getMolFraction

subroutine getChemicalPotential(i, value, ierr)
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
end subroutine getChemicalPotential

subroutine getElementPotential(i, value, ierr)
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
        write(*,*) 'Element idx', k,' ',cElementName(k)
     enddo

  else
     value=dElementPotential(i)*dTemperature*dIdealConstant
  endif

  return

end subroutine getElementPotential

subroutine SetReinitRequested(iRequested)

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

end subroutine SetReinitRequested

subroutine getReinitDataSizes(mElements,mSpecies)
  USE ModuleThermo, ONLY: nElements, nSpecies
  implicit none

  integer, intent(out) :: mElements, mSpecies

  mElements = nElements
  mSpecies = nSpecies

  return

end subroutine getReinitDataSizes

subroutine getReinitData(mAssemblage,mMolesPhase,mElementPotential, &
              mChemicalPotential,mMolFraction,mElementsUsed,mReinitAvailable)
  USE ModuleReinit
  USE ModuleThermoIO
  USE ModuleThermo, ONLY: nElements, nSpecies
  implicit none

  integer, intent(out)                           :: mReinitAvailable
  integer, intent(out), dimension(nElements)     :: mAssemblage
  real(8), intent(out), dimension(nElements)     :: mMolesPhase, mElementPotential
  real(8), intent(out), dimension(nSpecies)      :: mChemicalPotential, mMolFraction
  integer, intent(out), dimension(0:168) :: mElementsUsed


  if (lReinitAvailable) then
    mAssemblage = iAssemblage_Old
    mMolesPhase = dMolesPhase_Old
    mElementPotential = dElementPotential_Old
    mChemicalPotential = dChemicalPotential_Old
    mMolFraction =  dMolFraction_Old
    mElementsUsed = iElementsUsed_Old
    mReinitAvailable = 1
  else
    mReinitAvailable = 0
  end if

  return

end subroutine getReinitData

subroutine setReinitData(mElements,mSpecies,mAssemblage,mMolesPhase, &
              mElementPotential,mChemicalPotential,mMolFraction,mElementsUsed)
  USE ModuleReinit
  USE ModuleThermoIO
  USE ModuleThermo
  implicit none

  integer, intent(in)                            :: mElements, mSpecies
  integer, intent(in), dimension(mElements)      :: mAssemblage
  real(8), intent(in), dimension(mElements)      :: mMolesPhase, mElementPotential
  real(8), intent(in), dimension(mSpecies)       :: mChemicalPotential, mMolFraction
  integer, intent(in), dimension(0:168)  :: mElementsUsed

  allocate(dMolesPhase_Old(mElements),dChemicalPotential_Old(mSpecies),dElementPotential_Old(mElements),&
  dMolFraction_Old(mSpecies))
  allocate(iAssemblage_Old(mElements))
  ! allocate(iElementsUsed_Old(0:168))

  iAssemblage_Old = mAssemblage
  dMolesPhase_Old = mMolesPhase
  dElementPotential_Old = mElementPotential
  dChemicalPotential_Old = mChemicalPotential
  dMolFraction_Old = mMolFraction
  iElementsUsed_Old = mElementsUsed
  lReinitAvailable = .TRUE.

  return

end subroutine setReinitData

subroutine GetMolesPhase(mMolesPhase)
  USE ModuleThermo, ONLY: nElements, dMolesPhase
  implicit none

  real(8), intent(out), dimension(nElements) :: mMolesPhase

  mMolesPhase = dMolesPhase

  return

end subroutine GetMolesPhase

subroutine GibbsEnergyOfReinitData(mGibbsEnergyOut)
    USE ModuleReinit
    USE ModuleThermo, ONLY: nElements
    implicit none

    real(8), intent(out) :: mGibbsEnergyOut
    integer              :: nSolnPhasesReinit, nConPhasesReinit, i

    nSolnPhasesReinit = 0
    nConPhasesReinit = 0
    do i = 1, nElements
        if (iAssemblage_Old(i) > 0) then
            nConPhasesReinit = nConPhasesReinit + 1
        else if (iAssemblage_Old(i) < 0) then
            nSolnPhasesReinit = nSolnPhasesReinit + 1
        end if
    end do

    call GibbsEnergy(nConPhasesReinit, nSolnPhasesReinit, iAssemblage_Old, &
                     dMolesPhase_Old, dMolFraction_Old, mGibbsEnergyOut)

    return

end subroutine GibbsEnergyOfReinitData

subroutine GetAssemblage(mAssemblage)
  USE ModuleThermo, ONLY: nElements, iAssemblage
  implicit none

  integer, intent(out), dimension(nElements)     :: mAssemblage

  mAssemblage = iAssemblage

  return

end subroutine GetAssemblage

subroutine GetAllElementPotential(mElementPotential)
  USE ModuleThermo, ONLY: nElements, dElementPotential
  implicit none

  real(8), intent(out), dimension(nElements)     :: mElementPotential

  mElementPotential = dElementPotential

  return

end subroutine GetAllElementPotential

subroutine GetElementFraction(iAtom, dFrac)

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

end subroutine GetElementFraction

subroutine SetHeatCapacityEnthalpyEntropyRequested(iRequested)
  ! Toggles whether heat capacity, entropy, and enthalpy calculation is requested
  USE ModuleThermoIO, ONLY: lHeatCapacityEntropyEnthalpy

  implicit none

  ! passing bool/logical was sketchy so just going with an int here
  integer, intent(in)::  iRequested
  if (iRequested == 0) then
    lHeatCapacityEntropyEnthalpy = .FALSE.
  else
    lHeatCapacityEntropyEnthalpy = .TRUE.
  end if

  return

end subroutine SetHeatCapacityEnthalpyEntropyRequested

subroutine GetHeatCapacityEnthalpyEntropy(dHeatCapacityOut, dEnthalpyOut, dEntropyOut)
  USE ModuleThermoIO, ONLY: dHeatCapacity, dEnthalpy, dEntropy

  implicit none
  real(8), intent(out):: dHeatCapacityOut, dEnthalpyOut, dEntropyOut

  dHeatCapacityOut = dHeatCapacity
  dEnthalpyOut     = dEnthalpy
  dEntropyOut      = dEntropy

  return
end subroutine GetHeatCapacityEnthalpyEntropy

subroutine SetFuzzyStoich(lFuzzyStoichIn)
  USE ModuleThermoIO, ONLY: lFuzzyStoich

  implicit none
  logical, intent(in) :: lFuzzyStoichIn

  lFuzzyStoich = lFuzzyStoichIn

  return
end subroutine SetFuzzyStoich

subroutine SetFuzzyMagnitude(dFuzzMagIn)
  USE ModuleThermoIO, ONLY: dFuzzMag

  implicit none
  real(8), intent(in) :: dFuzzMagIn

  dFuzzMag = dFuzzMagIn

  return
end subroutine SetFuzzyMagnitude

subroutine SetGibbsMinCheck(lGibbsMinCheckIn)
  USE ModuleThermoIO, ONLY: lGibbsMinCheck

  implicit none
  logical, intent(in) :: lGibbsMinCheckIn

  lGibbsMinCheck = lGibbsMinCheckIn

  return
end subroutine SetGibbsMinCheck
