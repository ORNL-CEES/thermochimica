subroutine getSpeciesMass(iSpecies, dMass, ierr)
  USE ModuleThermo
  implicit none

  integer, intent(in)  :: iSpecies
  integer, intent(out) :: ierr
  real(8), intent(out) :: dMass
  integer :: i

  ierr = 0
  dMass = 0D0
  if (iSpecies < 1 .OR. iSpecies > nSpecies) then
     ierr = 1
  else
    do i = 1, nElements
        dMass = dMass + dStoichSpecies(iSpecies,i) * dAtomicMass(i)
    end do
    dMass = dMass * dMolesSpecies(iSpecies)
  endif

  return
end subroutine getSpeciesMass

subroutine getPhaseMass(iPhaseID, dMass, ierr)
  USE ModuleThermo
  implicit none

  integer, intent(in)  :: iPhaseID
  integer, intent(out) :: ierr
  real(8), intent(out) :: dMass
  integer :: i
  real(8) :: dSpeciesMass

  ierr = 0
  dMass = 0D0

  if ((iPhaseID < 1) .OR. (iPhaseID > nSolnPhasesSys)) then
    ierr = 1
  else
    do i = nSpeciesPhase(iPhaseID - 1) + 1, nSpeciesPhase(iPhaseID)
        call getSpeciesMass(i, dSpeciesMass, ierr)
        dMass = dMass + dSpeciesMass
    end do
  end if

  return
end subroutine getPhaseMass

subroutine getTotalMass(dMass, ierr)
  USE ModuleThermo
  implicit none

  integer, intent(out) :: ierr
  real(8), intent(out) :: dMass
  real(8) :: dPhaseMass
  integer :: i, k

  ierr = 0
  dMass = 0D0

  do i = 1, nSolnPhases
      k = -iAssemblage(nElements + 1 - i)
      call getPhaseMass(k, dPhaseMass, ierr)
      dMass = dMass + dPhaseMass
  end do

  do i = 1, nConPhases
      k = iAssemblage(i)
      call getSpeciesMass(k, dPhaseMass, ierr)
      dMass = dMass + dPhaseMass
  end do

  return
end subroutine getTotalMass

subroutine getMassFraction(iSpecies, dFraction, ierr)
  USE ModuleThermo
  implicit none

  integer, intent(in)  ::  iSpecies
  integer, intent(out) :: ierr
  real(8), intent(out) :: dFraction

  ierr = 0
  dFraction = 0D0
  if (iSpecies < 1 .OR. iSpecies > nSpecies) then
     ierr = 1
  else
     dFraction = dMolFraction(iSpecies)
  endif

  return
end subroutine getMassFraction
