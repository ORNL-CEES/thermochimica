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

subroutine getSpeciesMassFraction(iSpecies, dFraction, ierr)
  USE ModuleThermo
  implicit none

  integer, intent(in)  :: iSpecies
  integer, intent(out) :: ierr
  real(8), intent(out) :: dFraction
  real(8) :: dSpeciesMass, dPhaseMass
  integer :: i

  ierr = 0
  dFraction = 0D0
  dSpeciesMass = 0D0
  dPhaseMass = 0D0
  if (iSpecies < 1 .OR. iSpecies > nSpecies) then
     ierr = 1
  else
      LOOP_CheckPhases: do i = 1, nSolnPhasesSys
          if (iSpecies <= nSpeciesPhase(i)) then
              call getPhaseMass(i, dPhaseMass, ierr)
              exit LOOP_CheckPhases
          end if
      end do LOOP_CheckPhases
      if (dPhaseMass > 0D0) then
          call getSpeciesMass(iSpecies, dSpeciesMass, ierr)
          dFraction = dSpeciesMass / dPhaseMass
      end if
  endif

  return
end subroutine getSpeciesMassFraction

subroutine getPhaseMassFraction(iPhaseID, dFraction, ierr)
  USE ModuleThermo
  implicit none

  integer, intent(in)  :: iPhaseID
  integer, intent(out) :: ierr
  real(8), intent(out) :: dFraction
  real(8) :: dTotalMass, dPhaseMass

  ierr = 0
  dFraction = 0D0
  dTotalMass = 0D0
  dPhaseMass = 0D0
  if (iPhaseID < 1 .OR. iPhaseID > nSpecies) then
     ierr = 1
  else
     call getTotalMass(dTotalMass, ierr)
     if (iPhaseID <= nSolnPhasesSys) call getPhaseMass(iPhaseID, dPhaseMass, ierr)
     if (iPhaseID > nSolnPhasesSys)  call getSpeciesMass(iPhaseID, dPhaseMass, ierr)
     dFraction = dPhaseMass / dTotalMass
  endif

  return
end subroutine getPhaseMassFraction
