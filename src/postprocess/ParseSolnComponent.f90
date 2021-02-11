subroutine ParseSolnComponent(iElem, cCompName, iSLen, dMolSum)
  ! return dMolSum, number of moles for element iElem that is in
  ! solution species cCompName
  !
  ! iElem      - index of element for which total number of moles is requested
  ! cCompName  - character string of one solution species to match
  ! iLCompName - length of cCompName, fortran origin, no chopping of null string
  ! dMolSum    - number of moles for element with index i in species that match names

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in):: iElem
    real(8), intent(out):: dMolSum

    character (len=*) :: cCompName
    character(25)     :: cChopName
    integer :: iSLen

    integer                               :: i, j, k

    cChopName=''
    cChopName = cCompName(1:iSLen)
    cChopName=trim(cChopName)

!    write(*,*) ' iElem ', iElem, '--', cChopName, '--', cCompName(1:iSLen), iSLen, dMolSum

    dMolSum=0D0
    do j = 1, nSolnPhases
       i = nElements - j + 1
       k = -iAssemblage(nElements - j + 1)
! debug
!       write(*,*) 'cChopName ', cChopName, ' cCompName ', cCompName(1:iSLen), ' solname ', cSolnPhaseName(k)
       if( cSolnPhaseName(k) == cChopName )then
          dMolSum = dMolesPhase(i) * dEffStoichSolnPhase(k,iElem)
! debug
!          write(*,*) '     Found ', cChopName, ' dMolSum ', dMolSum, ' = ', dMolesPhase(i), ' * ', dEffStoichSolnPhase(k,iElem)
       end if
    end do


    return
end subroutine ParseSolnComponent
