subroutine FindSolnSpecies(cCompName, iLCompName, iElem, dMolSum)
  ! return dMolSum, number of moles for element iElem that is in
  ! solution species cCompName
  !
  ! cCompName  - character string of comma separated names of solution species to match
  ! iLCompName - length of cCompName, 0 - fortran origin, no chopping of null string
  !                                   N - c string size, must be result of strlen
  ! iElem      - index of element for which total number of moles is requested
  ! dMolSum    - number of moles for element with index i in species that match names

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    character (len=*) :: cCompName
    integer, intent(in) :: iLCompName
    integer, intent(in):: iElem
    real(8), intent(out):: dMolSum

    character(25)     :: cPhase(50)  ! need to fix max n -> matchdict
    integer :: nwords
    integer :: lenCompName

    integer           :: i, j, k, imatch, nmatch

    if(iLCompName <= 0)then
       lenCompName=len(cCompName)
    else
       lenCompName=iLCompName
    end if

    call tokenize(cCompName(1:lenCompName),",",cPhase,25,nwords)

    if(nwords <= 0)then
       write(*,"(A,i5)") "FindSolnSpecies: no words found ",nwords
       stop
    end if

    nmatch=0
    dMolSum=0D0
    do j = 1, nSolnPhases
       i = nElements - j + 1
       k = -iAssemblage(nElements - j + 1)
! debug
!       write(*,*) cSolnPhaseName(k)

       call matchdict(cSolnPhaseName(k), cPhase,nwords,25,imatch)

       if( imatch > 0 )then
          nmatch = nmatch + 1
          dMolSum = dMolSum + dMolesPhase(i) * dEffStoichSolnPhase(k,iElem)
! debug
!          write(*,"(3A,2e13.6,A,e13.6,5x)", ADVANCE="NO") &
!               '  Found ', trim(cSolnPhaseName(k)), ' dMolSum ', &
!               dMolSum, dMolesPhase(i), ' * ', dEffStoichSolnPhase(k,iElem)
       end if
    end do
    if(nmatch > 0)then
!       write(*,"(A,i5)") " Matched ",nmatch
       ! stop
    end if

    return
end subroutine FindSolnSpecies
