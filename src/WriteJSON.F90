subroutine WriteJSON

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModulePrintUtils

    implicit none

    integer :: i

    iPrintResultsMode     = 2

    ! Return if the print results mode is zero.
    if (iPrintResultsMode == 0) return

    ! Only proceed for a successful calculation:
    IF_PASS: if (INFOThermo == 0) then
      ! inquire(file= DATA_DIRECTORY // '../thermoout.json', exist=exist)
      ! if (exist) then
      !   open(1, file= DATA_DIRECTORY // '../thermoout.json', status='old', &
      !         position='append', action='write')
      ! else
        ! open(1, file= DATA_DIRECTORY // '../thermoout.json', status='new', action='write')
      ! end if
        open(1, file= DATA_DIRECTORY // '../thermoout.json', status='REPLACE', &
             action='write')

        write(1,*) '{'

        ! Print the results for solution phases:
        call WriteJSONSolnPhase

        ! Print the results for pure condensed phases:
        call WriteJSONPureConPhase

        write(1,*) '  "elements": {'
        do i = 1, nElements
            write(1,*) '    "', TRIM(cElementName(i)), '": {'
            write(1,*) '      "moles": ', dMolesElement(i), ','
            write(1,*) '      "element potential": ', dElementPotential(i) * dIdealConstant * dTemperature
            if (i < nElements) then
              write(1,*) '    },'
            else
              write(1,*) '    }'
            end if
        end do
        write(1,*) '  },'

        write(1,*) '  "temperature": ', dTemperature, ','
        write(1,*) '  "pressure": ', dPressure, ','
        write(1,*) '  "Integral Gibbs energy": ', dGibbsEnergySys, ','
        write(1,*) '  "Functional norm": ', dGEMFunctionNorm, ','
        write(1,*) '  "solution phases": ', nSolnPhases, ','
        write(1,*) '  "pure condensed phases": ', nConPhases

        write(1,*) '}'

        close (1)
    end if IF_PASS

    return

end subroutine WriteJSON

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

subroutine WriteJSONSolnPhase

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModulePrintUtils

    implicit none

    integer                                 :: c, i, j, k, l, s, iFirst, iLast, iChargedPhaseID, nMax, nCutOff
    integer,    dimension(:),   allocatable :: iTempVec, iTempSpecies
    real(8)                                 :: dCutOff, Tcritical, B, StructureFactor
    real(8),    dimension(:),   allocatable :: dTempVec, dTempSpecies
    character(2)                            :: cDummy
    character(30)                           :: cDummyB


    ! Initiate variables:
    cDummy  = '  '
    dCutOff = 1D-70

    ! Allocate arrays to sort solution phases:
    if (allocated(iTempVec)) deallocate(iTempVec)
    if (allocated(dTempVec)) deallocate(dTempVec)

    allocate(dTempVec(nSolnPhases), iTempVec(nSolnPhases))
    do i = 1, nSolnPhases
        j = nElements - i + 1
        dTempVec(i) = dMolesPhase(j)
    end do

    ! Sort solution phases:
    call SortPick(nSolnPhases, dTempVec, iTempVec)

    write(1,*) '  "solution phases": {'

    ! Loop through solution phases:
    LOOP_SolnStable: do j = 1, nSolnPhases

        ! Relative and absolute solution phase indices, respectively:
        k      = nElements - iTempVec(j) + 1
        l      = -iAssemblage(k)

        ! First and last solution species indices, respectively:
        iFirst = nSpeciesPhase(l-1) + 1
        iLast  = nSpeciesPhase(l)

        ! Solution phase precharacter:
        if (j > 1) cDummy = '+ '
        nMax              = 1

        ! Print solution phase name:
        write(1,*) '    "', TRIM(cSolnPhaseName(l)), '": {'
        write(1,*) '      "moles": ', dMolesPhase(k), ','

        if ((cSolnPhaseType(l) == 'SUBLM') .OR. (cSolnPhaseType(l) == 'RKMPM')) then
            Tcritical = 0D0
            B = 0D0
            call CompMagneticTemperatureMoment(l,Tcritical,B)
            StructureFactor = dCoeffGibbsMagnetic(iFirst,3)
            if (Tcritical < 0D0) then
                Tcritical = -Tcritical * StructureFactor
                write(1,*) '      "Neel temperature": ', Tcritical, ','
            else
                write(1,*) '      "Curie temperature": ', Tcritical, ','
            end if
            if (B < 0D0) then
                B         = -B * StructureFactor
            end if
            write(1,*) '      "Magnetic moment per atom": ', B, ','
        end if

        ! if ((cSolnPhaseType(l)) == 'SUBG' .OR. (cSolnPhaseType(l) == 'SUBQ')) then
        !     call CalculateCompositionSUBG(iSolnIndex=l,dMolesPairs=dMolesPairs,lPrint=.TRUE.)
        ! end if

        if (allocated(iTempSpecies)) deallocate(iTempSpecies)
        if (allocated(dTempSpecies)) deallocate(dTempSpecies)

        write(1,*) '      "species": {'
        k = iLast - iFirst + 1
        allocate(iTempSpecies(k), dTempSpecies(k))
        nCutOff = k
        select case (cSolnPhaseType(l))
            case ('IDMX', 'RKMP', 'RKMPM', 'QKTO')

                ! Initialize temporary variables:
                dTempSpecies(1:k) = dmolFraction(iFirst:iLast)

                ! Sort species in phase:
                call SortPick(k, dTempSpecies, iTempSpecies)

                LOOP_CutOffX: do i = 1, k
                    ! Convert relative species index to an absolute index:
                    c = iTempSpecies(i) + iFirst - 1

                    if (dMolFraction(c) < dCutOff) then
                        nCutOff = i - 1
                        exit LOOP_CutOffX
                    end if
                end do LOOP_CutOffX
            case default
                ! The species in this phase will not be sorted.
                do i = 1, k
                    iTempSpecies(i) = i
                end do
        end select

        ! The minimum number of species that will be printed is 2:
        nCutOff = MAX(nCutOff, 2)

        k    = LEN_TRIM(cSpeciesName(iFirst)) - 1
        nMax = MAX(k, nMax)

        ! Print species:
        do i = iFirst, iFirst + nCutOff - 1
            c = iTempSpecies(i-iFirst+1) + iFirst - 1
            if (i < iFirst + nCutOff - 1) then
              write(1,*) '        "', TRIM(cSpeciesName(c)), '": ', dMolFraction(c), ","
            else
              write(1,*) '        "', TRIM(cSpeciesName(c)), '": ', dMolFraction(c)
            end if
        end do
        write(1,*) '      }'

        if (j < nSolnPhases) then
          write(1,*) '    },'
        else
          write(1,*) '    }'
        end if

        ! Check if this phase is represented by the Compound Energy Formalism:
        IF_SUBL: if ((csolnPhaseType(l) == 'SUBL').OR.(csolnPhaseType(l) == 'SUBLM')) then

            ! Store the index # of the charged phase:
            iChargedPhaseID = iPhaseSublattice(l)

            ! Loop through sublattices:
            LOOP_Sub: do s = 1, nSublatticePhase(iChargedPhaseID)

                write(1,'(A18,I1,A30,F6.3)') 'Sublattice ', s, '; stoichiometric coefficient: ', &
                    dStoichSublattice(iChargedPhaseID,s)

                cDummy  = '{ '
                cDummyB = ' '

                ! Loop through constituents in sublattice s:
                LOOP_SubCon: do c = 1, nConstituentSublattice(iChargedPhaseID,s)

                    if (c == nConstituentSublattice(iChargedPhaseID,s)) cDummyB = ' }'

                    write(1,'(A20,A8,F9.4,A4,A2)') cDummy, cConstituentNameSUB(iChargedPhaseID,s,c), &
                            dSiteFraction(iChargedPhaseID,s,c), ' ', cDummyB
                end do LOOP_SubCon

                if (s /= nSublatticePhase(iChargedPhaseID)) write(1,*)
            end do LOOP_Sub

        end if IF_SUBL

    end do LOOP_SolnStable

    write(1,*) '  },'

    ! Deallocate allocatable arrays:
    deallocate(dTempVec,iTempVec)

end subroutine WriteJSONSolnPhase

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

subroutine WriteJSONPureConPhase

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModulePrintUtils

    implicit none

    integer                                 :: i, j
    integer,    dimension(:),   allocatable :: iTempVec
    real(8),    dimension(:),   allocatable :: dTempVec

    ! Allocate arrays to sort pure condensed phases:
    allocate(dTempVec(nConPhases), iTempVec(nConPhases))

    ! Initialize variables:
    dTempVec(1:nConPhases) = dMolesPhase(1:nConPhases)
    iTempVec(1:nConPhases) = iAssemblage(1:nConPhases)

    write(1,*) '  "pure condensed phases": {'

    ! Sort pure condensed phases:
    call SortPick(nConPhases, dTempVec, iTempVec)

    ! Loop through stable pure condensed phases:
    LOOP_PureStable: do i = 1, nConPhases
        j = iTempVec(i)
        if (i < nConPhases) then
          write(1,*) '    "', TRIM(cSpeciesName(iAssemblage(j))), '": ', dMolesPhase(j), ","
        else
          write(1,*) '    "', TRIM(cSpeciesName(iAssemblage(j))), '": ', dMolesPhase(j)
        end if
    end do LOOP_PureStable

    write(1,*) '  },'

    ! Deallocate allocatable arrays:
    deallocate(iTempVec,dTempVec)

end subroutine WriteJSONPureConPhase
