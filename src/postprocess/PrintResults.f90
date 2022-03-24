
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    PrintResults.f90
    !> \brief   Print results to screen in a style similar to FactSage.
    !> \author  M.H.A. Piro
    !> \date    January 14, 2013
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/14/2013      M.H.A. Piro         Original code.
    !   02/05/2013      M.H.A. Piro         Sort solution species by the mole fraction and apply a cut off
    !                                        value.
    !   02/05/2013      M.H.A. Piro         Apply variable formatting to pure condensed phases to constrain
    !                                        output to 5 significant figures.
    !
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to print the results to screen in a style similar to FactSage.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nConPhases        An integer scalar representing the number of stable pure condensed phases.
    ! nSolnPhases       An integer scalar representing the number of stable solution phases.
    ! iAssemblage       An integer vector representing the absolute indices of stable phases.
    ! dMolesPhase       A double real vector representing the number of moles for each stable phase in the system.
    ! dMolFraction      A double real vector representing the mole fraction of each species in the system.
    ! dGibbsEnergySys   A double real scalar representing the integral Gibbs energy of the system.
    ! dGEMFunctionNorm  A double real scalar representing the functional norm applied in Gibbs energy.
    ! dCutOff           A double real scalar representing the minimum mole fraction to be cut-off from output.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine PrintResults

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i


    ! Return if the print results mode is zero.
    if (iPrintResultsMode == 0) return

    ! Only proceed for a successful calculation:
    IF_PASS: if (INFOThermo == 0) then

        print *
        print *, '==========================================================='
        print *, '|                  THERMOCHIMICA RESULTS                  |'
        print *, '==========================================================='
        print *

        ! Print the results for solution phases:
        call PrintResultsSolnPhase

        ! Print the results for pure condensed phases:
        call PrintResultsPureConPhase

        print *, '==========================================================='
        print *, '|                    System properties                    |'
        print *, '==========================================================='
        print *

        if ((dPressure < 1D3).AND.(dPressure > 1D-1)) then
            print '(A16,F10.2,A4)', 'Temperature = ', dTemperature, ' [K]'
            print '(A16,F10.4,A6)', 'Pressure    = ', dPressure, ' [atm]'
        else
            print '(A16,F11.2,A4)', 'Temperature = ', dTemperature, ' [K]'
            print '(A16,ES11.3,A6)', 'Pressure    = ', dPressure, ' [atm]'
        end if
        print *
        print *, ' System Component | Mass [mol] | Chemical potential [J/mol]'
        print *, '------------------|------------|---------------------------'
        do i = 1, nElements
            print '(A18,A3,ES10.4,A3, ES14.6)', cElementName(i), ' | ', dMolesElement(i), ' | ', &
                dElementPotential(i) * dIdealConstant * dTemperature
        end do
        print *
        print '(A26,ES12.5,A4)', ' Integral Gibbs energy = ', dGibbsEnergySys, ' [J]'
        if (lHeatCapacityEntropyEnthalpy) then
            print '(A26,ES12.5,A6)', ' Entropy               = ', dEntropy, ' [J/K]'
            print '(A26,ES12.5,A4)', ' Enthalpy              = ', dEnthalpy, ' [J]'
            print '(A26,ES12.5,A6)', ' Heat capacity         = ', dHeatCapacity, ' [J/K]'
        end if

        if (iPrintResultsMode == 2) then
            print '(A26, ES12.5,A11)', ' Functional norm       = ', dGEMFunctionNorm, ' [unitless]'
        end if

        ! Provide additional details if in advanced post-processing mode:
        if (iPrintResultsMode == 2) then
            print *
            print '(A38,I3)', ' # of stable pure condensed phases = ', nConPhases
            print '(A38,I3)', ' # of stable solution phases       = ', nSolnPhases

        end if

        print *
        print *, '==========================================================='
        print *

    else
        ! Do nothing, let the debugger take over.

    end if IF_PASS

    return

end subroutine PrintResults

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to print the results from solution
    ! phases.
    !
    !---------------------------------------------------------------------------

subroutine PrintResultsSolnPhase

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                                 :: c, i, j, k, l, s, iFirst, iLast, iChargedPhaseID, nMax
    integer,    dimension(:),   allocatable :: iTempVec, iTempSpecies
    real(8)                                 :: dCutOff, dTemp, Tcritical, B, StructureFactor, dMolesPairs
    real(8),    dimension(:),   allocatable :: dTempVec, dTempSpecies
    character(2)                            :: cDummy
    character(30)                           :: cDummyB, FMTA, FMTB
    character(30), dimension(:),allocatable :: cPair
    real(8), dimension(:), allocatable      :: dPair


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
        !if ((dMolesPhase(k) >= 1D4).OR.(dMolesPhase(k) <= 1D-1)) then
        if ((dMolesPhase(k) >= 999.95).OR.(dMolesPhase(k) <= 1D-1)) then
            print '(A3,ES10.4,A5,A12)', cDummy, dMolesPhase(k), ' mol ', cSolnPhaseName(l)
        elseif ((dMolesPhase(k) < 0.99949).AND.(dMolesPhase(k) > 1D-1)) then
            ! Format the output so that there are 5 significant figures (k = 6 because it includes ".").
            i = 7
            write (FMTA, *) i
            dTemp = DLOG10(dMolesPhase(k))
            i = i - INT(dTemp) - 2
            write (FMTB, *) i
            dTemp = dMolesPhase(k)
            write (FMTA, "(F" // ADJUSTL(FMTA) // "." // ADJUSTL(FMTB) // ")") dTemp
            print '(A6,A7,A5,A15)', cDummy, FMTA, ' mol ', cSolnPhaseName(l)
        else
            i = 6
            write (FMTA, *) i
            dTemp = DLOG10(dMolesPhase(k))
            i = i - INT(dTemp) - 2
            write (FMTB, *) i
            dTemp = dMolesPhase(k)
            write (FMTA, "(F" // ADJUSTL(FMTA) // "." // ADJUSTL(FMTB) // ")") dTemp
            print '(A7,A6,A5,A15)', cDummy, FMTA, ' mol ', cSolnPhaseName(l)
        end if

        if ((cSolnPhaseType(l) == 'SUBLM') .OR. (cSolnPhaseType(l) == 'RKMPM')) then
            Tcritical = 0D0
            B = 0D0
            call CompMagneticTemperatureMoment(l,Tcritical,B)
            StructureFactor = dCoeffGibbsMagnetic(iFirst,3)
            if (Tcritical < 0D0) then
                Tcritical = -Tcritical * StructureFactor
                print '(A27,F10.2,A4)', 'Neel temperature = ', Tcritical, ' [K]'
            else
                print '(A28,F10.2,A4)', 'Curie temperature = ', Tcritical, ' [K]'
            end if
            if (B < 0D0) then
                B         = -B * StructureFactor
            end if
            print '(A35,F10.5)', 'Magnetic moment per atom = ', B
        end if

        if ((cSolnPhaseType(l)) == 'SUBG' .OR. (cSolnPhaseType(l) == 'SUBQ')) then
            allocate(cPair(nPairsSRO(iPhaseSublattice(l),1)))
            allocate(dPair(nPairsSRO(iPhaseSublattice(l),1)))
            call CalculateCompositionSUBG(l,dMolesPairs,.TRUE.,cPair,dPair)
            deallocate(cPair,dPair)
            print *, '   Quadruplet fractions:'
        end if

        if (allocated(iTempSpecies)) deallocate(iTempSpecies)
        if (allocated(dTempSpecies)) deallocate(dTempSpecies)

        k = iLast - iFirst + 1
        allocate(iTempSpecies(k), dTempSpecies(k))
        select case (cSolnPhaseType(l))

            case ('IDMX', 'RKMP', 'RKMPM', 'QKTO')
                ! Initialize temporary variables:
                dTempSpecies(1:k) = dmolFraction(iFirst:iLast)

                ! Sort species in phase:
                call SortPick(k, dTempSpecies, iTempSpecies)

                LOOP_CutOffX: do i = 1, k
                    ! Convert relative species index to an absolute index:
                    c = iTempSpecies(i) + iFirst - 1
                end do LOOP_CutOffX

            case default
                ! The species in this phase will not be sorted.
                do i = 1, k
                    iTempSpecies(i) = i
                end do
        end select

        ! Print middle species:
        do i = iFirst, iLast
            c = iTempSpecies(i-iFirst+1) + iFirst - 1
            cDummyB  = TRIM(cSpeciesName(c))
            k    = LEN_TRIM(cSpeciesName(c)) - 1
            nMax = MAX(k, nMax)
            if (i == iFirst) then
                cDummy = '{ '
            else
                cDummy = '+ '
            end if
            if (i == iLast) cDummyB(nMax+2:nMax+3) = '}'
            if (dmolFraction(c) >= 1D-1) then
                print '(A20,F7.5,A3,A35)', cDummy, dMolFraction(c), ' ', cDummyB
            else
                print '(A20,ES10.4,A35)', cDummy, dMolFraction(c), cDummyB
            end if
        end do
        print *

        ! Check if this phase is represented by the Compound Energy Formalism:
        IF_SUBL: if ((csolnPhaseType(l) == 'SUBL') .OR. (csolnPhaseType(l) == 'SUBLM') .OR. &
        (csolnPhaseType(l) == 'SUBI')) then

            ! Store the index # of the charged phase:
            iChargedPhaseID = iPhaseSublattice(l)

            print *, '      ------------------------------------------------'

            ! Loop through sublattices:
            LOOP_Sub: do s = 1, nSublatticePhase(iChargedPhaseID)

                print '(A18,I1,A30,F6.3)', 'Sublattice ', s, '; stoichiometric coefficient: ', &
                    dStoichSublattice(iChargedPhaseID,s)

                cDummy  = '{ '
                cDummyB = ' '

                ! Loop through constituents in sublattice s:
                LOOP_SubCon: do c = 1, nConstituentSublattice(iChargedPhaseID,s)

                    if (c == nConstituentSublattice(iChargedPhaseID,s)) cDummyB = ' }'

                    if (dSiteFraction(iChargedPhaseID,s,c) >= 0.999949D0) then
                        print '(A20,A8,F9.4,A4,A2)', cDummy, cConstituentNameSUB(iChargedPhaseID,s,c), &
                            dSiteFraction(iChargedPhaseID,s,c), ' ', cDummyB
                    elseif (dSiteFraction(iChargedPhaseID,s,c) >= 1D-1) then
                        print '(A20,A8,F10.5,A3,A2)', cDummy, cConstituentNameSUB(iChargedPhaseID,s,c), &
                            dSiteFraction(iChargedPhaseID,s,c), ' ', cDummyB
                    else
                        print '(A20,A8,ES13.4,A2)', cDummy, cConstituentNameSUB(iChargedPhaseID,s,c), &
                            dSiteFraction(iChargedPhaseID,s,c), cDummyB
                    end if

                    cDummy = '+ '


                end do LOOP_SubCon

                if (s /= nSublatticePhase(iChargedPhaseID)) print *
            end do LOOP_Sub

            print *, '      ------------------------------------------------'
            print *

        end if IF_SUBL

    end do LOOP_SolnStable

    ! Deallocate allocatable arrays:
    deallocate(dTempVec,iTempVec)

end subroutine PrintResultsSolnPhase

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to print the results for pure condensed
    ! phases.
    !
    !---------------------------------------------------------------------------


subroutine PrintResultsPureConPhase

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer                                 :: i, j, k
    integer,    dimension(:),   allocatable :: iTempVec
    real(8)                                 :: dTemp
    real(8),    dimension(:),   allocatable :: dTempVec
    character(2)                            :: cDummy
    character(20)                           :: FMTA, FMTB

    ! Allocate arrays to sort pure condensed phases:
    allocate(dTempVec(nConPhases), iTempVec(nConPhases))

    ! Initialize variables:
    dTempVec(1:nConPhases) = dMolesPhase(1:nConPhases)
    iTempVec(1:nConPhases) = iAssemblage(1:nConPhases)

    ! Initialize prefix:
    if (nSolnPhases > 0) then
        cDummy = '+ '
    else
        cDummy = ' '
    end if

    ! Sort pure condensed phases:
    call SortPick(nConPhases, dTempVec, iTempVec)

    ! Loop through stable pure condensed phases:
    LOOP_PureStable: do i = 1, nConPhases
        j = iTempVec(i)

        if ((dMolesPhase(j) >= 1D4).OR.(dMolesPhase(j) <= 1D-1)) then
            print '(A3,ES10.4,A5,A15)', cDummy, dMolesPhase(j), ' mol ', cSpeciesName(iAssemblage(j))
        elseif ((dMolesPhase(j) < 0.99949).AND.(dMolesPhase(j) > 1D-1)) then
            ! Format the output so that there are 5 significant figures (k = 6 because it includes ".").
            k = 7
            write (FMTA, *) k
            dTemp = DLOG10(dMolesPhase(j))
            k = k - INT(dTemp) - 2
            write (FMTB, *) k
            dTemp = dMolesPhase(j)
            write (FMTA, "(F" // ADJUSTL(FMTA) // "." // ADJUSTL(FMTB) // ")") dTemp
            print '(A6,A7,A5,A15)', cDummy, FMTA, ' mol ', cSpeciesName(iAssemblage(j))
        else
            ! Format the output so that there are 5 significant figures (k = 6 because it includes ".").
            k = 6
            write (FMTA, *) k
            dTemp = DLOG10(dMolesPhase(j))
            k = k - INT(dTemp) - 2
            write (FMTB, *) k
            dTemp = dMolesPhase(j)
            write (FMTA, "(F" // ADJUSTL(FMTA) // "." // ADJUSTL(FMTB) // ")") dTemp
            print '(A7,A6,A5,A15)', cDummy, FMTA, ' mol ', cSpeciesName(iAssemblage(j))
        end if
        print *

        cDummy = '+ '

    end do LOOP_PureStable

    ! Deallocate allocatable arrays:
    deallocate(iTempVec,dTempVec)

end subroutine PrintResultsPureConPhase

    !---------------------------------------------------------------------------
    !                      END - PrintResults.f90
    !---------------------------------------------------------------------------
