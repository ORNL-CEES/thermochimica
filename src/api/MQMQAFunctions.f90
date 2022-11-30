
subroutine GetMqmqaMolesPairs(cPhase, dMolesPairsOut, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    character(*),  intent(in)                :: cPhase
    integer,       intent(out)               :: INFO
    integer                                  :: i, j
    real(8),       intent(out)               :: dMolesPairsOut
    character(25)                            :: cSearchPhase
    character(30), dimension(:), allocatable :: cPair
    real(8), dimension(:), allocatable       :: dPair
    ! integer,  intent(in)                     :: lcPhase

    ! It seems that a single string argument works only without length supplied,
    ! but if there are multiple strings, they require lengths.
    cSearchPhase = cPhase
    ! cSearchPhase = TRIM(ADJUSTL(cSearchPhase))

    ! Initialize variables:
    INFO            = 0
    dMolesPairsOut   = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Remove trailing blanks:
        cSearchPhase = TRIM(ADJUSTL(cSearchPhase))

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_SOLN: do i = 1, nSolnPhasesSys
            if (cSearchPhase == cSolnPhaseName(i)) then
                ! Solution phase found.  Record integer index and exit loop.
                j = i
                exit LOOP_SOLN
            end if

        end do LOOP_SOLN

        ! Check to make sure that the solution phase was found:
        IF_SOLN: if (j > 0) then
            allocate(cPair(nPairsSRO(iPhaseSublattice(j),1)))
            allocate(dPair(nPairsSRO(iPhaseSublattice(j),1)))
            ! Call function to get moles of pairs
            call CalculateCompositionSUBG(j,dMolesPairsOut,.FALSE.,cPair,dPair)
            deallocate(cPair,dPair)
        else
           ! This solution phase was not found.  Report an error:
           INFO = 1
        end if IF_SOLN

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetMqmqaMolesPairs

subroutine GetMqmqaPairMolFraction(cPhase, lcPhase, cPairIn, lcPairIn, dMolesPairOut, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    character(*),  intent(in)                :: cPhase, cPairIn
    integer,       intent(out)               :: INFO
    integer                                  :: i, j ,k
    real(8),       intent(out)               :: dMolesPairOut
    real(8)                                  :: dMolesPairs
    character(25)                            :: cSearchPhase, cSearchPair
    character(30), dimension(:), allocatable :: cPair
    real(8), dimension(:), allocatable       :: dPair
    integer,  intent(in)                     :: lcPhase, lcPairIn

    ! It seems that a single string argument works only without length supplied,
    ! but if there are multiple strings, they require lengths.
    cSearchPhase = cPhase
    cSearchPair  = cPairIn
    cSearchPhase = TRIM(cSearchPhase(1:min(30,lcPhase)))
    cSearchPair = TRIM(cSearchPair(1:min(30,lcPairIn)))

    ! Initialize variables:
    INFO            = 0
    dMolesPairOut   = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Remove trailing blanks:
        cSearchPhase = TRIM(ADJUSTL(cSearchPhase))
        cSearchPair  = TRIM(ADJUSTL(cSearchPair))

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_SOLN: do i = 1, nSolnPhasesSys
            if (cSearchPhase == cSolnPhaseName(i)) then
                ! Solution phase found.  Record integer index and exit loop.
                j = i
                exit LOOP_SOLN
            end if

        end do LOOP_SOLN

        ! Check to make sure that the solution phase was found:
        IF_SOLN: if (j > 0) then
            allocate(cPair(nPairsSRO(iPhaseSublattice(j),1)))
            allocate(dPair(nPairsSRO(iPhaseSublattice(j),1)))
            ! Call function to get moles of pairs
            call CalculateCompositionSUBG(j,dMolesPairs,.FALSE.,cPair,dPair)
            ! Find selected pair
            INFO = 2
            LOOP_pairs: do k = 1, nPairsSRO(iPhaseSublattice(j),1)
                if (cSearchPair == cPair(k)) then
                    INFO = 0
                    dMolesPairOut = dPair(k)
                    exit LOOP_pairs
                end if
            end do LOOP_pairs
            deallocate(cPair,dPair)
        else
           ! This solution phase was not found.  Report an error:
           INFO = 1
        end if IF_SOLN

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetMqmqaPairMolFraction

subroutine GetMqmqaNumberPairsQuads(cPhase, nPairs, nQuads, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    character(*),  intent(in)                :: cPhase
    integer,       intent(out)               :: INFO, nPairs, nQuads
    integer                                  :: i, j
    character(25)                            :: cSearchPhase
    ! integer,  intent(in)                     :: lcPhase

    ! GetMqmqaNumberPairsQuads was given a mismatched number of arguments between Fortran and C
    ! The last argument in C is the string length.
    cSearchPhase = cPhase
    ! cSearchPhase = TRIM(cSearchPhase(1:min(30,lcPhase)))

    ! Initialize variables:
    INFO            = 0
    nPairs          = 0
    nQuads          = 0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Remove trailing blanks:
        cSearchPhase = TRIM(ADJUSTL(cSearchPhase))

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_SOLN: do i = 1, nSolnPhasesSys
            if (cSearchPhase == cSolnPhaseName(i)) then
                ! Solution phase found.  Record integer index and exit loop.
                j = i
                exit LOOP_SOLN
            end if

        end do LOOP_SOLN

        ! Check to make sure that the solution phase was found:
        IF_SOLN: if (j > 0) then
            ! Get numbers of pairs and quadruplets
            nPairs = nPairsSRO(iPhaseSublattice(j),1)
            nQuads = nSpeciesPhase(j) - nSpeciesPhase(j - 1)
        else
           ! This solution phase was not found.  Report an error:
           INFO = 1
        end if IF_SOLN

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetMqmqaNumberPairsQuads

subroutine GetMqmqaConstituentFraction(cPhase, iSublattice, cConstituent, dConstituentFractionOut, INFO)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    character(*),  intent(in)  :: cPhase, cConstituent
    character(25)              :: cPhaseSearch, cConstituentSearch
    integer,       intent(in)  :: iSublattice
    integer,       intent(out) :: INFO
    integer                    :: i, j, k, l
    integer                    :: iSolnIndex, iConIndex, iSPI
    integer                    :: iFirst, iLast, nSub
    real(8),       intent(out) :: dConstituentFractionOut
    real(8)                    :: dConstituentFraction, dSum, dZa, dZb, dZx, dZy

    ! Initialize output variables:
    INFO                    = 0
    dConstituentFractionOut = 0D0

    ! There are only 2 sublattices in MQMQA, index cannot be otherwise
    if (.NOT. ((iSublattice == 1) .OR. (iSublattice == 2))) then
        INFO = 3
        return
    end if

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Remove trailing blanks:
        cPhaseSearch       = TRIM(ADJUSTL(cPhase))
        cConstituentSearch = TRIM(ADJUSTL(cConstituent))

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        iSolnIndex = 0
        LOOP_SOLN: do i = 1, nSolnPhasesSys
            if (cPhaseSearch == cSolnPhaseName(i)) then
                ! Solution phase found.  Record integer index and exit loop.
                iSolnIndex = i
                exit LOOP_SOLN
            end if
        end do LOOP_SOLN

        ! Check to make sure that the solution phase was found:
        if (iSolnIndex > 0) then
            ! Only proceed if the correct phase type is selected:
            if (.NOT. (cSolnPhaseType(iSolnIndex) == 'SUBG' .OR. cSolnPhaseType(iSolnIndex) == 'SUBQ')) then
                INFO = 2
                return
            end if
        else
            ! This solution phase was not found.  Report an error:
            INFO = 1
            return
        end if
        
        iSPI = iPhaseSublattice(iSolnIndex)
        nSub = nConstituentSublattice(iSPI,iSublattice)

        ! Look for requested constituent
        iConIndex = 0
        LOOP_CON: do i = 1, nSub
            if (cConstituentSearch == cConstituentNameSUB(iSPI,iSublattice,i)) then
                ! Solution phase found.  Record integer index and exit loop.
                iConIndex = i
                exit LOOP_CON
            end if
        end do LOOP_CON

        ! Check to make sure that the constituent was found:
        if (iConIndex <= 0) then
            ! This constituent phase was not found.  Report an error:
            INFO = 4
            return
        end if

        ! Define temporary variables for sake of convenience:
        iFirst = nSpeciesPhase(iSolnIndex-1) + 1
        iLast  = nSpeciesPhase(iSolnIndex)

        dSum = 0D0
        dConstituentFraction = 0D0
        if (iSublattice == 1) then
            ! Cations:
            do i = 1, nSub
                do k = 1, nPairsSRO(iSPI,2)
                    l = iFirst + k - 1
                    dZa = dCoordinationNumber(iSPI,k,1)
                    dZb = dCoordinationNumber(iSPI,k,2)
                    if (i == iPairID(iSPI,k,1))  then
                        dSum = dSum + (dMolFraction(l) / dZa)
                        if (i == iConIndex) dConstituentFraction = dConstituentFraction + (dMolFraction(l) / dZa)
                    end if
                    if (i == iPairID(iSPI,k,2))  then
                        dSum = dSum + (dMolFraction(l) / dZb)
                        if (i == iConIndex) dConstituentFraction = dConstituentFraction + (dMolFraction(l) / dZb)
                    end if
                end do
            end do
        else if (iSublattice == 2) then
            ! Anions:
            do i = 1, nSub
                j = i + nConstituentSublattice(iSPI,iSublattice)
                do k = 1, nPairsSRO(iSPI,2)
                    l = iFirst + k - 1
                    dZx = dCoordinationNumber(iSPI,k,3)
                    dZy = dCoordinationNumber(iSPI,k,4)
                    if (j == iPairID(iSPI,k,3))  then
                        dSum = dSum + (dMolFraction(l) / dZx)
                        if (i == iConIndex) dConstituentFraction = dConstituentFraction + (dMolFraction(l) / dZx)
                    end if
                    if (j == iPairID(iSPI,k,4))  then
                        dSum = dSum + (dMolFraction(l) / dZy)
                        if (i == iConIndex) dConstituentFraction = dConstituentFraction + (dMolFraction(l) / dZy)
                    end if
                end do
            end do
        end if

        dConstituentFractionOut = dConstituentFraction / dSum

    else
        ! Record an error with INFO if INFOThermo /= 0.
        INFO = -1
    end if

    return

end subroutine GetMqmqaConstituentFraction