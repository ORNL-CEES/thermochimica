
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
    ! cSearchPhase = TRIM(ADJUSTL(cSearchPhase))

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