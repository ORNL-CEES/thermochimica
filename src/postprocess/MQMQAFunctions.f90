
subroutine GetMqmqaMolesPairs(cPhaseName, lcPhaseName, dMolesPairsOut, INFO) &
    bind(C, name="TCAPI_getMqmqaMolesPairs")

    USE ModuleThermo
    USE ModuleThermoIO
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer,       intent(out)               :: INFO
    real(8),       intent(out)               :: dMolesPairsOut
    character(kind=c_char,len=1), target, intent(in) :: cPhaseName(*)
    integer(c_size_t), intent(in), value             :: lcPhaseName
    character(kind=c_char,len=lcPhaseName), pointer  :: fPhaseName
    character(30), dimension(:), allocatable :: cPair
    real(8), dimension(:), allocatable       :: dPair
    integer                                  :: i, j

    call c_f_pointer(cptr=c_loc(cPhaseName), fptr=fPhaseName)

    ! Initialize variables:
    INFO            = 0
    dMolesPairsOut   = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_SOLN: do i = 1, nSolnPhasesSys
            if (fPhaseName == cSolnPhaseName(i)) then
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


subroutine GetMqmqaPairMolFraction(cPhaseName, lcPhaseName, cPairIn, lcPairIn, dMolesPairOut, INFO) &
    bind(C, name="TCAPI_getMqmqaPairMolFraction")

    USE ModuleThermo
    USE ModuleThermoIO
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer,       intent(out)               :: INFO
    real(8),       intent(out)               :: dMolesPairOut
    character(kind=c_char,len=1), target, intent(in) :: cPhaseName(*), cPairIn(*)
    integer(c_size_t), intent(in), value             :: lcPhaseName, lcPairIn
    character(kind=c_char,len=lcPhaseName), pointer  :: fPhaseName
    character(kind=c_char,len=lcPairIn), pointer     :: fPairIn
    integer                                  :: i, j ,k
    real(8)                                  :: dMolesPairs
    character(30), dimension(:), allocatable :: cPair
    real(8), dimension(:), allocatable       :: dPair


    call c_f_pointer(cptr=c_loc(cPhaseName), fptr=fPhaseName)
    call c_f_pointer(cptr=c_loc(cPairIn), fptr=fPairIn)

    ! Initialize variables:
    INFO            = 0
    dMolesPairOut   = 0D0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_SOLN: do i = 1, nSolnPhasesSys
            if (fPhaseName == cSolnPhaseName(i)) then
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
                if (fPairIn == cPair(k)) then
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

subroutine GetMqmqaNumberPairsQuads(cPhaseName,lcPhaseName, nPairs, nQuads, INFO) &
    bind(C, name="TCAPI_getMqmqaNumberPairsQuads")

    USE ModuleThermo
    USE ModuleThermoIO
    USE,INTRINSIC :: ISO_C_BINDING

    implicit none

    integer,       intent(out)               :: INFO, nPairs, nQuads
    character(kind=c_char,len=1), target, intent(in) :: cPhaseName(*)
    integer(c_size_t), intent(in), value             :: lcPhaseName
    character(kind=c_char,len=lcPhaseName), pointer  :: fPhaseName
    integer                                  :: i, j

    call c_f_pointer(cptr=c_loc(cPhaseName), fptr=fPhaseName)

    ! Initialize variables:
    INFO            = 0
    nPairs          = 0
    nQuads          = 0

    ! Only proceed if Thermochimica solved successfully:
    if (INFOThermo == 0) then

        ! Loop through stable soluton phases to find the one corresponding to the
        ! solution phase being requested:
        j = 0
        LOOP_SOLN: do i = 1, nSolnPhasesSys
            if (fPhaseName == cSolnPhaseName(i)) then
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
