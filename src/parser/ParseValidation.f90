    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ParseValidation.F90
    !> \brief   Parse Validation data and run phase transition tests.
    !> \author  A.E.F. Fitzsimmons
    !
    !> \todo    Multiple tolerance check
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    06/17/2024    A.E.F. Fitzsimmons   Original Code
    !
    ! Purpose:
    ! ========
    !> \details 
    !! 
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine ParseValidation(cTestCSV,lPass)
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleTesting

    implicit none

    character(*), intent(in)::  cTestCSV

    integer :: INFO, nElementTest, nParamTest, nTestNumber, iTransitions, i, j, iPos, iPos2, iTransCount
    integer, allocatable :: iElementIndex(:)
    real(8) :: dTempMax, dTempMin, dTempTolerance
    real(8), dimension(nMaxPhaseTransition) :: dPhaseTransitionTemp, dTestTransitionTemp, dTestTolerance 
    real(8), dimension(2 * nMaxPhaseTransition) :: dTransition
    real(8), allocatable :: dMass(:)
    character(1024) :: cTestCSVLen
    character(250) :: cLines
    character(3), dimension(0:nElementsPT):: cElementNamePT
    character(3), allocatable :: cTestElement(:)
    logical, intent(out) :: lPass

    ! Initialize variables:
    INFOThermo   = 0
    INFO         = 0
    nTestNumber  = 0
    iTransitions = 0
    iTransCount  = 0
    lPass        = .FALSE.
    dPhaseTransitionTemp = 0D0
    dTestTransitionTemp  = 0D0
    dTestTolerance       = 0D0
    dTransition          = 0D0
    dTempTolerance       = 1D0
    cTestCSVLen  = cTestCSV(1:min(1024,len(cTestCSV)))
    
    ! Attempt to open the ChemSage datafile:
    open (UNIT = 2, FILE = cTestCSVLen, STATUS = 'old', ACTION = 'read', IOSTAT = INFO)

    ! Record an error if there are issues opening the data-file:
    if (INFO /= 0) then
        INFOThermo = 60 
        return
    end if

    !Skip Line 1
    read(2, '(A)', IOSTAT = INFO) cLines

    !Read line 2 Get Number of elements and parameters
    read(2, *, IOSTAT = INFO) nElementTest, nParamTest
    if (INFO /= 0) then
        INFOThermo = 61
    end if
    allocate(cTestElement(nElementTest), iElementIndex(nElementTest), dMass(nElementTest)) 
    nElements = nElementTest

    !Init Arrays
    iElementIndex = 0
    dMass = 0D0
    cTestElement = " "
    call GetElementName(cElementNamePT)

    !Reading Elements required
    read(2, *, IOSTAT = INFO) cTestElement(1:nElementTest)
    print *, cTestElement
    if (cTestElement(nElementTest) == " ") then
        INFOThermo = 67
        return
    end if
    if (INFO /= 0) then
        INFOThermo = 62
        return
    end if

    !Getting Periodic table index
    do j = 1, nElementTest
        do i = 1, nElementsPT
            if(trim(cTestElement(j)) == trim(cElementNamePT(i))) iElementIndex(j) = i
        end do 
    end do

    !Testing loop
    do i = 1, nParamTest

        read (2, *, IOSTAT = INFO) nTestNumber, iTransitions, dMass, dTempMin, dTempMax
        if (INFO /= 0) then
            INFOThermo = 63
            return
        end if
        
        backspace(2, IOSTAT = INFO)
        read (2, *, IOSTAT = INFO) nTestnumber, iTransitions, dMass, dTempMin, dTempMax, dTransition(1:2 * iTransitions)
        if (INFO /= 0) then
            INFOThermo = 64
            return
        end if

        !Minimum temperature is less than 0
        if (dTempMin < 0 .OR. dTempMin > dTempMax) then
            INFOThermo = 66
            return
        end if

        !Set Element mass
        do j = 1, nElementTest
            dElementMass(iElementIndex(j)) = dMass(j)
        end do
        
        !Seperate dTransition get Temperature and Tolerance
        iPos = 1
        iPos2 = 1
        do j = 1, SIZE(dTransition)
            if (MOD(j, 2) == 1) then
                dTestTransitionTemp(iPos) = dTransition(j)
                iPos = iPos + 1
            else
                dTestTolerance(iPos2) = dTransition(j)
                iPos2 = iPos2 + 1
            end if
        end do

        call PhaseTransition(dTempMin, dTempMax, dTempTolerance, dPhaseTransitionTemp, iTransCount)
        call checkTransitionTest(dPhaseTransitionTemp, dTestTransitionTemp, 1D0, lPass)
        
        if(.NOT. lPass) then
            print *, "Test Number: ", nTestNumber, " failed."
            EXIT
        end if

        !Reset test variables
        iTransCount          = 0
        iPos                 = 1
        iPos2                = 1
        dPhaseTransitionTemp = 0D0
        dTestTransitionTemp  = 0D0
        dTransition          = 0D0

    end do 
    
    deallocate(cTestElement, iElementIndex, dMass)
    close (INFO)

end subroutine ParseValidation        



