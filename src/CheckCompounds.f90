subroutine CheckCompounds

    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer                                 :: i, j, k, lwork, nElementsTemp
    integer, dimension(:), allocatable      :: iElementSystemTemp
    real(8)                                 :: dCompoundUnitMass
    real(8), dimension(:), allocatable      :: dStoichSpeciesTemp
    real(8), dimension(:), allocatable      :: work
    real(8), dimension(:,:), allocatable    :: dCompoundStoichSmall, dStoichSpeciesCompounds, A

    ! Allocate LAPACK workspace
    lwork = 32 * nCompounds
    allocate(work(lwork))

    allocate(iElementSystemTemp(nCompounds))
    allocate(dStoichSpeciesTemp(nElementsCS))
    allocate(dStoichSpeciesCompounds(nSpeciesCS,nCompounds))
    allocate(dCompoundStoichSmall(nCompounds,nElementsCS),A(nElementsCS,nCompounds))
    dCompoundStoichSmall = 0D0

    ! Check elements in compounds
    nElementsTemp = 0
    do i = 1, nCompounds
        LOOP_AllElements: do j = 1, nElementsPT
            if (dCompoundStoich(i,j) > 0d0) then
                do k = 1, nElementsCS
                    if (iElementSystem(k) == j) then
                        dCompoundStoichSmall(i,k) = dCompoundStoich(i,j)
                        nElementsTemp = nElementsTemp + 1
                        cycle LOOP_AllElements
                    end if
                end do
                ! There is an element in the compound that is not in the dat file
                INFOThermo = 40
                return
            end if
        end do LOOP_AllElements
    end do

    ! Perform a mass conversion:
    do i = 1, nCompounds
        dCompoundUnitMass = 0
        do j = 1, nElementsCS
            dCompoundUnitMass = dCompoundUnitMass + dCompoundStoich(i,iElementSystem(j)) * dAtomicMass(j)
        end do
        select case (cInputUnitMass)
            case ('mass fraction','kilograms','grams','pounds')
                ! Convert mass unit to moles:
                dCompoundMass(i) = dCompoundMass(i) / dCompoundUnitMass
            case ('mole fraction','atom fraction','atoms','moles','gram-atoms')
                ! Do nothing
            case default
                ! The character string representing input units is not recognized.
                INFOThermo = 4
                return
        end select
        iElementSystemTemp(i) = i + nElementsPT
    end do

    ! Convert stoichiometry of all species to compounds
    do i = 1, nSpeciesCS
        do j = 1, nElementsCS
            dStoichSpeciesTemp(j) = dStoichSpeciesCS(i,j)
        end do
        do j = 1, nCompounds
            do k = 1, nElementsCS
                A(k,j) = dCompoundStoichSmall(j,k)
            end do
        end do
        call DGELS( 'N', nElementsCS, nCompounds, 1, A, nElementsCS, &
                         dStoichSpeciesTemp, nElementsCS, work, lwork, INFO )
        if (INFO > 0) then
            INFOThermo = 41
            return
        end if
        do j = 1, nCompounds
            dStoichSpeciesCompounds(i,j) = NINT(dStoichSpeciesTemp(j))
        end do
    end do

    ! Use compounds instead of elements everwhere
    deallocate(iElementSystem)
    allocate(iElementSystem(nCompounds))
    do i = 1, nCompounds
        iElementSystem(i) = iElementSystemTemp(i)
    end do

    deallocate(dStoichSpeciesCS)
    allocate(dStoichSpeciesCS(nSpeciesCS,nCompounds))
    do i = 1, nSpeciesCS
        do j = 1, nCompounds
            dStoichSpeciesCS(i,j) = dStoichSpeciesCompounds(i,j)
        end do
    end do
    nElementsCS = nCompounds

    dElementMass = 0D0
    do i = 1, nCompounds
        dElementMass(i + nElementsPT) = dCompoundMass(i)
    end do

    deallocate(dCompoundStoichSmall,dStoichSpeciesTemp,dStoichSpeciesCompounds,A)

end subroutine CheckCompounds
