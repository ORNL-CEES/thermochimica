
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckCompounds.f90
    !> \brief   Calculate stoichiometry and masses in terms of compounds rather than elements.
    !> \author  M. Poschmann
    !> \date    Apr. 19, 2019
    !> \sa      CheckSystem.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/19/2019      M. Poschmann        File creation.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details To avoid a situation in which the Gibbs tangent plane is underdetermined (or just for convenience)
    !! input mass may be specified in terms of compounds rather than pure elements. If input is given in terms of
    !! compounds, all following calculations will progress in terms of the amounts of those compounds. In other
    !! words, compounds will replace elements everywhere. To do this, the stoichiometry matrix has to be
    !! recalculated in therms of compounds, and the element masses and names replaced by compound masses and names.
    !! Some index recalculation is required as well. All these changes are made here (called shortly after the
    !! start of CheckSystem) to simplify this process (except the name changes, see CheckSystem).
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nCompounds                        A scalar integer specifying the number of compounds.
    ! dCompoundStoich                   A double array supplied by the user which specifies the stoichiometry of
    !                                    each of the compounds to be used.
    !
    ! dElementMass                      Total mass of each element, where the coefficient corresponds to the
    !                                    atomic number (e.g., dMolesElement(92) refers to uranium). This gets
    !                                    by the masses of compounds (stored at the end of the array).
    ! dStoichSpecies                    A double array specifying the stoichiometry of each species. Usually this
    !                                    would be in terms of elements, but is converted here to compounds.
    ! iElementSystem                    An integer array that points to the index of each element in the system.
    !                                    Here the element indices are replaced by compound indices.
    ! INFOThermo                        A scalar integer that indicates a successful exit or identifies an error.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine CheckCompounds

    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer                                 :: i, j, k, lwork
    integer, dimension(:), allocatable      :: iElementSystemTemp
    real(8)                                 :: dCompoundUnitMass, dSum
    real(8), dimension(:), allocatable      :: dStoichSpeciesTemp
    real(8), dimension(:), allocatable      :: work
    real(8), dimension(:,:), allocatable    :: dCompoundStoichSmall, dStoichSpeciesCompounds, A

    ! Allocate LAPACK workspace
    lwork = 32 * nCompounds
    allocate(work(lwork))

    ! Allocate temporary arrays.
    allocate(iElementSystemTemp(nCompounds))
    allocate(dStoichSpeciesTemp(nElementsCS))
    allocate(dStoichSpeciesCompounds(nSpeciesCS,nCompounds))
    allocate(dCompoundStoichSmall(nCompounds,nElementsCS),A(nElementsCS,nCompounds))
    dCompoundStoichSmall = 0D0

    ! Check elements in compounds
    do i = 1, nCompounds
        LOOP_AllElements: do j = 1, nElementsPT
            if (dCompoundStoich(i,j) > 0d0) then
                do k = 1, nElementsCS
                    if (iElementSystem(k) == j) then
                        ! Create stoichiometry matrix nElementsCS in length
                        dCompoundStoichSmall(i,k) = dCompoundStoich(i,j)
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
    ! Do not re-calculate this stoichiometry matrix if this is a reinited calculation
    if (.NOT. lCompoundStoichCalculated) then
        do i = 1, nSpeciesCS
            do j = 1, nElementsCS
                ! Make vector to pass stoichiometry to LAPACK
                dStoichSpeciesTemp(j) = dStoichSpeciesCS(i,j)
            end do
            do j = 1, nCompounds
                do k = 1, nElementsCS
                    ! Make matrix for LAPACK
                    A(k,j) = dCompoundStoichSmall(j,k)
                end do
            end do
            ! Linear solve on rectangular (overdetermined) system
            call DGELS( 'N', nElementsCS, nCompounds, 1, A, nElementsCS, &
                             dStoichSpeciesTemp, nElementsCS, work, lwork, INFO )
            if (INFO > 0) then
                INFOThermo = 41
                return
            end if
            do j = 1, nCompounds
                ! Copy LAPACK solution vector to new stoichiometry
                dStoichSpeciesCompounds(i,j) = dStoichSpeciesTemp(j)
            end do
        end do
    end if

    ! Use compounds instead of elements everwhere (i.e. write permanent variables)
    deallocate(iElementSystem)
    allocate(iElementSystem(nElementsCS))
    iElementSystem = 0
    do i = 1, nCompounds
        iElementSystem(i) = iElementSystemTemp(i)
    end do

    ! Make sure that the sums of compound stoichiometries makes sense
    ! Do not re-calculate this stoichiometry matrix if this is a reinited calculation
    if (.NOT. lCompoundStoichCalculated) then
        LOOP_checkCompoundStoich: do i = 1, nSpeciesCS
            do j = 1, nElementsCS
                dSum = 0D0
                do k = 1, nCompounds
                    dSum = dSum + (dCompoundStoichSmall(k,j) * dStoichSpeciesCompounds(i,k))
                end do
                if (ABS(dSum - dStoichSpeciesCS(i,j)) > 1D-8) then
                    dStoichSpeciesCompounds(i,1:nCompounds) = 0D0
                    cycle LOOP_checkCompoundStoich
                end if
            end do
        end do LOOP_checkCompoundStoich

        ! Save new stoichiometry matrix in terms of compounds
        deallocate(dStoichSpeciesCS)
        allocate(dStoichSpeciesCS(nSpeciesCS,nElementsCS))
        LOOP_compoundStoich: do i = 1, nSpeciesCS
            do j = 1, nCompounds
                if (dStoichSpeciesCompounds(i,j) > 1D-8) then
                    dStoichSpeciesCS(i,j) = dStoichSpeciesCompounds(i,j)
                elseif (dStoichSpeciesCompounds(i,j) < -1D-8) then
                    dStoichSpeciesCS(i,1:nCompounds) = 0D0
                    cycle LOOP_compoundStoich
                else
                    dStoichSpeciesCS(i,j) = 0D0
                end if
            end do
        end do LOOP_compoundStoich
        lCompoundStoichCalculated = .TRUE.
    end if
    nElemOrComp = nCompounds

    dElementMass = 0D0
    do i = 1, nCompounds
        dElementMass(i + nElementsPT) = dCompoundMass(i)
    end do

    ! Clear temporary arrays
    deallocate(dCompoundStoichSmall,dStoichSpeciesTemp,dStoichSpeciesCompounds,A,work,iElementSystemTemp)

end subroutine CheckCompounds
