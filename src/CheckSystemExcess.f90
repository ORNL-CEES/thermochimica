
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckSystemExcess.f90
    !> \brief   Check the excess terms for the system.
    !> \author  M.H.A. Piro
    !> \date    Mar. 17, 2018
    !> \sa      Thermochimica.f90
    !> \sa      CheckSystem.f90
    !> \todo    This subroutine could be smarter...it is uncessary to perform this exhaustive exercise if the
    !!           system has not changed.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/14/2013      M.H.A. Piro         Original code (relocated from CheckSystem.f90).
    !   01/14/2013      M.H.A. Piro         Add check for SUBL and SUBLM phases.
    !   02/05/2013      M.H.A. Piro         Fix bugs in applying the number of constituents per sublattice
    !                                        and the constituent names in SUBL phases.
    !   09/06/2015      M.H.A. Piro         Fixed bug in checking mixing terms in SUBL when fewer phases
    !                                        are included in the system then from the database.  Specifically,
    !                                        the vector iPhaseSublattice was previously used when
    !                                        iPhaseSublatticeCS should be used instead in checking constituents.
    !   03/17/2018      M.H.A. Piro         Added capability to handle SUBG phases.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check if all mixing parameters should be considered
    !! or if some should be neglected if some of the constituents are not part of the system.
    !! The supported solution model types are represented by cSolnPhaseTypeSupport (ModuleParseCS).
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nParamPhaseCS     An integer vector representing the number of parameters for a phase.
    ! iParamPassCS      An integer vector representing whether a parameter should be considered (1) or not (0).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckSystemExcess

    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver, ONLY: lMiscibility

    implicit none

    integer::  c, i, j, k, l, m, n, s, nCounter, pa, pb, px, py, pe, pf, nRemove, p, iIndex, iCon1, iCon2
    integer, dimension(nElementsCS**2) :: iRemove


    ! Initialize variables:
    nCounter        = 0
    nCountSublatticeCS = 0
    nCountSublattice   = 0
    if (allocated(iPhaseSublattice)) iPhaseSublattice   = 0

    ! Check the mixing parameters:
    LOOP_SolnPhases: do i = 1, nSolnPhasesSysCS

        j = nSpeciesPhaseCS(i-1) + 1
        k = nSpeciesPhaseCS(i)
        l = MAXVAL(iSpeciesPass(j:k))

        if ((cSolnPhaseTypeCS(i) == 'SUBL').OR.(cSolnPhaseTypeCS(i) == 'SUBLM').OR. &
             (cSolnPhaseTypeCS(i) == 'SUBG').OR.(cSolnPhaseTypeCS(i) == 'SUBQ').OR. &
             (cSolnPhaseTypeCS(i) == 'SUBI')) then
            nCountSublatticeCS = nCountSublatticeCS + 1
        end if

        if (l > 0) then
            ! This solution phase will be included in the system.
            nCounter = nCounter + 1
            if (nCounter > nSolnPhasesSys) exit LOOP_SolnPhases
            cSolnPhaseName(nCounter) = cSolnPhaseNameCS(i)
            cSolnPhaseType(nCounter) = cSolnPhaseTypeCS(i)
        else
            ! This solution phase will no longer be part of the system.
            cycle LOOP_SolnPhases
        end if

        select case (cSolnPhaseTypeCS(i))
            case ('IDMX')
                ! Ideal mixture, do nothing.
            case ('QKTO', 'RKMP', 'RKMPM')
                ! Note that this is just checking whether the parameter should be considered.  This format
                ! is consistent amoungst the above list of phase types.

                ! Proceed if there are any mixing parameters for this phase:
                IF_Param: if (nParamPhaseCS(i) /= nParamPhaseCS(i-1)) then
                    ! Loop through mixing parameters:
                    LOOP_Param: do j = nParamPhaseCS(i-1) + 1, nParamPhaseCS(i)
                        if (iRegularParamCS(j,1) == 2) then
                            ! Binary term
                            k = iRegularParamCS(j,2) + nSpeciesPhaseCS(i-1)
                            l = iRegularParamCS(j,3) + nSpeciesPhaseCS(i-1)
                            if ((iSpeciesPass(k) > 0).AND.(iSpeciesPass(l) > 0)) then
                                nParam = nParam + 1
                                iParamPassCS(j) = 1
                            end if
                        elseif (iRegularParamCS(j,1) == 3) then
                            ! Ternary term
                            k = iRegularParamCS(j,2) + nSpeciesPhaseCS(i-1)
                            l = iRegularParamCS(j,3) + nSpeciesPhaseCS(i-1)
                            m = iRegularParamCS(j,4) + nSpeciesPhaseCS(i-1)
                            if ((iSpeciesPass(k) > 0).AND.(iSpeciesPass(l) > 0).AND.(iSpeciesPass(m) > 0)) then
                                nParam = nParam + 1
                                iParamPassCS(j) = 1
                            end if
                        elseif (iRegularParamCS(j,1) == 4) then
                            ! Quaternary term
                            k = iRegularParamCS(j,2) + nSpeciesPhaseCS(i-1)
                            l = iRegularParamCS(j,3) + nSpeciesPhaseCS(i-1)
                            m = iRegularParamCS(j,4) + nSpeciesPhaseCS(i-1)
                            n = iRegularParamCS(j,4) + nSpeciesPhaseCS(i-1)
                            if ((iSpeciesPass(k) > 0).AND.(iSpeciesPass(l) > 0).AND.(iSpeciesPass(m) > 0).AND.&
                                (iSpeciesPass(n) > 0)) then
                                nParam = nParam + 1
                                iParamPassCS(j) = 1
                            end if
                        else
                            ! An unsupported number of mixing terms.  Report an error:
                            INFOThermo = 32
                            return
                        end if
                    end do LOOP_Param
                end if IF_Param

                nParamPhase(nCounter) = nParam

            case ('SUBL', 'SUBLM')

                ! Check if the constituents pass for a phase with a sublattice:
                nCountSublattice                 = nCountSublattice + 1
                iPhaseSublattice(nCounter)       = nCountSublattice

                nSublatticePhase(nCountSublattice)  = nSublatticePhaseCS(nCountSublatticeCS)
                j = SIZE(nConstituentSublattice,DIM=2)
                n = nSublatticePhase(nCountSublattice)
                dStoichSublattice(nCountSublattice,1:n) = dStoichSublatticeCS(nCountSublatticeCS,1:n)
                k = SIZE(iConstituentSublattice, DIM=3)
                ! iConstituentSublattice(nCountSublattice,1:n,1:k) = iConstituentSublatticeCS(nCountSublatticeCS,1:n,1:k)
                k = iPhaseSublatticeCS(i)

                ! Loop through species in phase to determine which constituents are stable:
                do j = nSpeciesPhaseCS(i-1) + 1, nSpeciesPhaseCS(i)

                    if (iSpeciesPass(j) > 0) then
                        ! This species has passed.
                        n = j - nSpeciesPhaseCS(i-1)

                        ! Loop through sublattices per phase:
                        do s = 1, nSublatticePhaseCS(k)
                            m = iConstituentSublatticeCS(iPhaseSublatticeCS(i), s, n)
                            iConstituentPass(k, s, m) = 1
                        end do
                    else
                        ! This phase did not pass.
                    end if
                end do

                ! Count the number of constituents on each sublattice:
                do s = 1, nSublatticePhaseCS(k)

                    ! Loop through constituents on sublattice s:
                    j = 0
                    do c = 1, nConstituentSublatticeCS(k,s)
                        nConstituentSublattice(nCountSublattice,s) = nConstituentSublattice(nCountSublattice,s) &
                            + iConstituentPass(k, s, c)

                        ! Store the correct constituent name:
                        if (iConstituentPass(k,s,c) /= 0) then
                            iConstituentPass(k,s,c) = nConstituentSublattice(nCountSublattice,s)
                            j = j + 1
                            cConstituentNameSUB(nCountSublattice,s,j) = cConstituentNameSUBCS(nCountSublatticeCS,s,c)
                            ! iConstituentSublattice(nCountSublattice,s,j) = iConstituentSublatticeCS(nCountSublatticeCS,s,c)
                        end if
                    end do
                end do

                ! Get the constituents that passed on each sublattice
                j = 0
                LOOP_findPassed: do c = 1, SIZE(iConstituentSublatticeCS, DIM=3)
                    do s = 1, nSublatticePhaseCS(k)
                        if (iConstituentSublatticeCS(k,s,c) < 1) cycle LOOP_findPassed
                        if (iConstituentPass(k,s,iConstituentSublatticeCS(k,s,c)) == 0) cycle LOOP_findPassed
                    end do
                    ! If not cycled above, then all constituents passed
                    j = j + 1
                    do s = 1, nSublatticePhaseCS(k)
                        iConstituentSublattice(nCountSublattice,s,j) = iConstituentSublatticeCS(nCountSublatticeCS,s,c)
                    end do
                end do LOOP_findPassed

                ! Now reduce indices to account for removed constituents
                do s = 1, nSublatticePhaseCS(nCountSublatticeCS)
                    nRemove = 0
                    iRemove = 0
                    do c = nConstituentSublatticeCS(k,s), 1, -1
                        if (iConstituentPass(k,s,c) == 0) then
                            nRemove = nRemove + 1
                            iRemove(nRemove) = c
                        end if
                    end do
                    do c = 1, nRemove
                        do l = SIZE(iConstituentSublattice,3), 1, -1
                            if (iConstituentSublattice(nCountSublattice,s,l) > iRemove(c)) then
                                iConstituentSublattice(nCountSublattice,s,l) = iConstituentSublattice(nCountSublattice,s,l) - 1
                            end if
                        end do
                    end do
                end do

                ! Proceed if there are any mixing parameters for this phase:
                IF_Param_SUBL: if (nParamPhaseCS(i) /= nParamPhaseCS(i-1)) then

                    ! Loop through all mixing parameters for this phase:
                    LOOP_Param_SUBL: do j = nParamPhaseCS(i-1) + 1, nParamPhaseCS(i)

                        ! Loop through constituents associated with this parameter:
                        do k = 1, iRegularParamCS(j,1)

                            ! Store the constituent index to memory:
                            l = iRegularParamCS(j,1+k)

                            ! Loop through sublattices associated with this phase:
                            LOOP_C: do m = 1, nSublatticePhaseCS(iPhaseSublatticeCS(i))

                                ! Store the number of constituents for this sublattice:
                                n = nConstituentSublatticeCS(iPhaseSublatticeCS(i),m)

                                if (l <= n) then
                                    ! l is the constituent index on sublattice m.

                                    if (iConstituentPass(iPhaseSublatticeCS(i),m,l) == 0) then

                                        ! If any of the constituents associated with this parameter did not pass, then
                                        ! the mixing parameter will not be used.
                                        cycle LOOP_Param_SUBL
                                    else
                                        exit LOOP_C
                                    end if
                                else
                                    l = l - n
                                    cycle LOOP_C
                                end if
                            end do LOOP_C

                        end do

                        ! The parameter will be considered in the system.
                        nParam          = nParam + 1
                        iParamPassCS(j) = 1

                    end do LOOP_Param_SUBL

                end if IF_Param_SUBL

                nParamPhase(nCounter) = nParam

            case ('SUBI')

                ! Check if the constituents pass for a phase with a sublattice:
                nCountSublattice                 = nCountSublattice + 1
                iPhaseSublattice(nCounter)       = nCountSublattice

                nSublatticePhase(nCountSublattice)  = nSublatticePhaseCS(nCountSublatticeCS)
                j = SIZE(nConstituentSublattice,DIM=2)
                n = nSublatticePhase(nCountSublattice)
                k = SIZE(iConstituentSublattice, DIM=3)
                ! iConstituentSublattice(nCountSublattice,1:n,1:k) = iConstituentSublatticeCS(nCountSublatticeCS,1:n,1:k)
                k = iPhaseSublatticeCS(i)

                ! Loop through species in phase to determine which constituents are stable:
                ! Identical to SUBL doesn't work because a passed species might have 0 stoichiometry of a constituent
                do j = nSpeciesPhaseCS(i-1) + 1, nSpeciesPhaseCS(i)

                    if (iSpeciesPass(j) > 0) then
                        ! This species has passed.
                        n = j - nSpeciesPhaseCS(i-1)

                        ! On first sublattice check if second sublattice has neutral, don't pass if so
                        if (dSublatticeChargeCS(nCountSublatticeCS,2,n) /= 0D0) then
                            m = iConstituentSublatticeCS(k, 1, n)
                            iConstituentPass(k, 1, m) = 1
                        end if
                        m = iConstituentSublatticeCS(k, 2, n)
                        iConstituentPass(k, 2, m) = 1
                    else
                        ! This phase did not pass.
                    end if
                end do

                ! Count the number of constituents on each sublattice:
                do s = 1, nSublatticePhaseCS(k)

                    ! Loop through constituents on sublattice s:
                    j = 0
                    do c = 1, nConstituentSublatticeCS(k,s)
                        nConstituentSublattice(nCountSublattice,s) = nConstituentSublattice(nCountSublattice,s) &
                            + iConstituentPass(k, s, c)

                        ! Store the correct constituent name:
                        if (iConstituentPass(k,s,c) /= 0) then
                            iConstituentPass(k,s,c) = nConstituentSublattice(nCountSublattice,s)
                            j = j + 1
                            cConstituentNameSUB(nCountSublattice,s,j) = cConstituentNameSUBCS(nCountSublatticeCS,s,c)
                            dSublatticeCharge(nCountSublattice,s,j) = dSublatticeChargeCS(nCountSublatticeCS,s,c)
                            ! iConstituentSublattice(nCountSublattice,s,j) = iConstituentSublatticeCS(nCountSublatticeCS,s,c)
                        end if
                    end do
                end do

                ! Get the constituents that passed on each sublattice
                j = 0
                LOOP_findPassed2: do c = 1, SIZE(iConstituentSublatticeCS, DIM=3)
                    do s = 1, nSublatticePhaseCS(k)
                        if (iConstituentSublatticeCS(k,s,c) < 1) cycle LOOP_findPassed2
                    end do
                    if ((dSublatticeChargeCS(nCountSublatticeCS,2,iConstituentSublatticeCS(k,2,c)) /= 0D0) .AND. &
                       (iConstituentPass(k,1,iConstituentSublatticeCS(k,1,c)) == 0)) cycle LOOP_findPassed2
                    if (iConstituentPass(k,2,iConstituentSublatticeCS(k,2,c)) == 0)  cycle LOOP_findPassed2
                    ! If not cycled above, then all constituents passed
                    j = j + 1
                    do s = 1, nSublatticePhaseCS(k)
                        iConstituentSublattice(nCountSublattice,s,j) = iConstituentSublatticeCS(nCountSublatticeCS,s,c)
                    end do
                end do LOOP_findPassed2

                ! Now reduce indices to account for removed constituents
                do s = 1, nSublatticePhaseCS(nCountSublatticeCS)
                    nRemove = 0
                    iRemove = 0
                    do c = nConstituentSublatticeCS(k,s), 1, -1
                        if (iConstituentPass(k,s,c) == 0) then
                            nRemove = nRemove + 1
                            iRemove(nRemove) = c
                        end if
                    end do
                    do c = 1, nRemove
                        do l = SIZE(iConstituentSublattice,3), 1, -1
                            if (iConstituentSublattice(nCountSublattice,s,l) == iRemove(c)) then
                                iConstituentSublattice(nCountSublattice,s,l) = -1
                            else if (iConstituentSublattice(nCountSublattice,s,l) > iRemove(c)) then
                                iConstituentSublattice(nCountSublattice,s,l) = iConstituentSublattice(nCountSublattice,s,l) - 1
                            end if
                        end do
                    end do
                end do

                ! Proceed if there are any mixing parameters for this phase:
                if (nParamPhaseCS(i) /= nParamPhaseCS(i-1)) then

                    ! Loop through all mixing parameters for this phase:
                    LOOP_Param_SUBI: do j = nParamPhaseCS(i-1) + 1, nParamPhaseCS(i)

                        ! Loop through constituents associated with this parameter:
                        do k = 1, iRegularParamCS(j,1)

                            ! Store the constituent index to memory:
                            l = iRegularParamCS(j,1+k)

                            ! Loop through sublattices associated with this phase:
                            LOOP_C2: do m = 1, nSublatticePhaseCS(iPhaseSublatticeCS(i))

                                ! Store the number of constituents for this sublattice:
                                n = nConstituentSublatticeCS(iPhaseSublatticeCS(i),m)

                                if (l <= n) then
                                    ! l is the constituent index on sublattice m.

                                    if (iConstituentPass(iPhaseSublatticeCS(i),m,l) == 0) then

                                        ! If any of the constituents associated with this parameter did not pass, then
                                        ! the mixing parameter will not be used.
                                        cycle LOOP_Param_SUBI
                                    else
                                        exit LOOP_C2
                                    end if
                                else
                                    l = l - n
                                    cycle LOOP_C2
                                end if
                            end do LOOP_C2

                        end do

                        ! The parameter will be considered in the system.
                        nParam          = nParam + 1
                        iParamPassCS(j) = 1

                    end do LOOP_Param_SUBI

                end if

                nParamPhase(nCounter) = nParam

            case ('SUBG', 'SUBQ')
                ! Check if the constituents pass for a phase with a sublattice:
                nCountSublattice                 = nCountSublattice + 1
                iPhaseSublattice(nCounter)       = nCountSublattice

                nSublatticePhase(nCountSublattice)  = nSublatticePhaseCS(nCountSublatticeCS)
                do j = 1, nSublatticePhase(nCountSublattice)
                    m = 0
                    do k = 1, nSublatticeElementsCS(nCountSublatticeCS,j)
                        if (iConstituentPass(nCountSublatticeCS,j,k) /= 0) then
                            m = m + 1
                            dSublatticeCharge(nCountSublattice,j,m) = dSublatticeChargeCS(nCountSublatticeCS,j,k)
                            iChemicalGroup(nCountSublattice,j,m) = iChemicalGroupCS(nCountSublatticeCS,j,k)
                        end if
                    end do
                    nSublatticeElements(nCountSublattice,j) = m
                end do

                m = 0
                LOOP_iConstitSubl: do k = 1, SIZE(iConstituentSublatticeCS, DIM=3)
                    iCon1 = iConstituentSublatticeCS(nCountSublatticeCS,1,k)
                    iCon2 = iConstituentSublatticeCS(nCountSublatticeCS,2,k)
                    if (iCon1 == 0) cycle LOOP_iConstitSubl
                    if (iCon2 == 0) cycle LOOP_iConstitSubl
                    if (iConstituentPass(nCountSublatticeCS,1,iCon1) == 0) cycle LOOP_iConstitSubl
                    if (iConstituentPass(nCountSublatticeCS,2,iCon2) == 0) cycle LOOP_iConstitSubl
                    m = m + 1
                    iConstituentSublattice(nCountSublattice,1,m) = iCon1
                    iConstituentSublattice(nCountSublattice,2,m) = iCon2
                    cConstituentNameSUB(nCountSublattice,1,m) = cConstituentNameSUBCS(nCountSublatticeCS,1,iCon1)
                    cConstituentNameSUB(nCountSublattice,2,m) = cConstituentNameSUBCS(nCountSublatticeCS,2,iCon2)
                end do LOOP_iConstitSubl

                ! Save quadruplet data corresponding to the quadruplets remaining in the system
                k = SIZE(iPairID,DIM = 2)
                do k = 1, nPairsSROCS(nCountSublatticeCS,2)
                    ! Find indices of constituents in quadruplet
                    pa = iPairIDCS(nCountSublattice,k,1)
                    pb = iPairIDCS(nCountSublattice,k,2)
                    px = iPairIDCS(nCountSublattice,k,3) - nSublatticeElementsCS(nCountSublatticeCS,1)
                    py = iPairIDCS(nCountSublattice,k,4) - nSublatticeElementsCS(nCountSublatticeCS,1)
                    ! Find index of quadruplet
                    if (px == py) then
                        p = (px - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) &
                                     * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
                    else
                        p = (nSublatticeElementsCS(nCountSublatticeCS,2) + (px - 1) + ((py-2)*(py-1)/2)) &
                          * (nSublatticeElementsCS(nCountSublatticeCS,1) * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
                    end if
                    if (pa == pb) then
                        l = pa
                    else
                        l = nSublatticeElementsCS(nCountSublatticeCS,1) + pa + ((pb-2)*(pb-1)/2)
                    end if
                    iIndex = p + l + nSpeciesPhaseCS(i - 1)
                    ! Make sure this quadruplet is still part of the system, if so save data
                    if (iSpeciesPass(iIndex) > 0) then
                        nPairsSRO(nCountSublattice,2) = nPairsSRO(nCountSublattice,2) + 1
                        n = nPairsSRO(nCountSublattice,2)
                        iPairID(nCountSublattice,n,1:4) = iPairIDCS(nCountSublatticeCS,k,1:4)
                        dCoordinationNumber(nCountSublattice,n,1:4) = dCoordinationNumberCS(nCountSublatticeCS,k,1:4)
                    end if
                end do

                ! Must remove unused elements from iPairID
                nRemove = 0
                iRemove = 0
                do j = nSublatticePhaseCS(nCountSublatticeCS), 1, -1
                    do l = nSublatticeElementsCS(nCountSublatticeCS,j), 1, -1
                        if (iConstituentPass(nCountSublatticeCS,j,l) == 0) then
                            nRemove = nRemove + 1
                            iRemove(nRemove) = l + ((j - 1) * nSublatticeElementsCS(nCountSublatticeCS,1))
                        end if
                    end do
                end do

                do k = 1, nRemove
                    do j = 1, nPairsSRO(nCountSublattice,2)
                        do l = 1, 4
                            if (iPairID(nCountSublattice,j,l) > iRemove(k)) then
                                iPairID(nCountSublattice,j,l) = iPairID(nCountSublattice,j,l) - 1
                            end if
                        end do
                    end do
                end do

                ! Also remove from iConstituentSublattice
                do j = nSublatticePhaseCS(nCountSublatticeCS), 1, -1
                    nRemove = 0
                    iRemove = 0
                    do l = nSublatticeElementsCS(nCountSublatticeCS,j), 1, -1
                        if (iConstituentPass(nCountSublatticeCS,j,l) == 0) then
                            nRemove = nRemove + 1
                            iRemove(nRemove) = l
                        end if
                    end do
                    do k = 1, nRemove
                        do l = SIZE(iConstituentSublattice,3), 1, -1
                            if (iConstituentSublattice(nCountSublattice,j,l) > iRemove(k)) then
                                iConstituentSublattice(nCountSublattice,j,l) = iConstituentSublattice(nCountSublattice,j,l) - 1
                            end if
                        end do
                    end do
                end do

                ! Loop through excess parameters (same checks as for quadruplets):
                ! Except for ternary terms, must also check the additional species
                LOOP_excess: do j = nParamPhaseCS(i-1) + 1, nParamPhaseCS(i)
                    pa = iRegularParamCS(j,2)
                    pb = iRegularParamCS(j,3)
                    px = iRegularParamCS(j,4) - nSublatticeElementsCS(nCountSublatticeCS,1)
                    py = iRegularParamCS(j,5) - nSublatticeElementsCS(nCountSublatticeCS,1)
                    pe = iRegularParamCS(j,10)
                    pf = iRegularParamCS(j,11)

                    ! If this is a ternary and the 3rd element got removed, then skip parameter
                    ! lRemoved = .TRUE.
                    if (pe > 0) then
                        if (iConstituentPass(nCountSublatticeCS,1,pe) == 0) cycle LOOP_excess
                    end if
                    ! If this is a ternary and the 3rd element got removed, then skip parameter
                    ! Making a guess at how to handle this for the other entry
                    ! lRemoved = .TRUE.
                    if (pf > 0) then
                        pf = pf - nSublatticeElementsCS(nCountSublatticeCS,1)
                        if (iConstituentPass(nCountSublatticeCS,2,pf) == 0) cycle LOOP_excess
                    end if

                    if (px == py) then
                        p = (px - 1) * (nSublatticeElementsCS(nCountSublatticeCS,1) &
                                     * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
                    else
                        p = (nSublatticeElementsCS(nCountSublatticeCS,2) + (px - 1) + ((py-2)*(py-1)/2)) &
                          * (nSublatticeElementsCS(nCountSublatticeCS,1) * (nSublatticeElementsCS(nCountSublatticeCS,1) + 1) / 2)
                    end if
                    if (pa == pb) then
                        l = pa
                    else
                        l = nSublatticeElementsCS(nCountSublatticeCS,1) + pa + ((pb-2)*(pb-1)/2)
                    end if
                    iIndex = p + l + nSpeciesPhaseCS(i - 1)
                    if (iSpeciesPass(iIndex) > 0) then
                        nParam          = nParam + 1
                        iParamPassCS(j) = 1
                    end if
                end do LOOP_excess

                nParamPhase(nCounter) = nParam

            case default
                ! The character string representing input units is not recognized.
                INFOThermo = 17
                return
        end select

        ! Magnetic mixing terms:
        if (cSolnPhaseTypeCS(i) == 'RKMPM') then
            ! Proceed if there are any mixing parameters for this phase:
            if (nMagParamPhaseCS(i) /= nMagParamPhaseCS(i-1)) then
                ! Loop through mixing parameters:
                do j = nMagParamPhaseCS(i-1) + 1, nMagParamPhaseCS(i)
                    if (iMagneticParamCS(j,1) == 2) then
                        ! Binary term
                        k = iMagneticParamCS(j,2) + nSpeciesPhaseCS(i-1)
                        l = iMagneticParamCS(j,3) + nSpeciesPhaseCS(i-1)
                        if ((iSpeciesPass(k) > 0).AND.(iSpeciesPass(l) > 0)) then
                            nMagParam = nMagParam + 1
                            iMagParamPassCS(j) = 1
                        end if
                    elseif (iMagneticParamCS(j,1) == 3) then
                        ! Ternary term
                        k = iMagneticParamCS(j,2) + nSpeciesPhaseCS(i-1)
                        l = iMagneticParamCS(j,3) + nSpeciesPhaseCS(i-1)
                        m = iMagneticParamCS(j,4) + nSpeciesPhaseCS(i-1)
                        if ((iSpeciesPass(k) > 0).AND.(iSpeciesPass(l) > 0).AND.(iSpeciesPass(m) > 0)) then
                            nMagParam = nMagParam + 1
                            iMagParamPassCS(j) = 1
                        end if
                    elseif (iMagneticParamCS(j,1) == 4) then
                        ! Quaternary term
                        k = iMagneticParamCS(j,2) + nSpeciesPhaseCS(i-1)
                        l = iMagneticParamCS(j,3) + nSpeciesPhaseCS(i-1)
                        m = iMagneticParamCS(j,4) + nSpeciesPhaseCS(i-1)
                        n = iMagneticParamCS(j,4) + nSpeciesPhaseCS(i-1)
                        if ((iSpeciesPass(k) > 0).AND.(iSpeciesPass(l) > 0).AND.(iSpeciesPass(m) > 0).AND.&
                            (iSpeciesPass(n) > 0)) then
                            nMagParam = nMagParam + 1
                            iMagParamPassCS(j) = 1
                        end if
                    else
                        ! An unsupported number of mixing terms.  Report an error:
                        INFOThermo = 43
                        return
                    end if
                end do
            end if
        else if (cSolnPhaseTypeCS(i) == 'SUBLM') then
            ! Proceed if there are any mixing parameters for this phase:
            if (nMagParamPhaseCS(i) /= nMagParamPhaseCS(i-1)) then

                ! Loop through all mixing parameters for this phase:
                LOOP_MagParam_SUBL: do j = nMagParamPhaseCS(i-1) + 1, nMagParamPhaseCS(i)

                    ! Loop through constituents associated with this parameter:
                    do k = 1, iMagneticParamCS(j,1)

                        ! Store the constituent index to memory:
                        l = iMagneticParamCS(j,1+k)

                        ! Loop through sublattices associated with this phase:
                        LOOP_Cmag: do m = 1, nSublatticePhaseCS(iPhaseSublatticeCS(i))

                            ! Store the number of constituents for this sublattice:
                            n = nConstituentSublatticeCS(iPhaseSublatticeCS(i),m)

                            if (l <= n) then
                                ! l is the constituent index on sublattice m.

                                if (iConstituentPass(iPhaseSublatticeCS(i),m,l) == 0) then

                                    ! If any of the constituents associated with this parameter did not pass, then
                                    ! the mixing parameter will not be used.
                                    cycle LOOP_MagParam_SUBL
                                else
                                    exit LOOP_Cmag
                                end if
                            else
                                l = l - n
                                cycle LOOP_Cmag
                            end if
                        end do LOOP_Cmag

                    end do

                    ! The parameter will be considered in the system.
                    nMagParam          = nMagParam + 1
                    iMagParamPassCS(j) = 1

                end do LOOP_MagParam_SUBL
            end if
        end if
        nMagParamPhase(nCounter) = nMagParam

    end do LOOP_SolnPhases

    ! Check to see if the mixing terms need to be reallocated:
    if (allocated(dExcessGibbsParam)) then
        i = SIZE(dExcessGibbsParam)
        if (i /= nParam) then
            deallocate(iRegularParam,dExcessGibbsParam,cRegularParam,iSUBLParamData, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            allocate(iRegularParam(nParam,nParamMax*2+3),dExcessGibbsParam(nParam),cRegularParam(nParam),iSUBLParamData(nParam,7))
        end if
    else
        ! Allocate memory for excess parameters:
        allocate(iRegularParam(nParam,nParamMax*2+3),dExcessGibbsParam(nParam),cRegularParam(nParam),iSUBLParamData(nParam,7))
    end if

    ! Same as above, for magnetic mixing:
    if (allocated(dMagneticParam)) then
        i = SIZE(dMagneticParam)
        if (i /= nMagParam) then
            deallocate(iMagneticParam,dMagneticParam, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            allocate(iMagneticParam(nMagParam,nParamMax*2+3),dMagneticParam(nMagParam,2))
        end if
    else
        ! Allocate memory for excess parameters:
        allocate(iMagneticParam(nMagParam,nParamMax*2+3),dMagneticParam(nMagParam,2))
    end if

    ! Determine whether a solution phase is miscibile.  This flag will be used by the main solver.
    do i = 2, nSolnPhasesSys
        if (cSolnPhaseName(i) == cSolnPhaseName(i-1)) then
            lMiscibility(i)   = .TRUE.
        end if
    end do

    ! Initialize variables:
    dExcessGibbsParam = 0D0
    iRegularParam     = 0

    return

end subroutine CheckSystemExcess
