
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

    integer::  c, i, j, k, l, m, n, s, nCounter, pa, pb, px, py
    integer, dimension(4) :: iPairMatch


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
             (cSolnPhaseTypeCS(i) == 'SUBG').OR.(cSolnPhaseTypeCS(i) == 'SUBQ')) then
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
                iConstituentSublattice(nCountSublattice,1:n,1:k) = iConstituentSublatticeCS(nCountSublatticeCS,1:n,1:k)
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
                        end if
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

            case ('SUBG', 'SUBQ')
                ! Check if the constituents pass for a phase with a sublattice:
                nCountSublattice                 = nCountSublattice + 1
                iPhaseSublattice(nCounter)       = nCountSublattice

                nSublatticePhase(nCountSublattice)  = nSublatticePhaseCS(nCountSublatticeCS)
                do j = 1, nSublatticePhase(nCountSublattice)
                    m = 0
                    do k = 1, nElemOrComp
                        if ((iSublatticeElementsCS(nCountSublatticeCS,j,k) > 0)) then
                            if (iElementSystem(iSublatticeElementsCS(nCountSublatticeCS,j,k)) > 0) then
                                m = m + 1
                                iSublatticeElements(nCountSublattice,j,m) = iSublatticeElementsCS(nCountSublatticeCS,j,k)
                                dSublatticeCharge(nCountSublattice,j,m) = dSublatticeChargeCS(nCountSublatticeCS,j,k)
                            end if
                        end if
                    end do
                    nSublatticeElements(nCountSublattice,j) = m
                end do

                do j = 1, nSublatticePhase(nCountSublattice)
                    m = 0
                    LOOP_iConstitSubl: do k = 1, SIZE(iConstituentSublatticeCS, DIM=3)
                        if (iConstituentSublatticeCS(nCountSublatticeCS,j,k) == 0) cycle LOOP_iConstitSubl
                        n = iSublatticeElements(nCountSublattice,j,iConstituentSublatticeCS(nCountSublatticeCS,j,k))
                        if (iElementSystem(n) > 0) then
                            m = m + 1
                            iConstituentSublattice(nCountSublattice,j,m) = iConstituentSublatticeCS(nCountSublatticeCS,j,k)
                        end if
                    end do LOOP_iConstitSubl
                end do

                k = SIZE(iPairID,DIM = 2)
                do k = 1, nPairsSROCS(nCountSublatticeCS,2)
                    iPairMatch = 0
                    pa = iSublatticeElementsCS(nCountSublatticeCS,1,iPairIDCS(nCountSublattice,k,1))
                    pb = iSublatticeElementsCS(nCountSublatticeCS,1,iPairIDCS(nCountSublattice,k,2))
                    px = iSublatticeElementsCS(nCountSublatticeCS,2,iPairIDCS(nCountSublattice,k,3) &
                       - nSublatticeElementsCS(nCountSublatticeCS,1))
                    py = iSublatticeElementsCS(nCountSublatticeCS,2,iPairIDCS(nCountSublattice,k,4) &
                       - nSublatticeElementsCS(nCountSublatticeCS,1))
                    do l = 1, nSublatticeElements(nCountSublattice,1)
                        if (pa == iSublatticeElements(nCountSublattice,1,l)) iPairMatch(1) = 1
                        if (pb == iSublatticeElements(nCountSublattice,1,l)) iPairMatch(2) = 1
                    end do
                    do l = 1, nSublatticeElements(nCountSublattice,2)
                        if (px == iSublatticeElements(nCountSublattice,2,l)) iPairMatch(3) = 1
                        if (py == iSublatticeElements(nCountSublattice,2,l)) iPairMatch(4) = 1
                    end do
                    if (SUM(iPairMatch) == 4) then
                        nPairsSRO(nCountSublattice,2) = nPairsSRO(nCountSublattice,2) + 1
                        n = nPairsSRO(nCountSublattice,2)
                        iPairID(nCountSublattice,n,1:4) = iPairIDCS(nCountSublatticeCS,k,1:4)
                        dCoordinationNumber(nCountSublattice,n,1:4) = dCoordinationNumberCS(nCountSublatticeCS,k,1:4)
                        dZetaSpecies(nCountSublattice,n) = dZetaSpeciesCS(nCountSublatticeCS,k)
                    end if
                end do

                ! Loop through excess parameters:
                do j = nParamPhaseCS(i-1) + 1, nParamPhaseCS(i)
                    iPairMatch = 0
                    pa = iRegularParamCS(j,2)
                    pa = iSublatticeElementsCS(nCountSublatticeCS,1,pa)
                    pb = iRegularParamCS(j,3)
                    pb = iSublatticeElementsCS(nCountSublatticeCS,1,pb)
                    px = iRegularParamCS(j,4) - nSublatticeElementsCS(nCountSublatticeCS,1)
                    px = iSublatticeElementsCS(nCountSublatticeCS,2,px)
                    py = iRegularParamCS(j,5) - nSublatticeElementsCS(nCountSublatticeCS,1)
                    py = iSublatticeElementsCS(nCountSublatticeCS,2,py)
                    do l = 1, nSublatticeElements(nCountSublattice,1)
                        if (pa == iSublatticeElements(nCountSublattice,1,l)) iPairMatch(1) = 1
                        if (pb == iSublatticeElements(nCountSublattice,1,l)) iPairMatch(2) = 1
                    end do
                    do l = 1, nSublatticeElements(nCountSublattice,2)
                        if (px == iSublatticeElements(nCountSublattice,2,l)) iPairMatch(3) = 1
                        if (py == iSublatticeElements(nCountSublattice,2,l)) iPairMatch(4) = 1
                    end do
                    if (SUM(iPairMatch) == 4) then
                        nParam          = nParam + 1
                        iParamPassCS(j) = 1
                    end if
                end do

                nParamPhase(nCounter) = nParam

            case default
                ! The character string representing input units is not recognized.
                INFOThermo = 17
                return
        end select

    end do LOOP_SolnPhases

    ! Check to see if the mixing terms need to be reallocated:
    if (allocated(dExcessGibbsParam)) then
        i = SIZE(dExcessGibbsParam)
        if (i /= nParam) then
            deallocate(iRegularParam,dExcessGibbsParam, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            allocate(iRegularParam(nParam,nParamMax*2+1),dExcessGibbsParam(nParam))
        end if
    else
        ! Allocate memory for excess parameters:
        allocate(iRegularParam(nParam,nParamMax*2+1),dExcessGibbsParam(nParam))
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
