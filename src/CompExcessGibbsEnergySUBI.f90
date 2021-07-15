
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergySUBI.f90
    !> \brief   Compute the partial molar excess Gibbs energy of mixing of solution phase constituents in a
    !!          SUBI solution phase.
    !> \author  M.H.A. Piro
    !> \date    January 19, 2021
    !> \sa      CompExcessGibbsEnergy.f90
    !> \sa      CompExcessGibbsEnergyRKMP.f90
    !> \sa      CompExcessGibbsEnergyQKTO.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   01/19/2021      M. Poschmann     Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the partial molar excess Gibbs energy of mixing
    !! (dPartialExcessGibbs) of all constituents in a non-ideal solution phase designated as 'SUBI'
    !! (Ionic Liquid Model).
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] iSolnIndex    Absolute index of a solution phase
    !
    ! nSpeciesPhase             An integec vector representing the last index of a species in a particular
    !                            solution phase
    ! nParamPhase               An integer vector representing the number of parameters for a particular phase.
    ! nSublattice               An integer scalar representing the number of sublattices associated with this
    !                            particular phase (used for convenience).
    ! nSublatticePhase          An integer vector representing the number of sublattices for each charged phase.
    ! iSPI           An integer scalar representing the relative index of the charged phase.
    ! iParam                    An integer scalar representing the index number of a parameter.
    ! iPhaseSublattice          An integer vector representing the relative index of a charged phase corresponding
    !                            to each solution phase in the system.  Thus, it is equal to zero for a phase
    !                            that does not contain any sublattices.
    ! iRegularParam             An integer matrix representing the relative indices of species mixing for a
    !                            particular parmater.
    ! iFirst             An integer scalar representing the absolute index of the first species in a phase.
    ! iLast              An integer scalar representing the absolute index of the last species in a phase.
    ! dChemicalPotential        A double real vector representing the chemical potential for every species in the
    !                            system.
    ! dMolFraction              A double real vector representing hte mole fraction of all species in the system.
    ! dSiteFraction             A double real array representing the site fraction of each constituent on each
    !                            sublattice.  The first dimension corresponds to the charged phase index, the
    !                            second dimension corresponds to the sublattice index and the third dimension
    !                            corresponds to the constituent index.
    ! dExcessGibbsParam         A double real vector representing the molar excess Gibbs energy of mixing for
    !                            each subsystem.
    ! dPartialExcessGibbs       Partial molar excess Gibbs energy of mixing of species.
    ! dMolFraction              Current estimated mole fraction.
    ! cSolnPhaseType            A character vector representing the solution phase type.
    ! iConstituentSublattice    An integer array representing the constituent index for a particular sublattice
    !                            phase.  The first dimension refers to the relative phase index for charged phases,
    !                            the second dimension represents the sublattice index and the third dimension
    !                            corresponds to the relative component index (not absolute).
    ! iFirstParam               An integer scalar representing the constituent index of the first constituent
    !                            that is being mixed for a mixing parameter.
    ! iSeconParam               An integer scalar representing the constituent index of the second constituent
    !                            that is being mixed for a mixing parameter.
    ! iSubParam                 An integer scalar representing the sublattice index corresponding to the mixing
    !                            parameter.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompExcessGibbsEnergySUBI(iSolnIndex)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, l, k1, l1, k2, l2, m, n, c, d, k, s
    integer :: iCi, iCj, iCk, iBi, iBj, iBk, iAi, iAj, iDi, iDj
    integer :: iSolnIndex, nSublattice, iSPI, iExponent
    integer :: iFirst, iLast
    real(8) :: dSub1Total, dSub2Total, dydn, dTemp1, dTemp2
    real(8) :: dSum, p, q, kc1, kc2, lc1, lc2, cc1, gref, gideal, gexcess, natom, yva, dMol, dMolAtoms, dTest
    real(8), dimension(:), allocatable :: dgdc1, dgdc2, dMolDerivatives
    real(8) :: dPreFactor, v, f, chargeCi, chargeCj, chargeCk
    real(8) :: yCi, yCj, yCk, yAi, yAj, yBi, yBj, yBk, yDi, yDj, gex

 ! print*,""
 ! print*,""
 ! print*,"                          CompExcessGibbsEnergySUBI.f90"
 ! print*,""
 ! print*,""

    ! Only proceed if the correct phase type is selected:
    IF_SUBL: if (cSolnPhaseType(iSolnIndex) == 'SUBI') then

        ! Define temporary variables for sake of convenience:
        iSPI = iPhaseSublattice(iSolnIndex)
        nSublattice     = nSublatticePhase(iSPI)
        iFirst          = nSpeciesPhase(iSolnIndex-1) + 1
        iLast           = nSpeciesPhase(iSolnIndex)

        if (allocated(dgdc1)) deallocate(dgdc1)
        if (allocated(dgdc2)) deallocate(dgdc2)
        if (allocated(dMolDerivatives)) deallocate(dMolDerivatives)
        allocate(dgdc1(nConstituentSublattice(iSPI,1)),dgdc2(nConstituentSublattice(iSPI,2)))
        allocate(dMolDerivatives(iLast - iFirst + 1))

        ! Initialize variables:
        dSiteFraction(iSPI,1:nSublattice,1:nMaxConstituentSys) = 0D0
        dChemicalPotential(iFirst:iLast)                       = 0D0
        dPartialExcessGibbs(iFirst:iLast)                      = 0D0
        dgdc1                                                  = 0D0
        dgdc2                                                  = 0D0
        dMolDerivatives                                        = 0D0

        ! Compute site fractions on first sublattice:
        do i = iFirst, iLast
            ! Relative component index:
            m = i - iFirst + 1
            c = iConstituentSublattice(iSPI,1,m)
            d = iConstituentSublattice(iSPI,2,m)

            if (cConstituentNameSUB(iSPI,2,d) == 'Va') then
                lc2 = 1D0
            else
                lc2 = -dSublatticeCharge(iSPI,2,d)
            end if
            if (c > 0) then
                dSiteFraction(iSPI,1,c) = dSiteFraction(iSPI,1,c) + dMolFraction(i) * lc2
            end if
        end do

        ! Correct first site fraction
        dSub1Total = 0D0
        do i = 1, nConstituentSublattice(iSPI,1)
            dSub1Total = dSub1Total + dSiteFraction(iSPI,1,i)
        end do
        ! print*,"dSub1Total",dSub1Total
        do i = 1, nConstituentSublattice(iSPI,1)
            dSiteFraction(iSPI,1,i) = dSiteFraction(iSPI,1,i) / dSub1Total
        end do

        ! print*,"cConstituentNameSUB(iSPI,1,i):   ",cConstituentNameSUB(iSPI,1,:)
        ! print*,"dSiteFraction(iSPI,1,:):   ",dSiteFraction(iSPI,1,:)
        !dSiteFraction(iSPI,1,1) = 0.5013D0
        !dSiteFraction(iSPI,1,2) = 0.4987D0

        ! Compute Q
        q = 0D0
        do i = 1, nConstituentSublattice(iSPI,1)
            q = q + dSublatticeCharge(iSPI,1,i) * dSiteFraction(iSPI,1,i)
        end do

        ! Compute site fractions on second sublattice:
        do i = iFirst, iLast
            ! Relative component index:
            m = i - iFirst + 1
            c = iConstituentSublattice(iSPI,2,m)
            if (cConstituentNameSUB(iSPI,2,c) == 'Va') then
                ! Vacancy gets scaled by Q
                d = iConstituentSublattice(iSPI,1,m)
                dSiteFraction(iSPI,2,c) = dSiteFraction(iSPI,2,c) + dMolFraction(i) * dSublatticeCharge(iSPI,1,d) / q
            else if (dSublatticeCharge(iSPI,2,c) == 0D0) then
                ! Neutral counts as 1
                dSiteFraction(iSPI,2,c) = dSiteFraction(iSPI,2,c) + dMolFraction(i)
            else
                d = iConstituentSublattice(iSPI,1,m)
                dSiteFraction(iSPI,2,c) = dSiteFraction(iSPI,2,c) + dMolFraction(i) * dSublatticeCharge(iSPI,1,d)
            end if
        end do

        ! Correct second site fraction
        dSub2Total = 0D0
        yva = 0D0
        do i = 1, nConstituentSublattice(iSPI,2)
            dSub2Total = dSub2Total + dSiteFraction(iSPI,2,i)
        end do
        ! print*,"dSub2Total",dSub2Total
        do i = 1, nConstituentSublattice(iSPI,2)
            dSiteFraction(iSPI,2,i) = dSiteFraction(iSPI,2,i) / dSub2Total
            ! find site fraction of vacancies
            if (cConstituentNameSUB(iSPI,2,i) == 'Va') yva = dSiteFraction(iSPI,2,i)
        end do

        ! print*,"cConstituentNameSUB(iSPI,2,i):   ",cConstituentNameSUB(iSPI,2,:)
        ! print*,"dSiteFraction(iSPI,2,:):   ",dSiteFraction(iSPI,2,:)
        !dSiteFraction(iSPI,2,1) = 0.40820D0
        !dSiteFraction(iSPI,2,2) = 0.58966D0
        !dSiteFraction(iSPI,2,3) = 4.8623D-12
        !dSiteFraction(iSPI,2,4) = 2.1448D-3

        ! Compute P
        p = 0D0
        do i = 1, nConstituentSublattice(iSPI,2)
            if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                ! Use Q as charge if this constituent is vacancy
                p = p + q * dSiteFraction(iSPI,2,i)
            else
                p = p - dSublatticeCharge(iSPI,2,i) * dSiteFraction(iSPI,2,i)
            end if
        end do

        ! Compute number of moles and its derivatives
        dMol = (p + (q * (1D0 - yva)))
        ! print*,"dMol",dMol
        do j = iFirst, iLast

            dTemp1 = 0D0
            dTemp2 = 0D0

            ! Relative species index:
            n = j - iFirst + 1

            ! Store constituent indices:
            l1 = iConstituentSublattice(iSPI,1,n)
            l2 = iConstituentSublattice(iSPI,2,n)

            ! Allocating the correct constituent charges
            lc1 = 1D0
            if (l1 > 0) lc1 = dSublatticeCharge(iSPI,1,l1)
            if ((cConstituentNameSUB(iSPI,2,l2) == 'Va') .OR. &
                (dSublatticeCharge(iSPI,2,l2) == 0D0)) then
                lc2 = 1D0
            else
                lc2 = -dSublatticeCharge(iSPI,2,l2)
            end if

            if (cConstituentNameSUB(iSPI,2,l2) == 'Va') then
                ! cation / vacancy
                dMolDerivatives(n) = (q-lc1)*(1D0+(p-q*yva)*yva/q)/(dSub1Total*dMol**2)
                dMolDerivatives(n) = dMolDerivatives(n) + (p-q*yva)*lc1/(dSub2Total*q*dMol**2)
            else if (dSublatticeCharge(iSPI,2,l2) == 0D0) then
                ! neutral
                dMolDerivatives(n) = (p-q*yva)/(dSub2Total*dMol**2)
            else
                ! cation / anion
                dMolDerivatives(n) = (q-lc1)*(1D0+(p-q*yva)*yva/q)*lc2/(dSub1Total*dMol**2)
                dMolDerivatives(n) = dMolDerivatives(n) + (p-q*yva-lc2)*lc1/(dSub2Total*dMol**2)
            end if
        end do
        ! print*,""
        ! print*,"dMolDerivatives(:)",dMolDerivatives(:)
        ! print*,""

        ! Correct the mole fractions of phase components by the site fractions of the constituents
        dSum = 0D0
        LOOP_CorrectX: do i = iFirst, iLast
            dMolFraction(i) = 1D0
            m = i - iFirst + 1

            c = iConstituentSublattice(iSPI,2,m)
            d = iConstituentSublattice(iSPI,1,m)

            if ((cConstituentNameSUB(iSPI,2,c) == 'Va') .OR. &
                (dSublatticeCharge(iSPI,2,c) == 0D0)) then
                if (d > 0) dMolFraction(i) = dMolFraction(i) * dSiteFraction(iSPI,1,d)

                if (c > 0) dMolFraction(i) = dMolFraction(i) * dSiteFraction(iSPI,2,c)

            else
                if (d > 0) dMolFraction(i) = dMolFraction(i) * dSiteFraction(iSPI,1,d) / dSublatticeCharge(iSPI,1,d)

                if (c > 0) dMolFraction(i) = dMolFraction(i) * dSiteFraction(iSPI,2,c)

            end if
            dSum = dSum + dMolFraction(i)
            ! print*,"dSum",dSum
        end do LOOP_CorrectX

        ! Normalize mole fractions and compute number of mole atoms per mole (yes that makes sense, don't think about it)
        dMolAtoms = 0D0
        do i = iFirst, iLast

            dMolFraction(i)  = dMolFraction(i) / dSum

            m = i - iFirst + 1
            d = iConstituentSublattice(iSPI,1,m)
            c = iConstituentSublattice(iSPI,2,m)
            if ((cConstituentNameSUB(iSPI,2,c) == 'Va') .OR. &
                (dSublatticeCharge(iSPI,2,c) == 0D0)) then
                dMolAtoms = dMolAtoms + dMolFraction(i)
            else
                dMolAtoms = dMolAtoms + dMolFraction(i) * (dSublatticeCharge(iSPI,1,d) - dSublatticeCharge(iSPI,2,c))
            end if
        end do

        ! print*,""
        ! print*,"dMolFraction(:):   ",dMolFraction(:)
        ! print*,""

        dStoichSublattice(iSPI,1) = p
        dStoichSublattice(iSPI,2) = q

        !---------------------------------------------------------------
        !             EXCESS TERMS - With Derived Excess Terms
        ! --------------------------------------------------------------
        gexcess = 0D0
        ! Loop through parameters:
        LOOP_Param: do l = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)
            ! Return control to the parent subroutine if there aren't any interaction parameters for this phase:
            if (nParamPhase(iSolnIndex) - nParamPhase(iSolnIndex-1) == 0) exit LOOP_Param

            ! Reinitialize temporary variable:
            dPreFactor = 1D0

            f = 0D0
            v = 0D0

            chargeCi = 0D0
            chargeCj = 0D0
            chargeCk = 0D0

            iAi = 0
            iAj = 0
            iBi = 0
            iBj = 0
            iBk = 0
            iCi = 0
            iCj = 0
            iCk = 0
            iDi = 0
            iDj = 0

            yCi = 0D0
            yCj = 0D0
            yCk = 0D0
            yAi = 0D0
            yAj = 0D0
            yBi = 0D0
            yBj = 0D0
            yBk = 0D0
            yDi = 0D0
            yDj = 0D0

            gex = 0D0

            ! Store the number of constituents involved in this parameter:
            n = iRegularParam(l,1)

            ! Case 1: Determine the mixing parameter type: L_Ci,Cj:Ak
            if (iSUBIMixType(l) == 1) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    ! Store the first and second site fractions:
                    if (k == 2) then
                        ! cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                        chargeCi = dSublatticeCharge(iSPI,s,c)

                    else if (k == 3) then
                        ! cation - Cj
                        yCj = dSiteFraction(iSPI,s,c)
                        iCj = c
                        chargeCj = dSublatticeCharge(iSPI,s,c)

                    else if (k == 4) then
                        ! anion - Ak
                        yAi =  dSiteFraction(iSPI,s,c)
                        iAi = c

                    end if
                end do

                ! Multiply prefactor term by excess Gibbs energy parameter:
                iExponent = iRegularParam(l,n+2)
                ! Excess Gibbs energy equation for L_Ci,Cj:Ak case
                gex = dPreFactor * dExcessGibbsParam(l) * (yCi - yCj)**(iExponent)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                ! Avoiding a situation with 0^(-1)
                if ((yCi - yCj) /= 0) then

                    ! First sublattice - Cation:vacancy contributions
                    do i = 1, nConstituentSublattice(iSPI,1)

                        ! Derivative with respect to Ci
                        if (i == iCi) then
                            dgdc1(i) = dgdc1(i) + gex / yCi

                            dgdc1(i) = dgdc1(i) + gex * iExponent / (yCi - yCj)

                        ! Derivative with respect to Cj
                        else if (i == iCj) then
                            dgdc1(i) = dgdc1(i) + gex / yCj

                            dgdc1(i) = dgdc1(i) + gex * iExponent * (-1) / (yCi - yCj)

                        end if
                    end do
                end if

                ! Second sublattice - Cation:vacancy contributions
                do i = 1, nConstituentSublattice(iSPI,2)
                    if (i == iAi) then
                        dgdc2(i) = dgdc2(i) + gex / yAi

                    end if
                end do

            ! Case 2: Determine the mixing parameter type: L_Ci:Aj,Dk
            else if (iSUBIMixType(l) == 2) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    ! Establishing the site fractions for the L_Ci:Aj,Dk cases
                    if (k == 2) then
                        ! cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c

                    else if (k == 3) then
                        ! anion - Aj
                        yAi = dSiteFraction(iSPI,s,c)
                        iAi = c

                    else if (k == 4) then
                        ! Dl - anion, vacancy or neutral
                        yDi = dSiteFraction(iSPI,s,c)
                        iDi = c

                    end if
                end do

                ! Multiply prefactor term by excess Gibbs energy parameter:
                iExponent = iRegularParam(l,n+2)
                ! Excess Gibbs energy equation for L_Ci:Aj,Dk case
                gex =  dPreFactor * dExcessGibbsParam(l) * (yAi - yDi)**(iExponent)
                ! Total Excess Gibbs
                gexcess = gexcess + gex

                ! Derived Excess Equations for Mixing Case 1: L_Ci:Aj,Dk
                ! Chemical potential with respect to the first sublattice
                do i = 1, nConstituentSublattice(iSPI,1)
                    if (i == iCi) then
                        dgdc1(i) = dgdc1(i) + gex / yCi

                    end if
                end do

                ! Avoiding a situation with 0^(-1)
                if ((yAi - yDi) /= 0D0) then
                    ! Chemical potential with respect to the second sublattice
                    do i = 1, nConstituentSublattice(iSPI,2)
                        if (i == iAi) then
                            ! cation / anion
                            dgdc2(i) = dgdc2(i) + gex / yAi

                            dgdc2(i) = dgdc2(i) + gex * iExponent / (yAi - yDi)

                        else if (i == iDi) then
                            ! Dl - anion, vacancy or neutral
                            dgdc2(i) = dgdc2(i) + gex / yDi

                            dgdc2(i) = dgdc2(i) + gex * iExponent * (-1) / (yAi - yDi)

                        end if
                    end do
                end if

            ! Case 3: Determine the mixing parameter type: L_Ci,Cj:Va
            else if (iSUBIMixType(l) == 3) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Store the first and second site fractions:
                    if (k == 2) then
                        ! cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                        chargeCi = dSublatticeCharge(iSPI,s,c)
                        ! Compute prefactor term:
                        dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    else if (k == 3) then
                        ! cation - Cj
                        yCj = dSiteFraction(iSPI,s,c)
                        iCj = c
                        chargeCj = dSublatticeCharge(iSPI,s,c)
                        ! Compute prefactor term:
                        dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    else if (k == 4) then
                        ! Vacancy
                        yva = dSiteFraction(iSPI,s,c)
                        ! Compute prefactor term:
                        dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)**2

                    end if
                end do

                ! Multiply prefactor term by excess Gibbs energy parameter:
                iExponent = iRegularParam(l,n+2)
                ! Excess Gibbs energy equation for L_Ci,Cj:Va case
                gex = q * dPreFactor * dExcessGibbsParam(l) * (yCi - yCj)**(iExponent)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                ! Derived Excess Equations for Mixing Case 2: L_Ci,Cj:Va
                ! Avoiding a situation with 0^(-1)
                if ((yCi - yCj) /= 0) then
                    ! Chemical potential with respect to the first sublattice
                    do i = 1, nConstituentSublattice(iSPI,1)

                        ! Derivative with respect to Ci
                        if (i == iCi) then
                            dgdc1(i) = dgdc1(i) + gex / yCi

                            dgdc1(i) = dgdc1(i) + gex * chargeCi / q

                            dgdc1(i) = dgdc1(i) + gex * iExponent / (yCi - yCj)

                        ! Derivative with respect to Cj
                        else if (i == iCj) then
                            dgdc1(i) = dgdc1(i) + gex / yCj

                            dgdc1(i) = dgdc1(i) + gex * chargeCj / q

                            dgdc1(i) = dgdc1(i) + gex * iExponent * (-1) / (yCi - yCj)

                       else
                            !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                            ! Not yet tested
                            ! If there are additional cations that are not i or j
                            dgdc1(i) = dgdc1(i) + gex * dSublatticeCharge(iSPI,1,i) / q

                        end if
                    end do
                end if

                ! Second sublattice - Cation:vacancy contributions
                do i = 1, nConstituentSublattice(iSPI,2)
                    if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                        dgdc2(i) = dgdc2(i) + 2 * gex / yva

                    end if
                end do

            ! Case 4: Determine the mixing parameter type: L_Ci:Va,Bj
            else if (iSUBIMixType(l) == 4) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    ! Store the first and second site fractions:
                    if (k == 2) then
                        ! cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                        chargeCi = dSublatticeCharge(iSPI,s,c)

                    else if (k == 3) then
                        ! vacancy - Va
                        yva = dSiteFraction(iSPI,s,c)

                    else if (k == 4) then
                        ! nuetral - Bj
                        yBi = dSiteFraction(iSPI,s,c)
                        iBi = c

                    end if
                end do

                ! Multiply prefactor term by excess Gibbs energy parameter:
                iExponent = iRegularParam(l,n+2)
                ! Excess Gibbs energy equation for L_Ci:Va,Bj case
                gex = q * dPreFactor * dExcessGibbsParam(l) * (yCi * yva - yBi)**(iExponent)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                ! Avoiding a situation with 0^(-1)
                if ((yCi * yva - yBi) /= 0) then

                    ! First sublattice
                    do i = 1, nConstituentSublattice(iSPI,1)
                        if (i == iCi) then
                            ! cation - Ci
                            dgdc1(i) = dgdc1(i) + gex / yCi

                            dgdc1(i) = dgdc1(i) + gex * chargeCi / q

                            dgdc1(i) = dgdc1(i) + gex * iExponent * yva / (yCi * yva - yBi)

                        else
                            !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                            ! Not yet tested
                            ! If there are additional cations that are not i
                            dgdc1(i) = dgdc1(i) + gex * dSublatticeCharge(iSPI,1,i) / q

                        end if
                    end do

                    ! Second sublattice
                    do i = 1, nConstituentSublattice(iSPI,2)
                        if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                            ! vacancy contributions
                            dgdc2(i) = dgdc2(i) + gex / yva

                            dgdc2(i) = dgdc2(i) + gex * iExponent * yCi / (yCi * yva - yBi)

                        else if (i == iBi) then
                            ! neutral contributions
                            dgdc2(i) = dgdc2(i) + gex / yBi

                            dgdc2(i) = dgdc2(i) + gex * iExponent * (-1) / (yCi * yva - yBi)

                        end if
                    end do
                end if

            ! Begining of ternary mixing cases
            ! Case 5: Determine the mixing parameter type: L_Ci:Bj,Bk
            else if (iSUBIMixType(l) == 5) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Store the first and second site fractions:
                    if (k == 2) then
                        ! cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                        chargeCi = dSublatticeCharge(iSPI,s,c)

                    else if (k == 3) then
                        ! nuetral - Bj
                        yBi = dSiteFraction(iSPI,s,c)
                        iBi = c
                        ! Compute prefactor term:
                        dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    else if (k == 4) then
                        ! nuetral - Bk
                        yBj = dSiteFraction(iSPI,s,c)
                        iBj = c
                        ! Compute prefactor term:
                        dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    end if
                end do

                ! Multiply prefactor term by excess Gibbs energy parameter:
                iExponent = iRegularParam(l,n+2)
                ! Excess Gibbs energy equation for L_Ci:Bj,Bk case
                gex = q * dPreFactor * dExcessGibbsParam(l) * (yCi*yBi - yBj)**(iExponent)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                ! Derived Excess Equations for Mixing Case 1: L_Ci:Aj,Dk
                ! Chemical potential with respect to the first sublattice
                do i = 1, nConstituentSublattice(iSPI,1)
                    dgdc1(i) = dgdc1(i) + gex * dSublatticeCharge(iSPI,1,i) / q

                end do

                ! Avoiding a situation with 0^(-1)
                if ((yBi - yBj) /= 0D0) then
                    do i = 1, nConstituentSublattice(iSPI,2)
                        if (i == iBi) then
                            ! cation / anion
                            dgdc2(i) = dgdc2(i) + gex / yBi

                            dgdc2(i) = dgdc2(i) + gex * iExponent / (yBi - yBj)

                        else if (i == iBj) then
                            ! Dl - anion, vacancy or neutral
                            dgdc2(i) = dgdc2(i) + gex / yBj

                            dgdc2(i) = dgdc2(i) + gex * iExponent * (-1) / (yBi - yBj)

                        end if
                    end do
                end if

            ! Case 6: Determine the mixing parameter type: L_Ci,Cj,Ck:Al
            else if (iSUBIMixType(l) == 6) then
                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    if (k == 2) then
                        ! Cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                    else if (k == 3) then
                        ! Cation - Cj
                        yCj = dSiteFraction(iSPI,s,c)
                        iCj = c
                    else if (k == 4) then
                        ! Cation - Ck
                        yCk = dSiteFraction(iSPI,s,c)
                        iCk = c
                    else if (k == 5) then
                        ! Anion - Al
                        yAi = dSiteFraction(iSPI,s,c)
                        iAi = c
                    end if
                end do

            ! Multiply prefactor term by excess Gibbs energy parameter:
            iExponent = iRegularParam(l,n+2)
            ! Calculate f and v
            f = (1D0 - yCi - yCj - yCk)/3D0
            if (iRegularParam(l,n+2) == 0) then
                v = yCi + f
            else if (iRegularParam(l,n+2) == 1) then
                v = yCj + f
            else if (iRegularParam(l,n+2) == 2) then
                v = yCk + f
            end if

            ! Excess Gibbs energy equation for L_Ci,Cj,Ck:Al case
            gex = dPreFactor * v * dExcessGibbsParam(l)
            ! Total Excess Gibbs Energy
            gexcess = gexcess + gex

            do i = 1, nConstituentSublattice(iSPI, 1)
                ! Derivative with respect to C
                if (i == iCi) then
                    ! prefactor part
                    dgdc1(i) = dgdc1(i) + gex / yCi
                    ! f part
                    dgdc1(i) = dgdc1(i) - gex / (3D0 * v)
                    ! other v part
                    if (iRegularParam(l,n+2) == 0) dgdc1(i) = dgdc1(i) + gex / v
                ! Derivative with respect to Cj
                else if (i == iCj) then
                    ! prefactor part
                    dgdc1(i) = dgdc1(i) + gex / yCj
                    ! f part
                    dgdc1(i) = dgdc1(i) - gex / (3D0 * v)
                    ! other v part
                    if (iRegularParam(l,n+2) == 1) dgdc1(i) = dgdc1(i) + gex / v
                else if (i == iCk) then
                    ! prefactor part
                    dgdc1(i) = dgdc1(i) + gex / yCk
                    ! f part
                    dgdc1(i) = dgdc1(i) - gex / (3D0 * v)
                    ! other v part
                    if (iRegularParam(l,n+2) == 2) dgdc1(i) = dgdc1(i) + gex / v
                end if
            end do

            do i = 1, nConstituentSublattice(iSPI,2)
              ! Derivative with respect to Bk
                if (i == iAi) then
                    ! prefactor part
                    dgdc2(i) = dgdc2(i) + gex / yAi
                end if
            end do

            ! Case 7: Determine the mixing parameter type: L_Ci:Aj,Dk,Dl
            else if (iSUBIMixType(l) == 7) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    if (k == 2) then
                        ! Cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                    else if (k == 3) then
                        ! Anion - Aj
                        yAi = dSiteFraction(iSPI,s,c)
                        iAi = c
                    else if (k == 4) then
                        ! Dk
                        yDi = dSiteFraction(iSPI,s,c)
                        iDi = c
                    else if (k == 5) then
                        ! Dl
                        yDj = dSiteFraction(iSPI,s,c)
                        iDj = c
                    end if
                end do

                ! Multiply prefactor term by excess Gibbs energy parameter:
                iExponent = iRegularParam(l,n+2)
                ! Calculate f and v
                f = (1D0 - yAi - yDi - yDj)/3D0
                if (iRegularParam(l,n+2) == 0) then
                    v = yAi + f
                else if (iRegularParam(l,n+2) == 1) then
                    v = yDi + f
                else if (iRegularParam(l,n+2) == 2) then
                    v = yDj + f
                end if

                ! Excess Gibbs energy equation for L_Ci:Aj,Dk,Dl case
                gex = dPreFactor * v * dExcessGibbsParam(l)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex
                ! First Sublattice
                do i = 1, nConstituentSublattice(iSPI, 1)
                    ! Derivative with respect to C
                    if (i == iCi) then
                        ! prefactor part
                        dgdc1(i) = dgdc1(i) + gex / yCi
                    end if
                end do

                ! Second Sublattice
                do i = 1, nConstituentSublattice(iSPI,2)
                  ! Derivative with respect to Bk
                    if (i == iAi) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yAi
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 0) dgdc2(i) = dgdc2(i) + gex / v
                    else if (i == iDi) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yDi
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 1) dgdc2(i) = dgdc2(i) + gex / v
                    else if (i == iDj) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yDj
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 2) dgdc2(i) = dgdc2(i) + gex / v
                    end if
                end do

            ! Case 8: Determine the mixing parameter type: L_Ci,Cj,Ck:Va
            else if (iSUBIMixType(l) == 8) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    if (k == 2) then
                        ! Cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                    else if (k == 3) then
                        ! Cation - Cj
                        yCj = dSiteFraction(iSPI,s,c)
                        iCj = c
                    else if (k == 4) then
                        ! Cation - Ck
                        yCk = dSiteFraction(iSPI,s,c)
                        iCk = c
                    end if
                end do
                ! Compute prefactor term:
                dPreFactor = yCi * yCj * yCk * yva**3
                ! Multiply prefactor term by excess Gibbs energy parameter:
                iExponent = iRegularParam(l,n+2)
                ! Calculate f and v
                f = (1D0 - yva * yCi - yva * yCj - yva * yCk)/3D0
                if (iRegularParam(l,n+2) == 0) then
                    v = yva * yCi + f
                else if (iRegularParam(l,n+2) == 1) then
                    v = yva * yCj + f
                else if (iRegularParam(l,n+2) == 2) then
                    v = yva * yCk + f
                end if

                ! Excess Gibbs energy equation for L_Ci,Cj,Ck:Al case
                gex = q * dPreFactor * v * dExcessGibbsParam(l)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                do i = 1, nConstituentSublattice(iSPI, 1)
                    ! Q part:
                    dgdc1(i) = dgdc1(i) + gex * dSublatticeCharge(iSPI,1,i) / q
                    ! Derivative with respect to C
                    if (i == iCi) then
                        ! prefactor part
                        dgdc1(i) = dgdc1(i) + gex / yCi
                        ! f part
                        dgdc1(i) = dgdc1(i) - yva * gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 0) dgdc1(i) = dgdc1(i) + yva * gex / v
                    ! Derivative with respect to Cj
                    else if (i == iCj) then
                        ! prefactor part
                        dgdc1(i) = dgdc1(i) + gex / yCj
                        ! f part
                        dgdc1(i) = dgdc1(i) - yva * gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 1) dgdc1(i) = dgdc1(i) + yva * gex / v
                    else if (i == iCk) then
                        ! prefactor part
                        dgdc1(i) = dgdc1(i) + gex / yCk
                        ! f part
                        dgdc1(i) = dgdc1(i) - yva * gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 2) dgdc1(i) = dgdc1(i) + yva * gex / v
                    end if
                end do

                do i = 1, nConstituentSublattice(iSPI,2)
                    ! Derivative with respect to Va
                    if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex * 3D0 / yva
                        ! f part
                        dgdc2(i) = dgdc2(i) - (yCi + yCj + yCk) * gex / (3D0 * v)
                        ! other v parts
                        if (iRegularParam(l,n+2) == 0) dgdc2(i) = dgdc2(i) + yCi * gex / v
                        if (iRegularParam(l,n+2) == 1) dgdc2(i) = dgdc2(i) + yCj * gex / v
                        if (iRegularParam(l,n+2) == 2) dgdc2(i) = dgdc2(i) + yCk * gex / v
                    end if
                end do

            ! Case 9: Determine the mixing parameter type: L_Ci:Va,Bj,Bk
            else if (iSUBIMixType(l) == 9) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    if (k == 2) then
                        ! Cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                    else if (k == 3) then
                        ! vacancy
                        yva = dSiteFraction(iSPI,s,c)
                    else if (k == 4) then
                        ! neutral - Bi
                        yBi = dSiteFraction(iSPI,s,c)
                        iBi = c
                    else if (k == 5) then
                        ! Neutral - Bj
                        yBj = dSiteFraction(iSPI,s,c)
                        iBj = c
                    end if
                end do

                ! Compute prefactor term:
                dPreFactor = yCi * yva * yBi * yBj
                ! Calculate f and v
                f = (1D0 - yCi * yva - yBi - yBj)/3D0
                if (iRegularParam(l,n+2) == 0) then
                    v = yCi * yva + f
                else if (iRegularParam(l,n+2) == 1) then
                    v = yBi + f
                else if (iRegularParam(l,n+2) == 2) then
                    v = yBj + f
                end if
                ! Excess Gibbs energy equation for L_Ci:Va,Bj,Bk case
                gex = q * dPreFactor * v * dExcessGibbsParam(l)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                do i = 1, nConstituentSublattice(iSPI, 1)
                    ! Q part:
                    dgdc1(i) = dgdc1(i) + gex * dSublatticeCharge(iSPI,1,i) / q
                    ! Derivative with respect to Ci
                    if (i == iCi) then
                        ! prefactor part
                        dgdc1(i) = dgdc1(i) + gex / yCi
                        ! f part
                        dgdc1(i) = dgdc1(i) - yva * gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 0) dgdc1(i) = dgdc1(i) + yva * gex / v
                    end if
                end do

                do i = 1, nConstituentSublattice(iSPI,2)
                    ! Derivative with respect to Bk
                    if (i == iBi) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yBi
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 1) dgdc2(i) = dgdc2(i) + gex / v
                    ! Derivative with respect to Bl
                    else if (i == iBj) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yBj
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 2) dgdc2(i) = dgdc2(i) + gex / v
                    ! Derivative with respect to Va
                    else if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex * 2 / yva
                        ! f part
                        dgdc2(i) = dgdc2(i) - yCi * gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 0) dgdc2(i) = dgdc2(i) + yCi * gex / v
                    end if
                end do

            ! Case 10: Determine the mixing parameter type: L_Ci:Bj,Bk,Bl
          else if (iSUBIMixType(l) == 10) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    if (k == 3) then
                        ! neutral - Bj
                        yBi = dSiteFraction(iSPI,s,c)
                        iBi = c
                    else if (k == 4) then
                        ! neutral - Bk
                        yBj = dSiteFraction(iSPI,s,c)
                        iBj = c
                    else if (k == 5) then
                        ! neutral - Bl
                        yBk = dSiteFraction(iSPI,s,c)
                        iBk = c
                    end if
                end do

                ! Compute prefactor term:
                dPreFactor = yBi * yBj * yBk
                ! Calculate f and v
                f = (1D0 - yBi - yBj - yBk)/3D0
                if (iRegularParam(l,n+2) == 0) then
                    v = yBi + f
                else if (iRegularParam(l,n+2) == 1) then
                    v = yBj + f
                else if (iRegularParam(l,n+2) == 2) then
                    v = yBk + f
                end if
                ! Excess Gibbs energy equation for L_Ci:Bj,Bk,Bl case
                gex = q * dPreFactor * v * dExcessGibbsParam(l)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                do i = 1, nConstituentSublattice(iSPI, 1)
                    ! Q part:
                    dgdc1(i) = dgdc1(i) + gex * dSublatticeCharge(iSPI,1,i) / q
                end do

                do i = 1, nConstituentSublattice(iSPI,2)
                    ! Derivative with respect to Bk
                    if (i == iBi) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yBi
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 0) dgdc2(i) = dgdc2(i) + gex / v
                    ! Derivative with respect to Bl
                    else if (i == iBj) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yBj
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 1) dgdc2(i) = dgdc2(i) + gex / v
                    ! Derivative with respect to Va
                    else if (i == iBk) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yBk
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 2) dgdc2(i) = dgdc2(i) + gex / v
                    end if
                end do

            ! Case 11: Determine the mixing parameter type: L_Ci,Cj:Ak,Dl
            else if (iSUBIMixType(l) == 11) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    ! Store the first and second site fractions:
                    if (k == 2) then
                        ! cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                    else if (k == 3) then
                        ! cation - Cj
                        yCj = dSiteFraction(iSPI,s,c)
                        iCj = c
                    else if (k == 4) then
                        ! anion - Ak
                        yAi = dSiteFraction(iSPI,s,c)
                        iAi = c
                    else if (k == 5) then
                        ! Dl
                        yDi = dSiteFraction(iSPI,s,c)
                        iDi = c
                    end if
                end do

                ! Multiply prefactor term by excess Gibbs energy parameter:
                iExponent = iRegularParam(l,n+2)
                ! Excess Gibbs energy equation for L_Ci,Cj:Ak,Dl case
                ! Dealing with the reciprocals
                ! Part 1 of equation:
                if (MOD(iExponent,2) == 0) then
                    gex = dPreFactor * dExcessGibbsParam(l) * (yCi - yCj)**(iExponent)
                end if
                ! Part 2 of equation:
                if (MOD(iExponent,2) == 1) then
                    gex = dPreFactor * dExcessGibbsParam(l) * (yAi - yDi)**(iExponent)
                end if

                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                ! Avoiding a situation with 0^(-1)
                if ((yCi - yCj) /= 0) then
                    do i = 1, nConstituentSublattice(iSPI,1)
                        ! Derivative with respect to Ci
                        if (i == iCi) then
                            dgdc1(i) = dgdc1(i) + gex / yCi
                            ! Dealing with the reciprocals
                            if (MOD(iExponent,2) == 0) then
                                dgdc1(i) = dgdc1(i) + gex * iExponent / (yCi - yCj)
                            end if

                        ! Derivative with respect to Cj
                        else if (i == iCj) then
                            dgdc1(i) = dgdc1(i) + gex / yCj
                            ! Dealing with the reciprocals
                            if (MOD(iExponent,2) == 0) then
                                dgdc1(i) = dgdc1(i) + gex * iExponent * (-1) / (yCi - yCj)
                            end if
                        end if
                    end do
                end if

                ! Avoiding a situation with 0^(-1)
                if ((yAi - yDi) /= 0) then

                    ! Chemical potential with respect to the second sublattice
                    do i = 1, nConstituentSublattice(iSPI,2)
                        !! When Dl = Al or Bl - Cases untested... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if (i == iAi) then
                            ! anion - Ak
                            dgdc2(i) = dgdc2(i) + gex / yAi
                            ! Dealing with the reciprocals
                            if (MOD(iExponent,2) == 1) then
                                dgdc2(i) = dgdc2(i) + gex * iExponent / (yAi - yDi)
                            end if

                        else if (i == iDi) then
                            ! Dl
                            dgdc2(i) = dgdc2(i) + gex / yDi
                            ! Dealing with the reciprocals
                            if (MOD(iExponent,2) == 1) then
                                dgdc2(i) = dgdc2(i) + gex * iExponent * (-1) / (yAi - yDi)
                            end if
                        end if
                    end do
                end if

            ! Case 12: Determine the mixing parameter type: L_Ci,Cj:Va,Bk
            else if (iSUBIMixType(l) == 12) then

                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
                    ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    if (k == 2) then
                        ! Cation - Ci
                        yCi = dSiteFraction(iSPI,s,c)
                        iCi = c
                    else if (k == 3) then
                        ! Cation - Cj
                        yCj = dSiteFraction(iSPI,s,c)
                        iCj = c
                    else if (k == 4) then
                        ! Vacancy
                        yva = dSiteFraction(iSPI,s,c)
                    else if (k == 5) then
                        ! Neutral - Bk
                        yBi = dSiteFraction(iSPI,s,c)
                        iBi = c
                    end if
                end do

                ! Compute prefactor term:
                dPreFactor = yCi * yCj * yBi * yva**2
                ! Calculate f and v
                f = (1D0 - yCi * yva - yCj * yva - yBi)/3D0
                if (iRegularParam(l,n+2) == 0) then
                    v = yCi * yva + f
                else if (iRegularParam(l,n+2) == 1) then
                    v = yCj * yva + f
                else if (iRegularParam(l,n+2) == 2) then
                    v = yBi + f
                end if
                ! Excess Gibbs energy equation for L_Ci,Cj:Va,Bk case
                gex = q * dPreFactor * v * dExcessGibbsParam(l)
                ! Total Excess Gibbs Energy
                gexcess = gexcess + gex

                do i = 1, nConstituentSublattice(iSPI, 1)
                    ! Q part:
                    dgdc1(i) = dgdc1(i) + gex * dSublatticeCharge(iSPI,1,i) / q
                    ! Derivative with respect to C
                    if (i == iCi) then
                        ! prefactor part
                        dgdc1(i) = dgdc1(i) + gex / yCi
                        ! f part
                        dgdc1(i) = dgdc1(i) - yva * gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 0) dgdc1(i) = dgdc1(i) + yva * gex / v
                    ! Derivative with respect to Cj
                    else if (i == iCj) then
                        ! prefactor part
                        dgdc1(i) = dgdc1(i) + gex / yCj
                        ! f part
                        dgdc1(i) = dgdc1(i) - yva * gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 1) dgdc1(i) = dgdc1(i) + yva * gex / v
                    end if
                end do

                do i = 1, nConstituentSublattice(iSPI,2)
                  ! Derivative with respect to Bk
                    if (i == iBi) then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex / yBi
                        ! f part
                        dgdc2(i) = dgdc2(i) - gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 2) dgdc2(i) = dgdc2(i) + gex / v
                    ! Derivative with respect to Va
                    else if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                        ! prefactor part
                        dgdc2(i) = dgdc2(i) + gex * 2 / yva
                        ! f part
                        dgdc2(i) = dgdc2(i) - (yCi + yCj) * gex / (3D0 * v)
                        ! other v part
                        if (iRegularParam(l,n+2) == 0) dgdc2(i) = dgdc2(i) + yCi * gex / v
                        if (iRegularParam(l,n+2) == 1) dgdc2(i) = dgdc2(i) + yCj * gex / v
                    end if
                end do


            else
                print *, 'Unrecognized excess mixing term in SUBI phase ', cSolnPhaseName(iSolnIndex)
                INFOThermo = 36
                return

            end if

            !print*,"gexcess:", (gexcess)*dIdealConstant * dTemperature
            !print*,"l",l
            !print*,"gextest(l)",gextest(l)
            !print*,"-------------------------"
            !print*,""
        end do LOOP_Param
        !print*,"gextest(:)",gextest(:)*dIdealConstant * dTemperature
        !print*,"Sum(gextest(:))",Sum(gextest)*dIdealConstant * dTemperature

        !print*,"cConstituentNameSUB(iSPI,s,c)",cConstituentNameSUB(iSPI,s,c)
        !print*,"dSiteFraction(iSPI,s,c)",dSiteFraction(iSPI,s,c)


        ! REFERENCE GIBBS ENERGY AND IDEAL MIXING
        ! ---------------------------------------
        gref = 0D0
        do j = iFirst, iLast
            ! Relative species index:
            n = j - iFirst + 1

            ! Store constituent indices:
            l1 = iConstituentSublattice(iSPI,1,n)
            l2 = iConstituentSublattice(iSPI,2,n)

            if (cConstituentNameSUB(iSPI,2,l2) == 'Va') then
                ! cation / vacancy
                gref = gref + q * dSiteFraction(iSPI,1,l1) * dSiteFraction(iSPI,2,l2) * dStdGibbsEnergy(j)
            else if (dSublatticeCharge(iSPI,2,l2) == 0D0) then
                ! neutral
                gref = gref + q * dSiteFraction(iSPI,2,l2) * dStdGibbsEnergy(j)
            else
                ! cation / anion
                gref = gref + dSiteFraction(iSPI,1,l1) * dSiteFraction(iSPI,2,l2) * dStdGibbsEnergy(j)
            end if
        end do

        gideal = 0D0
        do i = 1, nConstituentSublattice(iSPI,1)
            gideal = gideal + p * dSiteFraction(iSPI,1,i) * DLOG(dSiteFraction(iSPI,1,i))
        end do

        do i = 1, nConstituentSublattice(iSPI,2)
            gideal = gideal + q * dSiteFraction(iSPI,2,i) * DLOG(dSiteFraction(iSPI,2,i))
        end do

        ! print*,"p",p,"   q",q
        ! print*,"gref",gref*dIdealConstant * dTemperature
        ! print*,"gideal",gideal*dIdealConstant * dTemperature
        ! print*,"gref+gideal:", (gref+gideal)*dIdealConstant * dTemperature
        ! print*,"gref+gideal+gexcess:", (gref+gideal+gexcess)*dIdealConstant * dTemperature
        ! print*,""

        ! For Sublattice Number 1
        do i = 1, nConstituentSublattice(iSPI,1)
            lc1 = dSublatticeCharge(iSPI,1,i)
            ! Reference
            do j = iFirst, iLast
                ! Relative species index:
                n = j - iFirst + 1
                ! Store constituent indices:
                l1 = iConstituentSublattice(iSPI,1,n)
                l2 = iConstituentSublattice(iSPI,2,n)

                if (cConstituentNameSUB(iSPI,2,l2) == 'Va') then
                    ! cation / vacancy
                    if (i == l1) dgdc1(i) = dgdc1(i) + q * dSiteFraction(iSPI,2,l2) * dStdGibbsEnergy(j)
                    dgdc1(i) = dgdc1(i) + lc1 * dSiteFraction(iSPI,1,l1) * dSiteFraction(iSPI,2,l2) * dStdGibbsEnergy(j)
                else if (dSublatticeCharge(iSPI,2,l2) == 0D0) then
                    ! neutral
                    dgdc1(i) = dgdc1(i) + lc1 * dSiteFraction(iSPI,2,l2) * dStdGibbsEnergy(j)
                else
                    ! cation / anion
                    if (i == l1) dgdc1(i) = dgdc1(i) + dSiteFraction(iSPI,2,l2) * dStdGibbsEnergy(j)
                end if

            end do

            ! Entropy
            dgdc1(i) = dgdc1(i) + (1 + DLOG(dSiteFraction(iSPI,1,i))) * p
            do j = 1, nConstituentSublattice(iSPI,1)
                dgdc1(i) = dgdc1(i) + dSublatticeCharge(iSPI,1,i) * yva &
                                    * dSiteFraction(iSPI,1,j)*DLOG(dSiteFraction(iSPI,1,j))
            end do
            do j = 1, nConstituentSublattice(iSPI,2)
                dgdc1(i) = dgdc1(i) + dSublatticeCharge(iSPI,1,i) * dSiteFraction(iSPI,2,j) * DLOG(dSiteFraction(iSPI,2,j))
            end do
        end do

        ! For Sublattice Number 2
        do i = 1, nConstituentSublattice(iSPI,2)
            do j = iFirst, iLast
                ! Relative species index:
                n = j - iFirst + 1

                ! Store constituent indices:
                l1 = iConstituentSublattice(iSPI,1,n)
                l2 = iConstituentSublattice(iSPI,2,n)

                if (cConstituentNameSUB(iSPI,2,l2) == 'Va') then
                    ! cation / vacancy
                    if (i == l2) dgdc2(i) = dgdc2(i) + q * dSiteFraction(iSPI,1,l1) * dStdGibbsEnergy(j)
                else if (dSublatticeCharge(iSPI,2,l2) == 0D0) then
                    ! neutral
                    if (i == l2) dgdc2(i) = dgdc2(i) + q * dStdGibbsEnergy(j)
                else
                    ! cation / anion
                    if (i == l2) dgdc2(i) = dgdc2(i) + dSiteFraction(iSPI,1,l1) * dStdGibbsEnergy(j)
                end if
            end do
            ! Entropy
            dgdc2(i) = dgdc2(i) + (1 + DLOG(dSiteFraction(iSPI,2,i))) * q
            do j = 1, nConstituentSublattice(iSPI,1)
                if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                    ! cation / vacancy
                    dgdc2(i) = dgdc2(i) + q * dSiteFraction(iSPI,1,j) * DLOG(dSiteFraction(iSPI,1,j))
                else if (dSublatticeCharge(iSPI,2,i) == 0D0) then
                    ! neutral
                else
                    ! cation / anion
                    dgdc2(i) = dgdc2(i) + (-dSublatticeCharge(iSPI,2,i)) * dSiteFraction(iSPI,1,j) * DLOG(dSiteFraction(iSPI,1,j))
                end if
            end do
        end do
        ! print*,"After All"
        ! print*,"dgdc1(:)",dgdc1(:)*dIdealConstant * dTemperature
        ! print*,"dgdc2(:)",dgdc2(:)*dIdealConstant * dTemperature
        ! print*,""

        cc1 = 0D0

        ! Compute the chemical potential for each phase component assuming ideal mixing:
        LOOP_Ideal: do i = iFirst, iLast
            ! print*,""
            ! print*,"i",i
            ! print*,""
            ! Relative species index:
            m = i - iFirst + 1
            k1 = iConstituentSublattice(iSPI,1,m)
            k2 = iConstituentSublattice(iSPI,2,m)
            !cc1 = dSublatticeCharge(iSPI,1,k1)
            !print*,"k1",k1

            kc1 = 1D0
            ! cation / vacancy
            if (cConstituentNameSUB(iSPI,2,k2) == 'Va') then
                kc2 = 1D0
                cc1 = dSublatticeCharge(iSPI,1,k1)
                natom = 1D0 / dMol + dMolDerivatives(m) * dMolAtoms
            ! neutral
            else if (dSublatticeCharge(iSPI,2,k2) == 0D0) then
                kc2 = 1D0
                cc1 = 1D0
                natom = 1D0 / dMol + dMolDerivatives(m) * dMolAtoms
            ! cation / anion
            else
                if (k1 > 0) kc1 = dSublatticeCharge(iSPI,1,k1)
                kc2 = -dSublatticeCharge(iSPI,2,k2)
                cc1 = dSublatticeCharge(iSPI,1,k1)
                natom = (kc1 + kc2) / dMol + dMolDerivatives(m) * dMolAtoms
            end if

            dChemicalPotential(i) = (gref + gideal + gexcess) * natom

            do j = 1, nConstituentSublattice(iSPI,1)
                ! cation / (anion or vacancy)
                if (k1 > 0) then
                    dydn = -kc2 * dSiteFraction(iSPI,1,j) / dSub1Total
                    ! print*,"dydn - Cation Sublattice - 1:   ",dydn
                    if (j == k1) dydn = dydn + kc2 / dSub1Total
                    dChemicalPotential(i) = dChemicalPotential(i) + dydn * dgdc1(j) * dMolAtoms / dMol
                    ! print*,"dydn - Cation Sublattice - Final:   ",dydn
                    ! print*,"dydn * dgdc1(j) * dMolAtoms / dMol",dydn * dgdc1(j) * dMolAtoms / dMol
                    ! print*,""
                end if
            end do
            ! print*,""

            dTest = 0D0

            do j = 1, nConstituentSublattice(iSPI,2)
                ! print*,"j",j
                ! cation / vacancy
                if (cConstituentNameSUB(iSPI,2,k2) == 'Va') then
                    dydn = -dSiteFraction(iSPI,2,j)*(dSub1Total*cc1+(q-cc1)*dSub2Total*yva)/(dSub1Total*dSub2Total*q)
                    if (cConstituentNameSUB(iSPI,2,j) == 'Va') then
                        ! vacancy
                        dydn = dydn + (dSub1Total*cc1+(q-cc1)*dSub2Total*yva)/(dSub1Total*dSub2Total*q)
                    else
                        ! neutral or anion
                    end if
                ! neutral
                else if (dSublatticeCharge(iSPI,2,k2) == 0D0) then
                    dydn = -dSiteFraction(iSPI,2,j) / dSub2Total
                    if (cConstituentNameSUB(iSPI,2,j) == 'Va') then
                        ! vacancy
                    else
                        ! neutral or anion
                        if (j == k2) dydn = dydn + 1 / dSub2Total
                    end if
                ! cation / anion
                else
                    dydn = -dSiteFraction(iSPI,2,j) * (cc1 / dSub2Total + kc2*yva*(q-cc1)/(dSub1Total*q))
                    if (cConstituentNameSUB(iSPI,2,j) == 'Va') then
                        dydn = dydn + dSiteFraction(iSPI,2,j) * kc2*(q-cc1)/(dSub1Total*q)
                    else if (dSublatticeCharge(iSPI,2,j) == 0D0) then
                        ! neutral
                    else
                        ! cation / anion
                        if (j == k2) dydn = dydn + cc1 / dSub2Total
                    end if
                end if
                ! print*,"dydn - Anion Sublattice:   ",dydn
                dChemicalPotential(i) = dChemicalPotential(i) + dydn * dgdc2(j) * dMolAtoms / dMol
            end do
            ! print*,""
            ! print*,"-------------------------------------"
            ! print*,""

        end do LOOP_Ideal

        ! print*,"_______________________________________________________________________________"
        ! print *, ' System Component 2', ' Mass [mol]  ', 'Chemical potential [J/mol]'
        ! print *, ' ---------------- ', ' ----------  ', '--------------------------'
        ! do i = 1, nSpecies
        !     print '(A14,A1,ES15.4,A1, ES14.6)', cSpeciesName(i), ' ', dMolesSpecies(i), ' ', &
        !     dChemicalPotential(i) * dIdealConstant * dTemperature
        !
        ! end do
        ! print*,"_______________________________________________________________________________"

        deallocate(dgdc1,dgdc2)
        deallocate(dMolDerivatives)
    end if IF_SUBL

    return

end subroutine CompExcessGibbsEnergySUBI
