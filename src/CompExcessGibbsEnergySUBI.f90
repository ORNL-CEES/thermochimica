
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
    integer :: iSolnIndex, nSublattice, iSPI, iExponent
    integer :: iMixTypeVa , iMixTypeBi, iMixTypeAni, iMixTypeCati, iMixTypeCatj
    integer :: iFirst, iLast, iFirstParam, iSecondParam, iSubParam, iFirstParam2, iSecondParam2, iSubParam2
    integer, dimension(:), allocatable:: iMixType
    real(8) :: dSub1Total, dSub2Total, dydn
    real(8) :: dSum, p, q, kc1, kc2, lc1, lc2, gref, gideal, gexcess, natom, yva, dMol, dMolAtoms
    real(8), dimension(:), allocatable :: dgdc1, dgdc2, dMolDerivatives
    real(8) :: dTemp, dPreFactor, dFirstParam, dSecondParam, dFirstParam2, dSecondParam2, dThirdParam
    real(8) :: y1, y2, y3
    character(8) :: cDummyBi, cDummyCi, cDummyCj, cDummyAi, cDummyAj

print*,""
print*,""
print*,"                          CompExcessGibbsEnergySUBI.f90"
print*,""
print*,""

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
        if (allocated(iMixType)) deallocate(iMixType)
        allocate(dgdc1(nConstituentSublattice(iSPI,1)),dgdc2(nConstituentSublattice(iSPI,2)))
        allocate(dMolDerivatives(iLast - iFirst + 1))
        allocate(iMixType(nParamPhase(iSolnIndex)))

        ! Initialize variables:
        dSiteFraction(iSPI,1:nSublattice,1:nMaxConstituentSys) = 0D0
        dChemicalPotential(iFirst:iLast)                       = 0D0
        dPartialExcessGibbs(iFirst:iLast)                      = 0D0
        dgdc1                                                  = 0D0
        dgdc2                                                  = 0D0
        dMolDerivatives                                        = 0D0
        iMixType                                               = 0D0

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
        do i = 1, nConstituentSublattice(iSPI,1)
            dSiteFraction(iSPI,1,i) = dSiteFraction(iSPI,1,i) / dSub1Total
        end do

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
            if ((cConstituentNameSUB(iSPI,2,c) == 'Va') .OR. &
                (dSublatticeCharge(iSPI,2,c) == 0D0)) then
                ! Vacancy or neutral get scaled by Q
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
        do i = 1, nConstituentSublattice(iSPI,2)
            dSiteFraction(iSPI,2,i) = dSiteFraction(iSPI,2,i) / dSub2Total
            ! find site fraction of vacancies
            if (cConstituentNameSUB(iSPI,2,i) == 'Va') yva = dSiteFraction(iSPI,2,i)
        end do

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
        do j = iFirst, iLast
            ! Relative species index:
            n = j - iFirst + 1

            ! Store constituent indices:
            l1 = iConstituentSublattice(iSPI,1,n)
            l2 = iConstituentSublattice(iSPI,2,n)

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
                dMolDerivatives(n) = -dSub2Total * lc1
                do i = 1, nConstituentSublattice(iSPI,1)
                    dMolDerivatives(n) = dMolDerivatives(n) + dSub2Total*dSiteFraction(iSPI,1,i)*dSublatticeCharge(iSPI,1,i)
                end do
                do i = 1, nConstituentSublattice(iSPI,2)
                    ! Only include anions (not neutrals or vacancies)
                    if (.NOT.(cConstituentNameSUB(iSPI,2,i) == 'Va') .OR. &
                        (dSublatticeCharge(iSPI,2,i) == 0D0)) &
                    dMolDerivatives(n) = dMolDerivatives(n) + dSub1Total*dSiteFraction(iSPI,2,i)*(-dSublatticeCharge(iSPI,2,i))
                end do
            else if (dSublatticeCharge(iSPI,2,l2) == 0D0) then
                ! neutral
                do i = 1, nConstituentSublattice(iSPI,2)
                    ! Only include anions (not neutrals or vacancies)
                    if (.NOT.(cConstituentNameSUB(iSPI,2,i) == 'Va') .OR. &
                        (dSublatticeCharge(iSPI,2,i) == 0D0)) &
                    dMolDerivatives(n) = dMolDerivatives(n) + dSub1Total*dSiteFraction(iSPI,2,i)*(-dSublatticeCharge(iSPI,2,i))
                end do
            else
                ! cation / anion
                dMolDerivatives(n) = -(dSub1Total + dSub2Total) * lc1 * lc2
                do i = 1, nConstituentSublattice(iSPI,1)
                    dMolDerivatives(n) = dMolDerivatives(n)+dSub2Total*lc2*dSiteFraction(iSPI,1,i)*dSublatticeCharge(iSPI,1,i)
                end do
                do i = 1, nConstituentSublattice(iSPI,2)
                    ! Only include anions (not neutrals or vacancies)
                    if (.NOT.((cConstituentNameSUB(iSPI,2,i) == 'Va') .OR. &
                    (dSublatticeCharge(iSPI,2,i) == 0D0))) &
                    dMolDerivatives(n) = dMolDerivatives(n)+dSub1Total*lc1*dSiteFraction(iSPI,2,i)*(-dSublatticeCharge(iSPI,2,i))
                end do
            end if
            dMolDerivatives(n) = dMolDerivatives(n) / (dSub1Total*dSub2Total*dMol**2)
        end do

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
!print*,"dMolFraction(i) - d",dMolFraction(i)
                if (c > 0) dMolFraction(i) = dMolFraction(i) * dSiteFraction(iSPI,2,c)
!print*,"dMolFraction(i) - c",dMolFraction(i)
            else
                if (d > 0) dMolFraction(i) = dMolFraction(i) * dSiteFraction(iSPI,1,d) / dSublatticeCharge(iSPI,1,d)
!print*,"dMolFraction(i) - d",dMolFraction(i)
                if (c > 0) dMolFraction(i) = dMolFraction(i) * dSiteFraction(iSPI,2,c)
!print*,"dMolFraction(i) - c",dMolFraction(i)
            end if
!print*,""
            dSum = dSum + dMolFraction(i)
        end do LOOP_CorrectX

        ! Normalize mole fractions and compute number of mole atoms per mole (yes that makes sense, don't think about it)
        dMolAtoms = 0D0
        do i = iFirst, iLast
            dMolFraction(i)  = dMolFraction(i) / dSum
!print*,"i",i
!print*,"dMolFraction(i) - normalized",dMolFraction(i)
!print*,""
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
!print*,"dMolAtoms",dMolAtoms
!print*,""
        dStoichSublattice(iSPI,1) = p
        dStoichSublattice(iSPI,2) = q

        ! EXCESS TERMS
        ! ------------
        gexcess = 0D0

        ! Loop through parameters:
        LOOP_Param: do l = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)
          ! Return control to the parent subroutine if there aren't any interaction parameters for this phase:
          if (nParamPhase(iSolnIndex) - nParamPhase(iSolnIndex-1) == 0) exit LOOP_Param
          ! Reinitialize temporary variable:
            dPreFactor = 1D0
            iFirstParam = 0
            iSecondParam = 0
            iFirstParam2 = 0
            iSecondParam2 = 0
            dFirstParam = 0
            dSecondParam = 0
            dThirdParam = 0

          ! Store the number of constituents involved in this parameter:
          n = iRegularParam(l,1)

          ! Determine the mixing parameter type: L_Ci:Aj,Dk
          if ((iSUBIParamData(l,1) == 1) .AND. &
              (iSUBIParamData(l,2) == 3) .AND. &
              (iSUBIParamData(l,3) == 2)) then
                ! Defines the mixing case type: L_Ci:Aj,Dk
                iMixType(l) = 1

            ! Loop through constituents associated with this parameter:
            do k = 2, n + 1

              ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    ! Store the first and second site fractions:
                    if (k == iSUBIParamData(l,2)) then
                        dFirstParam = dSiteFraction(iSPI,s,c)

                    else if (k == iSUBIParamData(l,2) + 1) then
                        dSecondParam = dSiteFraction(iSPI,s,c)

                    end if
            end do

            ! Multiply prefactor term by excess Gibbs energy parameter:
            iExponent = iRegularParam(l,n+2)
            ! Excess Gibbs energy equation for L_Ci:Aj,Dk case
            gexcess = gexcess + dPreFactor * dExcessGibbsParam(l) * (dFirstParam - dSecondParam)**(iExponent)


            ! Determine the mixing parameter type: L_Ci,Cj:(Va or Ak)
          else if ((iSUBIParamData(l,1) == 1) .AND. &
                   (iSUBIParamData(l,2) == 2) .AND. &
                   (iSUBIParamData(l,3) == 2)) then

                   ! Loop through constituents associated with this parameter:
                   do k = 2, n + 1

                     ! Determine constituent and sublattice indices:
                           c = MOD(iRegularParam(l,k), 10000)
                           s = iRegularParam(l,k) - c
                           s = s / 10000

                           ! Store the first and second site fractions:
                           if (k-1 == iSUBIParamData(l,1)) then
                               dFirstParam = dSiteFraction(iSPI,s,c)
                               ! Compute prefactor term:
                               dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                           else if (k == iSUBIParamData(l,2)+1) then
                               dSecondParam = dSiteFraction(iSPI,s,c)
                               ! Compute prefactor term:
                               dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                           ! Determine if the last mixing constituent is a vacancy or anion
                           else if ((cConstituentNameSUB(iSPI,s,c) == 'Va') .AND. &
                               (k == n + 1)) then
                               ! Defines the mixing case type: L_Ci,Cj:Va
                               iMixType(l) = 2

                               dThirdParam = dSiteFraction(iSPI,s,c)
                               ! Compute prefactor term:
                               dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)**2

                           else if ((dSublatticeCharge(iSPI,s,c) < 0D0) .AND. &
                                   (cConstituentNameSUB(iSPI,s,c) /= 'Va') .AND. &
                                   (k == n + 1)) then
                              ! Defines the mixing case type: L_Ci,Cj:Ak
                              iMixType(l) = 3
                              ! Compute prefactor term:
                              dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                            else

                              print *, 'Unrecognized excess mixing term in SUBI phase ', cSolnPhaseName(iSolnIndex)
                              INFOThermo = 36

                           end if
                   end do

                   ! Multiply prefactor term by excess Gibbs energy parameter:
                   iExponent = iRegularParam(l,n+2)

                   ! Excess Gibbs energy equation for L_Ci,Cj:Va case
                   if (iMixType(l) == 2) then

                     gexcess = gexcess + q * dPreFactor * dExcessGibbsParam(l) * &
                              (dFirstParam - dSecondParam)**(iExponent)

                   ! Excess Gibbs energy equation for L_Ci,Cj:Ak case
                   else if (iMixType(l) == 3) then

                     gexcess = gexcess + dPreFactor * dExcessGibbsParam(l) * &
                              (dFirstParam - dSecondParam)**(iExponent)

                   end if

          else

            print *, 'Unrecognized excess mixing term in SUBI phase ', cSolnPhaseName(iSolnIndex)
            INFOThermo = 36
            return

          end if

!          print*,"iMixType(l)",iMixType(l)
!          print*,"-------------------------"
!          print*,""
        end do LOOP_Param





        ! Derived Excess Mixing Terms
        ! ---------------------------

        ! Loop through parameters:
        LOOP_Param_Derivatives: do l = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)
          ! Return control to the parent subroutine if there aren't any interaction parameters for this phase:
          if (nParamPhase(iSolnIndex) - nParamPhase(iSolnIndex-1) == 0) exit LOOP_Param_Derivatives

          ! Reinitialize temporary variable:
            dPreFactor = 1D0
            iFirstParam = 0
            iSecondParam = 0
            iFirstParam2 = 0
            iSecondParam2 = 0

            y1 = 0D0
            y2 = 0D0
            y3 = 0D0

            iMixTypeVa = 0
            iMixTypeBi = 0
            iMixTypeAni = 0
            iMixTypeCati = 0
            iMixTypeCatj = 0

            cDummyAi = ' '
            cDummyAj = ' '
            cDummyBi = ' '
            cDummyCi = ' '
            cDummyCj = ' '



          ! Store the number of constituents involved in this parameter:
          n = iRegularParam(l,1)

          ! Determine the mixing parameter type
          ! This is Case -> L_Ci:Aj,Dk where Dk is a vacancy, neutral or anion
          if (iMixType(l) == 1) then
              !(iSUBIParamData(l,1) == 1) .AND. &
                !  (iSUBIParamData(l,2) == 3) .AND. &
                !  (iSUBIParamData(l,3) == 2)
            ! Loop through constituents associated with this parameter:
            do k = 2, n + 1

              ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Establishing the site fractions for this scenario
                    if (cConstituentNameSUB(iSPI,s,c) == 'Va') then
                        ! vacancy
                        y3 = dSiteFraction(iSPI,s,c)
                        iMixTypeVa = 1
                    else if (dSublatticeCharge(iSPI,s,c) == 0D0) then
                        ! neutral
                        y3 = dSiteFraction(iSPI,s,c)
                        iMixTypeBi = 1
                        cDummyBi = cConstituentNameSUB(iSPI,s,c)
                    else if ((dSublatticeCharge(iSPI,s,c) < 0D0) .AND. &
                             (cConstituentNameSUB(iSPI,s,c) /= 'Va')) then
                      if (y2 /= 0D0) then
                        y3 = dSiteFraction(iSPI,s,c)
                        iMixTypeAni = 1
                        cDummyAj = cConstituentNameSUB(iSPI,s,c)
                      end if
                        ! anion
                        y2 = dSiteFraction(iSPI,s,c)
                        cDummyAi = cConstituentNameSUB(iSPI,s,c)
                    else
                        ! cation
                        y1 = dSiteFraction(iSPI,s,c)
                        cDummyCi = cConstituentNameSUB(iSPI,s,c)
                    end if
                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iSPI,s,c)

                    ! Store the first and second site fractions:
                    if (k == iSUBIParamData(l,2)) then
                        dFirstParam = dSiteFraction(iSPI,s,c)
                    else if (k == iSUBIParamData(l,2) + 1) then
                        dSecondParam = dSiteFraction(iSPI,s,c)
                    end if
                    iExponent = iRegularParam(l,n+2)
            end do

                    ! Chemical potential with respect to the first sublattice
                    do i = 1, nConstituentSublattice(iSPI,1)
                      if (cDummyCi == cConstituentNameSUB(iSPI,1,i)) then
!print*,"cConstituentNameSUB(iSPI,1,i)  ",cConstituentNameSUB(iSPI,1,i)
!print*,"dgdc1(i)",dgdc1(i)*dIdealConstant * dTemperature
                        dgdc1(i) = dgdc1(i) + y2 * y3 * dExcessGibbsParam(l) * (dFirstParam - dSecondParam)**(iExponent)
!print*,"dgdc1(i)",dgdc1(i)*dIdealConstant * dTemperature
!print*,""
                      end if
                    end do

                    ! Chemical potential with respect to the second sublattice
                    do i = 1, nConstituentSublattice(iSPI,2)
                      do j = iFirst, iLast
                         ! Relative species index:
                         n = j - iFirst + 1

!print*,""
!print*,"i", i,"    n",n,"    j",j

                         ! Store constituent indices:
                         l1 = iConstituentSublattice(iSPI,1,n)
                         l2 = iConstituentSublattice(iSPI,2,n)
!print*,"l1",l1,"     l2",l2
!print*,"cConstituentNameSUB(iSPI,2,l2)",cConstituentNameSUB(iSPI,2,l2)
!print*,"cConstituentNameSUB(iSPI,1,n)",cConstituentNameSUB(iSPI,1,n)

                      if ((cConstituentNameSUB(iSPI,2,l2) == 'Va') .AND. &
                          (cConstituentNameSUB(iSPI,1,n) == cDummyCi) .AND. &
                            (iMixTypeVa == 1) .AND. &
                            (i == l2)) then
                            ! cation / vacancy
!print*,"Va 1 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                              dgdc2(i) = dgdc2(i) + y1 * y2 * dExcessGibbsParam(l) * (dFirstParam - dSecondParam)**(iExponent)
!print*,"Va 2 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                              dgdc2(i) = dgdc2(i) + y1 * y2 * y3 * dExcessGibbsParam(l) * &
                                        iExponent * (dFirstParam - dSecondParam)**(iExponent-1) * (-1)
!print*,"Va 3 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature


                        else if ((dSublatticeCharge(iSPI,2,l2) == 0D0) .AND. &
                                (iMixTypeBi == 1) .AND. &
                                (cConstituentNameSUB(iSPI,2,l2) == cDummyBi) .AND. &
                                (i == l2))then
!print*,"Bi 1 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                            ! neutral
                              dgdc2(i) = dgdc2(i) + y1 * y2 * dExcessGibbsParam(l) * (dFirstParam - dSecondParam)**(iExponent)
!print*,"Bi 2 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                              dgdc2(i) = dgdc2(i) + y1 * y2 * y3 * dExcessGibbsParam(l) * &
                                         iExponent * (dFirstParam - dSecondParam)**(iExponent-1) * (-1)
!print*,"Bi 3 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature

                        else if ((dSublatticeCharge(iSPI,2,l2) < 0D0) .AND. &
                                (cConstituentNameSUB(iSPI,2,l2) /= 'Va') .AND. &
                                (i == l2)) then
!print*,"An 1 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                              ! cation / anion
                              if (cConstituentNameSUB(iSPI,1,l1) == cDummyCi) then

                               dgdc2(i) = dgdc2(i) + y1 * y3 * dExcessGibbsParam(l) * (dFirstParam - dSecondParam)**(iExponent)
!print*,"An 2 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                               dgdc2(i) = dgdc2(i) + y1 * y2 * y3 * dExcessGibbsParam(l) * &
                                        iExponent * (dFirstParam - dSecondParam)**(iExponent-1) * (1)
!print*,"An 3 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                                 if ((iMixTypeAni == 1)) then
                                   ! If D_k is an anion
                                   dgdc2(i) = dgdc2(i) + y1 * y2 * dExcessGibbsParam(l) * (dFirstParam - dSecondParam)**(iExponent)
                                   dgdc2(i) = dgdc2(i) + y1 * y2 * y3 * dExcessGibbsParam(l) * &
                                            iExponent * (dFirstParam - dSecondParam)**(iExponent-1) * (-1)
                                 end if
                               end if
                        end if
                      end do
                  end do

          ! Determine the mixing parameter type
          ! This is Case -> L_Ci,Cj:Va
          else if (iMixType(l) == 2) then

            ! Loop through constituents associated with this parameter:
            do k = 2, n + 1

              ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Establishing the site fractions for this case
                    if (k-1 == iSUBIParamData(l,1)) then
                        dFirstParam = dSiteFraction(iSPI,s,c)
                        ! Cation - Ci
                        y1 = dSiteFraction(iSPI,s,c)
                        cDummyCi = cConstituentNameSUB(iSPI,s,c)

                    else if (k == iSUBIParamData(l,2) + 1) then
                        dSecondParam = dSiteFraction(iSPI,s,c)
                        ! Cation - Cj
                        y2 = dSiteFraction(iSPI,s,c)
                        cDummyCj = cConstituentNameSUB(iSPI,s,c)

                    else if (cConstituentNameSUB(iSPI,s,c) == 'Va') then
                      ! vacancy
                      y3 = dSiteFraction(iSPI,s,c)

                    end if
                    iExponent = iRegularParam(l,n+2)
            end do

            ! First sublattice - Cation:vacancy contributions
            do i = 1, nConstituentSublattice(iSPI,1)
              do j = iFirst, iLast

                ! Relative species index
                n = j - iFirst + 1

                ! Store constituent indices:
                l1 = iConstituentSublattice(iSPI,1,n)
                l2 = iConstituentSublattice(iSPI,2,n)

                ! Derivative with respect to Ci
                if ((i == l1) .AND. &
                    (cDummyCi == cConstituentNameSUB(iSPI,1,n)) .AND. &
                    (cConstituentNameSUB(iSPI,2,l2) == 'Va')) then

                    dgdc1(i) = dgdc1(i) + q * y2 * y3**2 * dExcessGibbsParam(l) * &
                              (dFirstParam - dSecondParam)**(iExponent)

                    dgdc1(i) = dgdc1(i) + dSublatticeCharge(iSPI,1,1) * y1 * y2 * y3**2 * dExcessGibbsParam(l) * &
                              (dFirstParam - dSecondParam)**(iExponent)

                    dgdc1(i) = dgdc1(i) + q * y1 * y2 * y3**2 * dExcessGibbsParam(l) * &
                               iExponent * (dFirstParam - dSecondParam)**(iExponent - 1) * (1)

                ! Derivative with respect to Cj
                else if ((i == l1) .AND. &
                    (cDummyCj == cConstituentNameSUB(iSPI,1,n)) .AND. &
                    (cConstituentNameSUB(iSPI,2,l2) == 'Va')) then

                    dgdc1(i) = dgdc1(i) + q * y1 * y3**2 * dExcessGibbsParam(l) * &
                              (dFirstParam - dSecondParam)**(iExponent)

                    dgdc1(i) = dgdc1(i) + dSublatticeCharge(iSPI,1,2) * y1 * y2 * y3**2 * dExcessGibbsParam(l) * &
                              (dFirstParam - dSecondParam)**(iExponent)

                    dgdc1(i) = dgdc1(i) + q * y1 * y2 * y3**2 * dExcessGibbsParam(l) * &
                               iExponent * (dFirstParam - dSecondParam)**(iExponent - 1) * (-1)

                end if
              end do
            end do

            ! Second sublattice - Cation:vacancy contributions
            do i = 1, nConstituentSublattice(iSPI,2)
                if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                  dgdc2(i) = dgdc2(i) + q * y1 * y2 * 2 * y3 * dExcessGibbsParam(l) * &
                            (dFirstParam - dSecondParam)**(iExponent)
                end if
            end do

          else if (iMixType(l) == 3) then

            ! Loop through constituents associated with this parameter:
            do k = 2, n + 1

              ! Determine constituent and sublattice indices:
                    c = MOD(iRegularParam(l,k), 10000)
                    s = iRegularParam(l,k) - c
                    s = s / 10000

                    ! Establishing the site fractions for this case
                    if (k-1 == iSUBIParamData(l,1)) then
                        dFirstParam = dSiteFraction(iSPI,s,c)
                        ! Cation - Ci
                        y1 = dSiteFraction(iSPI,s,c)
                        cDummyCi = cConstituentNameSUB(iSPI,s,c)

                    else if (k == iSUBIParamData(l,2) + 1) then
                        dSecondParam = dSiteFraction(iSPI,s,c)
                        ! Cation - Cj
                        y2 = dSiteFraction(iSPI,s,c)
                        cDummyCj = cConstituentNameSUB(iSPI,s,c)

                    else if (dSublatticeCharge(iSPI,s,c) < 0D0) then
                      ! Anion
                      y3 = dSiteFraction(iSPI,s,c)

                    end if
                    iExponent = iRegularParam(l,n+2)
            end do

            ! First sublattice - Cation:vacancy contributions
            do i = 1, nConstituentSublattice(iSPI,1)
              do j = iFirst, iLast

                ! Relative species index
                n = j - iFirst + 1

                ! Store constituent indices:
                l1 = iConstituentSublattice(iSPI,1,n)
                l2 = iConstituentSublattice(iSPI,2,n)
!print*,"cDummyCi",cDummyCi
!print*,"cDummyCj",cDummyCj
!print*,"cConstituentNameSUB(iSPI,1,n)",cConstituentNameSUB(iSPI,1,n)
!print*,"i",i,"       n",n
                ! Derivative with respect to Ci
                if ((i == l1) .AND. &
                    (cDummyCi == cConstituentNameSUB(iSPI,1,n)) .AND. &
                    (dSublatticeCharge(iSPI,2,l2) < 0D0)) then
!print*,"Ci 1 - dgdc1(i)",dgdc1(i)*dIdealConstant * dTemperature
                    dgdc1(i) = dgdc1(i) + y2 * y3 * dExcessGibbsParam(l) * &
                              (dFirstParam - dSecondParam)**(iExponent)
!print*,"Ci 2 - dgdc1(i)",dgdc1(i)*dIdealConstant * dTemperature
                    dgdc1(i) = dgdc1(i) + y1 * y2 * y3 * dExcessGibbsParam(l) * &
                               iExponent * (dFirstParam - dSecondParam)**(iExponent - 1) * (1)
!print*,"Ci 3 - dgdc1(i)",dgdc1(i)*dIdealConstant * dTemperature
                ! Derivative with respect to Cj
                else if ((i == l1) .AND. &
                    (cDummyCj == cConstituentNameSUB(iSPI,1,n)) .AND. &
                    (dSublatticeCharge(iSPI,2,l2) < 0D0)) then
!print*,"Cj 1 - dgdc1(i)",dgdc1(i)*dIdealConstant * dTemperature
                    dgdc1(i) = dgdc1(i) + y1 * y3 * dExcessGibbsParam(l) * &
                              (dFirstParam - dSecondParam)**(iExponent)
!print*,"Cj 2 - dgdc1(i)",dgdc1(i)*dIdealConstant * dTemperature
                    dgdc1(i) = dgdc1(i) + y1 * y2 * y3 * dExcessGibbsParam(l) * &
                               iExponent * (dFirstParam - dSecondParam)**(iExponent - 1) * (-1)
!print*,"Cj 3 - dgdc1(i)",dgdc1(i)*dIdealConstant * dTemperature
                end if
              end do
            end do

            ! Second sublattice - Cation:vacancy contributions
            do i = 1, nConstituentSublattice(iSPI,2)
                if ((dSublatticeCharge(iSPI,2,i) < 0D0) .AND. &
                    (cConstituentNameSUB(iSPI,2,i) /= 'Va')) then
!print*,"Ak 1 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                  dgdc2(i) = dgdc2(i) + y1 * y2 * dExcessGibbsParam(l) * &
                            (dFirstParam - dSecondParam)**(iExponent)
!print*,"Ak 2 - dgdc2(i)",dgdc2(i)*dIdealConstant * dTemperature
                end if
            end do

          else

            print *, 'Unrecognized excess mixing term in SUBI phase ', cSolnPhaseName(iSolnIndex)
            INFOThermo = 36
            return

          end if

!print*,"y1",y1
!print*,"y2",y2
!print*,"y3",y3
!print*,"dFirstParam",dFirstParam
!print*,"dSecondParam",dSecondParam
!print*,"iMixTypeVa",iMixTypeVa
!print*,"iMixTypeBi",iMixTypeBi
!print*,"iMixTypeAni",iMixTypeAni
!print*,"iMixTypeCati",iMixTypeCati
!print*,"iMixTypeCatj",iMixTypeCatj
!print*,""

!print*,"nParamPhase(iSolnIndex)",nParamPhase(iSolnIndex)

end do LOOP_Param_Derivatives

        ! REFERENCE GIBBS ENERGY AND IDEAL MIXING
        ! ---------------------------------------
        gref = 0D0
        do j = iFirst, iLast
            ! Relative species index:
            n = j - iFirst + 1

            ! Store constituent indices:
            l1 = iConstituentSublattice(iSPI,1,n)
            l2 = iConstituentSublattice(iSPI,2,n)

            lc1 = 1D0
            if ((cConstituentNameSUB(iSPI,2,l2) == 'Va') .OR. &
                (dSublatticeCharge(iSPI,2,l2) == 0D0)) then
                lc2 = 1D0
            else
                if (l1 > 0) lc1 = dSublatticeCharge(iSPI,1,l1)
                lc2 = -dSublatticeCharge(iSPI,2,l2)
            end if

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

print*,"gref+gideal:", (gref+gideal)*dIdealConstant * dTemperature/p*2
print*,"gref+gideal+gexcess:", (gref+gideal+gexcess)*dIdealConstant * dTemperature/p*2
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

                if ((cConstituentNameSUB(iSPI,2,l2) == 'Va') .OR. &
                    (dSublatticeCharge(iSPI,2,l2) == 0D0)) then
                    lc2 = 1D0
                else
                    lc2 = -dSublatticeCharge(iSPI,2,l2)
                end if

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
            do j = 1, nConstituentSublattice(iSPI,1)
                dgdc1(i) = dgdc1(i) + dSublatticeCharge(iSPI,1,i) * yva &
                                    * dSiteFraction(iSPI,1,j)*DLOG(dSiteFraction(iSPI,1,j))
            end do
            do j = 1, nConstituentSublattice(iSPI,2)
                dgdc1(i) = dgdc1(i) + dSublatticeCharge(iSPI,1,i) * dSiteFraction(iSPI,2,j) * DLOG(dSiteFraction(iSPI,2,j))
                if (cConstituentNameSUB(iSPI,2,j) == 'Va') then
                    ! cation / vacancy
                    dgdc1(i) = dgdc1(i) + (1 + DLOG(dSiteFraction(iSPI,1,i))) * q * dSiteFraction(iSPI,2,j)
                else if (dSublatticeCharge(iSPI,2,j) == 0D0) then
                    ! neutral
                else
                    ! cation / anion
                    dgdc1(i) = dgdc1(i) + (1+DLOG(dSiteFraction(iSPI,1,i)))*(-dSublatticeCharge(iSPI,2,j))*dSiteFraction(iSPI,2,j)
                end if
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

                lc1 = 1D0
                if ((cConstituentNameSUB(iSPI,2,l2) == 'Va') .OR. &
                    (dSublatticeCharge(iSPI,2,l2) == 0D0)) then
                    lc2 = 1D0
                else
                    if (l1 > 0) lc1 = dSublatticeCharge(iSPI,1,l1)
                    lc2 = -dSublatticeCharge(iSPI,2,l2)
                end if

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
!print*,"i",i
!print*,"l1: ",l1,"  l2:",l2
!print*,"dgdc2(i)",dgdc2(i)
!print*
            end do
            ! Entropy
            dgdc2(i) = dgdc2(i) + (1 + DLOG(dSiteFraction(iSPI,2,i))) * q
            do j = 1, nConstituentSublattice(iSPI,1)
                if (cConstituentNameSUB(iSPI,2,i) == 'Va') then
                    ! cation / vacancy
                    dgdc2(i) = dgdc2(i) + q * dSiteFraction(iSPI,1,j) * DLOG(dSiteFraction(iSPI,1,j))
                else if (dSublatticeCharge(iSPI,2,j) == 0D0) then
                    ! neutral
                else
                    ! cation / anion
                    dgdc2(i) = dgdc2(i) + (-dSublatticeCharge(iSPI,2,i)) * dSiteFraction(iSPI,1,j) * DLOG(dSiteFraction(iSPI,1,j))
                end if
            end do
        end do
print*,"After All"
print*,"dgdc1(:)",dgdc1(:)*dIdealConstant * dTemperature
print*,"dgdc2(:)",dgdc2(:)*dIdealConstant * dTemperature
print*,""
        ! Compute the chemical potential for each phase component assuming ideal mixing:
        LOOP_Ideal: do i = iFirst, iLast
            ! Relative species index:
            m = i - iFirst + 1
            k1 = iConstituentSublattice(iSPI,1,m)
            k2 = iConstituentSublattice(iSPI,2,m)

            kc1 = 1D0
            ! cation / vacancy
            if (cConstituentNameSUB(iSPI,2,k2) == 'Va') then
                kc2 = 1D0
                natom = 1D0 / dMol + dMolDerivatives(m) * dMolAtoms
            ! neutral
            else if (dSublatticeCharge(iSPI,2,k2) == 0D0) then
                kc2 = 1D0
                natom = 1D0 / dMol + dMolDerivatives(m) * dMolAtoms
            ! cation / anion
            else
                if (k1 > 0) kc1 = dSublatticeCharge(iSPI,1,k1)
                kc2 = -dSublatticeCharge(iSPI,2,k2)
                natom = (kc1 + kc2) / dMol + dMolDerivatives(m) * dMolAtoms
            end if

            dChemicalPotential(i) = (gref + gideal + gexcess) * natom

            do j = 1, nConstituentSublattice(iSPI,1)
                ! cation / (anion or vacancy)
                if (k1 > 0) then
                    dydn = -kc2 * dSiteFraction(iSPI,1,j) / dSub1Total
                    if (j == k1) dydn = dydn + kc2 / dSub1Total
                    dChemicalPotential(i) = dChemicalPotential(i) + dydn * dgdc1(j) * dMolAtoms / dMol
                end if
            end do

            do j = 1, nConstituentSublattice(iSPI,2)
                dydn = -kc1 * dSiteFraction(iSPI,2,j) / dSub2Total
                if (j == k2) dydn = dydn + kc1 / dSub2Total
                dChemicalPotential(i) = dChemicalPotential(i) + dydn * dgdc2(j) * dMolAtoms / dMol
            end do
        end do LOOP_Ideal

        deallocate(dgdc1,dgdc2)
        deallocate(dMolDerivatives)
        deallocate(iMixType)
    end if IF_SUBL

    return

end subroutine CompExcessGibbsEnergySUBI
