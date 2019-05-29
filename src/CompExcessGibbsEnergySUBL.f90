
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergySUBL.f90
    !> \brief   Compute the partial molar excess Gibbs energy of mixing of solution phase constituents in a SUBL
    !!           or SUBLM solution phase.
    !> \author  M.H.A. Piro
    !> \date    January 17, 2013
    !> \sa      CompExcessGibbsEnergy.f90
    !> \sa      CompExcessGibbsEnergyRKMP.f90
    !> \sa      CompExcessGibbsEnergyQKTO.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/17/2013      M.H.A. Piro         Original code.
    !   02/11/2013      M.H.A. Piro         Added capability of handling non-ideal terms.
    !   02/14/2013      M.H.A. Piro         Corrected handling of non-ideal terms with higher order
    !                                        mixing parameters (happy valentine's day).
    !   02/15/2013      M.H.A. Piro         Fix bug in handling higher order terms from yesterday.
    !   02/27/2013      M.H.A. Piro         The sum of site fractions may not necessarily equal unity. Apply
    !                                        an appropriate correction.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the partial molar excess Gibbs energy of mixing
    !! (dPartialExcessGibbs) of all constituents in a non-ideal solution phase designated as 'SUBL'
    !! (Compound Energy Formalism) or 'SUBLM' (SUBL-Magnetic) phases.  For more information on the
    !! SUBL model and the derivation of equations, the reader is referred to the following literature:
    !!
    !!      B. Sundman and J. Agren, "A Regular Solution Model for Phases with Several Components
    !!      and Sublattices, Suitable for Computer Applications," Journal of Physical Chemistry of
    !!      Solids, 42 (1981) 297-301.
    !!
    !!      M. Hillert, "The Compound Energy Formalism," Journal of Alloys and Compounds, 320 (2001)
    !!      161-176.
    !!
    !! An important distinction between a phase represented by the SUBL model and any other model
    !! is that the phase components and the constituents are not the same entity (refer to above
    !! literature for a theoretical discussion).  As a result of this, the contribution to the chemical
    !! potential term from the reference molar Gibbs energy is NOT \f$ g_i^{\circ} \f$, but rather a
    !! more complicated expression involved all components in this phase.
    !!
    !! The chemical potential of a component in a SUBL phase is defined by the following equation:
    !!
    !! 	\f$ \mu_{i(\lambda)} = \sum_{j=1}^{N_{\lambda}} x_{j(\lambda)} g_{j(\lambda)}^{\circ}
    !!      \left(1 - N_s + \sum_{s=1}^{N_s} \frac{\delta_{i,j}}{y_{j(s)}} \right)
    !!      + RT \sum_{s=1}^{N_s} a_s \mathrm{ln} (y_{i(s)}) + g_{i(\lambda)}^{ex} \f$
    !!
    !! The partial molar excess Gibbs energy of mixing of a component of a SUBL phase is:
    !!
    !! 	\f$ g_{i(\lambda)}^{ex} = \sum_{p=1}^{N_p} \left( \prod_{m=1} y_{m(s)}  \right)
    !!      \sum_{z=0} {^zL_{j,k}} \left(  1 - (N_s + z) + \sum_{s=1}^{N_s}
    !!      \frac{\delta _{i,p}}{y_{i(s)}}   \right) \f$
    !!
    !! The mole fraction of any component in a SUBL phase is related to the site fractions of the
    !! constituents through the following multiplicative relationship:
    !!
    !! \f$ x_{i(\lambda)} = \prod_{s=1}^{N_s} y_{i(s)} \f$
    !!
    !! Similarly, the site fraction of any constituent is related to the mole fraction through the following
    !! summation:
    !!
    !! \f$ y_{c(s)} = \sum_{i=1}^{N_{\lambda}} x_{i(\lambda)} \delta_{i,c(s)} \f$
    !!
    !! Since the molar Gibbs energy terms are defined by the phase components and NOT the constituents, the
    !! numerics necessarily works with mole fractions (corresponding to components) rather than site fractions
    !! (corresponding to constituents).  Therefore, the mole fractions are modified at any iteration and the
    !! site fractions are computed from the mole fractions.  Note that the computation of the site fractions
    !! is a summation operation and it is thus insensitive to components with relatively small mole fractions.
    !! However, the calculation of a mole fraction from the site fractions is sensitive because it is a
    !! multiplicative function. A correction is applied to the mole fractions of all components after the site
    !! fractions are computed by defining the mole fraction as a function of the site fractions.
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
    ! iChargedPhaseID           An integer scalar representing the relative index of the charged phase.
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


subroutine CompExcessGibbsEnergySUBL(iSolnIndex)

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, l, m, n, s, c, d, KD
    integer :: iSolnIndex, nSublattice, iChargedPhaseID
    integer :: iFirst, iLast, iFirstParam, iSecondParam, iSubParam
    real(8) :: dTemp, dPreFactor, dFirstParam, dSecondParam
    real(8), dimension(nMaxSublatticeSys):: dTempVec


    ! Only proceed if the correct phase type is selected:
    IF_SUBL: if ((cSolnPhaseType(iSolnIndex) == 'SUBL').OR.(cSolnPhaseType(iSolnIndex) == 'SUBLM')) then

        ! Define temporary variables for sake of convenience:
        iChargedPhaseID = iPhaseSublattice(iSolnIndex)
        nSublattice     = nSublatticePhase(iChargedPhaseID)
        iFirst          = nSpeciesPhase(iSolnIndex-1) + 1
        iLast           = nSpeciesPhase(iSolnIndex)

        ! Initialize variables:
        dSiteFraction(iChargedPhaseID,1:nSublattice,1:nMaxConstituentSys) = 0D0
        dChemicalPotential(iFirst:iLast)                                  = 0D0
        dPartialExcessGibbs(iFirst:iLast)                                 = 0D0
        dTempVec                                                          = 0D0

        ! Compute site fractions on each sublattice:
        LOOP_SiteFraction: do i = nSpeciesPhase(iSolnIndex-1) + 1, nSpeciesPhase(iSolnIndex)

            ! Relative component index:
            m = i - iFirst + 1

            ! Loop through sublattices:
            do s = 1, nSublattice
                ! Store constituent index on sublattice s:
                c = iConstituentSublattice(iChargedPhaseID,s,m)
                dSiteFraction(iChargedPhaseID,s,c) = dSiteFraction(iChargedPhaseID,s,c) + dMolFraction(i)
            end do

        end do LOOP_SiteFraction

        ! Compute sum of site fractions on each sublattice:
        do s = 1, nSublattice
            do c = 1, nMaxConstituentSys
                dTempVec(s) = dTempVec(s) + dSiteFraction(iChargedPhaseID,s,c)
            end do
            dTempVec(s) = 1D0 / dTempVec(s)
        end do

        ! Correct the mole fractions of phase components by the site fractions of the constituents
        ! (see top for details):
        LOOP_CorrectX: do i = iFirst, iLast
            dTemp = 1D0
            m = i - iFirst + 1

            ! Loop through sublattices:
            do s = 1, nSublattice
                c     = iConstituentSublattice(iChargedPhaseID,s,m)
                dTemp = dTemp * dSiteFraction(iChargedPhaseID,s,c) * dTempVec(s)
            end do

            dMolFraction(i)  = dTemp

        end do LOOP_CorrectX

        ! Correct the number of moles for each species if the solution phase is stable:
        if (lSolnPhases(iSolnIndex) .EQV. .TRUE.) then
            l = 0
            ! Determine the relative phase index:
            LOOP_C: do j = 1, nSolnPhases
                l = nElements - j + 1
                k = -iAssemblage(l)
                if (k == iSolnIndex) exit LOOP_C
            end do LOOP_C
            ! Correct molar quantities:
            do i = iFirst, iLast
                !dMolesSpecies(i) = dMolesPhase(l) * dMolFraction(i)
            end do
        end if

        ! REFERENCE GIBBS ENERGY AND IDEAL MIXING
        ! ---------------------------------------

        ! Compute the chemical potential for each phase component assuming ideal mixing:
        LOOP_Ideal: do i = iFirst, iLast
            ! Relative species index:
            m = i - iFirst + 1

            ! Loop through components:
            LOOP_Ideal_Components: do j = iFirst, iLast
                ! Relative species index:
                n = j - iFirst + 1

                ! Compute pre-factor term:
                dTemp = 1D0 - DFLOAT(nSublattice)

                ! Loop through sublattices:
                do s = 1, nSublattice
                    ! Store constituent indices:
                    k = iConstituentSublattice(iChargedPhaseID,s,m)
                    l = iConstituentSublattice(iChargedPhaseID,s,n)

                    ! Effectively apply Kronecker-Delta term to pre-factor:
                    if (k == l)  dTemp = dTemp + 1D0 / dSiteFraction(iChargedPhaseID,s,k)
                end do

                ! Update the reference molar Gibbs energy:
                dChemicalPotential(i) = dChemicalPotential(i) + dTemp * dMolFraction(j) * dStdGibbsEnergy(j)

            end do LOOP_Ideal_Components

            ! Add ideal mixing contribution:
            do s = 1, nSublattice
                c = iConstituentSublattice(iChargedPhaseID,s,m)
                dChemicalPotential(i) = dChemicalPotential(i) + dStoichSublattice(iChargedPhaseID,s) &
                    * DLOG(dSiteFraction(iChargedPhaseID,s,c))
            end do

        end do LOOP_Ideal


        ! MAGNETIC TERMS
        ! --------------

        ! Compute magnetic ordering terms if this phase is magnetic:
        if (cSolnPhaseType(iSolnIndex) == 'SUBLM') call CompGibbsMagneticSoln(iSolnIndex)


        ! EXCESS TERMS
        ! ------------

        ! Return control to the parent subroutine if there aren't any interaction parameters for this phase:
        if (nParamPhase(iSolnIndex) - nParamPhase(iSolnIndex-1) == 0) return

        ! Loop through parameters:
        LOOP_Param: do l = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)

            ! Reinitialize temporary variable:
            dPreFactor = 1D0

            ! Store the number of constituents involved in this parameter:
            n = iRegularParam(l,1)

            ! Loop through constituents associated with this parameter:
            do k = 2, n + 1

                ! Determine constituent and sublattice indices:
                c = MOD(iRegularParam(l,k), 10000)
                s = iRegularParam(l,k) - c
                s = s / 10000

                ! Compute prefactor term:
                dPreFactor = dPreFactor * dSiteFraction(iChargedPhaseID,s,c)

                ! Store the first and second site fractions:
                ! This assumes that the constituents that are mixing are the first two listed:
                if (k == 2) then
                    dFirstParam = dSiteFraction(iChargedPhaseID,s,c)
                    iFirstParam = c
                    iSubParam   = s
                elseif (k == 3) then
                    dSecondParam = dSiteFraction(iChargedPhaseID,s,c)
                    iSecondParam = c
                end if
            end do

            ! Multiply prefactor term by excess Gibbs energy parameter:
            dPreFactor = dPreFactor * dExcessGibbsParam(l) * (dFirstParam - dSecondParam)**(iRegularParam(l,n+2))

            ! Loop through species in phase:
            LOOP_Param_Species: do i = iFirst, iLast

                ! Reinitialize variables:
                KD    = 0
                m     = i - iFirst + 1
                dTemp = -DFLOAT(nSublattice + iRegularParam(l,n+2))

                ! Loop through sublattices associated with this phase:
                LOOP_Param_Sub: do s = 1, nSublattice
                    ! Store constituent index corresponding to component i on sublattice s:
                    c = iConstituentSublattice(iChargedPhaseID,s,m)

                    ! Assign Kronecker-Delta term if this component contains the constituent corresponding
                    ! to the mixing parameter:
                    if (s == iSubParam) then
                        if (c == iFirstParam) then
                            ! This is the first mixing constituent:
                            KD = 1
                        elseif (c == iSecondParam) then
                            ! This is the second mixing constituent:
                            KD = -1
                        else
                            ! The constituents don't match:
                            KD = 0
                        end if
                    end if

                    ! Loop through constituents involved in this mixing parameter:
                    LOOP_Param_Const: do j = 2, n + 1
                        k = MOD(iRegularParam(l,j), 10000)  ! Constituent index corresponding to parameter.
                        d = iRegularParam(l,j) - k
                        d = d / 10000                       ! Sublattice index corresponding to parameter.

                        ! Cycle if they are on different sublattices:
                        if (d /= s) cycle LOOP_Param_Const

                        ! Cycle if they are different constituents:
                        if (k /= c) cycle LOOP_Param_Const

                        ! Include contribution from this site fraction:
                        dTemp = dTemp + 1D0 / dSiteFraction(iChargedPhaseID,s,c)

                    end do LOOP_Param_Const ! j
                end do LOOP_Param_Sub       ! s

                ! Apply higher order terms (only if dFirstParam and dSecondParam are not the same):
                if (dFirstParam /= dSecondParam) then
                    dTemp = dTemp + DFLOAT(KD * iRegularParam(l,n+2)) / (dFirstParam - dSecondParam)
                end if

                ! Apply partial molar excess Gibbs energy of mixing:
                dPartialExcessGibbs(i) = dPartialExcessGibbs(i) + dPreFactor * dTemp

            end do LOOP_Param_Species       ! i

        end do LOOP_Param                   ! l

    end if IF_SUBL

    return

end subroutine CompExcessGibbsEnergySUBL
