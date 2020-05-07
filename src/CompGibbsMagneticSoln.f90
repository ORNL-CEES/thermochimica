
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompGibbsMagneticSoln.f90
    !> \brief   Compute magnetic contributions to the Gibbs energy terms for a solution phase.
    !> \author  M.H.A. Piro
    !> \date    March 7, 2013
    !> \sa      CompExcessGibbsEnergy.f90
    !
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    03/07/2013    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the magnetic contribution to the standard molar
    !! Gibbs energy of a solution species.  This contribution is given by
    !! \f$ \Delta g_{mag} = RT ln(B_o + 1) g(\tau ) \f$, where \f$ B_o \f$ is the average magnetic moment per
    !! atom, \f$ \tau \f$  is the absolute temperature divided by the critical temperature (i.e., the Curie
    !! temperature for ferromagnetic materials or the Neel temperature for antiferromagnetic materials) and
    !! \f$ g \f$ is a function of \f$ \tau \f$, given by:
    !! \f$ g(\tau ) = 1 - \left( \frac{79\tau^-1}{140p} + \frac{474}{497}(\frac{1}{p} - 1)(\frac{\tau ^3}{6}
    !!   + \frac{\tau ^9}{135} \frac{\tau ^15}{600}   )  \right) /D, \tau \leq 1 \f$
    !! and \f$ g(\tau ) = - \left( \frac{\tau ^{-5}}{10} + \frac{\tau ^{-15}}{315} + \frac{\tau ^{-25}}{1500}
    !!   \right) /D, \tau > 1 \f$, where \f$ D = \frac{518}{1125} + \frac{11692}{15975} \left( \frac{1}{p} -1
    !!   \right) \f$.
    !!
    !
    ! References:
    ! ===========
    !
    !> \details The following references explain the magnetic contribution to the Gibbs energy term that is
    !! used in this subroutine:
    !!
    !!   M. Hillert and M. Jarl, "A Model for Alloying Effects in Ferromagnetic Alloys," CALPHAD, 2, 3
    !!   (1978) 227-238.
    !!
    !!   A.T. Dinsdale, "SGTE Data for Pure Elements," CALPHAD, 15, 4 (1991) 317-425.
    !!
    !!   H.L. Lukas, S.G. Fries and B. Sundman, "Computational Thermodynamics: The Calphad Method," Cambridge
    !!   University Press, New York (2007).
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iSolnPhaseIndex   An integer scalar representing the absolute solution phase index.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompGibbsMagneticSoln(iSolnPhaseIndex)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer :: i, iFirst, iLast, iSolnPhaseIndex, iParam, iExponent, n, c, s, iChargedPhaseID, k
    real(8) :: B, D, p, invpmone, tau, Tcritical, g, StructureFactor
    real(8) :: dTemp, dTempA, dTempB, dTempC, dTempD
    real(8) :: x1, x2, xprod, dx, dPreFactor


    ! Initialize variables:
    Tcritical = 0D0
    B         = 0D0

    ! Store the first and last species indices:
    iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnPhaseIndex)

    ! Note: the structure factor and p factor are the same for all constituents in the phase.
    StructureFactor = dCoeffGibbsMagnetic(iFirst,3)
    p               = dCoeffGibbsMagnetic(iFirst,4)
    invpmone        = 1D0/p - 1D0

    ! The critical temperature can either be the Curie temperature for ferromagnetic materials or the Neel
    ! temperature for antiferromagnetic materials.  For a solution phase, this is a linear function of the
    ! mole fractions of solution phase constituents:
    do i = iFirst, iLast
        Tcritical = Tcritical + dMolFraction(i) * dCoeffGibbsMagnetic(i,1)
        B         = B + dMolFraction(i) * dCoeffGibbsMagnetic(i,2)
    end do
    ! Tcritical = Tcritical + dMolFraction(iFirst) * dMolFraction(iLast) * (-3605.0D0)
    ! B         = B + dMolFraction(iFirst) * dMolFraction(iLast) * (-1.91D0)

    LOOP_Param: do iParam = nMagParamPhase(iSolnPhaseIndex-1)+1, nMagParamPhase(iSolnPhaseIndex)
        if (cSolnPhaseType(iSolnPhaseIndex) == 'RKMPM') then
            ! Compute temporary variables for sake of convenience:
            if (iMagneticParam(iParam,1) == 2) then
                ! Binary parameter:
                x1    = dMolFraction(iFirst + iMagneticParam(iParam,2) - 1)
                x2    = dMolFraction(iFirst + iMagneticParam(iParam,3) - 1)
                iExponent = iMagneticParam(iParam,4)
                xprod = x1 * x2
                dx    = x1 - x2
                ! Cycle if dx = 0 to prevent calculating either an INF or a NAN:
                if (dx == 0D0) cycle LOOP_Param
                Tcritical = Tcritical + xprod * dx**iExponent * dMagneticParam(iParam,1)
                B         = B         + xprod * dx**iExponent * dMagneticParam(iParam,2)
            else
                ! The parameter index is not supported/recognized.  Report an error and exit.
                INFOThermo = 43
                exit LOOP_Param
            end if
        else if (cSolnPhaseType(iSolnPhaseIndex) == 'SUBLM') then
            iChargedPhaseID = iPhaseSublattice(iSolnPhaseIndex)
            ! Reinitialize temporary variable:
            dPreFactor = 1D0

            ! Store the number of constituents involved in this parameter:
            n = iMagneticParam(iParam,1)

            ! Loop through constituents associated with this parameter:
            do k = 2, n + 1

                ! Determine constituent and sublattice indices:
                c = MOD(iMagneticParam(iParam,k), 10000)
                s = iMagneticParam(iParam,k) - c
                s = s / 10000

                ! Compute prefactor term:
                dPreFactor = dPreFactor * dSiteFraction(iChargedPhaseID,s,c)

                ! Store the first and second site fractions:
                ! This assumes that the constituents that are mixing are the first two listed:
                if (k == 2) then
                    x1 = dSiteFraction(iChargedPhaseID,s,c)
                elseif (k == 3) then
                    x2 = dSiteFraction(iChargedPhaseID,s,c)
                end if
            end do

            ! Multiply prefactor term by excess Gibbs energy parameter:
            dPreFactor = dPreFactor * (x1 - x2)**(iMagneticParam(iParam,n+2))
            Tcritical = Tcritical + dPreFactor * dMagneticParam(iParam,1)
            B         = B         + dPreFactor * dMagneticParam(iParam,2)
        end if

    end do LOOP_Param

    ! print *, cSolnPhaseName(iSolnPhaseIndex), Tcritical, B

    ! ChemSage files store the critical temperature for antiferromagnetic materials
    ! (i.e., the Neel temperature) as a negative real value divided by the structure factor.
    if (Tcritical < 0D0) then
        Tcritical = -Tcritical * StructureFactor
        B         = -B * StructureFactor
    end if
    ! print *, cSolnPhaseName(iSolnPhaseIndex), Tcritical, B

    ! Only proceed if the critical temperature is not zero:
    IF_Proceed: if (Tcritical /= 0D0) then

        ! Compute model parameters:
        tau = dTemperature / Tcritical
        D   = (518D0/1125D0) + (11692D0/15975D0) * invpmone

        ! The magnetic model of Hillert and Jarl is empirical and depends on tau:
        IF_Tau: if (tau > 1D0) then
            dTempA = tau**(-5)                  ! tau^(-5)
            dTempB = dTempA**(3)                ! tau^(-15)
            dTempC = dTempA * dTempA * dTempB   ! tau^(-25)
            dTempD = (1D0 / (D * Tcritical)) * (dTempA / 2D0 + dTempB / 21D0 + dTempC / 60D0)
            g      = -(dTempA/10D0 + dTempB/315D0 + dTempC/1500D0) / D

            ! Loop through species in this phase and update the chemical potential:
            do i = iFirst, iLast
                dTemp = (Tcritical - dCoeffGibbsMagnetic(i,1)) * dTempD
                dTemp = g * ((dCoeffGibbsMagnetic(i,2) - B) / (1D0 + B)) + DLOG(1D0 + B) * (dTemp + g)
                dMagGibbsEnergy(i) = dTemp
            end do
        else
            dTempA = tau**(3)                   ! tau^(3)
            dTempB = dTempA**(3)                ! tau^(9)
            dTempC = dTempA * dTempA * dTempB   ! tau^(15)
            dTempD = (-79D0 / (140D0 * p * dTemperature))
            dTempD = dTempD + (474D0/(497D0*Tcritical)) * invpmone * (dTempA / 2D0 + dTempB / 15D0 + dTempC / 40D0)
            dTempD = dTempD / D
            g      = 1D0 - (79D0/(140D0*p*tau) + (474D0/497D0)*invpmone*(dTempA/6D0 + dTempB/135D0 + dTempC/600D0)) / D

            ! Loop through species in this phase and update the chemical potential:
            do i = iFirst, iLast
                dTemp = (dCoeffGibbsMagnetic(i,1) - Tcritical) * dTempD
                dTemp = g * ((dCoeffGibbsMagnetic(i,2) - B) / (1D0 + B)) + DLOG(1D0 + B) * (dTemp + g)
                ! print *, dTemp
                dMagGibbsEnergy(i) = dTemp
            end do

        end if IF_Tau

        ! if (cSolnPhaseName(iSolnPhaseIndex) == 'FCC_A1') print *, g, Tcritical, dTempD, B

    end if IF_Proceed

    return

end subroutine CompGibbsMagneticSoln
