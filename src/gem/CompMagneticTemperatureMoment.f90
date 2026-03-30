subroutine CompMagneticTemperatureMoment(iSolnPhaseIndex,Tcritical,B)

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer, intent(in)  :: iSolnPhaseIndex
    real(8), intent(out) :: B, Tcritical
    integer :: i, iFirst, iLast, iParam, iExponent, n, c, s, iChargedPhaseID, k
    real(8) :: x1, x2, xprod, dx, dPreFactor


    ! Initialize variables:
    Tcritical = 0D0
    B         = 0D0

    ! Store the first and last species indices:
    iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnPhaseIndex)

    ! The critical temperature can either be the Curie temperature for ferromagnetic materials or the Neel
    ! temperature for antiferromagnetic materials.  For a solution phase, this is a linear function of the
    ! mole fractions of solution phase constituents:
    do i = iFirst, iLast
        Tcritical = Tcritical + dMolFraction(i) * dCoeffGibbsMagnetic(i,1)
        B         = B + dMolFraction(i) * dCoeffGibbsMagnetic(i,2)
    end do

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
            x1 = 0D0
            x2 = 0D0
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

end subroutine CompMagneticTemperatureMoment
