subroutine CompGibbsMagneticSolnInit(j)

    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer, intent(in) :: j
    real(8) :: B, D, p, invpmone, tau, Tcritical, g, StructureFactor


    ! Assign the critical temperature can either be the Curie temperature for ferromagnetic materials
    ! or the Neel temperature for antiferromagnetic materials:
    Tcritical       = dCoeffGibbsMagnetic(j,1)
    B               = dCoeffGibbsMagnetic(j,2)
    StructureFactor = dCoeffGibbsMagnetic(j,3)
    p               = dCoeffGibbsMagnetic(j,4)
    invpmone        = 1D0/p - 1D0

    ! ChemSage files store the critical (i.e., Neel) temperature for antiferromagnetic materials
    ! as a negative real value divided by the structure factor.  Correct Tcritical and B:
    if (Tcritical < 0D0) then
        Tcritical = -Tcritical * StructureFactor
        B         = -B * StructureFactor
    end if

    tau = dTemperature / Tcritical
    D   = 518D0/1125D0 + (11692D0/15975D0) * invpmone

    if (tau > 1D0) then
        g = -((tau**(-5))/10D0 + (tau**(-15))/315D0 + (tau**(-25))/1500D0) / D
    else
        g = 1D0 - (79D0/(140D0*p*tau) + (474D0/497D0)*invpmone*((tau**3)/6D0 + (tau**9)/135D0 + (tau**15)/600D0)) / D
    end if

    ! Add the magnetic contribution to the Gibbs energy term:
    dMagGibbsEnergy(j) = DLOG(B + 1D0) * g
    ! print *, cSpeciesName(j), dMagGibbsEnergy(j)

    return

end subroutine CompGibbsMagneticSolnInit
