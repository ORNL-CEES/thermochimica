
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompGibbsMagnetic.f90
    !> \brief   Compute magnetic contributions to the Gibbs energy terms.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      CompThermoData.f90
    !> \param[in] i Gibbs energy equation coefficient
    !> \param[in] j Species index
    !> \param[in,out] dChemicalPotential A double real vector representing the chemical potential.
    ! 
    !
    !
    ! Revisions:
    ! ==========
    ! 
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    10/20/2011    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the magnetic contribution to the standard molar 
    !! Gibbs energy of a pure species.  This contribution is given by 
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
    !!   A.T. Dinsdale, "SGTE Data for Pure Elements," CALPHAD, 15, 4 (1991) 317-425.
    !!
    !!   H.L. Lukas, S.G. Fries and B. Sundman, "Computational Thermodynamics: The Calphad Method," Cambridge
    !!   University Press, New York (2007).
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   i   An integer scalar corresponding to the Gibbs energy index.
    !> \param[in]   j   An integer scalar representing the species index.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine CompGibbsMagnetic(i,j)

    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none
    
    integer :: i, j
    real(8) :: B, D, p, invpmone, tau, Tcritical, g, StructureFactor
    
    
    ! Assign the critical temperature can either be the Curie temperature for ferromagnetic materials
    ! or the Neel temperature for antiferromagnetic materials:
    Tcritical       = dGibbsMagneticCS(i,1)
    B               = dGibbsMagneticCS(i,2)
    StructureFactor = dGibbsMagneticCS(i,3)
    p               = dGibbsMagneticCS(i,4)
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
    dChemicalPotential(j) = dChemicalPotential(j) + dIdealConstant * dTemperature * DLOG(B + 1D0) * g

    return

end subroutine CompGibbsMagnetic
