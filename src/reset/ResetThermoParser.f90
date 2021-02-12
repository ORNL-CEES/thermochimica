
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ResetThermoParser.f90
    !> \brief   Deallocate allocatable variables used by the ModuleParseCS.f90 module.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      ModuleParseCS.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer         Description of change
    !    ----          ----------         ---------------------
    !    02/17/2012    M.H.A. Piro        Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to deallocate all allocatable arrays from the
    !! ModuleParseCS.f90 module, which are allocated in the ParseCS*.f90 subroutines and then used in several
    !! subroutines by Thermochimica.  A value of INFOThermo = 18 is returned if an error has occured during
    !! deallocation.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                  An error is returned if deallocation is unsuccessful.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                       it encounters an error.  A description for each error is given in ThermoDebug.f90.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine ResetThermoParser

    USE ModuleParseCS
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer :: i

    ! Initialize variables:
    i    = 0
    INFO = 0

    ! Deallocate integer arrays from ModuleParseCS:
    if (allocated(nSpeciesPhaseCS)) deallocate(nSpeciesPhaseCS,STAT = INFO)
    i = i + INFO
    if (allocated(nGibbsEqSpecies)) deallocate(nGibbsEqSpecies,STAT = INFO)
    i = i + INFO
    if (allocated(iPhaseCS)) deallocate(iPhaseCS,STAT = INFO)
    i = i + INFO
    if (allocated(iParticlesPerMoleCS)) deallocate(iParticlesPerMoleCS,STAT = INFO)
    i = i + INFO
    if (allocated(nParamPhaseCS)) deallocate(nParamPhaseCS,STAT = INFO)
    i = i + INFO
    if (allocated(iParamPassCS)) deallocate(iParamPassCS,STAT = INFO)
    i = i + INFO
    if (allocated(dStoichSpeciesCS)) deallocate(dStoichSpeciesCS,STAT = INFO)
    i = i + INFO
    if (allocated(iRegularParamCS)) deallocate(iRegularParamCS,STAT = INFO)
    i = i + INFO
    if (allocated(cRegularParamCS)) deallocate(cRegularParamCS,STAT = INFO)
    i = i + INFO
    if (allocated(nMagParamPhaseCS)) deallocate(nMagParamPhaseCS,STAT = INFO)
    i = i + INFO
    if (allocated(iMagneticParamCS)) deallocate(iMagneticParamCS,STAT = INFO)
    i = i + INFO
    if (allocated(dMagneticParamCS)) deallocate(dMagneticParamCS,STAT = INFO)
    i = i + INFO
    if (allocated(iMagParamPassCS)) deallocate(iMagParamPassCS,STAT = INFO)
    i = i + INFO
    if (allocated(dAtomicMassCS)) deallocate(dAtomicMassCS,STAT = INFO)
    i = i + INFO
    if (allocated(dGibbsCoeffSpeciesTemp)) deallocate(dGibbsCoeffSpeciesTemp,STAT = INFO)
    i = i + INFO
    if (allocated(dRegularParamCS)) deallocate(dRegularParamCS,STAT = INFO)
    i = i + INFO
    if (allocated(dGibbsMagneticCS)) deallocate(dGibbsMagneticCS,STAT = INFO)
    i = i + INFO
    if (allocated(cElementNameCS)) deallocate(cElementNameCS,STAT = INFO)
    i = i + INFO
    if (allocated(cSolnPhaseTypeCS)) deallocate(cSolnPhaseTypeCS,STAT = INFO)
    i = i + INFO
    if (allocated(cSolnPhaseNameCS)) deallocate(cSolnPhaseNameCS,STAT = INFO)
    i = i + INFO
    if (allocated(cSpeciesNameCS)) deallocate(cSpeciesNameCS,STAT = INFO)
    i = i + INFO
    if (allocated(cConstituentNameSUBCS)) deallocate(cConstituentNameSUBCS,STAT = INFO)
    i = i + INFO
    if (allocated(iPhaseSublatticeCS)) deallocate(iPhaseSublatticeCS,STAT = INFO)
    i = i + INFO
    if (allocated(nConstituentSublatticeCS)) deallocate(nConstituentSublatticeCS,STAT = INFO)
    i = i + INFO
    if (allocated(iConstituentSublatticeCS)) deallocate(iConstituentSublatticeCS,STAT = INFO)
    i = i + INFO
    if (allocated(dStoichSublatticeCS)) deallocate(dStoichSublatticeCS,STAT = INFO)
    i = i + INFO
    if (allocated(nSublatticePhaseCS)) deallocate(nSublatticePhaseCS,STAT = INFO)
    i = i + INFO
    if (allocated(dZetaSpeciesCS)) deallocate(dZetaSpeciesCS,STAT = INFO)
    i = i + INFO
    if (allocated(dConstituentCoefficientsCS)) deallocate(dConstituentCoefficientsCS,STAT = INFO)
    i = i + INFO
    if (allocated(dSublatticeChargeCS)) deallocate(dSublatticeChargeCS,STAT = INFO)
    i = i + INFO
    if (allocated(dStoichPairsCS)) deallocate(dStoichPairsCS,STAT = INFO)
    i = i + INFO
    if (allocated(iChemicalGroupCS)) deallocate(iChemicalGroupCS,STAT = INFO)
    i = i + INFO
    if (allocated(cPairNameCS)) deallocate(cPairNameCS,STAT = INFO)
    i = i + INFO
    if (allocated(dGibbsMagneticCS)) deallocate(dGibbsMagneticCS, STAT = INFO)
    i = i + INFO
    if (allocated(nPairsSROCS)) deallocate(nPairsSROCS, STAT = INFO)
    i = i + INFO
    if (allocated(iPairIDCS)) deallocate(iPairIDCS, STAT = INFO)
    i = i + INFO
    if (allocated(dCoordinationNumberCS)) deallocate(dCoordinationNumberCS, STAT = INFO)
    i = i + INFO

    ! Return an INFOThermo if deallocation of any of the allocatable variables failed:
    if (i > 0) then
        INFOThermo = 18
    end if

    return

end subroutine ResetThermoParser
