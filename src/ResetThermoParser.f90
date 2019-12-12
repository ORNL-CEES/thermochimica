
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

    integer::   i


    ! Initialize variables:
    i    = 0
    INFO = 0

    ! Deallocate integer arrays from ModuleParseCS:
    if (allocated(nSpeciesPhaseCS)) then
        deallocate (nSpeciesPhaseCS, nGibbsEqSpecies, iPhaseCS, iParticlesPerMoleCS, nParamPhaseCS, &
            iParamPassCS, dStoichSpeciesCS, iRegularParamCS, cRegularParamCS, STAT = INFO)
        i = i + INFO
    end if

    ! Deallocate real arrays from ModuleParseCS:
    if (allocated(dAtomicMass)) then
        deallocate (dAtomicMass, dGibbsCoeffSpeciesTemp, dRegularParamCS, dGibbsMagneticCS, STAT = INFO)
        i = i + INFO
    end if

    ! Deallocate character arrays from ModuleParseCS:
    if (allocated(cElementNameCS)) then
        deallocate (cElementNameCS,cSolnPhaseTypeCS,cSolnPhaseNameCS,cSpeciesNameCS, STAT = INFO)
        i = i + INFO
    end if

    ! Deallocate arrays used for sublattice phases:
    if (allocated(cConstituentNameSUBCS)) then
        deallocate(cConstituentNameSUBCS,iPhaseSublatticeCS,nConstituentSublatticeCS, &
            iConstituentSublatticeCS,dStoichSublatticeCS,nSublatticePhaseCS, &
            iSublatticeElementsCS, nSublatticeElementsCS, dZetaSpeciesCS, &
            dSublatticeChargeCS, dStoichPairsCS, iChemicalGroupCS, STAT = INFO)
        i = i + INFO
    end if

    if (allocated(dGibbsMagneticCS)) then
        deallocate(dGibbsMagneticCS,STAT = INFO)
        i = i + INFO
    end if

    if (allocated(nPairsSROCS)) then
        deallocate(nPairsSROCS,iPairIDCS,dCoordinationNumberCS, STAT = INFO)
        i = i + INFO
    end if

    ! Return an INFOThermo if deallocation of any of the allocatable variables failed:
    if (i > 0) then
        INFOThermo = 18
    end if

    return

end subroutine ResetThermoParser
