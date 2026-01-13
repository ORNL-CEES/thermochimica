
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TestThermo70.F90
    !> \brief   Spot test - 1500K with 70% Sn, 30% O.
    !> \author  M.H.A. Piro, B.A.T. Breeden
    !
    ! DISCLAIMER
    ! ==========
    ! All of the programming herein is original unless otherwise specified.  Details of contributions to the
    ! programming are given below.
    !
    ! Revisions:
    ! ==========
    !    Date          Programmer          Description of change
    !    ----          ----------          ---------------------
    !    05/14/2013    M.H.A. Piro         Original code
    !    01/10/2022    B.A.T. Breeden      Modification to use Dupin's Zirc Data base with SUBI
    !    08/27/2024    A.E.F. Fitzsimmons  Remodel
    !
    ! Purpose:
    ! ========
    !> \details The purpose of this application test is to ensure that Thermochimica computes the correct
    !!  results for the Zirc Data file at 1500K with 0.7 mols of Sn and 0.3 of O. It also tests mixing term Case #2
    !!  an #4 of the SUBI phase, with the presence of a miscibility gap. Permission was granted from N. Dupin
    !!  to make use of the Zirc DAT file.
    !
    !-------------------------------------------------------------------------------------------------------------

program TestThermo45

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleTesting

    implicit none

    !Init variables
    logical :: lPass
    real(8) :: dGibbsCheck, dHeatCapacityCheck
    integer :: nSpeciesTest
    integer, allocatable :: iSpeciesIndexTest(:)
    real(8), allocatable :: dMolFractionTest(:)

    ! Specify units:
    cInputUnitTemperature  = 'K'
    cInputUnitPressure     = 'atm'
    cInputUnitMass         = 'moles'
    cThermoFileName        = DATA_DIRECTORY // 'ZIRC_no_liq.dat'

    ! Specify values:
    dPressure              = 1D0
    dTemperature           = 1500D0
    dElementMass(8)        = 0.3D0          ! O
    dElementMass(50)       = 0.7D0          ! Sn

    !Init test values
    dGibbsCheck            = -1.80962D05
    dHeatCapacityCheck     = -1.198592D07
    nSpeciesTest           = 4
    iSpeciesIndexTest      = [2, 7, 11, 14] !O2, Sn2, O2Sn, Sn:O
    dMolFractionTest       = [3.0860D-10, 1.5042D-08, 3.5945D-26, 8.0039D-03]
    lPass                  = .FALSE.

    ! Parse the ChemSage data-file:
    call ParseCSDataFile(cThermoFileName)

    ! Call Thermochimica:
    call Thermochimica
    call HeatCapacity

    !Execute the test for mole fractions, gibbs energy and heat capacity
    call testProperties(dGibbsCheck, dHeatCapacityCheck, nSpeciesTest, iSpeciesIndexTest, dMolFractionTest, lPass)

    ! Deallocation
    deallocate(iSpeciesIndexTest, dMolFractionTest)
    
    if (lPass) then
        ! The test passed:
        print *, 'TestThermo45: PASS'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(0)
    else
        ! The test failed.
        print *, 'TestThermo45: FAIL <---'
        ! Reset Thermochimica:
        call ResetThermo
        call EXIT(1)
    end if

end program TestThermo45
